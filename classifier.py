
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import re

pd.set_option("future.no_silent_downcasting", True)

import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

from constants import *


class Block:
    """
    A class to represent a portion of a dataframe.

    Attributes
    ----------
    ind : list/Index object
        indicices of where the block is located
    df : DataFrame
        contents of the block
    """
    def __init__(self, ind, df):
        self.ind = ind
        self.df = df


class Classifier:
    """
    Classifies LIS plankton as mixotrophs based on Mixoplankton Database.

    Parameters
    ----------
    csv_name : String
        name of the csv to classify
        * csv should have skipped line between different phylums (ex. empty row above Diatom, Dinoflagellate, etc.)
    confirmed_genus_before : Array of String
        names of confirmed mixoplankton genuses
    confirmed_species_before : Array of String
        names of confirmed mixoplankton genuses

    Attributes
    ----------
    lis_cleaned : DataFrame
        cleaned version of lis data
    all_classified : DataFrame
        all classified plankton as mixoplankton (yes, unsure, or no) based on hard coded rules
    mixoplankton : DataFrame
        only the mixoplankton (yes or unsure)
    mixoplankton_with_header : DataFrame
        mixoplankton with header
    totals : DataFrame
        totals of each genus
    pretty : DataFrame
        classified mixoplankton w/ totals after each genus and header
        
    ----------
    HARD CODED RULES FOR CLASSIFICATION
    Assumed all {Ochromonas, } are mixotrophs.

    1. assume everything after "Unknown flagellates" is irrelevant (to be deleted)
    2. diatoms are NOT mixotrophs
    3. remove all "[name]-like" (without genus specified)
    4. remove all "[genus name] spp." AND "[genus name] sp."
    5. check "cysts of"

    Status Key--  
    Mixoplankton := explicitly in the Mixotroph Database  
    Unsure (sp. in mdb) := genus in Mixotroph Database lists "[genus name] sp." (ex. Ochromonas sp. for Ochromonas danica)  
    Unsure (inexact name):= LIS name is in a longer Mixotroph Database name or vice versa (ex. Chattonella marina in Chattonella marina var. ovata)   
    Zooplankton := manually confirmed as Zooplankton species in "confirmed_zoo_genus_before" or "confirmed_zoo_species_before"
    Phytoplankton := not in Mixoplankton Database or priorly confirmed as Zooplankton
    """
    
    def __init__(self, csv_name, confirmed_mixo_genus_before=["Ochromonas"], confirmed_mixo_species_before=['Chattonella marina'], confirmed_zoo_genus_before=["Protoperidinium"], confirmed_zoo_species_before=[]):
        self.csv_name = csv_name
        self.year = re.search(r'\b(\d{4})\b', self.csv_name).group(1) if re.search(r'\b(\d{4})\b', self.csv_name) else ""
        self.confirmed_mixo_genus_before = confirmed_mixo_genus_before
        self.confirmed_mixo_species_before = confirmed_mixo_species_before
        self.confirmed_zoo_genus_before = confirmed_zoo_genus_before
        self.confirmed_zoo_species_before = confirmed_zoo_species_before

        # import and clean mixoplankton database
        mdb = pd.read_csv(MDB_PATH)
        mdb.columns = mdb.iloc[1]
        mdb = mdb.drop([0, 1]).reset_index(drop=True)
        mdb['Species Name'] = mdb['Species Name'].str.replace(r'sp$', 'sp.', regex=True) # edit so that species ending in "sp" now end in "sp."
        self.mdb = mdb

        # calculate average volumes for each plankton size class
        mdb_volume = mdb.loc[:, ['size class', 'L (μm)', 'W (μm) or diameter (μm)']]
        mdb_volume[['L (μm)', 'W (μm) or diameter (μm)']] = mdb_volume[['L (μm)', 'W (μm) or diameter (μm)']].replace(r'[~]?(\d+)-(\d+)', lambda x: str((float(x.group(1)) + float(x.group(2))) / 2), regex=True)
        mdb_volume[['L (μm)', 'W (μm) or diameter (μm)']] = mdb_volume[['L (μm)', 'W (μm) or diameter (μm)']].replace(r'[~≤]', '', regex=True)
        mdb_volume[['L (μm)', 'W (μm) or diameter (μm)']] = mdb_volume[['L (μm)', 'W (μm) or diameter (μm)']].apply(pd.to_numeric, errors='coerce')
        mdb_volume = mdb_volume.dropna(subset=['L (μm)', 'W (μm) or diameter (μm)'])
        mdb_volume['Volume'] = np.where(mdb_volume['L (μm)'] == mdb_volume['W (μm) or diameter (μm)'], (4/3) * np.pi * (mdb_volume['W (μm) or diameter (μm)'] / 2) ** 3, np.nan)
        mdb_volume['Volume'] = np.where((mdb_volume['L (μm)'].notna()) & (mdb_volume['W (μm) or diameter (μm)'].notna()) & (mdb_volume['L (μm)'] != mdb_volume['W (μm) or diameter (μm)']), (4/3) * np.pi * ((mdb_volume['W (μm) or diameter (μm)'] / 2) ** 2) * (mdb_volume['L (μm)'] / 2), mdb_volume['Volume'])
        mdb_volume = mdb_volume.drop(['L (μm)', 'W (μm) or diameter (μm)'], axis=1).groupby('size class', as_index=False).mean()
        self.mdb_volume = mdb_volume

        # import and clean LIS data and save original header
        lis = pd.read_csv(INPUTS_PATH + csv_name)
        self.orig_header = lis.columns  # save original column headers

        self.lis_cleaned = self.clean_lis(lis)
        self.all_classified = self.classify_lis(self.lis_cleaned)
        self.mixoplankton = self.only_mixotrophs(self.all_classified)
        self.mixoplankton_with_header = self.add_multiheader(self.mixoplankton)

        self.totals = self.calc_totals(self.all_classified)

        self.pretty = self.make_pretty()
        

    def clean_lis(self, lis):
        phylum_ind = lis[lis.iloc[:, 0] == "Phylum"].index[0]
        lis.columns = lis.iloc[phylum_ind]  # reset column headers
        lis = lis.iloc[phylum_ind+2:].reset_index(drop=True)

        # remove rows after unknown flagellates
        unknown_flagellates_ind = lis[lis["Phylum"] == "Unknown flagellates"].index[0] 
        lis = lis.iloc[:unknown_flagellates_ind]
        lis = lis.iloc[:lis.last_valid_index()+1]  # remove trailing nan rows

        # remove rows that contain "TOTAL"
        lis = lis[~lis["Phylum"].str.contains("TOTAL", na=False)].reset_index(drop=True)  

        # construct correct phylum column
        actual_phylum_ind = lis[lis["Species"].isna() & lis["Phylum"].isna()].index + 1
        lis = lis.rename(columns={"Phylum": "Genus"}) # rename phylum column to genus
        lis.insert(0, 'Phylum', lis["Genus"].iloc[actual_phylum_ind])  # reconstruct phylum column
        lis['Phylum'] = lis['Phylum'].ffill()  # forwardfill phylum
        
        lis['Genus'] = lis['Species'].str.split().str[0]  # fill genus using first word of species name

        lis = lis.dropna(subset=['Species']).reset_index(drop=True) # delete rows with na in Species column

        # ensure numerical values are floats and not strings
        lis = lis.fillna(0)
        SPECIES_COL = lis.columns.get_loc("Species")
        lis.iloc[:, SPECIES_COL+1:] = lis.iloc[:, SPECIES_COL+1:].replace(",| ", "", regex=True).replace("", 0).astype(float).astype(int)
        
        # add totals for each row
        lis['Totals'] = lis.loc[:, ~lis.columns.isin(['Status', 'Phylum', 'Genus', 'Species'])].sum(axis=1)
        lis = pd.concat([lis.iloc[:, :3], lis.iloc[:, -1:], lis.iloc[:, 3:-1]], axis=1)

        return lis


    def classify_lis(self, lis):

        # add Status column
        lis = lis.copy()
        lis.insert(0, 'Status', None)

        # store blocks of known mixotroph genuses and species (given beforehand) 
        known_mixo_blocks = []
        for genus in self.confirmed_mixo_genus_before:
            ind = lis[lis["Species"].str.contains(genus)].index
            df = lis.iloc[ind]
            known_mixo_blocks.append(Block(ind, df))
        for species in self.confirmed_mixo_species_before:
            ind = lis[lis["Species"].str.contains(species)].index
            df = lis.iloc[ind]
            known_mixo_blocks.append(Block(ind, df))

        # store blocks of known zooplankton genuses and species (given beforehand)
        known_zoo_blocks = []
        for genus in self.confirmed_zoo_genus_before:
            ind = lis[lis["Species"].str.contains(genus)].index
            df = lis.iloc[ind]
            known_zoo_blocks.append(Block(ind, df))
        for species in self.confirmed_zoo_species_before:
            ind = lis[lis["Species"].str.contains(species)].index
            df = lis.iloc[ind]
            known_zoo_blocks.append(Block(ind, df))

        # remove based on hard coded rules (NOT RESETTING INDEX IN ORDER TO ADD CONFIRMED_BEFORE GENUSES BACK CORRECTLY)
        lis = lis[~lis["Species"].str.contains("unknown|other|cysts")]
        lis = lis[~lis["Species"].str.contains("-like")] # remove species ending with "-like"
        lis = lis[~lis["Species"].str.contains("sp.|spp.")]  # remove all sp. / spp.

        # add back stored blocks of known mixoplankton and mark as "Mixoplankton"
        for known_block in known_mixo_blocks:
            lis = pd.concat([lis, known_block.df]).sort_index().drop_duplicates()
            lis.loc[known_block.ind, "Status"] = "Mixoplankton"

        # add back stored blocks of known zooplankton and mark as "Zooplankton"
        for known_block in known_zoo_blocks:
            lis = pd.concat([lis, known_block.df]).sort_index().drop_duplicates()
            lis.loc[known_block.ind, "Status"] = "Zooplankton"

        # check if (in none status) direct match and mark all Trues as "Yes"
        filtered = lis[lis['Status'].isnull()]["Species"].isin(self.mdb['Species Name'])
        lis.loc[filtered[filtered].index, "Status"] = "Mixoplankton"
        
        # check (in remaining none status) if the genus has sp. and mark all Trues as "Unsure (sp. in mdb)"
        genus_to_check = lis[lis['Status'].isnull()]['Species'].str.split().str[0].drop_duplicates() + " sp."
        filtered = genus_to_check.isin(self.mdb['Species Name'])
        lis.loc[filtered[filtered].index, "Status"] = "Unsure (sp. mdb)"

        # check (in remaining none status) if the name is contained in the mdb and vice versa and mark all Trues as "Unsure (inexact name)"
        filtered = lis[lis['Status'].isnull()]["Species"].apply(lambda x: self.mdb["Species Name"].str.contains(x, regex=False).any())
        lis.loc[filtered[filtered].index, "Status"] = "Unsure (inexact name)"
        
        pattern = '|'.join(self.mdb['Species Name'])
        filtered = lis[lis['Status'].isnull()]["Species"].str.contains(pattern, regex=True)
        lis.loc[filtered[filtered].index, "Status"] = "Unsure (inexact name)"
        
        # replace None's in Status with No's 
        lis["Status"] = lis["Status"].replace(np.nan, 'Phytoplankton')

        return lis.reset_index(drop=True)
    

    def only_mixotrophs(self, lis):
        # drop all rows with Status != "Mixoplankton"
        lis = lis[lis["Status"] == "Mixoplankton"].reset_index(drop=True)

       # merge additional columns from mdb
        lis = pd.merge(lis, self.mdb[['Species Name', 'MFT', 'Evidence of mixoplankton activity', 'size class', 'L (μm)', 'W (μm) or diameter (μm)']], left_on='Species', right_on='Species Name', how='left').drop(columns=['Species Name']).reset_index(drop=True) 
        lis = pd.concat([lis.iloc[:, :4], lis.iloc[:, -5:], lis.iloc[:, 4:-5]], axis=1)
        lis[['MFT', 'Evidence of mixoplankton activity', 'size class', 'L (μm)', 'W (μm) or diameter (μm)']] = lis[['MFT', 'Evidence of mixoplankton activity', 'size class', 'L (μm)', 'W (μm) or diameter (μm)']].fillna("")

        # clean length/width columns and calculate volumes
        lis[['L (μm)', 'W (μm) or diameter (μm)']] = lis[['L (μm)', 'W (μm) or diameter (μm)']].replace(r'[^\d\-]', '', regex=True).apply(lambda col: col.str.split('-').apply(lambda x: (float(x[0]) + float(x[1])) / 2 if len(x) == 2 else float(x[0]) if x[0] else None))
        lis['Volume'] = np.where(lis['L (μm)'] == lis['W (μm) or diameter (μm)'], (4/3) * np.pi * (lis['W (μm) or diameter (μm)'] / 2) ** 3, np.nan)
        lis['Volume'] = np.where((lis['L (μm)'].notna()) & (lis['W (μm) or diameter (μm)'].notna()) & (lis['L (μm)'] != lis['W (μm) or diameter (μm)']), (4/3) * np.pi * ((lis['W (μm) or diameter (μm)'] / 2) ** 2) * (lis['L (μm)'] / 2), lis['Volume'])

        # manually add values for additional columns from mdb for confirmed_befores
        lis.loc[(lis['Genus'] == 'Ochromonas'), ['MFT', 'Evidence of mixoplankton activity', 'size class']] = ['CM', 'uptake of eubacteria', 'nano']
        lis.loc[(lis['Species'] == 'Chattonella marina'), ['MFT', 'Evidence of mixoplankton activity', 'size class']] = ['CM', 'uptake of eubacteria', 'micro']

        # fill unknown volumes with averages from mdb_volume based on size class and convert
        lis['Volume'] = lis['Volume'].fillna(lis.merge(self.mdb_volume, on='size class', how='left')['Volume_y']) 
        lis['Total Biomass (pgC)'] = (((lis['Volume'])**0.939) * 0.216) * lis['Totals']
        lis = lis.rename(columns={'Volume':'Volume (µm³/cell)'})

        # move volume and biomass column near front
        lis = pd.concat([lis.iloc[:, :9], lis.iloc[:, -2:], lis.iloc[:, 9:-2]], axis=1)

        # remove Status, length, and width columns
        lis = lis.drop(columns=["Status", 'L (μm)', 'W (μm) or diameter (μm)'])

        return lis

    def calc_totals(self, lis):
        totals = lis.groupby('Phylum', as_index=False, sort=False).sum()

        # empty text-containing columns
        totals = totals.drop(columns=["Status"], axis=1, errors='ignore')
        totals["Genus"] = ""
        totals["Species"] = ""
        totals["MFT"] = ""
        totals["Evidence of mixoplankton activity"] = ""
        totals["size class"] = ""
        totals['Volume (µm³/cell)'] = ""

        # rename to TOTAL "   "
        totals["Phylum"] = totals["Phylum"].str.upper().apply(lambda x: "TOTAL " + x + "S")
        return totals


    # add back multiheader with sorted dates and stations
    def add_multiheader(self, lis):
        lis = lis.copy()
        needed_cols = pd.Series(np.full(len(lis.columns) - len(self.orig_header), None))  # create series of "None"'s to be concatted
        original_headers = pd.concat([needed_cols, self.orig_header.to_series()], ignore_index=True)  # concat "None"'s so lines up correctly
        lis.columns = pd.MultiIndex.from_arrays([original_headers, lis.columns])

        # Isolate Station/Date columns
        species_columns = lis.columns.get_level_values(1).isin(['Status', 'Phylum', 'Genus', 'Species', 'MFT', 'Evidence of mixoplankton activity', 'size class', 'Totals', 'Volume (µm³/cell)', 'Total Biomass (pgC)'])
        removed_columns = lis.loc[:, species_columns].copy()

        lis = lis.loc[:, ~species_columns]
        lis.columns = pd.MultiIndex.from_tuples(lis.columns, names=['Station', 'Date'])

        # Extract the 'Date' level and convert to Series
        date_level = lis.columns.get_level_values('Date').to_series()

        # Extract month numbers and day numbers from date strings
        date_parts = date_level.str.extract(r'(\d{1,2})/(\d{1,2})/(\d{2})', expand=False)
        month_numbers = pd.to_numeric(date_parts[0], errors='coerce')
        day_numbers = pd.to_numeric(date_parts[1], errors='coerce')

        # Map for month names
        month_dict = {1: 'January', 2: 'February', 3: 'March', 4: 'April', 5: 'May', 6: 'June', 7: 'July', 8: 'August', 9: 'September', 10: 'October', 11: 'November', 12: 'December'}

        # Adjust month if day >= 26
        month_numbers = month_numbers + (day_numbers >= 26).astype(int)
        month_names = month_numbers.map(month_dict).fillna(date_level.str.extract(r'(\b\w+\b)')[0]).fillna('Unknown')

        # Normalize station names and create new MultiIndex
        station_level = lis.columns.get_level_values('Station').str.split(' ').str[0]
        lis.columns = pd.MultiIndex.from_arrays([month_names, station_level, date_level], names=['Month', 'Station', 'Date'])

        # Add back initially removed columns
        needed_cols = pd.Series(np.full(len(removed_columns.columns), None))
        removed_columns.columns = pd.MultiIndex.from_arrays([needed_cols, removed_columns.columns.get_level_values(0), removed_columns.columns.get_level_values(1)])
        lis = pd.concat([removed_columns, lis], axis=1)

        return lis


    # adds in totals w/ line skips and adds back multiheader
    def make_pretty(self):
        # add in line skips
        totals = self.calc_totals(self.mixoplankton).set_index(self.mixoplankton.groupby(['Phylum']).tail(1).index + 0.1)
        empty_df = pd.DataFrame("", index=self.mixoplankton.groupby(['Phylum']).tail(1).index+0.2, columns=totals.columns)
        totals = pd.concat([totals, empty_df]).sort_index()

        with_totals = pd.concat([self.mixoplankton, totals]).sort_index().reset_index(drop=True)  # add totals w/ line skips
        
        with_headers = self.add_multiheader(with_totals)

        return with_headers