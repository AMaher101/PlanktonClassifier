
import pandas as pd
import numpy as np

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
    confirmed_before : Array of String
        names of confirmed mixoplankton genuses

    Attributes
    ----------
    lis_cleaned : DataFrame
        cleaned version of lis data
    mixoplankton : DataFrame
        classified plankton as confirmed or unsure mixoplankton based on hard coded rules
    totals : DataFrame
        totals of each genus
    pretty : DataFrame
        classified mixoplankton w/ totals after each genus and original header


    ----------
    HARD CODED RULES FOR CLASSIFICATION
    Assumed all {Ochromonas, } are mixotrophs.

    1. assume everything after "Unknown flagellates" is irrelevant (to be deleted)
    2. diatoms are NOT mixotrophs
    3. remove all "[name]-like" (without genus specified)
    4. remove all "[genus name] spp." AND "[genus name] sp."
    5. check "cysts of"

    Status Key--  
    Confirmed := explicitly in the Mixotroph Database  
    Unsure (sp. in mdb) := genus in Mixotroph Database lists "[genus name] sp." (ex. Ochromonas sp. for Ochromonas danica)  
    Unsure (inexact name):= LIS name is in a longer Mixotroph Database name or vice versa (ex. Chattonella marina in Chattonella marina var. ovata)   
    """
    
    def __init__(self, csv_name, confirmed_before=["Ochromonas"]):
        self.csv_name = csv_name
        self.confirmed_before = confirmed_before

        # import and clean mixoplankton database
        mdb = pd.read_csv(MDB_PATH)
        mdb.columns = mdb.iloc[1]
        mdb = mdb.drop([0, 1]).reset_index(drop=True)
        mdb['Species Name'] = mdb['Species Name'].str.replace(r'sp$', 'sp.', regex=True) # edit so that species ending in "sp" now end in "sp."
        self.mdb = mdb

        # import and clean LIS data and save original header
        lis = pd.read_csv(INPUTS_PATH + csv_name)
        self.orig_header = lis.columns  # save original column headers

        self.lis_cleaned = self.clean_lis(lis)
        self.mixoplankton = self.classify_lis(self.lis_cleaned)
        self.totals = self.calc_totals(self.mixoplankton)
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
        SPECIES_COL = lis.columns.get_loc("Species")
        lis.iloc[:, SPECIES_COL+1:] = lis.iloc[:, SPECIES_COL+1:].replace(",| ", "", regex=True).replace("", np.nan).astype(float)
        lis = lis.fillna(0)

        return lis


    def classify_lis(self, lis):
        # add Status column
        lis = lis.copy()
        lis.insert(0, 'Status', None)

        # store blocks of known mixotroph genuses (given beforehand) 
        known_blocks = []
        for genus in self.confirmed_before:
            ind = lis[lis["Species"].str.contains(genus)].index
            df = lis.iloc[ind]
            known_blocks.append(Block(ind, df))

        # remove based on hard coded rules (NOT RESETTING INDEX IN ORDER TO ADD CONFIRMED_BEFORE GENUSES BACK CORRECTLY)
        lis = lis[lis["Phylum"] != "Diatom"] # remove all diatoms
        lis = lis[~lis["Species"].str.contains("-like")] # remove species ending with "-like"
        lis = lis[~lis["Species"].str.contains("sp.|spp.")]  # remove all sp. / spp.

        # check "cysts of"
        CYSTS_LEN = len("cysts of ")
        cysts_of = lis[lis["Species"].str.contains("cysts of", regex=False)]["Species"].str.slice(CYSTS_LEN)
        filtered = cysts_of.isin(self.mdb['Species Name'])
        lis.loc[filtered[filtered].index, "Status"] = "Confirmed"

        # add back stored blocks of known mixotrophs and mark as Confirmed
        for known_block in known_blocks:
            lis = pd.concat([lis, known_block.df]).sort_index().drop_duplicates()
            lis.loc[known_block.ind, "Status"] = "Confirmed"

        # check if (in none status) direct match and mark all Trues as "Confirmed"
        filtered = lis[lis['Status'].isnull()]["Species"].isin(self.mdb['Species Name'])
        lis.loc[filtered[filtered].index, "Status"] = "Confirmed"
        
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

        # drop all rows with Status = "None"
        lis = lis.dropna(subset=['Status']).reset_index(drop=True)

        return lis


    def calc_totals(self, lis):
        totals = lis.groupby('Phylum', as_index=False, sort=False).sum()

        # empty text-containing columns
        totals = totals.drop("Status", axis=1)
        totals.insert(0, 'Status', "") 
        totals["Genus"] = ""
        totals["Species"] = ""

        # rename to TOTAL "   "
        totals["Phylum"] = totals["Phylum"].str.upper().apply(lambda x: "TOTAL " + x + "S")
        return totals


    # add back multiheader
    def add_multiheader(self, df):
        needed_cols = pd.Series(np.full(len(df.columns) - len(self.orig_header), None))  # create series of "None"'s to be concatted
        original_headers = pd.concat([needed_cols, self.orig_header.to_series()], ignore_index=True)  # concat "None"'s so lines up correctly
    
        df.columns = pd.MultiIndex.from_arrays([original_headers, df.columns])
        return df


    # adds in totals w/ line skips and adds back multiheader
    def make_pretty(self):
        # add in line skips
        totals = self.totals.set_index(self.mixoplankton.groupby(['Phylum']).tail(1).index + 0.1)
        empty_df = pd.DataFrame("", index=self.mixoplankton.groupby(['Phylum']).tail(1).index+0.2, columns=totals.columns)
        totals = pd.concat([totals, empty_df]).sort_index()

        with_totals = pd.concat([self.mixoplankton, totals]).sort_index().reset_index(drop=True)  # add totals w/ line skips
        with_headers = self.add_multiheader(with_totals)

        return with_headers