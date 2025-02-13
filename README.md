# PlanktonClassifier
Constructs formatted dataframe with accurately classified plankton and attributes using plankton abundance and taxonomic composition data.

## Usage
### Preparing a dataset
PlanktonClassifier accepts CSVs with multiheader columns formatted as below:

<table>
    <th></th>
    <th></th>
    <th>A[any]S</th>
    <th>B[any]S</th>
    <th>C[any]S</th>
    <th>...</th>
    <th>A[any]B</th>
    <th>B[any]B</th>
    <th>C[any]B</th>
    <th>...</th>
  </tr>
    <th>Phylum</th>
    <th>Species</th>
    <th>(Date)</th>
    <th>(Date)</th>
    <th>(Date)</th>
    <th>...</th>
    <th>(Date)</th>
    <th>(Date)</th>
    <th>(Date)</th>
    <th>...</th>
  </tr>
    <th>Genus</th>
    <th></th>
    <th></th>
    <th></th>
    <th></th>
    <th></th>
    <th></th>
    <th></th>
    <th></th>
    <th></th>
  </tr>
</table>

The first header indicates site (or "Station") labels. (S = "Surface," B = "Bottom").

The second header indicates Phylum and Species distinctions in the first two columns, and the date in which samples were taken for the following columns. (All samples taken should be within the bounds of a single calendar year).

The third header is necessary for indication of genus. In this column, an initial empty row is reserved for the Phylum name with following rows specifiying Genus name—the adjacent column in the same row will specify Species name.

PlanktonClassifier requires a format of a skipped line between different phylums (ex. empty row above Diatom, Dinoflagellate, etc.)
Example datasets are provided in the inputs folder for reference.

### Classifying a plankton dataset
1. Place your appropriately formatted CSV file in the _inputs_ folder.
2. In _run.py_, ensure all inputs listed in _csvs_ are correct.
3. Run _run.py_.

Excel files of the classified databases are saved to `PlanktonClassifier/outputs`.

#### Useful variables to call from Classifier
- `all_classified`: Cleans dataset to specify genus, phylum, and total cell counts for every species, classifies plankton based on trophic strategy, and removes multiheader for simpler analysis.
- `mixoplankton`: Provides dataset with only mixoplankton, specifying phylum, genus, species, functional type, size class, evidence of mixoplankton activity and calculating total cell count, volume ((µm³/cell), and total biomass (pgC) for every species. 
- `mixoplankton_with_header`: Provides dataset with specifications as with 'mixoplankton' in addition to multiheader with stations and organization based on month.
- `totals`: Provides total cell counts for each phylum throughout an entire calendar year as well as at each station at every date.
- `pretty`: Provides dataset with specifications as with 'mixoplankton_with_header' in addition to total cell count and total biomass (pgC) for every phylum in a full calendar year and at each station on every date.

## Analysis
The folder labeled `Analysis` contains notebooks analyzing ecological communities of plankton as well as elements of their physicochemical environments throughout the Long Island Sound from 2014-2021. 

### Findings
Elevated cell counts of phytoplankton and mixoplankton and drops in dissolved oxygen levels observed in waters near urban environments indicate a connection to the development of acute hypoxic "dead zones" found in parts of the Long Island Sound. 

These results suggest that the release of excess nitrogen from wastewater treatment plants surrounding urban cities, such as New York City and New Haven, are likely providing an abundant amount of a key nutrient conducive to such heightened growth and proliferation of phytoplankton and mixoplankton detected, resulting in large algal blooms. Thus, when bacteria subsequently begin the process of decomposition, they presumably consume inordinate amounts of dissolved oxygen faster than it can be replenished, creating the hypoxic zones that have significantly disrupted ecosystem balance.   

This research emphasizes the urgency with which we must address these interconnected environmental challenges to protect the health and sustainability of the Long Island Sound's ecosystem. 

---
For questions, contact Aiden Maher at aidenhmaher@gmail.com 
