# PlanktonClassifier
Constructs database with accurately classified plankton and attributes using plankton abundance and taxonomic composition data.

## Usage
### Preparing a dataset
1. PlanktonClassifier accepts CSVs with multiheader columns formatted as below:

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

2. PlanktonClassifier requires a format of a skipped line between different phylums (ex. empty row above Diatom, Dinoflagellate, etc.)
3. Example datasets are provided in the inputs folder for reference.

### Classifying a plankton dataset
1. Place your appropriately formatted CSV file in the _inputs_ folder.
2. In _run.ipynb_, ensure all inputs listed in _csvs_ are correct.
3. Run _run.ipynb_.
