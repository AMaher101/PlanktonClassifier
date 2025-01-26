# PlanktonClassifier
Constructs database with accurately classified plankton and attributes using plankton abundance and taxonomic composition data.

## Usage
### Preparing a dataset
1. PlanktonClassifier accepts CSVs with multiheader columns formatted as below:
| (any) |(any)| A[any]S | B[any]S | C[any]S | ... | A[any]B | B[any]B | C[any]B | ... | 
| --- | --- | --- | --- | --- | --- | --- | --- | --- |---|
<br>
The first header indicates site (or "Station") labels.
The second header indicates the date in which samples were taken.

2. PlanktonClassifier requires a format of a skipped line between different phylums (ex. empty row above Diatom, Dinoflagellate, etc.)
