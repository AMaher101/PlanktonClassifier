# PlanktonClassifier
Constructs database with accurately classified plankton and attributes using plankton abundance and taxonomic composition data.

## Usage
### Preparing a dataset
1. PlanktonClassifier accepts CSVs with multiheader columns formatted as below:

<table>
  <tr>
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
  <tr>
    <th>(any)</th>
    <th>(any)</th>
    <td></td>
    <td></td>
    <td></td>
    <td></td>
    <td></td>
    <td></td>
    <td></td>
    <td></td>
  </tr>
</table>

The first header indicates site (or "Station") labels.
The second header indicates the date in which samples were taken.

2. PlanktonClassifier requires a format of a skipped line between different phylums (ex. empty row above Diatom, Dinoflagellate, etc.)
