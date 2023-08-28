## Protein Structure

The below scripts downloads and parses protein 3d structure from EMBL-EBI / Alpha Fold. The proteins found within `multi-modal_data` with valid protein sequences are automatically used as input to the `downloader.py` script. The structure information and errors from prediction are automatically downloaded into the `protein_3d_structure` and `prediction_errors` folder, respectively. These 3d structure are converted to a binary contact map, indicating pairs of amino acid residues within close contact with one another, represented as a .npy file.

To run the script, perform the following commands:
```bash
python downloader.py # download the files
(optional) python checker.py # redownload failed files, investigate missing structures
python pdb2cont_map.py # convert to binary map
```
