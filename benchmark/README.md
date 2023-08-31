# Benchmarking

## Setup
- Anaconda 3 (recommended): please follow the [guide](https://github.com/Yijia-Xiao/Know2BIO#getting-started)
- Pip requirements: `pip install -r requirements`


## Datasets
- Prepare the data.
  - Dataset folder (e.g. named `$D`) should be put under the ./data/ folder: ./data/$D/
  - Train, valid, test splits should be put under `$D`: ./data/$D/{train,valid,test}
  - In data files, each triple occupied 1 line, in the order of (head, rel, tail), separated by `tab`.

- Preprocess the datasets.
  - Make sure all the subfolders under ./data/ are correctly formatted. Below script will iterate through all the subfolders in it and preprocess the data.
  - `python datasets/process.py`

## Models

### Euclidean Space
- TransE, TransR
- RotE, RefE, AttE, MurE
- CTDecomp, DistMult

### Complex Space
- ComplEx
- RotatE

#### Hyperbolic Space
- RotH
- RefH
- AttH

## Experiments
### Running Examples
We provide `run.sh` as an example to train models on Know2BIO.

### Configurations for benchmarked models
The complete configuration files are available under the `./configs/` folder.

## Acknowledgement
The code is adapted from KGEmb. The original repository is available at: https://github.com/HazyResearch/KGEmb.