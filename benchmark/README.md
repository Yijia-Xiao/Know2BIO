# Benchmarking

## Setup
- Anaconda 3 (recommended): `conda env create -f env.yml` (please follow the [guide](https://github.com/Yijia-Xiao/Know2BIO#getting-started))
- Pip requirements: `pip install -r requirements`


## Datasets
- Prepare the data.
  - Dataset folder (e.g. named `$D`) should be put under the ./data/ folder: ./data/$D/
  - Train, valid, test splits should be put under `$D`: ./data/$D/{train,valid,test}
  - In data files, each triple occupied 1 line, in the order of (head, rel, tail), separated by `tab`.

- Preprocess the datasets.
  - Make sure all the subfolders under ./data/ are correctly formatted. Below script will iterate through all the subfolders in it and preprocess the data.
  - Preprocess the dataset: `DATA_PATH=data python datasets/process.py --dataset $D`

## Models Overview
The benchmarked models can be categorized into Euclidean, Complex, and Hyperbolic spaces.

**Euclidean Space** models can be further classified into distance similarity based: TransE, TransR, RotE, RefE, AttE, MurE; dot similarity based (i.e. semantic matching) CP (CTDecomp), DistMult.

**Complex Space** models embed entities and relations with complex space: ComplEx and RotatE models.

**Hyperbolic Space** models make use of Hyperbolic Space's properties (excels in the modeling of hierarchical knowledge graphs): RotH, RefH, AttH.

## Experiments
### Script usage
```
usage: run.py [-h] [--dataset {ontology,instance,whole,FB15K,WN,WN18RR,FB237,YAGO3-10}]
              [--model {TransE,TransR,DistMult,CP,MurE,RotE,RefE,AttE,RotH,RefH,AttH,ComplEx,RotatE}] [--regularizer {N3,F2}] [--reg REG]
              [--optimizer {Adagrad,Adam,SparseAdam}] [--max_epochs MAX_EPOCHS] [--patience PATIENCE] [--valid VALID] [--rank RANK] [--batch_size BATCH_SIZE] [--neg_sample_size NEG_SAMPLE_SIZE]
              [--init_size INIT_SIZE] [--learning_rate LEARNING_RATE]

Knowledge Graph Embedding

options:
  -h, --help            show this help message and exit
  --dataset {ontology,instance,whole,FB15K,WN,WN18RR,FB237,YAGO3-10}
                        Knowledge Graph dataset
  --model {TransE,TransR,DistMult,CP,MurE,RotE,RefE,AttE,RotH,RefH,AttH,ComplEx,RotatE}
                        Knowledge Graph embedding model
  --optimizer {Adagrad,Adam,SparseAdam}
                        Optimizer
  --max_epochs MAX_EPOCHS
                        Maximum number of epochs to train for
  --patience PATIENCE   Number of epochs before early stopping
  --valid VALID         Number of epochs before validation
  --rank RANK           Embedding dimension
  --batch_size BATCH_SIZE
                        Batch size
  --neg_sample_size NEG_SAMPLE_SIZE
                        Negative sample size, -1 to not use negative sampling
  --dropout DROPOUT     Dropout rate
  --init_size INIT_SIZE
                        Initial embeddings' scale
  --learning_rate LEARNING_RATE
                        Learning rate
```

### Running Examples
We provide `run.sh` as an example to train models on Know2BIO. For example:

```bash
CUDA_VISIBLE_DEVICES=0 python main.py --model TransE --dataset whole --valid 10 --patience 5 --rank 512 --neg_sample_size 150 --optimizer Adam --learning_rate 0.001
```

### Reproducibility for benchmarked models
The complete configuration files are available under the `./configs/` folder. Subfolders (ontology, instance, whole) hold configuration files for all the benckmarked models on the corresponding views.
```
.
├── instance
│   ├── AttE-config.json
│   ├── ComplEx-config.json
│   ...
│   ├── RotatE-config.json
│   └── TransE-config.json
├── ontology
│   ├── AttE-config.json
│   ├── ComplEx-config.json
│   ...
│   ├── RotatE-config.json
│   └── TransE-config.json
└── whole
    ├── AttE-config.json
    ├── ComplEx-config.json
    ...
    ├── RotatE-config.json
    └── TransE-config.json
```

## Acknowledgement
The code is adapted from KGEmb. The original repository is available at: https://github.com/HazyResearch/KGEmb. We would like to express our acknowledgement to the authors.