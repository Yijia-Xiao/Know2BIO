<p align="center">
  <img src="./assets/logo.svg">
</p>

------------------------------------------------

# Know2BIO

Know2BIO is a comprehensive biomedical knowledge graph benchmark harmonizing heterogeneous database sources.

## Getting Started
### Environment Setup
We recommend using Anaconda3 to manage the environment.
- Install Anaconda3.
- Edit `env.yaml`: set $USER_PATH to user's directory.
- Create `know2bio` environment using `conda env create -f env.yml`.

### Hardware Requirements
- Server: AMD EPYC 7542 Processor (124 cores), 1.73 TB RAM, and 8 NVIDIA A100-SXM4-40GB GPUs.
- Operating system: Ubuntu 20.04 LTS.


## Benchmarking
### Setup
- Python environment: follow the guide in `Environment Setup` Section.
- C++ extensions: compile C++ executable file (`Base.so`)
```bash
cd know2bio/
bash make.sh
```

### Experiments
- To run the experiments, please excute `main.py` script. Arguments are listed below.
```
usage: main.py [-h] [--model_name MODEL_NAME] [--split SPLIT] [--data_path DATA_PATH] [--ckpt_path CKPT_PATH] [--lr LR] [--batch_size BATCH_SIZE] [--train_epochs TRAIN_EPOCHS] [--opt_method OPT_METHOD]

options:
  -h, --help            show this help message and exit
  --model_name MODEL_NAME
                        Name of the model to be evaluated.
  --split SPLIT         The dataset/split to be used.
  --data_path DATA_PATH
                        Path to the datasets: GO, Know2BIO's ontology/instance/aggregate view, etc.
  --ckpt_path CKPT_PATH
                        Path to the model checkpoints storage.
  --lr LR               Learning rate in training process.
  --batch_size BATCH_SIZE
                        Batch size of for modeling training.
  --train_epochs TRAIN_EPOCHS
                        Number of training epochs.
  --opt_method OPT_METHOD
                        Optimization method for model training.
```

- Example: Train TransE model on Know2BIO's aggregate view
```bash
python main.py --model_name transe --split BMKG --lr 1e-5 --batch_size 1024 --train_epochs 3000 --opt_method adam
```

For the `split` argument: Know2BIO's ontology view, instance view and aggregate view are represented by `ONTO`, `INST` and `BMKG` arguments. GO's biological process, cellular component and molecular function are represented by `BP`, `CC` and `MF` arguments.

Code and README for the benchmarking Know2BIO can be found in [**benchmark**](./benchmark).

## Dataset Construction
### Dataset Schema
![Know2BIO Schema](./assets/Know2BIO_Schema.png)

### Data Source and Relationships
![Know2BIO Data Source](./assets/Know2BIO_Data_Source.jpg)

### Usage and Datasheet
Code and README for the construction of Know2BIO can be found in [**dataset**](./dataset).