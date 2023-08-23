# Benchmarking

## Setup
- Python environment: please follow the [guide](https://github.com/Yijia-Xiao/Know2BIO#getting-started)

- C++ extensions: compile C++ executable file (`Base.so`)
```bash
cd know2bio/
bash make.sh
```


## Data
There are five required files and one optional file to run the benchmark evaluation.
- `entity2id.txt`: This file contains a list of all entities and their corresponding IDs, with each entity and ID on a separate line. The first line of the file indicates the total number of entities.
- `relation2id.txt`: This file contains a list of all relations and their corresponding IDs, with each relation and ID on a separate line. The first line of the file indicates the total number of relations.
- `train2id.txt`: This file is used for training and contains triples in the following format: (e1, e2, rel), where e1 and e2 are entity IDs and rel is the relation ID. The first line of the file specifies the total number of triples for training. The IDs in this file correspond to the entities and relations listed in entity2id.txt and relation2id.txt, respectively.
- `valid2id.txt`: This file is used for validation and has the same format as test2id.txt. It contains triples in the format (e1, e2, rel), with the first line specifying the total number of triples for validation.
- `test2id.txt`: This file is used for testing and follows the same format as train2id.txt. It contains triples in the format (e1, e2, rel), with the first line specifying the total number of triples for testing.
- `type_constrain.txt` (optional): This file provides type constraints for relations. The first line indicates the total number of relations. Each subsequent line specifies the type constraints for a relation. The type constrain file can be generated using the n-n.py in data/ folder.


## Experiments
Follow the guide [here](../README.md#experiments).

## Quick Test
- Example: Train TransE model on Know2BIO's sampled aggregate view
```bash
python main.py --model_name transe --split K2BIO --lr 1e-5 --batch_size 1024 --train_epochs 3000 --opt_method adam
```
- Results:

```
metric:                  MRR             MR              hit@10          hit@3           hit@1
l(raw):                  0.038962        3380.391846     0.100278        0.047883        0.000000
r(raw):                  0.059106        2195.340332     0.157876        0.085010        0.000000
averaged(raw):           0.049034        2787.866211     0.129077        0.066447        0.000000

l(filter):               0.052115        3206.843506     0.138446        0.070784        0.000000
r(filter):               0.074869        2125.733398     0.200555        0.113116        0.000000
averaged(filter):        0.063492        2666.288574     0.169500        0.091950        0.000000
0.169500
mrr, mr, hit10, hit3, hit1 0.06349201500415802 2666.28857421875 0.16950035095214844 0.09195003658533096 0.0
```


## Acknowledgement
The code is adapted from OpenKE, an open toolkit for knowledge embedding. The original repository is available at https://github.com/thunlp/OpenKE/tree/OpenKE-PyTorch.
