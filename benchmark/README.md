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
- `type_constrain.txt` (optional): This file provides type constraints for relations. The first line indicates the total number of relations. Each subsequent line specifies the type constraints for a relation. For example, if a relation with ID 1200 has 4 types of head entities (3123, 1034, 58, and 5733) and 4 types of tail entities (12123, 4388, 11087, and 11088), they would be listed accordingly. This file can be generated using the n-n.py in data/ folder.


## Experiments
Follow the guide [here](../README.md#experiments).


## Acknowledgement
The code is adapted from OpenKE, an open toolkit for knowledge embedding. The original repository is available at https://github.com/thunlp/OpenKE/tree/OpenKE-PyTorch.
