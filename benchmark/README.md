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
no type constraint results:                                                                                               
metric:                  MRR             MR              hit@10          hit@3           hit@1                            
l(raw):                  0.146781        2710.333008     0.373404        0.271277        0.002128                         
r(raw):                  0.224690        2179.098877     0.506383        0.425532        0.006383                         
averaged(raw):           0.185735        2444.715820     0.439894        0.348404        0.004255                         
                                                                                                                          
l(filter):               0.163650        2524.646729     0.388298        0.293617        0.006383                         
r(filter):               0.261193        1975.074463     0.569149        0.478723        0.024468                         
averaged(filter):        0.212421        2249.860596     0.478723        0.386170        0.015426                         
0.478723                                                                                                                  
mrr, mr, hit10, hit3, hit1 0.21242141723632812 2249.860595703125 0.478723406791687 0.3861702084541321 0.015425531193614006
```


## Acknowledgement
The code is adapted from OpenKE, an open toolkit for knowledge embedding. The original repository is available at https://github.com/thunlp/OpenKE/tree/OpenKE-PyTorch.
