# Know2BIO Knowledge Graph Preparation

Scripts for constructing the knowledge graph for Know2BIO dataset. 
Run `prepare_kgs.py` to construct the knowledge graph after preparing the Know2BIO Dataset.
Set up the environment to run the script, using `pip install -r requirements.txt`

The script takes as input the directory where Know2BIO is downloaded, the output directory, and a list of all input files separated by which kg it is a part of. An example is included in this folder, in `ont_bridge_inst_lists.txt`


`python prepare_kgs.py <input_directory> <output_directory> <input_lists>`


`example: python prepare_kgs.py Know2BIO_files Know2BIO_kgs ont_bridge_inst_lists.txt`


Run `kg_sampler.py` to sample a subset of the knowledge graph for testing purposes.

