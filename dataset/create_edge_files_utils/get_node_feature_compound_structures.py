from rdkit import Chem
import networkx as nx
import json
import pickle


def map_mesh_id_to_compound_structure():
    smiles_to_mesh = json.load(open('output/compound2compound/smiles2mesh.json'))
    
    compound_graphs = dict()
    for smiles, mesh_ids in smiles_to_mesh.items(): 
        try:
            compound_rdkit = Chem.MolFromSmiles(smiles)
            compound_graph = convert_molecule_to_networkx_object(compound_rdkit)
        except:
            print(smiles,'did not convert to graph')
        for mesh_id in mesh_ids:
            compound_graphs[mesh_id] = compound_graph


    with open('output/node_features/structures/mesh_id_to_compound_structures.pkl','wb') as fout:
        pickle.dump(compound_graphs, fout)


    with open('output/node_features/structures/mesh_id_to_compound_structures.pkl','rb') as fin:
        compound_graphs = pickle.load(fin)
        
    return compound_graphs
    
    
def map_drugbank_id_to_compound_structure():
    smiles_to_drugbank = json.load(open('output/compound2compound/smiles2db.json'))

    # Make molecules RDKit objects
    compound_graphs = dict()
    for smiles, drugbank_ids in smiles_to_drugbank.items():    
        try:
            compound_rdkit = Chem.MolFromSmiles(smiles)
            compound_graph = convert_molecule_to_networkx_object(compound_rdkit)
        except:
            print(smiles,'did not convert to graph')        
        for drugbank_id in drugbank_ids:
            compound_graphs[drugbank_id] = compound_graph


    with open('output/node_features/structures/drugbank_id_to_compound_structures.pkl','wb') as fout:
        pickle.dump(compound_graphs, fout)


    with open('output/node_features/structures/drugbank_id_to_compound_structures.pkl','rb') as fin:
        compound_graphs = pickle.load(fin)
        
    return compound_graphs
    
    
if __name__ == '__main__':
    mesh_to_compound_graphs = map_mesh_id_to_compound_structure()
    drugbank_to_compound_graphs = map_drugbank_id_to_compound_structure()



'''
Example of graph embeddings of the compounds
import os
os.system('pip install karateclub')
from karateclub import Graph2Vec
import pandas as pd

# Make graph2vec embeddings
model = Graph2Vec()
model.fit(pd.Series(drugbank_to_compound_graphs.values()))
drugbank_compound_graph2vec = model.get_embedding()

drugbank_compound_graph2vec_df = pd.DataFrame(drugbank_compound_graph2vec, index=list(drugbank_to_compound_graphs.keys()))
drugbank_compound_graph2vec_dict = {'DrugBank_Compound:'+drug:list(row) for drug,row in drugbank_compound_graph2vec_df.iterrows()}
with open('output/node_features/structures/drugbank_compound_graph2vec_dict.json','w') as fout:
    json.dump(drugbank_compound_graph2vec_dict, fout)
'''
