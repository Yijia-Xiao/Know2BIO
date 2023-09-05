from rdkit import Chem
import networkx as nx
import json
import pickle


def convert_molecule_to_networkx_object(mol):
    G = nx.Graph()

    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(),
                   atomic_num=atom.GetAtomicNum(),
                   is_aromatic=atom.GetIsAromatic(),
                   atom_symbol=atom.GetSymbol())
        
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   bond_type=bond.GetBondType())
        
    return G


def map_mesh_id_to_compound_structure():
    smiles_to_mesh = json.load(open('output/compound2compound/smiles2mesh.json'))
    
    compound_graphs = dict()
    for smiles, mesh_ids in smiles_to_mesh.items(): 
        try:
            compound_rdkit = Chem.MolFromSmiles(smiles)
            compound_graph = convert_molecule_to_networkx_object(compound_rdkit)
        except:
            print(smiles,'did not convert to graph')
            continue
        for mesh_id in mesh_ids:
            compound_graphs['MeSH_Compound:'+mesh_id] = [compound_graph]


    #with open('output/node_features/structures/mesh_id_to_compound_structures.pkl','wb') as fout:
    #    pickle.dump(compound_graphs, fout)


    #with open('output/node_features/structures/mesh_id_to_compound_structures.pkl','rb') as fin:
    #    compound_graphs = pickle.load(fin)
        
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
            continue
        for drugbank_id in drugbank_ids:
            compound_graphs['DrugBank_Compound:'+drugbank_id] = [compound_graph]


    #with open('output/node_features/structures/drugbank_id_to_compound_structures.pkl','wb') as fout:
    #    pickle.dump(compound_graphs, fout)


    #with open('output/node_features/structures/drugbank_id_to_compound_structures.pkl','rb') as fin:
    #    compound_graphs = pickle.load(fin)
        
    return compound_graphs
    
    
if __name__ == '__main__':
    mesh_to_compound_graphs = map_mesh_id_to_compound_structure()
    drugbank_to_compound_graphs = map_drugbank_id_to_compound_structure()
    
    # Export
    compound_graphs = {}
    compound_graphs.update(drugbank_to_compound_graphs)
    print(len(compound_graphs))
    compound_graphs.update(mesh_to_compound_graphs)
    print(len(compound_graphs))
    print(len(mesh_to_compound_graphs))
    with open('output/node_features/structures/compound_id_to_compound_structures.pkl','wb') as fout:
        pickle.dump(compound_graphs, fout)

        
        
'''
Example of graph embeddings of the compounds
import os
os.system('pip install karateclub')
from karateclub import Graph2Vec
import pandas as pd

# Make graph2vec embeddings
model = Graph2Vec()
model.fit(pd.Series(compound_id_to_compound_graphs.values()))
compound_graph2vec = model.get_embedding()

compound_graph2vec_df = pd.DataFrame(compound_graph2vec, index=list(compound_to_compound_graphs.keys()))
compound_graph2vec_dict = {'DrugBank_Compound:'+drug:list(row) for drug,row in compound_graph2vec_df.iterrows()}
with open('output/node_features/structures/compound_id_to_compound_graph2vec_dict.json','w') as fout:
    json.dump(compound_graph2vec_dict, fout)
'''
