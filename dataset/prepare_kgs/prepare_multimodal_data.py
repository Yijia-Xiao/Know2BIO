import json
import pandas as pd
import numpy as np
import random
import os


def read_kg(kg_file):
    '''
    Reads and returns kg_file as Pandas DataFrame
    '''
    df = pd.read_csv(kg_file,sep="\t",header=None)
    df.columns = ['h','r','t']
    return df


def remap_protein2sequence(protein2sequence):
    '''
    Remaps protein2sequence to include UniProt: name space prefix
    '''
    return {"UniProt:"+k:v for k,v in protein2sequence.items()}
    
    
def remap_db2smiles(db2smiles):
    '''
    Remaps db2smiles to include DrugBank_Compound: name space prefix
    '''
    return {"DrugBank_Compound:"+k:v for k,v in db2smiles.items()}


def remap_mesh2smiles(mesh2smiles):
    '''
    Remaps mesh2smiles to include MeSH_Compound: name space prefix
    '''
    return {"MeSH_Compound:"+k:v for k,v in mesh2smiles.items()}


def remap_entrez2seq(entrez2seq):
    '''
    Remaps entrez2seq to include Entrez: name space prefix
    '''
    return {"Entrez:"+k:v for k,v in entrez2seq.items()}


def combine_compound(db2smiles,mesh2smiles):
    '''
    Combines db2smiles and mesh2smile dictionaries
    '''
    combined_dict = {**dict1, **dict2}
    return combined_dict


def remap_id2name_dicts(id2name_dicts, remap_namespaces=None):
    '''
    Remaps id2name_dicts to include namespace prefixes. 
    Also manually aligns namespace prefixes if included in remap_namespaces dict
    '''
    
    if remap_namespaces:
        print("Remapping namespaces")
        
        remapped_id2name_dicts = {}
        for namespace,nodes_to_desc in id2name_dicts.items():
            if namespace == 'GO':
                # split GO into 3 categories
                bp_dict = {k.split(":")[1]:v for k,v in nodes_to_desc.items() if "biological_process" in k}
                mf_dict = {k.split(":")[1]:v for k,v in nodes_to_desc.items() if "molecular_function" in k}
                cc_dict = {k.split(":")[1]:v for k,v in nodes_to_desc.items() if "cellular_component" in k}
                remapped_id2name_dicts['biological_process'] = bp_dict
                remapped_id2name_dicts['molecular_function'] = mf_dict
                remapped_id2name_dicts['cellular_component'] = cc_dict
            elif namespace == 'MeSH_ID':
                # split into MeSH_Disease and MeSH_Compound
                md_dict = {k:v for k,v in nodes_to_desc.items() if k not in remap_namespaces['MeSH_Compound']}
                mc_dict = {k:v for k,v in nodes_to_desc.items() if k in remap_namespaces['MeSH_Compound']}
                
                remapped_id2name_dicts[remap_namespaces[namespace]] = md_dict
                remapped_id2name_dicts['MeSH_Compound'] = mc_dict
            else:
                remapped_id2name_dicts[remap_namespaces[namespace]] = nodes_to_desc
        id2name_dicts = remapped_id2name_dicts
    
    # Collapse namespaces
    print("Collapsing namespaces")
    node2desc = {}
    for namespace, nodes_to_desc in id2name_dicts.items():
        
        for node, desc in nodes_to_desc.items():
            
            node2desc[namespace+":"+node]=desc
            
    print("%d node descriptions"%(len(node2desc)))
    
    return node2desc
    
    
def get_protein_contact_maps(directory, debug=False):
    
    file_mapping = {}
    for filename in os.listdir(directory):
        if filename.endswith('.cont_map.npy'):
            uid = filename.split('.')[0]
            full_path = os.path.join(directory, filename)
            file_mapping[uid] = full_path
            
    if debug:
        print("%d protein 3d structure contact map files"%len(file_mapping))
    return file_mapping


def get_multimodal_data_files(directory, contact_map_directory=None, node_type_to_nodes = None, debug=False):
    # For reformatting changed namespaces
    remap_namespaces_dict = {'DrugBank':'DrugBank_Compound',
                         'MeSH_ID':'MeSH_Disease', 
                         'MeSH_Tree':'MeSH_Tree_Disease', 
                         'Protein':'UniProt', 
                         'Gene':'Entrez', 
                         'GO':None, 
                         'SMPDB':'SMPDB', 
                         'KEGG':"KEGG_Pathway", 
                         'Reactome_Pathway':'Reactome_Pathway', 
                         'Reactome_Reaction':'Reactome_Reaction',
                         'MeSH_Compound':node_type_to_nodes['MeSH_Compound'] # this wasn't defined previously
    }
    
    # store parsed files in file_mapping
    file_mapping = {}
    
    # store multiple descriptions in id2desc
    id2desc = {}
    
    # read through files in multi-modal_data
    for filename in os.listdir(directory):
        
        if filename == "id2name_dicts.json":
            with open(os.path.join(directory,filename)) as json_file:
                id2name = json.load(json_file)
                file_mapping['id2name'] = id2name 
                
                if debug:
                    print("%d entity name categories"%(len(id2name)))
                    for category, entities in id2name.items():
                        print("%s has %d entities"%(category, len(entities)))
                        
        if filename == "protein2sequence.json":
            with open(os.path.join(directory,filename)) as json_file:
                protein2sequence = json.load(json_file)
                remapped_protein2sequence = remap_protein2sequence(protein2sequence)
                file_mapping['protein2sequence'] = remapped_protein2sequence
                if debug:
                    print("%d protein sequences"%(len(protein2sequence)))
      
        if filename == "db2smiles.json":
            with open(os.path.join(directory,filename)) as json_file:
                db2smiles = json.load(json_file)
                remapped_db2smiles = remap_db2smiles(db2smiles)
                file_mapping['db2smiles'] = remapped_db2smiles 
                
                if debug:
                    print("%d compound sequences from DrugBank"%(len(db2smiles)))
                    
        if filename == "mesh2smiles.json":
            with open(os.path.join(directory,filename)) as json_file:
                mesh2smiles = json.load(json_file)
                remapped_mesh2smiles = remap_mesh2smiles(mesh2smiles)
                file_mapping['mesh2smiles'] = remapped_mesh2smiles  
                
                if debug:
                    print("%d compound sequences from MeSH"%(len(mesh2smiles)))         
    
        if filename == "entrez2seq.json":
            with open(os.path.join(directory,filename)) as json_file:
                entrez2seq = json.load(json_file)
                remapped_entrez2seq = remap_entrez2seq(entrez2seq)
                file_mapping['entrez2seq'] = remapped_entrez2seq  
                
                if debug:
                    print("%d gene sequences from Entrez"%(len(entrez2seq)))    
    
        if "description" in filename:
            if filename == "smpdb_pathway2description.json":
                with open(os.path.join(directory,'smpdb_pathway2description.json')) as json_file:
                    smpdb2desc = json.load(json_file)
                    id2desc.update({"SPMDB_Pathway:"+k:v for k,v in smpdb2desc.items()})
            elif filename == 'go2text_description.json':
                with open(os.path.join(directory,'go2text_description.json')) as json_file:
                    go2desc = json.load(json_file)
                    id2desc.update(go2desc)

    # protein contact map
    if contact_map_directory:              
        protein2contactmap = get_protein_contact_maps(contact_map_directory, debug=True)
        remapped_protein2contactmap = remap_protein2sequence(protein2contactmap)
        file_mapping['protein2contactmap'] = remapped_protein2contactmap
                      
    file_mapping['id2desc'] = id2desc
    if debug:
        print("%d node to descriptions"%len(id2desc))
    
#     compound2smiles = combine_compound(remapped_db2smiles, remapped_mesh2smiles)
    compound2smiles = file_mapping['mesh2smiles'] # skip drugbank since it's not in the safe release
    file_mapping['compound2smiles'] = compound2smiles
                      
    node2desc = remap_id2name_dicts(id2name, remap_namespaces=remap_namespaces_dict)
    file_mapping['node2desc'] = node2desc     
    return file_mapping
  
    
def sample_node_features(node_features, n=100):
    
    # Randomly sample n keys
    keys = [k for k,v in node_features.items() if len(v) > 0]
    sampled_keys = random.sample(keys, n)  
    return {key: node_features[key] for key in sampled_keys}

                      
def create_node_features(multimodal_data, kg_df=None):
    '''
    Creates node_features.json dict, mapping node name to its node features.
    If kg_df is included, only include nodes which are found in the kg.
    '''
    
    protein2sequence = multimodal_data['protein2sequence']
    compound2smiles = multimodal_data['compound2smiles']
    node2desc = multimodal_data['node2desc']
    id2desc = multimodal_data['id2desc']
    entrez2seq = multimodal_data['entrez2seq']
    protein2contactmap = multimodal_data['protein2contactmap']
                      
    if kg_df is not None:
        nodes = set(kg_df['h']).union(set(kg_df['t']))
        print("Nodes extracted from kg: %s"%len(nodes))
    else:
        nodes = list(node2desc.keys())
        print("Node descriptions: %s"%len(nodes))
        
    node_features = {}
    node_feature_counters={'names':0, 'protein_sequence':0, 'compound_sequence':0, 
                           'description':0, 'gene_sequence':0, 'protein_3d_structure':0}
    for n in nodes:
        features = {}
        if n in node2desc:
            feat_name = 'names'
            features[feat_name] = node2desc[n]
            node_feature_counters[feat_name] += 1
        if n in protein2sequence:
            feat_name = 'protein_sequence'
            features[feat_name] = protein2sequence[n]
            node_feature_counters[feat_name] += 1
        if n in compound2smiles:
            feat_name = 'compound_sequence'
            features[feat_name] = compound2smiles[n]
            node_feature_counters[feat_name] += 1
        if n in id2desc:
            feat_name = 'description'
            features[feat_name] = id2desc[n]
            node_feature_counters[feat_name] += 1
        if n in entrez2seq:
            feat_name = 'gene_sequence'
            features[feat_name] = entrez2seq[n]
            node_feature_counters[feat_name] += 1
        if n in protein2contactmap:
            feat_name = 'protein_3d_structure'
            features[feat_name] = protein2contactmap[n]
            node_feature_counters[feat_name] += 1
        node_features[n] = features
        
    print("Done. Number of nodes: %d"%len(node_features))
    for node_feature,count in node_feature_counters.items():
        print(node_feature+":",count)
            
    return node_features

    

# input files
multimodal_data_folder = './multi-modal_data'
protein_contact_map_directory = './protein_structure/protein_contact_map/'
whole_kg_file = './sampled_know2bio_safe_release/whole_kg_sampled.txt'


# read whole_kg only output node features within this kg
kg_df = read_kg(whole_kg_file)
nodes = set(kg_df['h']).union(set(kg_df['t']))
node_types = set([n.split(":")[0] for n in nodes])
node_type_to_nodes = {n_t:[n.split(":")[1] for n in nodes if n_t in n] for n_t in node_types}

# parse multimodal data
multimodal_data = get_multimodal_data_files(multimodal_data_folder, 
                                            contact_map_directory = protein_contact_map_directory,
                                            node_type_to_nodes = node_type_to_nodes,
                                            debug=True)

# generate node_features
node_features = create_node_features(multimodal_data, kg_df=kg_df)

sampled_node_features = sample_node_features(node_features, n=100)

with open("node_features.json", "w") as json_file:
    json_file.write(json.dumps(node_features, indent=4))
    
with open("sampled_node_features.json", "w") as json_file:
    json_file.write(json.dumps(sampled_node_features, indent=4))
