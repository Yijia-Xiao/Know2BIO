from datetime import datetime
import urllib.request
import xml.etree.ElementTree as ET
import json
import os
from biomedkg_utils import output_edgefile_onerel_noweight, switch_dictset_to_dictlist
import pandas as pd
import obonet

def download_mesh_xml():
    url = f'https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/desc{year}.xml'
    dest = f'input/desc{year}.xml'
    urllib.request.urlretrieve(url, dest);

    
def parse_mesh_xml():
    year = datetime.now().year
    tree = ET.parse(f'input/desc{year}.xml')
    root = tree.getroot()
    
    return root


def align_anatomy_identifiers_within_mesh(root):
    name2id, id2name, id2tree, tree2id = dict(), dict(), dict(), dict()
    name2tree, tree2name = {}, {}
    all_tree_numbers = list()
    # Anatomy tree terms. Modify as you wish. Source: https://meshb.nlm.nih.gov/treeView
    anatomy_mesh_tree_prefix = ('A01','A02','A03','A04',
                                'A05','A06','A07','A08',
                                'A09','A10','A11','A12',
                                'A14','A15','A16','A17.815',
                                'A17.360.710','A17.360.421',
                                'A17.360.296')

    for ele in root:
        try:
            # MeSH Tree Number
            tree_numbers = ele.find('TreeNumberList').findall('TreeNumber')

            # If anatomy
            for tree_number in tree_numbers:
                if tree_number.text.startswith(anatomy_mesh_tree_prefix):

                    tree_number = tree_number.text
                    all_tree_numbers.append(tree_number)

                    # ID to Tree
                    try:
                        ID = ele.find('DescriptorUI').text
                        id2tree.setdefault(ID,set()).add(tree_number)
                        tree2id.setdefault(tree_number,set()).add(ID)
                    except:
                        pass

                    # ID to Name
                    try:
                        ID = ele.find('DescriptorUI').text
                        name = ele.find('DescriptorName').find('String').text
                        name2id.setdefault(name,set()).add(ID)
                        id2name.setdefault(ID,set()).add(name)
                    except:
                        pass
                    
                    # Name to Tree
                    try:
                        name = ele.find('DescriptorName').find('String').text
                        name2tree.setdefault(name, set()).add(tree_number)
                        tree2name.setdefault(tree_number, set()).add(name)
                    except:
                        pass
        except:
            continue        

    all_tree_numbers = sorted(all_tree_numbers)
    tree2id = dict(sorted(tree2id.items()))

    for k,v in name2id.copy().items():
        name2id[k] = list(name2id[k])
        
    # MeSH Term -[is]- MeSH ID
    with open('output/anatomy2anatomy/meshterm-IS-meshid.json','w') as fout:
        json.dump(name2id, fout)
        
    name2tree = switch_dictset_to_dictlist(name2tree)
    with open('output/anatomy2anatomy/name2tree.json','w') as fout:
        json.dump(name2tree, fout)
        
    tree2name = switch_dictset_to_dictlist(tree2name)
    with open('output/anatomy2anatomy/tree2name.json','w') as fout:
        json.dump(tree2name, fout)
        
    return id2tree, all_tree_numbers


def export_mesh_tree_and_id_alignment(id2tree):
    # MeSH Tree Number -[is]- MeSH ID
    file = 'Anatomy_(MeSH)_2_Anatomy_(MeSH_Tree).csv'
    outpath = os.path.join('output/anatomy2anatomy/',file)
    output_edgefile_onerel_noweight(outpath = outpath,
                                    columns = ['Anatomy (MeSH)','Anatomy (MeSH Tree)','Relationship'],
                                    dictionary = id2tree,
                                    rel = '-is-',
                                    prefix_col1 = 'MeSH_Anatomy:',
                                    prefix_col2 = 'MeSH_Tree_Anatomy:')

    tree_to_id = pd.read_csv(outpath).drop_duplicates()
    tree_to_id.to_csv(os.path.join('output/edges/', file), index=False)
    tree_to_id.to_csv(os.path.join('output/edges_to_use/', file), index=False)
    NUM_MAPPED_MESH_IDS = len(set(tree_to_id['Anatomy (MeSH)']))
    NUM_MAPPED_MESH_TREE_NUMS = len(set(tree_to_id['Anatomy (MeSH Tree)']))
    print(f'Aligned {NUM_MAPPED_MESH_IDS} MeSH IDs to '+\
          f'{NUM_MAPPED_MESH_TREE_NUMS} MeSH Tree Numbers '+\
          f'\nTotal of {len(tree_to_id)} relationships, averaging '+\
          f'{round(NUM_MAPPED_MESH_TREE_NUMS/NUM_MAPPED_MESH_IDS, 2)} Tree Numbers per ID')
    

def map_and_export_mesh_tree_number_hierarchy(all_tree_numbers):
    tree2tree = dict()

    # Tree Number
    for tree_num in all_tree_numbers:
        if '.' in tree_num:

            # Parent of Tree Number
            parent = ''
            for num in tree_num.split('.')[:len(tree_num.split('.'))-1]:
                parent += num+'.'
            parent = parent.strip('.')

            # Tree Number -[subclass of]-> Tree Number
            tree2tree[tree_num] = [parent]


    # MeSH Tree Number -[subclass of]-> MeSH Tree Number
    file = 'Anatomy_(MeSH_Tree)_2_Anatomy_(MeSH_Tree).csv'
    outpath = os.path.join('output/anatomy2anatomy/',file)
    output_edgefile_onerel_noweight(outpath = outpath,
                                    columns = ['Anatomy (MeSH Tree)','Anatomy (MeSH Tree)','Relationship'],
                                    dictionary = tree2tree,
                                    rel = '-subclass_of->',
                                    prefix_col1 = 'MeSH_Tree_Anatomy:',
                                    prefix_col2 = 'MeSH_Tree_Anatomy:')

    tree_df = pd.read_csv(outpath).drop_duplicates()
    tree_df.to_csv(os.path.join('output/edges/', file), index=False)
    tree_df.to_csv(os.path.join('output/edges_to_use/', file), index=False)
    print(f'{len(tree_df)} MeSH Tree Number intra relationships')
    
    
def download_uberon_anatomy():
    # Download Uberon obo
    url = 'http://purl.obolibrary.org/obo/uberon.obo'
    dest = 'input/uberon.obo'
    urllib.request.urlretrieve(url, dest);
    
def align_uberon_to_mesh_anatomy(dest='input/uberon.obo'):
    # Parse Uberon obo, align to MeSH
    obo_file = dest
    graph = obonet.read_obo(obo_file)

    uberon_to_mesh = {}
    for node_term, data in graph.nodes(data=True):
        if not node_term.startswith('UBERON'):
            continue    
        if 'xref' in data:
            for xref in data['xref']:
                if xref.startswith('MESH:'):
                    uberon_to_mesh.setdefault(node_term, []).append(xref)
                    
    # Export Uberon-to-MeSH alignment
    for k,v in uberon_to_mesh.items():
        assert len(v) == 1

    with open('output/anatomy2anatomy/uberon2mesh.json','w') as fout:
        json.dump(uberon_to_mesh, fout)                

if __name__ == '__main__':        
    download_mesh_xml()
    root = parse_mesh_xml()
    id2tree, all_tree_numbers = align_anatomy_identifiers_within_mesh(root)
    export_mesh_tree_and_id_alignment(id2tree)
    map_and_export_mesh_tree_number_hierarchy(all_tree_numbers)

    download_uberon_anatomy()
    uberon_to_mesh = align_uberon_to_mesh_anatomy()