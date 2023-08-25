import os
import pandas as pd
import networkx as nx
from biomedkg_utils import output_edgefile_onerel_noweight


def map_reactome_pathway_ontology():
    pw_prefix = 'Reactome_Pathway:'
    headers = ['Pathway (Reactome)','Pathway (Reactome)', 'Relationship', 'Score']
    os.system('wget -N -P input/ https://reactome.org/download/current/ReactomePathwaysRelation.txt')

    pw_tree_df = pd.read_table('input/ReactomePathwaysRelation.txt', header=None)
    pw_tree_df.columns = ['Parent', 'Child']
    human_pw_tree_df = pw_tree_df[pw_tree_df['Parent'].str.contains('-HSA-')].copy()
    human_pw_tree_df['Parent'] = [pw_prefix+pw for pw in human_pw_tree_df['Parent']]
    human_pw_tree_df['Child'] = [pw_prefix+pw for pw in human_pw_tree_df['Child']]
    relations = ['-pathway_is_parent_of->']*len(human_pw_tree_df)
    scores = [1.0]*len(human_pw_tree_df)
    human_pw_tree_df.insert(2, 'relation', relations)
    human_pw_tree_df.insert(3, 'scores', scores)
    human_pw_tree_df.columns = headers

    human_pw_tree_df.to_csv('output/pathway2pathway/edges_reactomePathwayHierarchy.csv', index = False)
    human_pw_tree_df.to_csv('output/edges_to_use/Pathway_(Reactome)_2_Pathway_(Reactome).csv', index = False)
    df = pd.read_csv('output/pathway2pathway/edges_reactomePathwayHierarchy.csv')

    

def parse_kegg_hierarchy(file_name):
    nodes = []
    edges = []
    curr_top_level = ""
    curr_mid_level = ""

    lines = [l.strip("\n") for l in open(file_name,"r").readlines() if len(l)>1]
    for l in lines:

        # determine pathway level
        level = 3 # by default leaf pathway
        if ". " in l: # top level pathways don't have a number after '.'
            level = 0
            curr_top_level = l
        elif "." in l: # mid level pathways follow format "X.X" or "X.XX"
            level = 1
            curr_mid_level = l
            edges += [(curr_top_level,curr_mid_level)]
        else: # leaf pathway
            edges += [(curr_mid_level,l)]
        nodes += [l]
    G = nx.Graph()
    G.add_edges_from(edges)
    

    return G


def map_human_pathways_to_kegg_pathway_id(human_pathways, kegg_pathway_ids):

    human_pathway_to_kegg_pathway_id = {}

    for hpw in human_pathways:
        h_key = int(hpw.strip("path:hsa"))
        for kpw in kegg_pathway_ids:
            k_key = int(kpw.split(" ")[0])
            if h_key == k_key:
                human_pathway_to_kegg_pathway_id[hpw] = kpw
                break
    return human_pathway_to_kegg_pathway_id


def parse_mapping_table(input_file):
    df = pd.read_csv(input_file,sep="\t", header=None)
    mapping = {x:y for x,y in zip(df[0],df[1])}
    return mapping

def map_kegg_pathways():
    os.system('curl https://rest.kegg.jp/list/pathway/hsa > input/kegg_human_pathways.tsv')

    kegg_pathway_mapping_file = "input/kegg_human_pathways.tsv"
    kegg_pathway_mapping = parse_mapping_table(kegg_pathway_mapping_file)

    from_human_pathway_mapping = set(int(k.strip("path:hsa")) for k in kegg_pathway_mapping.keys())
    from_kegg_pathway_hierarchy = set(int(k.split(" ")[0]) for k in kegg_pathway_G.nodes() if "." not in k)
    intersection = from_human_pathway_mapping.intersection(from_kegg_pathway_hierarchy)

    print("%d human_pathways.tsv; %d intersection; %d kegg_pathway_hierarchy.csv"%(len(from_human_pathway_mapping),len(intersection),len(from_kegg_pathway_hierarchy)))

    human_pathways =  kegg_pathway_mapping.keys()
    kegg_pathway_ids = [k for k in kegg_pathway_G.nodes() if "." not in k]
    human_pathway_to_kegg_pathway_id = map_human_pathways_to_kegg_pathway_id(human_pathways, kegg_pathway_ids)
    # output kegg pathway to human pathway
    with open("input/kegg_human_pathway_to_pathway_id.tsv","w") as out_file:
        out_file.write("\n".join(["\t".join([x,y]) for x,y in human_pathway_to_kegg_pathway_id.items()]))
    # output kegg pathway hierarchy
    with open("input/kegg_pathway_hierchy.tsv","w") as out_file:
        out_file.write("\n".join(["\t".join([x,y]) for x,y in kegg_pathway_G.edges()]))
    
    
def process_pathway_name(name):
    name = name.split(' ')[1:]
    try:
        name_starts_with_letter = name[1] in ('M','N','R')
    except:
        name_starts_with_letter = False
    if name_starts_with_letter:
        name = ' '.join(name[2:])
    elif name[0] in ('M','N','R'):
        name = ' '.join(name[1:])
    else:
        name = ' '.join(name)
    return name

    
def align_kegg_pathway_id_to_name():
    '''Align Pathway ID -is- Pathway Name'''
    pathway_id2name = dict()
    pathway_name2id = dict()

    for line in open('input/kegg_human_pathway_to_pathway_id.tsv'):
        line = line.strip().split('\t')

        # Pathway ID, Pathway Name
        pathway_id = line[0].replace('path:','path_')
        pathway_name = process_pathway_name(' '.join(line[1:]))

        # Pathway Name -is- Pathway ID
        pathway_id2name[pathway_id] = pathway_name
        pathway_name2id[pathway_name] = pathway_id
    return pathway_name2id
    
    
def check_if_pathway_category(pathway_name):
    '''
    FUNCTION:
    - Some 'pathway names' are categories that begin with a decimal number.
      Some pathway names begin with the pathway ID. We want to determine
      which 'pathway names' are actually categories. We will use these as categories.
      The other pathway names will be included if they are human pathways.
    
    PARAMS:
    - pathway_name: The alleged 'pathway name' which includes prefixes.
    '''
    first_word = pathway_name.split(' ')[0]
    
    # Check if the first word is a decimal number (e.g., 1.1 Global and overview maps)
    
    # Decimal in first word?
    if '.' not in first_word:
        return False
    else:
        
        # First word == section or ID?
        try:
            first_word_float = float(first_word)
            return True
        except:
            return False
    

def map_kegg_pathway_ontology(pathway_name2id):
    ''' Pathway -parent of-> Pathway '''
    parent2child_kegg_pathway = dict()

    for line in open('input/kegg_pathway_hierchy.tsv'):
        line = line.strip().split('\t')

        # Pathways' Names
        parent_pathway_name = line[0].replace('path:','path_')
        child_pathway_name  = line[1].strip().replace('path:','path_')

        # Is the parent pathway good (human or a category)?
        cleaned_parent_pathway_name = process_pathway_name(parent_pathway_name)
        if cleaned_parent_pathway_name in pathway_name2id: 
            parent_human_pathway = True
            parent_pathway_category = False
        else:
            parent_human_pathway = False
            parent_pathway_category = check_if_pathway_category(parent_pathway_name)

        if parent_pathway_category:
            parent_pathway_id = parent_pathway_name
        elif parent_human_pathway:
            parent_pathway_id = pathway_name2id[cleaned_parent_pathway_name]
        else:
            continue


        # Is the child pathway good (human or a category)?
        cleaned_child_pathway_name = process_pathway_name(child_pathway_name)
        if cleaned_child_pathway_name in pathway_name2id:
            child_human_pathway = True
            child_pathway_category = False
        else:
            child_human_pathway = False
            child_pathway_category = check_if_pathway_category(child_pathway_name)

        if child_pathway_category:
            child_pathway_id = child_pathway_name
        elif child_human_pathway:
            child_pathway_id = pathway_name2id[cleaned_child_pathway_name]
        else:
            continue

        # Parent Pathway - Child Pathway
        parent2child_kegg_pathway.setdefault(parent_pathway_id, set()).add(child_pathway_id)


    # Output edges    
    file = 'Pathway_(KEGG)_2_Pathway_(KEGG).csv'
    outpath = os.path.join('output/pathway2pathway/', file)
    output_edgefile_onerel_noweight(
        outpath = outpath,
        columns = ['Pathway (KEGG)','Pathway (KEGG)','Relationship'],
        dictionary = parent2child_kegg_pathway,
        rel = '-parent_of->',
        prefix_col1='KEGG_Pathway:',
        prefix_col2='KEGG_Pathway:'
    )
    df = pd.read_csv(outpath)
    df.to_csv(os.path.join('output/edges', file), index=False)
    df.to_csv(os.path.join('output/edges_to_use/', file), index=False)
    

if __name__ == '__main__':
    map_reactome_pathway_ontology()
    kegg_pathway_G = parse_kegg_hierarchy('input/kegg_pathway_hierarchy.csv')
    # INSTRUCTIONS: Obtain the KEGG pathway hierarchy, manually formatted from the website
    # https://www.kegg.jp/kegg/pathway.html or downloaded from the GitHub: kegg_pathway_hierarchy.csv
    map_kegg_pathways()
    pathway_name2id = align_kegg_pathway_id_to_name()
    map_kegg_pathway_ontology(pathway_name2id)