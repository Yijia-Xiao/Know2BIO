import os
import pandas as pd
import networkx as nx
from biomedkg_utils import output_edgefile_onerel_noweight
import requests
from bs4 import BeautifulSoup
import re
import csv

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

    
def get_kegg_pathway_hierarchy():
    url = 'https://www.genome.jp/kegg/pathway.html'
    response = requests.get(url)
    soup = BeautifulSoup(response.text, 'html.parser')
    text = soup.get_text()

    lines = []
    print_the_line = False
    for line in text.split('\n'):
        if print_the_line and line != '':
            line = re.sub(r'(\d)([a-zA-Z])', r'\1 \2', line)
            line = re.sub(r'(\d) ([A-Z]{2})', r'\1 \2', line)
            lines.append(line)
        if line.startswith('N - network'):
            print_the_line = True
    with open('input/KEGG/kegg_pathway_hierarchy.csv','w') as fout:
        for line in lines:
            fout.write(line+'\n')

            
def remove_letter_before_kegg_pathway(pathway_name):
    pathway_id = pathway_name[0]
    pathway_name = pathway_name[1:]
    
    # Remove the 'M', 'R', or 'N' that's on its own
    pathway_name_starts_with_one_letter_word = len(pathway_name[0]) == 1
    if pathway_name_starts_with_one_letter_word: 
        pathway_name = pathway_name[1:]
        
    # Remove the 'M', 'R', or 'N' that got merged with the pathway name
    letter_got_merged_with_name = pathway_name[0][:2].isupper()
    mrn = pathway_name[0].startswith(('M','R','N'))
    if letter_got_merged_with_name and mrn:
        pathway_name[0] = pathway_name[0][1:]
    
    pathway_name = ' '.join(pathway_name)
    return pathway_id, pathway_name
            
            
def map_kegg_pathway_to_pathway_hierarchy():            
    # Finish KEGG pathway to pathway
    file = 'Pathway_(KEGG)_to_Pathway_(KEGG).csv'
    with open(os.path.join('output/pathway2pathway/',file),'w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Pathway (KEGG)', 'Pathway (KEGG)', 'Relationship'])

        with open('input/KEGG/kegg_pathway_hierarchy.csv') as fin:
            for line in fin:
                line = line.split(' ')
                first_number = line[0]
                if first_number.endswith('.'):
                    line = ' '.join(line[1:])
                    top_level = line.strip().replace(' ','_')
                elif '.' in first_number:
                    line = ' '.join(line[1:])
                    middle_level = line.strip().replace(' ','_')
                    writer.writerow(['KEGG_Pathway:'+top_level, 
                                     'KEGG_Pathway:'+middle_level, 
                                     '-pathway_is_parent_of->'])
                else:
                    pathway_id, _ = remove_letter_before_kegg_pathway(line)
                    writer.writerow(['KEGG_Pathway:'+middle_level, 
                                     'KEGG_Pathway:path_hsa'+pathway_id,
                                     '-pathway_is_parent_of->'])

    kegg_pathway_df = pd.read_csv(os.path.join('output/pathway2pathway/',file))
    kegg_pathway_df.to_csv(os.path.join('output/edges/',file), index=False)
    kegg_pathway_df.to_csv(os.path.join('output/edges_to_use/',file), index=False)    

    

if __name__ == '__main__':
    map_reactome_pathway_ontology()
    get_kegg_pathway_hierarchy()
    map_kegg_pathway_to_pathway_hierarchy()
    
    
    
    