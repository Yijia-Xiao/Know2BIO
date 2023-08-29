import json
import requests
import pandas as pd
from biomed_apis import *

def download_human_proteome():
    # Human proteome
    url = 'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Cgene_primary%2Cgene_synonym%2Creviewed%2Ccc_tissue_specificity%2Corganism_id&format=tsv&query=%28proteome%3AUP000005640%29%20AND%20%28reviewed%3Atrue%29%20AND%20%28model_organism%3A9606%29'
    request = requests.get(url)

    lines = request.text.split('\n')
    line_entries = [line.split('\t') for line in lines]
    header = line_entries[0]
    df = pd.DataFrame(line_entries[1:], columns=header)
    human_protein_ids = set(df['Entry'])
    human_protein_ids.remove('')
    human_protein_ids = list(human_protein_ids)
    uniprot_ids = ['UniProt:'+id_ for id_ in human_protein_ids]

    json.dump(human_protein_ids, open('output/protein2protein/human_protein_ids_no_prefix','w'))
    json.dump(uniprot_ids, open('output/protein2protein/human_protein_ids_with_prefix','w'))
    
    
def get_human_proteins():

    protein_list = list(json.load(open('output/protein2protein/human_protein_ids_no_prefix')))
    protein_set = set()
    for protein_ids in protein_list:

        # Split each "protein" into multiple proteins if there are multiple
        protein_ids = protein_ids.split(', ')

        for protein_id in protein_ids:

            # Remove "UniProt:" prefix
            if 'UniProt:' in protein_id:
                protein_id = protein_id.split('UniProt:')[1]

            # Don't add numbers (genes), add UniProt IDs
            if not protein_id[0:].isnumeric():
                protein_set.add(protein_id)
    print(len(protein_set), 'human proteins')
                
    return protein_set

def align_protein_ids_uniprot_to_string(protein_set):
    # Call API
    job_id = submit_id_mapping_UniProtAPI(
                      from_db = 'UniProtKB_AC-ID',
                      to_db = 'STRING',
                      ids = protein_set)

    if check_id_mapping_results_ready_UniProtAPI(job_id):
        link = get_id_mapping_results_link_UniProtAPI(job_id)
        results = get_id_mapping_results_search_UniProtAPI(link)
                
    # Map API results
    string_is_uniprot, uniprot_is_string = dict(), dict()
    for entry in results['results']:
        string_id = entry['to']
        uniprot_id = entry['from']

        if '9606' in string_id:
            #uniprot_is_string[uniprot_id] = string_id
            string_is_uniprot[string_id] = uniprot_id

    print(len(string_is_uniprot),'/',len(protein_set), 'UniProt proteins aligned with STRING')
    return string_is_uniprot


def find_associated_proteins(string_is_uniprot, NUM_HOPS=5, CONF_SCORE=0.70):
    string_ids = list(string_is_uniprot.keys())
    string_ppi_df, counts = k_hop_interactors_STRINGAPI(
                                              string_ids, 
                                              k=NUM_HOPS, 
                                              score_thresh = CONF_SCORE,
                                              debug = True,
                                              return_counts = True,
    )
    return string_ppi_df
    

def align_string_to_uniprot(string_ppi_df, string_is_uniprot):
    ppi_string_ids = set(string_ppi_df['query_ensp']).union(set(string_ppi_df['partner_ensp']))

    ### Select the PPI's STRING IDs that still need UniProt IDs
    string_ids2 = set()
    for string_id in ppi_string_ids:
        if string_id not in string_is_uniprot.keys():
            string_ids2.add(string_id)

    string_ids2.add(list(ppi_string_ids)[0])
    
    print(len(string_is_uniprot),'/',len(protein_set), 'UniProt proteins aligned with STRING')
    
    
def export_protein_protein_interactions_from_string(string_is_uniprot, string_ppi_df):
    with open('output/protein2protein/edges_protein2protein_string.csv','w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Protein (UniProt)','Protein (UniProt)','Relationship','Weight'])
        string2string_ppi = dict()
        for i in range(0,len(string_ppi_df)):
            try:
                protein1 = string_is_uniprot[string_ppi_df['query_ensp'].iloc[i]]
                protein2 = string_is_uniprot[string_ppi_df['partner_ensp'].iloc[i]]
                score = str(string_ppi_df['combined_score'].iloc[i])
                writer.writerow(['UniProt:'+protein1, 'UniProt:'+protein2, '-ppi-',score])
                string2string_ppi.setdefault(protein1,[]).append(protein2)
            except:
                continue

        string_ppi = pd.read_csv('output/protein2protein/edges_protein2protein_string.csv')
        string_ppi = string_ppi[~string_ppi.apply(frozenset, axis=1).duplicated()] # remove duplicates
        string_ppi.to_csv('output/edges/edges_protein-INTERACTS-protein_string.csv', index=False)
    
# Reactome is nearly a proper subset of STRING
#def get_ppi_from_reactome():
#    os.system('wget -N -P input/ https://reactome.org/download/current/interactors/reactome.homo_sapiens.interactions.tab-delimited.txt')
#    reactome_ppi = pd.read_table('input/reactome.homo_sapiens.interactions.tab-delimited.txt')#

#    with open('output/protein2protein/edges_all_protein2protein_reactome.csv', 'w') as fout:
#        writer = csv.writer(fout)
#        writer.writerow(['Protein (UniProt)','Protein (UniProt)', 'Relationship'])
#        for prot1, prot2 in zip(reactome_ppi['# Interactor 1 uniprot id'], reactome_ppi['Interactor 2 uniprot id']):
#            if 'uniprot' in prot1.lower() and 'uniprot' in prot2.lower():
#                writer.writerow([prot1.replace('uniprotkb','UniProt'), prot2.replace('uniprotkb','UniProt'), '-ppi-'])
#    reactome_ppi = pd.read_csv('output/protein2protein/edges_all_protein2protein_reactome.csv')

    
def export_ppi_from_string():
    string_ppi = pd.read_csv('output/edges/edges_protein-INTERACTS-protein_string.csv')
    
    file = 'Protein_(UniProt)_2_Protein_(UniProt).csv'
    string_ppi.to_csv(os.path.join('output/protein2protein', file), index=False)
    string_ppi.to_csv(os.path.join('output/edges/', file), index=False)
    string_ppi.to_csv(os.path.join('output/edges_to_use/', file), index=False)
    
if __name__ == '__main__':
    download_human_proteome()
    protein_set = get_human_proteins()
    string_is_uniprot = align_protein_ids_uniprot_to_string(protein_set)
    string_ppi_df = find_associated_proteins(string_is_uniprot)
    align_string_to_uniprot(string_ppi_df, string_is_uniprot)
    export_protein_protein_interactions_from_string(string_is_uniprot, string_ppi_df)
    #get_ppi_from_reactome()
    export_ppi_from_string()