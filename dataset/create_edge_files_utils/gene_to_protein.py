import os
import requests
from biomed_apis import *


def map_gene_to_protein():
    os.system('wget -N -P input/ https://www.genenames.org/cgi-bin/download/custom?col=gd_app_name&col=gd_pub_acc_ids&col=md_prot_id&col=gd_pub_eg_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit')
    os.system('mv input/custom?col=gd_app_name input/custom_entrez_uniprot.txt')

    protein_id_2_gene_id, gene_id_2_protein_id = dict(), dict()
    for line in open('input/custom_entrez_uniprot.txt','r'):
        line = line.split('\t')
        protein_ids = line[2].strip() # UniProt Protein ID
        gene_id = line[3].strip()     # Entrez Gene ID

        if 'UniProt ID(supplied by UniProt)' not in protein_ids:
            protein_ids = protein_ids.replace('_','').split(',')
            for protein_id in protein_ids:
                if protein_id != '' and gene_id != '':
                    protein_id_2_gene_id.setdefault(protein_id.strip(), set()).add(gene_id.strip())
                    gene_id_2_protein_id.setdefault(gene_id.strip(),set()).add(protein_id.strip())

    protein_id_2_gene_id = switch_dictset_to_dictlist(protein_id_2_gene_id)
    gene_id_2_protein_id = switch_dictset_to_dictlist(gene_id_2_protein_id)
    
    return protein_id_2_gene_id, gene_id_2_protein_idd
    
    
def download_human_proteins():
    # Reviewed human proteome
    url = 'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names&format=json&query=%28reviewed%3Atrue%20AND%20proteome%3Aup000005640%29'
    res = requests.get(url)
    r = res.json()
    protein_set = set(['UniProt:'+entry['primaryAccession'] for entry in r['results']])
    print(len(protein_set))
    return protein_set
    

def save_human_proteins(protein_set):
    protein_list = list(protein_set)
    protein_set = set()
    for protein_ids in protein_list:

        # Split each "protein" into multiple proteins if there are multiple
        protein_ids = protein_ids.split(', ')

        for protein_id in protein_ids:

            # Remove "UniProt:" prefix
            if 'UniProt:' in protein_id:
                protein_id = protein_id.split('UniProt:')[1].strip().split('_')[0]

            # Don't add numbers (genes), add UniProt IDs
            if not protein_id[0:].isnumeric():
                protein_set.add(protein_id)

    json.dump(list(protein_set), open('output/protein2protein/human_proteins.json','w'))
    return protein_set

def map_protein_to_gene_via_uniprot_api(protein_id_2_gene_id, gene_id_2_protein_id, protein_set):
    job_id = submit_id_mapping_UniProtAPI(
                      from_db = 'UniProtKB_AC-ID',
                      to_db = 'GeneID', 
                      ids = protein_set)

    # This checks on the job until it is finished
    if check_id_mapping_results_ready_UniProtAPI(job_id):
        link = get_id_mapping_results_link_UniProtAPI(job_id)
        results = get_id_mapping_results_search_UniProtAPI(link)


    print('Before UniProt API')
    print('Entrez-is-UniProt', len(gene_id_2_protein_id), 
          'UniProt-is-Entrez', len(protein_id_2_gene_id))

    gene_id_2_protein_id = switch_dictlist_to_dictset(gene_id_2_protein_id)
    protein_id_2_gene_id = switch_dictlist_to_dictset(protein_id_2_gene_id)

    for protein_to_gene in results['results']:
        protein_id = protein_to_gene['from'].strip()
        gene_id = protein_to_gene['to'].strip()
        if gene_id != '' and protein_id != '':
            gene_id_2_protein_id.setdefault(gene_id, set()).add(protein_id)
            protein_id_2_gene_id.setdefault(protein_id, set()).add(gene_id)

    print('\nAfter UniProt API')
    print('Entrez-is-UniProt', len(gene_id_2_protein_id), 
          'UniProt-is-Entrez', len(protein_id_2_gene_id))    
    return gene_id_2_protein_id, protein_id_2_gene_id
    
    
def export_gene_to_protein_mappings(protein_id_2_gene_id, gene_id_2_protein_id):
    gene_id_2_protein_id = switch_dictset_to_dictlist(gene_id_2_protein_id)
    protein_id_2_gene_id = switch_dictset_to_dictlist(protein_id_2_gene_id)
    
    json.dump(gene_id_2_protein_id, open('output/protein2gene/all_entrez2uniprot.json','w'))
    json.dump(protein_id_2_gene_id, open('output/protein2gene/all_uniprot2entrez.json','w'))

    # Edges
    gene2prot_dict = json.load(open('output/protein2gene/all_entrez2uniprot.json'))
    with open('output/edges/edges_gene-ENCODES->protein.csv','w') as fout1:
        writerE = csv.writer(fout1)
        writerE.writerow(['Gene (Entrez)','Protein (UniProt)', 'Relationship'])

        with open('output/protein2gene/edges_gene-ENCODES->protein.csv','w') as fout:
            writer = csv.writer(fout)
            writer.writerow(['Gene (Entrez)','Protein (UniProt)', 'Relationship'])
            for gene, proteins in gene2prot_dict.items():
                for protein in proteins:
                    writer.writerow(['Entrez:'+gene, 'UniProt:'+protein, '-encodes->'])
                    writerE.writerow(['Entrez:'+gene, 'UniProt:'+protein, '-encodes->'])

    os.system('cp "output/edges/edges_gene-ENCODES->protein.csv" "output/edges/Gene_(Entrez)_2_Protein_(UniProt).csv"')
    os.system('cp "output/edges/edges_gene-ENCODES->protein.csv" "output/edges_to_use/Gene_(Entrez)_2_Protein_(UniProt).csv"')
    

if __name__ == '__main__':
    #map_gene_to_protein()
    protein_set = download_human_proteins()
    protein_set_no_prefix = save_human_proteins(protein_set)
    protein_id_2_gene_id, gene_id_2_protein_id = {}, {} # while HGNC is broken
    gene_id_2_protein_id, protein_id_2_gene_id = map_protein_to_gene_via_uniprot_api(protein_id_2_gene_id,
                                                                                     gene_id_2_protein_id, 
                                                                                     protein_set_no_prefix)
    export_gene_to_protein_mappings(protein_id_2_gene_id, gene_id_2_protein_id)