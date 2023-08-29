import requests
import sys
import json
from gene_to_anatomy import batch_iterator        
from multiprocessing import cpu_count
from joblib import Parallel, delayed
    
'''Genes'''
def map_gene_id_to_dna_sequence():
    ensembl_to_entrez = json.load(open('output/gene2gene/ensembl2entrez.json'))
    entrez_id_to_dna_sequence = {}

    # Use API
    url = 'https://rest.ensembl.org/sequence/id'
    headers = {'Content-Type':'application/json',
               'Accept' : 'application/json'}
    headers={"Content-Type" : "application/json", 
             "Accept" : "application/json"}
    BATCH_SIZE = 50
    ensembl_ids = list(ensembl_to_entrez.keys())

    for idx, gene_ids in enumerate(batch_iterator(ensembl_ids, BATCH_SIZE)):
        print(idx*BATCH_SIZE, '/', len(ensembl_ids), end='\r')

        data = {'ids': gene_ids}
        r = requests.post(url, headers=headers, json=data)

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        results = r.json()

        # Extract results
        for entry in results:
            ensembl_id = entry['query']
            try:
                entrez_id = ensembl_to_entrez[ensembl_id]
            except:
                continue
            dna_sequence = entry['seq']

            entrez_id_to_dna_sequence[entrez_id] = dna_sequence
    # Export 
    with open('output/node_features/sequences/gene_id_to_dna_sequences.json','w') as fout:
        json.dump(entrez_id_to_dna_sequence, fout)

        
        
'''Proteins'''
def extract_amino_acid_sequence(data):
    return ''.join(data.split(' ')[-1].split('\n')[1:-1])


def remove_ids_not_in_all_ids(my_dict, all_ids):
    return {key:list(value) for key,value in my_dict.items() if key in all_ids}


def batch_map_protein_id_to_sequence(protein_ids):
    protein_id_to_sequence = {}
    for protein_id in protein_ids:
        url = f'http://www.uniprot.org/uniprot/{protein_id}.fasta'
        response = requests.post(url)
        data = ''.join(response.text)
        sequence = extract_amino_acid_sequence(data)
        protein_id_to_sequence[protein_id] = sequence
    return protein_id_to_sequence


'''
def map_protein_id_to_sequence_in_batches(protein_ids):
    BATCH_SIZE = 70 if len(protein_ids) >= 70 else len(protein_ids)
    protein_id_to_sequence_batch = {}
    for idx, protein_batch in enumerate(batch_iterator(protein_ids, BATCH_SIZE)):

        # Query protein sequence using ID
        protein_batch = ','.join(protein_batch)
        url = f'https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=uniprotkb&id={protein_batch}&format=fasta&style=raw'
        res = requests.get(url).text

        # Extract protein IDs & sequences
        results = res.split('>')
        results.remove('')
        for idx in range(len(results)):
            protein_id = results[idx].split('|')[1]
            sequence = extract_amino_acid_sequence(results[idx])
            protein_id_to_sequence_batch[protein_id] = sequence
    
    return protein_id_to_sequence_batch
'''


def split_a_list_into_batches(input_list, BATCHES):
    '''BATCHES = number of batches, not size of batch'''
    batch_size = len(input_list) // BATCHES
    remainder = len(input_list) % BATCHES
    
    batched_lists = []
    start_index = 0
    
    for _ in range(BATCHES):
        if remainder > 0:
            end_index = start_index + batch_size + 1
            remainder -= 1
        else:
            end_index = start_index + batch_size
        
        batched_lists.append(input_list[start_index:end_index])
        start_index = end_index
    
    return batched_lists


def map_protein_id_to_sequence_via_uniprot_human_proteome_download():
    protein_id_to_sequence = {}

    # Call API
    url = 'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Csequence&format=json&query=%28*%29+AND+%28reviewed%3Atrue%29+AND+%28model_organism%3A9606%29'
    res = requests.get(url)
    results = res.json()['results']

    # Extract results from API
    for result in results:
        protein_id = result['primaryAccession']
        sequence = result['sequence']['value']
        protein_id_to_sequence[protein_id] = sequence
        
    return protein_id_to_sequence


def map_protein_id_to_amino_acid_sequence():
    # Get protein sequences from UniProt
    protein_id_to_sequence_uniprot = (
        map_protein_id_to_sequence_via_uniprot_human_proteome_download())

    # Get protein sequences from EBI that were missing from UniProt
    proteins_without_sequences = list(set(protein_ids).difference(protein_id_to_sequence_uniprot))
    protein_batches = split_a_list_into_batches(proteins_without_sequences, BATCHES)
    protein_to_sequence_dicts = (
        Parallel(n_jobs=-1)(delayed(batch_map_protein_id_to_sequence)(protein_batch) 
                            for protein_batch in protein_batches))

    # Combine results
    protein_to_sequence = {}
    for protein_to_sequence_batch in protein_to_sequence_dicts:
        protein_to_sequence.update(protein_to_sequence_batch)
    protein_to_sequence.update(protein_id_to_sequence_uniprot)

    all_ids = json.load(open('output/otherMappings/all_ids.json'))
    all_protein_ids = list({id_.split(':')[1] for id_ in all_ids if id_.startswith('UniProt')})
    protein_to_sequence = remove_ids_not_in_all_ids(protein_to_sequence, all_protein_ids)
    print(f'Mapped {len(protein_to_sequence)}/{len(protein_ids)} proteins to sequences')

    # Export results                                           
    with open('output/node_features/sequences/protein_id_to_sequences.json','w') as fout:
        json.dump(protein_to_sequence, fout)      
        
        
def map_compound_ids_to_chemical_smiles_sequence():
    all_ids = json.load(open('output/otherMappings/all_ids.json'))

    # DrugBank Compounds
    drugbank_to_smiles = json.load(open('output/compound2compound/db2smiles.json'))
    drugbank_to_sequence = {'DrugBank_Compound:'+drugbank_id:sequence for drugbank_id, sequence in drugbank_to_smiles.items()}
    all_drugbank_ids = list({id_ for id_ in all_ids if id_.startswith('DrugBank_Compound:')})
    print(len(drugbank_to_sequence), 'DrugBank compounds before filtering by compounds in the KG')
    drugbank_to_sequence = remove_ids_not_in_all_ids(drugbank_to_sequence, all_ids)
    print(f'Mapped {len(drugbank_to_sequence)}/{len(all_drugbank_ids)} DrugBank compounds to sequences')
    with open('output/node_features/sequences/drugbank_compounds_to_chemical_sequence.json','w') as fout:
        json.dump(drugbank_to_sequence, fout)

    # MeSH Compounds
    mesh_to_smiles = json.load(open('output/compound2compound/mesh2smiles.json'))
    mesh_to_sequence = {'MeSH_Compound:'+mesh_id:sequence for mesh_id, sequence in mesh_to_smiles.items()}
    all_mesh_ids = list({id_ for id_ in all_ids if id_.startswith('MeSH_Compound:')})
    mesh_to_sequence = remove_ids_not_in_all_ids(mesh_to_sequence, all_ids)
    print(len(mesh_to_sequence), 'MeSH compounds before filtering by compounds in the KG')
    print(f'Mapped {len(mesh_to_sequence)}/{len(all_mesh_ids)} MeSH Compounds to sequences')
    with open('output/node_features/sequences/mesh_compounds_to_chemical_sequence.json','w') as fout:
        json.dump(mesh_to_sequence, fout)

    mesh_sequences = set([seq[0] for seq in mesh_to_sequence.values()])
    drugbank_sequences = set([seq[0] for seq in drugbank_to_sequence.values()])

    drugbank_mesh_sequences = mesh_sequences.intersection(drugbank_sequences)
    print(f'{len(drugbank_mesh_sequences)} compounds are the same (have a DrugBank and a MeSH ID we know of)')

    

if __name__ == '__main__':
    map_gene_id_to_dna_sequence()
    map_protein_id_to_amino_acid_sequence()
    map_compound_ids_to_chemical_smiles_sequence()  # Compound sequences obtained through compound_to_compound_alignment
    