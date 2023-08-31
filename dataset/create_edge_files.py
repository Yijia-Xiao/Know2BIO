import os

os.chdir('create_edge_files_utils/')

file_names = [
    'compound_to_compound_alignment.py', # required for many compound files
    'gene_to_protein.py', # required for some gene/protein files
    'anatomy_to_anatomy.py',
    'disease_to_disease.py',
    'pathway_to_pathway.py',
    'protein_to_protein.py',
    'compound_to_compound_interactions.py',
    'compound_to_disease.py',
    'compound_to_drug_class.py',
    'compound_to_gene.py',
    'compound_to_protein.py',
    'compound_to_pathway.py',
    'compound_to_side_effect.py',
    'protein_and_compound_to_reaction.py',
    'protein_and_gene_to_pathway.py',
    'protein_to_gene_ie_transcription_factor_edges.py',
    'disease_to_pathway.py',
    'gene_to_anatomy.py',
    'gene_to_disease.py',
    'go_to_go.py',
    'go_to_protein.py',
    'reaction_to_pathway.py',
    'reaction_to_reaction.py',
    'disease_to_anatomy.py',
    'protein_and_compound_to_reaction.py',
    'protein_and_gene_to_pathway.py',
    'get_node_feature_natural_language_names.py',
    'get_node_feature_sequences.py',
    'get_node_feature_compound_structures.py',
]

with open('output/log_file.txt','w') as fout:
    for file_name in file_names:
        print('*'*50)
        print(f"{'*'*10} Running {file_name} {'*'*10}")
        fout.write(f"{'*'*10} Running {file_name} {'*'*10}")
        os.system(f'python create_edge_files_utils/{file_name}')
        print('*'*50)
