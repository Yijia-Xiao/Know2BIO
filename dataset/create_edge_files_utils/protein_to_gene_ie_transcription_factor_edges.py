import os
import mygene
import json
import pandas as pd
from biomedkg_utils import switch_dictset_to_dictlist, switch_dictlist_to_dictset, output_edgefile_onerel_noweight
def download_transcription_factor_data_from_grndb():
    grndb_link_prefix = "http://www.grndb.com/download/txt?condition="
    tissues = ['Heart_GTEx','Adult-Heart', 'Fetal-Heart','whole_NeonatalHeart',
    'Blood_Vessel_GTEx','Adipose_Tissue_GTEx', 
    'Blood_GTEx', 'Adrenal_Gland_GTEx', 'Breast_GTEx', 'Colon_GTEx', 
    'Esophagus_GTEx','Kidney_GTEx','Liver_GTEx','Lung_GTEx',
    'Muscle_GTEx','Esophagus_GTEx', 'Nerve_GTEx', 'Ovary_GTEx', 
    'Pancreas_GTEx','Pituitary_GTEx', 'Prostate_GTEx', 'Salivary_Gland_GTEx', 
    'Skin_GTEx','Small_Intestine_GTEx', 'Spleen_GTEx', 'Stomach_GTEx', 
    'Testis_GTEx', 'Thyroid_GTEx', 'Uterus_GTEx', 'Vagina_GTEx']

    if not os.path.exists('input/GRNdb'):
        os.mkdir('input/GRNdb')

    for tissue in tissues:
        url = grndb_link_prefix+tissue
        os.system(f'wget -N -P input/GRNdb/ {url}')      


def map_tf_name_to_target_name():
    # Collect gene names of all TFs and targets from each tissue
    target_gene_names = set()
    tf_gene_names = set()
    tf_gene_name_to_target_gene_name = dict()
    root = 'input/GRNdb/'

    for file in os.listdir(root):
        if 'condition' not in file:
            continue
        tissue_df = pd.read_table(os.path.join(root, file))
        try:
            high_conf_df = tissue_df[tissue_df['Confidence']=='High']
        except:
            print(f'Issue with {file}. Columns = {tissue_df.columns}')
            continue
        target_gene_names = target_gene_names.union(set(high_conf_df['gene']))
        tf_gene_names = tf_gene_names.union(set(high_conf_df['TF']))
        high_conf_dict = high_conf_df.groupby('TF')['gene'].apply(set).to_dict()

        # Add this file to the TF-to-Target dictionary
        for tf, target in high_conf_dict.items():
            tf_gene_name_to_target_gene_name.setdefault(tf, []).extend(target)

    gene_names = tf_gene_names.union(target_gene_names)
    tf_gene_name_to_target_gene_name = switch_dictset_to_dictlist(tf_gene_name_to_target_gene_name)
    json.dump(tf_gene_name_to_target_gene_name, open('output/gene2gene/tf_gene_name_to_target_gene_name.json','w'))
    json.dump(list(gene_names), open('output/gene2gene/gene_names_tf.json','w'))
    

def map_gene_name_to_id_via_mygeneinfo():
    # Use MyGeneInfo Python client to map gene names to IDs
    mg = mygene.MyGeneInfo()

    gene_names = json.load(open('output/gene2gene/gene_names_tf.json'))
    params = {'qterms':gene_names, 
              'species':9606,
              'scopes':'symbol',
              'entrezonly':True,
              'as_dataframe':True}
    mg_df = mg.querymany(**params)
    mg_df = mg_df[mg_df['symbol'].isin(gene_names)]
    gene_name_to_gene_id = dict(zip(mg_df['symbol'], [int(id_) for id_ in mg_df['_id']]))
    json.dump(gene_name_to_gene_id, open('output/gene2gene/gene_name_to_gene_id.json','w'))
    
    
def map_tf_to_target_with_ids():
    # Import
    with open('output/gene2gene/tf_gene_name_to_target_gene_name.json') as fin:
         tf_gene_name_to_target_gene_name = json.load(fin)
    with open('output/gene2gene/gene_name_to_gene_id.json') as fin:
        gene_name_to_gene_id = json.load(fin)
    with open('output/protein2gene/all_entrez2uniprot.json') as fin:
        gene_id_to_protein_id = json.load(fin)

    #print(len(tf_gene_name_to_target_gene_name))
    #print(len(gene_name_to_gene_id))
    #print(len(gene_id_to_protein_id))

    tf_gene_id_to_target_gene_id = {}
    tf_protein_id_to_target_gene_id = {}
    for tf_name, target_names in tf_gene_name_to_target_gene_name.items():

        # Map TF gene ID & protein ID -to- target gene ID

        # TF gene ID and protein ID
        try:
            tf_id = str(gene_name_to_gene_id[tf_name])
        except:
            continue
        try:
            tf_protein_ids = gene_id_to_protein_id[tf_id]
        except:
            tf_protein_ids = -1

        # Target gene ID
        for target_name in target_names:
            try:
                target_id = gene_name_to_gene_id[target_name]
            except:
                continue

            # TF -to- Target
            tf_gene_id_to_target_gene_id.setdefault(tf_id, set()).add(str(target_id))
            if tf_protein_ids != -1:
                for tf_protein_id in tf_protein_ids:
                    tf_protein_id_to_target_gene_id.setdefault(tf_protein_id, set()).add(str(target_id))    
    
    # Export
    tf_gene_id_to_target_gene_id = switch_dictset_to_dictlist(tf_gene_id_to_target_gene_id)
    with open('output/gene2gene/tf_gene_id_to_target_gene_id.json','w') as fout:
        json.dump(tf_gene_id_to_target_gene_id, fout)
        
    tf_protein_id_to_target_gene_id = switch_dictset_to_dictlist(tf_protein_id_to_target_gene_id)
    with open('output/protein2gene/tf_protein_id_to_target_gene_id.json','w') as fout:
        json.dump(tf_protein_id_to_target_gene_id, fout)
        
        
def export_transcription_factor_edges():
    tf_gene_id_to_target_gene_id = json.load(open('output/gene2gene/tf_gene_id_to_target_gene_id.json'))
    tf_protein_id_to_target_gene_id = json.load(open('output/protein2gene/tf_protein_id_to_target_gene_id.json'))

    output_edgefile_onerel_noweight(outpath='output/gene2gene/Gene_(Entrez)_targets_Gene_(Entrez).csv',
                                    columns=['Gene (Entrez)','Gene (Entrez)', 'Relationship'],
                                    dictionary=tf_gene_id_to_target_gene_id,
                                    rel='-genes_transcription_factor_targets->', 
                                    prefix_col1='Entrez:', 
                                    prefix_col2='Entrez:',
                                    edges_to_use_folder=False,
                                   )

    output_edgefile_onerel_noweight(outpath='output/protein2gene/Protein_(UniProt)_targets_Gene_(Entrez).csv',
                                    columns=['Protein (UniProt)','Gene (Entrez)', 'Relationship'],
                                    dictionary=tf_protein_id_to_target_gene_id,
                                    rel='-transcription_factor_targets->', 
                                    prefix_col1='UniProt:', 
                                    prefix_col2='Entrez:',
                                   )
    print('Exported transcription factor edges')
    
if __name__ == '__main__':    
    #download_transcription_factor_data_from_grndb()
    map_tf_name_to_target_name()
    map_gene_name_to_id_via_mygeneinfo()
    map_tf_to_target_with_ids()
    export_transcription_factor_edges()