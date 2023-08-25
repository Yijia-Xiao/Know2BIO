import os
import pandas as pd
import json
import requests as req
from biomedkg_utils import switch_dictset_to_dictlist, output_edgefile_onerel_noweight

def download_bgee_diff_expr_in_anatomy():
    urls = ['https://bgee.org/ftp/bgee_v13_2/download/calls/diff_expr_calls/Homo_sapiens_diffexpr-anatomy-simple.tsv.zip']

    for url in urls:
        filename = url.split('/')[-1] 
        #os.remove(filename[:-4])
        os.system(f'wget -N -P input/ {url}')       # Download file
        os.system(f'unzip input/{filename}')        # Unzip file
        os.system(f'mv filename[:-4] input/{filename[:-4]}') 


        
def align_gene_ids_entrez_to_ensembl():

    # Get Entrez Gene IDs
    #entrez_df = pd.read_csv('output/nodes/genes_nodes.csv')[['Gene (Entrez)']]
    #entrez_genes = sorted([gene.split(':')[1] for gene in list(entrez_df['Gene (Entrez)'])[1:]])
    #entrez_genes = sorted(list(json.load(open('output/protein2gene/all_entrez2uniprot.json')).keys()))[1:]
    os.system('wget -N -P input/ https://www.genenames.org/cgi-bin/download/custom?col=md_eg_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit')
    entrez_df = pd.read_csv('input/custom?col=md_eg_id')
    id_col = list(entrez_df.columns)[0]
    entrez_genes = [str(gene) for gene in sorted(list(entrez_df[id_col]))]
    entrez_genes = str(",".join(entrez_genes))
    
    
    # Use MyGene.Info to get Ensembl IDs
    headers = {'content-type': 'application/x-www-form-urlencoded'}
    params = 'q=%s&scopes=entrezgene&fields=ensembl&species=human'%entrez_ids
    res = req.post('http://mygene.info/v3/query', \
                    data=params, \
                    headers=headers).json()
    
    # Align/map Ensembl to Entrez
    ensembl_is_entrez = dict()

    for r in res:

        # Entrez ID
        try:entrez_id = r['_id']
        except:continue

        # Ensembl ID
        try:
            if type(r['ensembl']) == list:
                for entry in r['ensembl']:
                    ensembl = entry['gene']
                    ensembl_is_entrez[ensembl] = entrez_id # 1 to 1

            elif type(r['ensembl']) == dict:
                ensembl = r['ensembl']['gene']
                ensembl_is_entrez[ensembl] = entrez_id # 1 to 1
        except: continue
            
    json.dump(ensembl_is_entrez, open('output/gene2gene/ensembl2entrez.json','w'))    

    
def gene_expressed_in_anatomy(df, gene2anat, anat2gene, uberon2mesh, ensembl_is_entrez):
    '''
    FUNCTION:
    - Gene -[over/under]expressed_in- Anatomy relationships
    - Uses Bgee data to make a dictionary of those relationships^
    
    PARAMS:
    - df: dataframe of under or overexpressed genes to anatomy from Bgee
    - gene2anat: dictionary mappings expressed genes to anatomy 
    - anat2gene: dictionary mappings anatomy to expressed genes
    - uberon2mesh: dictionary mapping UBERON -is- MeSH
    - ensembl_is_entrez: dictionary mapping Ensembl -is- Entrez ID
    '''
    for i in range(len(df)):
        
        try:
            gene = df['Gene ID'].iloc[i]
            entrez_gene = ensembl_is_entrez[gene]
            
            anatomy = df['Anatomical entity ID'].iloc[i]
            mesh_anatomies = uberon2mesh[anatomy]
            
            # Gene -[over/under]expressed in- Anatomy
            for mesh_anatomy in mesh_anatomies:
                gene2anat.setdefault(entrez_gene, set()).add(mesh_anatomy)
                anat2gene.setdefault(mesh_anatomy, set()).add(entrez_gene)

        except:
            continue
                
    
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
    
    
def map_gene_to_anatomy():
    # Differential expression dataframes
    diffex_df = pd.read_table('input/Homo_sapiens_diffexpr-anatomy-simple.tsv.zip', compression='zip')
    diffex_df = diffex_df[diffex_df['Call quality']=='high quality']
    underexp_df = diffex_df[diffex_df['Differential expression']=='under-expression']
    overexp_df = diffex_df[diffex_df['Differential expression']=='over-expression']

    try:
        uberon2mesh = json.load(open('output/anatomy2anatomy/uberon2mesh.json'))
    except:
        download_uberon_anatomy()
        align_uberon_to_mesh_anatomy()
        uberon2mesh = json.load(open('output/anatomy2anatomy/uberon2mesh.json'))
        
    # Gene -overexpressed in- Anatomy
    underexp_gene2anat, underexp_anat2gene = dict(), dict()
    gene_expressed_in_anatomy(underexp_df, underexp_gene2anat, underexp_anat2gene, 
                              uberon2mesh, ensembl_is_entrez)

    # Gene -underexpressed in- Anatomy
    overexp_gene2anat, overexp_anat2gene  = dict(), dict()
    gene_expressed_in_anatomy(overexp_df, overexp_gene2anat, overexp_anat2gene, 
                              uberon2mesh, ensembl_is_entrez)

    print('Under', len(underexp_gene2anat), len(underexp_anat2gene))
    print('Over', len(overexp_gene2anat), len(overexp_anat2gene))

    # Prepare to export
    underexp_gene2anat = switch_dictset_to_dictlist(underexp_gene2anat)
    overexp_gene2anat = switch_dictset_to_dictlist(overexp_gene2anat)

    # Export
    json.dump(underexp_gene2anat, open('output/gene2anatomy/underexpressed_gene2anat.json','w'))
    json.dump(overexp_gene2anat, open('output/gene2anatomy/overexpressed_gene2anat.json','w'))

    # Output underexpressed genes
    output_edgefile_onerel_noweight(
        outpath = 'output/gene2anatomy/edges_gene-underexpressed_in-anatomy.csv', 
        columns = ['Gene (Entrez)', 'Anatomy (MeSH)','Relationship'], 
        dictionary = underexp_gene2anat, 
        rel = '-underexpressed_in->', 
        prefix_col1='Entrez:', 
        prefix_col2='MeSH_Anatomy:')

    # Output overexpressed genes
    output_edgefile_onerel_noweight(
        outpath = 'output/gene2anatomy/edges_gene-overexpressed_in-anatomy.csv', 
        columns = ['Gene (Entrez)', 'Anatomy (MeSH)','Relationship'], 
        dictionary = overexp_gene2anat, 
        rel = '-overexpressed_in->', 
        prefix_col1='Entrez:', 
        prefix_col2='MeSH_Anatomy:')
    
    
if __name__ == '__main__':
    download_bgee_diff_expr_in_anatomy()
    align_gene_ids_entrez_to_ensembl()
    ensembl_is_entrez = json.load(open('output/gene2gene/ensembl2entrez.json'))         
    map_gene_to_anatomy()