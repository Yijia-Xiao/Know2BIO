import os
import shutil
import requests
import gzip
import pandas as pd
import csv
import json
import zipfile
import mygene
from disease_to_disease import align_all_mondo_omim_umls
from biomedkg_utils import switch_dictset_to_dictlist, output_edgefile_onerel_noweight
from biomed_apis import *

def download_disgenet_disease_gene_associations():
    # Download and uncompress the gene-disease association file from DisGeNET
    file_url = 'https://www.disgenet.org/static/disgenet_ap1/files/downloads/all_gene_disease_associations.tsv.gz'
    compressed_file_path = 'input/all_gene_disease_associations.tsv.gz'
    uncompressed_file_path = 'input/all_gene_disease_associations.tsv'

    try:
        os.remove(compressed_file_path)
    except:
        pass

    try:
        os.remove(uncompressed_file_path)
    except:
        pass

    response = requests.get(file_url, stream=True)
    with open(compressed_file_path, 'wb') as file:
        shutil.copyfileobj(response.raw, file)

    with gzip.open(compressed_file_path, 'rb') as compressed_file:
        with open(uncompressed_file_path, 'wb') as uncompressed_file:
            shutil.copyfileobj(compressed_file, uncompressed_file)

    os.remove(compressed_file_path)
    
    
def map_gene_disease_assocation_from_disgenet():
    # Entrez Gene ID -associated_with-> MeSH Disease
    umls2mesh = json.load(open('output/otherMappings/umls2mesh.json'))
    entrez_gene2mesh_disease_fromdisgenet, mesh_disease2entrez_gene_fromdisgenet = dict(), dict()

    evidenceCutoff = 0.5
    score_thresh = 0.06 # https://think-lab.github.io/d/105/
    ev, noev = 0, 0
    with open('output/gene2disease/edges_gene-ASSOCIATED-WITH->disease_disgenet.csv','w',newline='') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Gene (Entrez)', 'Disease (MeSH)', 'Relationship','Weight'])

        with open('input/all_gene_disease_associations.tsv') as fin:
            for i, line in enumerate(fin):
                if i == 0:
                    continue
                line = line.strip().split('\t')

                # Disease
                umls_dis = line[4]
                try: 
                    mesh_dises = umls2mesh[umls_dis]
                except: 
                    continue

                # Gene
                entrez_id = line[0]

                # Evidence
                score = float(line[9])

                # Check evidence, write row
                try:
                    if score >= score_thresh:
                        for mesh_dis in mesh_dises:
                            writer.writerow(['Entrez:'+str(entrez_id), 
                                             'MeSH_Disease:'+mesh_dis, 
                                             '-associated_with-',
                                             score])
                        ev += 1
                except:
                    noev += 1

    print(ev, 'Gene-Disease Associations with Evidence')
    print(noev,'without Evidence')

    df = pd.read_csv('output/gene2disease/edges_gene-ASSOCIATED-WITH->disease_disgenet.csv')
    df.to_csv('output/edges/edges_gene-ASSOCIATED-WITH->disease_disgenet.csv', index=False)

    
def download_pharmgkb_gene_disease_data():
    zip_url = 'https://api.pharmgkb.org/v1/download/file/data/relationships.zip'
    zip_file_path = 'input/pharmgkb_relationships.zip'

    response = requests.get(zip_url)
    with open(zip_file_path, 'wb') as fout:
        fout.write(response.content)

    with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
        zip_ref.extractall('input')

    os.remove(zip_file_path)
    
    
def get_pharmgkb_df():
    pgkb_df = pd.read_table('input/relationships.tsv')
    pgkb_df = pd.read_table('input/relationships.tsv')
    pgkb_df = pgkb_df.loc[(pgkb_df['Entity1_type']=='Gene') \
                         &(pgkb_df['Entity2_type']=='Disease')]
    return pgkb_df    


def align_gene_names_to_ids_via_mygeneinfo(pgkb_df):
    # Use MyGeneInfo Python client to map gene names to IDs
    mg = mygene.MyGeneInfo()
    gene_names = list(set(pgkb_df['Entity1_name']))

    params = {'qterms':gene_names, 
              'species':9606,
              'scopes':'symbol',
              'entrezonly':True,
              'as_dataframe':True}
    mg_df = mg.querymany(**params)
    mg_df = mg_df[mg_df['symbol'].isin(gene_names)]
    gene_name_to_gene_id = dict(zip(mg_df['symbol'], [int(id_) for id_ in mg_df['_id']]))


# Needs to be modified. This is the most "manual" of the curation, and it's not getting any
# associations anymore. It only got about 450-660 before, so it was a very minor contribution. 
def map_disease_names_to_ids(pgkb_df):
    # Name --> UMLS ID
    name2umls = dict()
    with open('input/MRCONSO.RRF') as fin:
        for i,line in enumerate(fin):
            line = line.strip().split('|')
            if 'ENG' not in line:
                continue
            umls = line[0]
            name = line[14]
            name2umls.setdefault(name.lower(), set()).add(umls)

    print(len(name2umls))

    # All PharmGKB-provided disease names
    disease_names = [dis.lower() for dis in pgkb_df['Entity2_name']]
    print(len(disease_names))
    
    # Removing non-disease PharmGKB-provided 'disease names'
    not_diseases = ['Elderly Adult','Pregnancy','Recurrence','cessation','retreatment failure', 'time above therapeutic range',
     'time below therapeutic range', 'time in therapeutic range','time to achieve stable dose', 'time to delivery',
     'time to relapse', 'time to therapeutic inr', 'tolerance']
    not_diseases = [not_dis.lower() for not_dis in not_diseases]
    for notdis in not_diseases:
        disease_names.remove(notdis)


    # DiseaseName -> Disease UMLS ID
    dis_name_to_mesh = dict()
    for dis_name in disease_names:
        try:
            dis_umlses = name2umls[dis_name]
            dis_meshes = list()
            for dis_umls in dis_umlses:
                dis_meshes += umls2mesh[dis_umls]
            for dis_mesh in dis_meshes:
                dis_name_to_mesh.setdefault(dis_name, set()).add(dis_mesh)
        except:
            continue
    print(len(dis_name_to_mesh), 'diseases mapped from name to MeSH')

    dis_name_to_mesh_str = dict()
    for k,v in dis_name_to_mesh.items():
        dis_name_to_mesh_str[k] = str(v)
        

def finish_mapping_pharmgkb_gene_disease_associations():
    # Gene ID --> Disease ID

    # Filter by / preprocess disease
    pgkb_df['Entity2_name'] = pgkb_df['Entity2_name'].str.lower()
    pgkb_df['Entity2_name'].replace(dis_name_to_mesh_str, inplace=True)
    pgkb_df = pgkb_df.replace(np.nan, 'not a num')
    pgkb_df = pgkb_df.drop(pgkb_df[pgkb_df['Entity2_name'].str.contains('}') == False].index)
    pgkb_df['Entity2_name'] = pgkb_df['Entity2_name'].str.upper()

    # Filter by / preprocess gene
    pgkb_df['Entity1_name'].replace(gene_name_to_gene_id, inplace=True)

    pgkb_df = pgkb_df[['Entity1_name','Entity2_name','Association']]
    pgkb_df.columns = ['Gene', 'Disease', 'Relationship']

    # Turn into dictionaries
    pgkb_pos = pgkb_df[pgkb_df['Relationship'] == 'associated'].set_index('Gene')['Disease'].to_dict()
    pgkb_neg = pgkb_df[pgkb_df['Relationship'] == 'not associated'].set_index('Gene')['Disease'].to_dict()
    pgkb_pos_temp = dict()
    for gene, diseases in pgkb_pos.items():
        diseases = ast.literal_eval(diseases)
        for disease in diseases:
            pgkb_pos_temp.setdefault(gene, list()).append(disease)
    pgkb_neg_temp = dict()
    for gene, diseases in pgkb_neg.items():
        diseases = ast.literal_eval(diseases)
        for disease in diseases:
            pgkb_neg_temp.setdefault(gene, list()).append(disease)
            
            
    with open('output/gene2disease/edges_gene-ASSOCIATED_WITH->disease_pharmgkb.csv','w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Gene (Entrez)','Disease (MeSH)','Relationship'])
    
    for gene, diseases in pgkb_pos_temp.items():
        for disease in diseases:
            writer.writerow(['Entrez:'+str(gene),'MeSH_Disease:'+disease,'-associated_with-'])
    os.system('cp output/gene2disease/edges_gene-ASSOCIATED_WITH->disease_pharmgkb.csv output/edges/edges_gene-ASSOCIATED_WITH->disease_pharmgkb.csv')
        
        
def map_clinvar_gene_disease():
    os.system('wget -N -P input/ https://ftp.ncbi.nlm.nih.gov/pub/clinvar/gene_condition_source_id')
    try:
        # MONDO-is-OMIM-is-MeSH 
        mondo2omim = json.load(open('output/disease2disease/mondo2omim.json'))
        omim2mesh = json.load(open('output/disease2disease/omim2mesh.json'))
    except:
        align_all_mondo_omim_umls()
        # MONDO-is-OMIM-is-MeSH 
        mondo2omim = json.load(open('output/disease2disease/mondo2omim.json'))
        omim2mesh = json.load(open('output/disease2disease/omim2mesh.json'))
    for omim,mesh in omim2mesh.copy().items():
        #del omim2mesh[omim]
        omim = omim.replace('OMIM_','')
        omim2mesh[omim] = mesh


    entrez_gene2omim_disease, omim_disease2entrez_gene = dict(), dict()

    for i, line in enumerate(open('input/gene_condition_source_id')):
        line = line.split('\t')
        if i == 0:
            continue


        '''Entrez Gene'''
        entrez_gene_id = line[0]
        # If no gene ID, skip this line
        if entrez_gene_id == '':
            continue


        '''OMIM Disease ID'''
        disease_omim_ids = line[7]

        # If no OMIM
        if disease_omim_ids == '':
            # Try to get OMIM via MONDO2OMIM
            disease_other_id_type = line[5]        
            disease_other_ids = line[6]

            if disease_other_id_type == 'MONDO':
                mondo_id = disease_other_ids.replace(':','_')
                if mondo_id in mondo2omim:
                    disease_omim_ids = mondo2omim[mondo_id]
                # Otherwise, skip this line
                else:
                    continue


        ''' Entrez Gene - OMIM Disease '''
        if type(disease_omim_ids) != list:
            disease_omim_ids = [disease_omim_ids]

        for disease_omim_id in disease_omim_ids:
            entrez_gene2omim_disease.setdefault(entrez_gene_id, set()).add(disease_omim_id)
            omim_disease2entrez_gene.setdefault(disease_omim_id, set()).add(entrez_gene_id)

    print(len(omim_disease2entrez_gene), 'OMIM Diseases to Entrez Genes')
    print(len(entrez_gene2omim_disease), 'Entrez Genes to OMIM Diseases')


    mesh_disease2entrez_gene, entrez_gene2mesh_disease = dict(), dict()

    for omim, entrez_genes in omim_disease2entrez_gene.items():
        if omim in omim2mesh:

            # MeSH Disease
            mesh_diseases = omim2mesh[omim]
            for mesh_disease in mesh_diseases:

                # Gene
                for entrez_gene in entrez_genes:

                    # MeSH Disease - Gene
                    mesh_disease2entrez_gene.setdefault(mesh_disease, set()).add(entrez_gene)
                    entrez_gene2mesh_disease.setdefault(entrez_gene, set()).add(mesh_disease)

    print(len(mesh_disease2entrez_gene), 'MeSH Diseases to Entrez Genes')
    print(len(entrez_gene2mesh_disease), 'Entrez Genes to MeSH Diseases')

    json.dump(switch_dictset_to_dictlist(entrez_gene2mesh_disease), open('output/gene2disease/entrez_gene2mesh_disease_fromclinvar.json','w'))
    json.dump(switch_dictset_to_dictlist(mesh_disease2entrez_gene), open('output/gene2disease/mesh_disease2entrez_gene_fromclinvar.json','w'))
    json.dump(switch_dictset_to_dictlist(entrez_gene2omim_disease), open('output/gene2disease/entrez_gene2omim_disease_fromclinvar.json','w'))
    json.dump(switch_dictset_to_dictlist(omim_disease2entrez_gene), open('output/gene2disease/omim_disease2entrez_gene_fromclinvar.json','w'))

    output_edgefile_onerel_noweight('output/gene2disease/edges_gene-ASSOCIATED-WITH->disease_clinvar.csv',
                                   ['Gene (Entrez)', 'Disease (MeSH)', 'Relationship'],
                                    entrez_gene2mesh_disease,
                                    '-associated_with-',
                                   'Entrez:',
                                   'MeSH_Disease:',
                                   edges_to_use_folder=False,
                                   )


    gene2protein = json.load(open('output/protein2gene/all_entrez2uniprot.json'))

    protein2mesh_disease = dict()

    for gene, mesh_diseases in entrez_gene2mesh_disease.items():
        try: 
            proteins = gene2protein[gene]
            for protein in proteins:
                for mesh_disease in mesh_diseases:
                    protein2mesh_disease.setdefault(protein, set()).add(mesh_disease)
        except:
            continue
            
        
        
def download_clingen_gene_disease():
    os.system('wget -N -P input/ https://search.clinicalgenome.org/kb/gene-validity/download')   
    clingen_df = pd.read_csv('input/download',header=4)
    clingen_df.drop(clingen_df.index[0], inplace=True)
    clingen_df = clingen_df.reset_index()
    return clingen_df
    
# credit to Alexander for these ClinGen functions
def extract_mapping_table(results):
    mapping_table = {}
    for entry in results['results']:
        from_ = entry['from']
        to_ = entry['to']['primaryAccession']
        if from_ not in mapping_table:
            mapping_table[from_] = []
        mapping_table[from_] += [to_]
    return mapping_table


def uniprot_convert_hgnc_to_uniprot_gene(genes):
    # send job
    job_id = submit_id_mapping_UniProtAPI(
                      from_db = 'HGNC',
                      #to_db = 'Gene_Name', 
                      to_db = 'UniProtKB', 
                      ids = genes)

    # check job until it is finished
    if check_id_mapping_results_ready_UniProtAPI(job_id):
        link = get_id_mapping_results_link_UniProtAPI(job_id)
        results = get_id_mapping_results_search_UniProtAPI(link)
    
    mapping_table = extract_mapping_table(results)
    
    print("%d results. %d out of %d sucessfully mapped (%d failed)"%(len(results['results']),
                                                        len(mapping_table),
                                                        len(genes),
                                                        len(results['failedIds'])))
    return mapping_table


# need to prepare clingen to have h,r,t encoding UniProt gene -[ClinGen_Classification]-> MeSH Term
def prepare_clingen(hgnc_to_uniprot, mondo_to_mesh, clingen_df, debug=False):
    mapped_clingen_dict = {x:[] for x in ["head","relation","tail","weight","edge_type"]}
    
    for hgnc_gene, classification, mondo_disease in zip(clingen_df['GENE ID (HGNC)'],
                                                        clingen_df['CLASSIFICATION'],
                                                        clingen_df['DISEASE ID (MONDO)']):
        if hgnc_gene in hgnc_to_uniprot and mondo_disease in mondo_to_mesh:
            uniprot_list = hgnc_to_uniprot[hgnc_gene]
            mesh_terms = mondo_to_mesh[mondo_disease]
                        
            for uniprot in uniprot_list:
                for disease in mesh_terms:
                
                    h = uniprot
                    r = "ClinGen_Classification:%s"%(classification)
                    t = disease
                    mapped_clingen_dict['head'] += ['UniProt:'+h]
                    mapped_clingen_dict['relation'] += [r]
                    mapped_clingen_dict['tail'] += ["MeSH_Disease:"+t]
                    mapped_clingen_dict['weight'] += [1]
                    mapped_clingen_dict['edge_type'] += ["ClinGen"]
    
    mapped_clingen_df = pd.DataFrame(mapped_clingen_dict)
    
    if debug:
        # print statistics
        hgnc_genes = set(clingen_df['GENE ID (HGNC)'])
        mondo_diseases = set(clingen_df['DISEASE ID (MONDO)'])
        uniprot_proteins = set(mapped_clingen_df['head'])
        mesh_terms = set(mapped_clingen_df['tail'])
        print("Originally %d HGNC genes and %d MONDO diseases"%(len(hgnc_genes), len(mondo_diseases)))
        print("%d edges between %d UniProt proteins and %d MeSH terms"%(mapped_clingen_df.shape[0],
                                                                        len(uniprot_proteins),
                                                                        len(mesh_terms)))

    return mapped_clingen_df


def remove_low_confidence_clingen_associations(mapped_clingen_df):
    bad_indices = list()
    bad_confidences = ['ClinGen_Classification:Disputed','ClinGen_Classification:No Known Disease Relationship',\
                       'ClinGen_Classification:Refuted','ClinGen_Classification:Limited']

    for bad_confidence in bad_confidences:
        bad_indices += list(np.where(mapped_clingen_df['relation']==bad_confidence)[0])
    bad_indices = sorted(bad_indices)
    
    confident_clingen_df = mapped_clingen_df.drop(bad_indices)
    confident_clingen_df = mapped_clingen_df.drop(bad_indices)
    
    return confident_clingen_df


def export_clingen_gene_disease_associations(confident_clingen_df):

    prot_ids_to_gene_ids = json.load(open('output/protein2gene/all_uniprot2entrez.json','r'))

    with open('output/gene2disease/gene_to_disease_clingen.csv','w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Gene (Entrez)','Disease (MeSH)','Relationship'])
        for i in range(0,len(confident_clingen_df)):
            protein = confident_clingen_df['head'].iloc[i].split('UniProt:')[1]
            try: 
                genes = prot_ids_to_gene_ids[protein]
                disease = confident_clingen_df['tail'].iloc[i]
                for gene in genes:
                    writer.writerow(['Entrez:'+str(gene),
                                     disease,
                                     '-associated_with-'])
            except:
                continue

if __name__ == '__main__':
    ''' DisGeNET '''
    download_disgenet_disease_gene_associations()
    map_gene_disease_assocation_from_disgenet()    


    ''' PharmGKB ''' #(not working, but it only got a very small amount of edges when it did work)
    #download_pharmgkb_gene_disease_data()
    #pgkb_df = get_pharmgkb_df()
    #align_gene_names_to_ids_via_mygeneinfo(pgkb_df)
    #map_disease_names_to_ids(pgkb_df)


    ''' ClinVar '''
    map_clinvar_gene_disease()


    ''' ClinGen '''
    clingen_df = download_clingen_gene_disease()
    mondo_to_mesh = json.load(open("output/disease2disease/mondo2mesh.json","r"))
    hgnc_to_uniprot = uniprot_convert_hgnc_to_uniprot_gene(clingen_df['GENE ID (HGNC)'])
    mapped_clingen_df = prepare_clingen(hgnc_to_uniprot, mondo_to_mesh, clingen_df, debug=True)
    confident_clingen_df = remove_low_confidence_clingen_associations(mapped_clingen_df)
    export_clingen_gene_disease_associations(confident_clingen_df)