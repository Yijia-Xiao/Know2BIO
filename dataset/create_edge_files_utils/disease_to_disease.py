import requests
import os
import gzip
import shutil
import csv
import pandas as pd
import json
from datetime import datetime
import time
import xml.etree.ElementTree as ET
import urllib.request
from biomedkg_utils import output_edgefile_onerel_noweight, \
switch_dictset_to_dictlist, switch_dictlist_to_dictset



def download_disgenet_curated_disease_to_disease():
    url = 'https://www.disgenet.org/static/disgenet_ap1/files/downloads/disease_to_disease_CURATED.tsv.gz'
    folder = 'input/'

    # Download
    response = requests.get(url)
    file_path = os.path.join(folder, 'disease_to_disease_CURATED.tsv.gz')
    with open(file_path, 'wb') as fin:
        fin.write(response.content)

    # Extract 
    output_file_path = os.path.join(folder, 'disease_to_disease_CURATED.tsv')
    with gzip.open(file_path, 'rb') as gz_file:
        with open(output_file_path, 'wb') as fout:
            fout.write(gz_file.read())

    # Remove the gzipped file
    os.remove(file_path)

    os.system('wget -N -P input/ https://www.disgenet.org/static/disgenet_ap1/files/downloads/disease_to_disease_CURATED.tsv.gz')
    try:
        os.remove('input/disease_to_disease_CURATED.tsv')
    except:
        pass
    os.system('gunzip input/disease_to_disease_CURATED.tsv.gz')
    

def map_disease_shares_genes_with_umls_disease(dda_df):
    '''Shared genes'''
    dda_gene_df = dda_df.query('jaccard_genes != 0')

    # UMLS Disease IDs
    folder = 'output/disease2disease/'
    umls_edge_file = 'Disease_(UMLS)_shares-genes_Disease_(UMLS).csv'
    umls_edge_path = os.path.join(folder, umls_edge_file)
    with open(umls_edge_path,'w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Disease (UMLS)','Disease (UMLS)','Relationship','Weight'])

        for dis1, dis2, jac_gene in zip(dda_gene_df['diseaseId1'], 
                                        dda_gene_df['diseaseId2'], 
                                        dda_gene_df['jaccard_genes']):
            if jac_gene > 0:
                writer.writerow(['UMLS_Disease:'+dis1, 
                                 'UMLS_Disease:'+dis2, 
                                 '-diseases_share_genes-', jac_gene])

    dda_gene_df = pd.read_csv(umls_edge_path).drop_duplicates()
    dda_gene_df.to_csv('output/edges/'+umls_edge_file, index=False)

    num_disease_pairs = len(dda_gene_df)
    unique_diseases = len(set(dda_gene_df['Disease (UMLS)']).union(set(dda_gene_df['Disease (UMLS).1'])))
    print(num_disease_pairs, 'disease pairs with shared genes')
    print(unique_diseases, 'diseases that share genes with another disease')
    

def map_disease_shares_variants_with_umls_disease(dda_df):
    '''Shared variants'''
    dda_variant_df = dda_df.query('jaccard_variant != 0')

    folder = 'output/disease2disease/'
    umls_edge_file = 'Disease_(UMLS)_shares-variants_Disease_(UMLS).csv'
    umls_edge_path = os.path.join(folder, umls_edge_file)
    with open(umls_edge_path,'w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Disease (UMLS)','Disease (UMLS)','Relationship','Weight'])

        for dis1, dis2, jac_variant in zip(dda_variant_df['diseaseId1'], 
                                           dda_variant_df['diseaseId2'], 
                                           dda_variant_df['jaccard_variant']):
            if jac_variant > 0: 
                writer.writerow(['UMLS_Disease:'+dis1, 
                                 'UMLS_Disease:'+dis2, 
                                 '-diseases_share_variants-', 
                                 jac_variant])


    dda_variant_df = pd.read_csv(umls_edge_path).drop_duplicates()
    dda_variant_df.to_csv('output/edges/'+umls_edge_file, index=False)

    num_disease_pairs = len(dda_variant_df)
    unique_diseases = len(set(dda_variant_df['Disease (UMLS)']).union(set(dda_variant_df['Disease (UMLS).1'])))
    print(num_disease_pairs, 'disease pairs with shared variants')
    print(unique_diseases, 'diseases that share variants with another disease')
    
    
def map_disease_shares_genes_with_disease():
    # MeSH Disease IDs
    umls2mesh = json.load(open('output/otherMappings/umls2mesh.json'))
    file = 'Disease_(MeSH)_sharesgenes_Disease_(MeSH).csv'
    outpath = os.path.join('output/disease2disease/',file)
    with open(outpath, 'w', newline='') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Disease (MeSH)','Disease (MeSH)','Relationship','Weight'])
        with open('input/disease_to_disease_CURATED.tsv') as fin:
            for i,line in enumerate(fin):
                if i == 0:
                    continue
                line = line.strip().split('\t')
                try:
                    dis1s = umls2mesh[line[0]]
                    dis2s = umls2mesh[line[1]]
                    jac_gen = float(line[12])
                    if jac_gen <= 0:
                        continue
                except:
                    continue
                for dis1 in dis1s:
                    for dis2 in dis2s:
                        writer.writerow(['MeSH_Disease:'+dis1, 'MeSH_Disease:'+dis2, '-diseases_share_genes-', jac_gen])

    df = pd.read_csv(outpath)
    df.to_csv(os.path.join('output/edges', file), index=False)
    df.to_csv(os.path.join('output/edges_to_use', file), index=False)


def map_disease_shares_variants_with_disease():
    # MeSH Disease IDs
    umls2mesh = json.load(open('output/otherMappings/umls2mesh.json'))
    file = 'Disease_(MeSH)_sharesvariants_Disease_(MeSH).csv'
    outpath = os.path.join('output/disease2disease/',file)
    with open(outpath, 'w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Disease (MeSH)',
                         'Disease (MeSH)',
                         'Relationship',
                         'Weight'])

        with open('input/disease_to_disease_CURATED.tsv') as fin:
            for i,line in enumerate(fin):
                if i == 0:
                    continue
                line = line.strip().split('\t')
                try:
                    dis1s = umls2mesh[line[0]]
                    dis2s = umls2mesh[line[1]]
                    jacvar = float(line[12])
                    if jacvar <= 0:
                        continue
                except:
                    continue
                for dis1 in dis1s:
                    for dis2 in dis2s:
                        writer.writerow(['MeSH_Disease:'+dis1, 
                                         'MeSH_Disease:'+dis2, 
                                         '-diseases_share_variants-', 
                                         jacvar])

    df = pd.read_csv(outpath)
    df.to_csv(os.path.join('output/edges', file), index=False)
    df.to_csv(os.path.join('output/edges_to_use', file), index=False)
    
    
def download_mesh_xml():
    year = datetime.now().year
    url = f'https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/desc{year}.xml'
    dest = f'input/desc{year}.xml'
    urllib.request.urlretrieve(url, dest);

def parse_mesh_xml():
    year = datetime.now().year
    tree = ET.parse(f'input/desc{year}.xml')
    root = tree.getroot()   
    return root


def align_disease_identifiers_within_mesh(root):
    name2id, id2name, id2tree, tree2id = dict(), dict(), dict(), dict()
    all_tree_numbers = list()

    for ele in root:
        try:
            # MeSH Tree Number
            tree_numbers = ele.find('TreeNumberList').findall('TreeNumber')

            # If disease
            for tree_number in tree_numbers:
                if tree_number.text.startswith(('C','F03')): # exclude animal diseases C22?

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
        except:
            continue        

    all_tree_numbers = sorted(all_tree_numbers)
    tree2id = dict(sorted(tree2id.items()))

    for k,v in name2id.copy().items():
        name2id[k] = list(name2id[k])
        
    # MeSH Term -[is]- MeSH ID
    with open('output/disease2disease/meshterm-IS-meshid.json','w') as fout:
        json.dump(name2id, fout)
        
    return id2tree, all_tree_numbers


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
    file = 'Disease_(MeSH_Tree)_2_Disease_(MeSH_Tree).csv'
    outpath = os.path.join('output/disease2disease/',file)
    output_edgefile_onerel_noweight(outpath = outpath,
                                    columns = ['Disease (MeSH Tree)','Disease (MeSH Tree)','Relationship'],
                                    dictionary = tree2tree,
                                    rel = '-subclass_of->',
                                    prefix_col1 = 'MeSH_Tree_Disease:',
                                    prefix_col2 = 'MeSH_Tree_Disease:')

    tree_df = pd.read_csv(outpath).drop_duplicates()
    tree_df.to_csv(os.path.join('output/edges/', file), index=False)
    tree_df.to_csv(os.path.join('output/edges_to_use/', file), index=False)
    print(f'{len(tree_df)} MeSH Tree Number intra relationships')
    
    
def export_mesh_tree_and_id_alignment_disease(id2tree):
    # MeSH Tree Number -[is]- MeSH ID
    file = 'Disease_(MeSH)_2_Disease_(MeSH_Tree).csv'
    outpath = os.path.join('output/disease2disease/',file)
    output_edgefile_onerel_noweight(outpath = outpath,
                                    columns = ['Disease (MeSH)','Disease (MeSH Tree)','Relationship'],
                                    dictionary = id2tree,
                                    rel = '-is-',
                                    prefix_col1 = 'MeSH_Disease:',
                                    prefix_col2 = 'MeSH_Tree_Disease:')

    tree_to_id = pd.read_csv(outpath).drop_duplicates()
    tree_to_id.to_csv(os.path.join('output/edges/', file), index=False)
    tree_to_id.to_csv(os.path.join('output/edges_to_use/', file), index=False)
    NUM_MAPPED_MESH_IDS = len(set(tree_to_id['Disease (MeSH)']))
    NUM_MAPPED_MESH_TREE_NUMS = len(set(tree_to_id['Disease (MeSH Tree)']))
    print(f'Aligned {NUM_MAPPED_MESH_IDS} MeSH IDs to '+\
          f'{NUM_MAPPED_MESH_TREE_NUMS} MeSH Tree Numbers '+\
          f'\nTotal of {len(tree_to_id)} relationships, averaging '+\
          f'{round(NUM_MAPPED_MESH_TREE_NUMS/NUM_MAPPED_MESH_IDS, 2)} Tree Numbers per ID')
    
    
def align_disease_ontology_to_mesh():
    os.system('wget -N -P input/ https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/HumanDO.obo')   
    doid2mesh = dict()
    with open('input/HumanDO.obo') as fin:
        IDs = []
        meshes = []
        for i,line in enumerate(fin):
            if i < 29:
                continue

            # New DOID term
            if '[Term]' in line:
                if len(IDs) > 0:
                    meshes = set(meshes)
                    for ID in IDs:
                        for mesh in meshes:
                            doid2mesh.setdefault(ID,set()).add(mesh)
                    IDs = []
                    meshes = []

            # DOID IDs
            elif 'id' in line.split(':')[0]:
                IDs.append(line.split(':')[1].strip()+':'+line.split(':')[2].strip())

            # Get UMLS CUIs
            if 'xref: UMLS_CUI:' in line:
                umls = line.split(':')[2].strip()
                try: meshes += umls2mesh[umls]
                except: continue

            # Get MeSH via UMLS CUIs
            if 'xref: MESH:' in line:
                mesh = line.split(':')[2].strip()
                meshes.append(mesh)


    for k in doid2mesh.copy():
        doid2mesh[k] = list(doid2mesh[k])

    with open('output/disease2disease/doid2mesh.json','w') as fout:
        json.dump(doid2mesh, fout)
        

        
def get_meshes(temp_mesh, meshes):
    if type(temp_mesh) == list:
        temp_meshes = temp_mesh
    else:
        temp_meshes = [temp_mesh]
    
    for mesh in temp_meshes:
        meshes.append(mesh)
    
    return meshes


def align_all_mondo_omim_umls():
    os.system('wget -N -P input/ http://purl.obolibrary.org/obo/mondo/mondo-with-equivalents.json')
    mondo2others = json.load(open('input/mondo-with-equivalents.json'))

    omim2mesh, omim2umls, omim2doid, omim2mondo = dict(), dict(), dict(), dict()
    mondo2mesh, mondo2umls, mondo2doid, mondo2omim = dict(), dict(), dict(), dict()
    mesh2mondo = dict()
    mondos = set()

    entries = mondo2others['graphs'][0]['equivalentNodesSets']
    for i in range(len(entries)):
        node = mondo2others['graphs'][0]['equivalentNodesSets'][i]

        umls, doid, mesh, omim = '','','',''
        for orig_id_link in node['nodeIds']:
            id_link = orig_id_link.split('/')
            if '/umls/' in id_link:
                umls = id_link[-1]
            elif 'DOID' in id_link[-1]:
                doid = id_link[-1]
            elif 'MESH_' in id_link[-1]:
                mesh = id_link[-1]
            elif 'MONDO_' in id_link[-1]:
                mondo = id_link[-1]
            elif 'omim' in orig_id_link:
                omim = id_link[-1]

        if omim != '':
            if mesh != '':
                omim2mesh.setdefault(omim, set()).add(mesh)
            if umls != '':
                omim2umls.setdefault(omim, set()).add(umls)
            if doid != '':
                omim2doid.setdefault(omim, set()).add(doid)
            if mondo != '':
                omim2mondo.setdefault(omim, set()).add(mondo)

        if mondo != '':
            if mesh != '':
                mondo2mesh.setdefault(mondo, set()).add(mesh)
                mesh2mondo.setdefault(mesh, set()).add(mondo)
            if umls != '':
                mondo2umls.setdefault(mondo, set()).add(umls)
            if doid != '':
                mondo2doid.setdefault(mondo, set()).add(doid)
            if omim != '':
                mondo2omim.setdefault(mondo, set()).add(omim)
            mondos.add(mondo)

    mondo_list = mondos

    #omim2mesh = dict()
    mesh2omim = dict()
    for omim, mondos in omim2mondo.items():

        # MONDO
        for mondo in mondos:
            if mondo in mondo2mesh:

                # MeSH
                meshes = mondo2mesh[mondo]
                for mesh in meshes:
                    mesh = mesh.split('MESH_')[1]

                    # OMIM - MeSH
                    omim2mesh.setdefault(omim, set()).add(mesh)
                    mesh2omim.setdefault(mesh, set()).add(omim)

    print(len(omim2mesh), len(mesh2omim))
    json.dump(switch_dictset_to_dictlist(omim2mesh), open('output/disease2disease/omim2mesh.json','w'))
    json.dump(switch_dictset_to_dictlist(mondo2omim), open('output/disease2disease/mondo2omim.json','w'))

    mondos = {mondo.replace('_',':') for mondo in mondo_list}
    mondo_list = list(mondos)
    mesh2mondo = {mesh.replace('MESH_',''):{mondo.replace('_',':') for mondo in mondos} for mesh,mondos in mesh2mondo.items()}
    mondo2mesh = {mondo.replace('_',':'):{mesh.replace('MESH_','') for mesh in meshes} for mondo,meshes in mondo2mesh.items()}


    ''' OMIM - UMLS - MeSH '''
    omim2mesh = json.load(open('output/disease2disease/omim2mesh.json'))
    omim2mesh = switch_dictlist_to_dictset(omim2mesh)

    # UMLS - OMIM
    umls2omim = dict()
    with open('input/MRCONSO.RRF') as fin:
        for line in fin:
            line = line.strip().split('|')
            umls = line[0]
            ontology = line[11]

            if 'OMIM' in ontology:
                disease = (line[12] == 'PT')
                if disease:
                    omim_id = line[13]
                    umls2omim.setdefault(umls, set()).add(omim_id)


    # OMIM - UMLS - MeSH
    print(len(omim2mesh), 'OMIM-MeSH Before')
    print(len(mesh2omim), 'MeSH-OMIM Before')
    umls2mesh = json.load(open('output/otherMappings/umls2mesh.json'))

    # UMLS
    for umls, omims in umls2omim.items():
        if umls in umls2mesh:

            # MeSH
            meshes = umls2mesh[umls]

            # OMIM
            for omim in omims:

                # OMIM - MeSH
                mesh2omim.setdefault(mesh, set()).add(omim)
                omim2mesh.setdefault(omim, set()).add(mesh)

    print(len(omim2mesh), 'OMIM-MeSH After')
    print(len(mesh2omim), 'MeSH-OMIM After')





    BATCH_SIZE = 500
    for batch_i in range(0,int(len(mondo_list)/BATCH_SIZE)):
        start_ind = batch_i*500
        end_ind = (batch_i+1)*500
        if end_ind > len(mondo_list):
            end_ind = len(mondo_list)-1

        print(start_ind, end_ind)

        params = {'q': ','.join(mondo_list[start_ind:end_ind])}
        res = requests.post('http://mydisease.info/v1/query', params).json()
        time.sleep(3)

        for i in range(0,len(res)):

            ''' MONDO '''
            mondo = res[i]['query']


            ''' OMIM '''
            try:
                omim = res[i]['disease_ontology']['xrefs']['omim']
            except:
                try:
                    omim = res[i]['mondo']['xrefs']['omim']
                except:
                    try:
                        omim = res[i]['disgenet']['xrefs']['omim'] 
                    except:
                        try:
                            omim = res[i]['hpo']['omim']
                        except:
                            continue


            ''' OMIM - MONDO '''
            if type(omim) == list:
                omims = omim
            else:
                omims = [omim]

            for omim in omims:
                omim2mondo.setdefault(omim, set()).add(mondo)
                mondo2omim.setdefault(mondo, set()).add(omim)


            ''' MeSH '''
            meshes = list()
            try: meshes = get_meshes(res[i]['umls']['mesh']['preferred'], meshes)
            except: pass

            try: meshes = get_meshes(res[i]['mondo']['xrefs']['mesh'], meshes)
            except: pass

            try: meshes = get_meshes(res[i]['ctd']['mesh'], meshes)
            except: pass

            try: meshes = get_meshes(res[i]['disease_ontology']['xrefs']['mesh'], meshes)
            except: pass

            try: meshes = get_meshes(res[i]['disgenet']['xrefs']['mesh'], meshes)
            except: pass

            try: meshes = get_meshes(res[i]['umls']['mesh']['non-preferred'], meshes)
            except: pass

            ''' MeSH - MONDO '''
            for mesh in meshes:
                mesh2mondo.setdefault(mesh, set()).add(mondo)
                mondo2mesh.setdefault(mondo, set()).add(mesh)


    '''OMIM MeSH'''
    # MONDO
    for mondo, meshes in mondo2mesh.items():
        if mondo in mondo2omim:

            # OMIM
            omims = mondo2omim[mondo]

            # MeSH
            for mesh in meshes:

                # OMIM - MeSH
                omim2mesh.setdefault(omim, set()).add(mesh)
                mesh2omim.setdefault(mesh, set()).add(omim)

    print(len(omim2mesh), len(mesh2omim))
    print(len(mondo2mesh), len(mesh2mondo))

    json.dump(switch_dictset_to_dictlist(omim2mesh), open('output/disease2disease/omim2mesh.json','w'))
    json.dump(switch_dictset_to_dictlist(mesh2omim), open('output/disease2disease/mesh2omim.json','w'))
    json.dump(switch_dictset_to_dictlist(mondo2mesh), open('output/disease2disease/mondo2mesh.json','w'))
    json.dump(switch_dictset_to_dictlist(mesh2mondo), open('output/disease2disease/mesh2mondo.json','w'))
    
    
    
if __name__ == '__main__':
    # DisGeNET: diseases share genes and variants
    download_disgenet_curated_disease_to_disease()
    dda_df = pd.read_table('input/disease_to_disease_CURATED.tsv')[['diseaseId1','diseaseId2','jaccard_genes','jaccard_variant']]
    map_disease_shares_genes_with_umls_disease(dda_df)
    map_disease_shares_variants_with_umls_disease(dda_df)
    map_disease_shares_genes_with_disease()
    map_disease_shares_variants_with_disease()
    umls2mesh = json.load(open('output/otherMappings/umls2mesh.json'))


    # MeSH Tree and IDs
    download_mesh_xml()
    root = parse_mesh_xml()
    id2tree, all_tree_numbers = align_disease_identifiers_within_mesh(root)
    map_and_export_mesh_tree_number_hierarchy(all_tree_numbers)
    export_mesh_tree_and_id_alignment_disease(id2tree)

    # DOID to MeSH
    align_disease_ontology_to_mesh()

    # MONDO, OMIM, MeSH
    align_all_mondo_omim_umls()