import requests 
import json 
import pandas as pd 
import numpy as np 
import os 
import time
from datetime import datetime
import zipfile
import xml.etree.ElementTree as ET 
from parse_xml import *
from biomedkg_utils import *
from biomed_apis import *
direct = 'output/compound2compound/'



def create_empty_folders_for_data():
    '''
    Make the folders if they're not there
    '''
    paths = ['input','output','input/KEGG','input/SMPDB']
    for path in paths:
        if not os.path.exists(path):
            os.mkdir(path)

    folders = ['anatomy2anatomy','compound2drugclass','compound2compound',
               'compound2disease','compound2gene','compound2pathway',
               'compound2protein','compound2reaction','disease2anatomy',
               'disease2disease','disease2pathway','gene2anatomy','gene2disease',
               'gene2gene','go2go','otherMappings','pathway2gene','pathway2pathway',
               'pathway2protein','pathway2reaction','protein2gene','protein2go',
               'protein2protein','protein2reaction','reaction2reaction',
               'edges','edges_to_use', 'node_features', 'node_features/sequences',
               'node_features/natural_language_names','node_features/structures']

    for folder in folders:
        path = os.path.join('output', folder)
        if not os.path.exists(path):
            print(f'Making folder {path}')
            os.mkdir(path)

            
def extract_mr_conso_file():
    print('Extracting MRCONSO')
    year = datetime.now().year
    zip_path = f'input/umls-{year}AA-mrconso.zip'
    dest = 'input/'
    
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(dest)
        
    
def align_db_to_mesh_via_x_to_mesh_set(db2x, x2mesh, db2mesh, mesh2db):
    for db, xes in db2x.items():
        try:
            for x in xes:
                meshes = x2mesh[x]
                for mesh in meshes:
                    db2mesh.setdefault(db, set()).add(mesh)
                    mesh2db.setdefault(mesh, set()).add(db)
        except:
            continue        
    db2mesh = switch_dictset_to_dictlist(db2mesh)
    mesh2db = switch_dictset_to_dictlist(mesh2db)
    json.dump(db2mesh, open('output/compound2compound/db2mesh.json','w'))
    json.dump(mesh2db, open('output/compound2compound/mesh2db.json','w'))        
        
        
def align_db_to_mesh_via_x_to_mesh_list(db2x, x2mesh, db2mesh, mesh2db):
    for db, xes in db2x.items():
        try:
            for x in xes:
                meshes = x2mesh[x]
                for mesh in meshes:
                    db2mesh.setdefault(db, []).append(mesh)
                    mesh2db.setdefault(mesh, []).append(db)
        except:
            continue        
    json.dump(db2mesh, open('output/compound2compound/db2mesh.json','w'))
    json.dump(mesh2db, open('output/compound2compound/mesh2db.json','w'))

        
def extract_mesh_umls_drugbank_from_mr_conso(direct='output/compound2compound/',
                                             verbose=True):
    print('Extracting the MeSH, UMLS, and DrugBank interconnections/mappings')
    umls2mesh, mesh2umls = dict(), dict()
    umls2db, db2umls = dict(), dict()
    db2mesh, mesh2db = dict(), dict()
    atc2umls, umls2atc = dict(), dict()

    with open('input/MRCONSO.RRF') as fin:
        for line in fin:
            line = line.strip().split('|')
            umls = line[0]
            ontology = line[11]

            if 'MSH' == ontology:
                mesh = line[10]
                umls2mesh.setdefault(umls,set()).add(mesh)
                mesh2umls.setdefault(mesh,set()).add(umls)

            if 'DRUGBANK' in line:    
                db_ind = (line.index('DRUGBANK'))-2
                db = line[db_ind]
                db2umls.setdefault(db,set()).add(umls)    
                umls2db.setdefault(umls,set()).add(db)
                
            if 'ATC' in line[11]:
                atc = line[13]
                umls2atc.setdefault(umls, set()).add(atc)
                atc2umls.setdefault(atc, set()).add(umls)
      
    
    align_db_to_mesh_via_x_to_mesh_set(db2umls, umls2mesh, db2mesh, mesh2db)

    '''Export'''
    db2umls_json = switch_dictset_to_dictlist(db2umls)
    umls2db_json = switch_dictset_to_dictlist(umls2db)
    umls2mesh_json = switch_dictset_to_dictlist(umls2mesh)
    mesh2umls_json = switch_dictset_to_dictlist(mesh2umls)
    db2mesh_json = switch_dictset_to_dictlist(db2mesh)
    mesh2db_json = switch_dictset_to_dictlist(mesh2db)
    
    umls2atc_json = switch_dictset_to_dictlist(umls2atc)
    atc2umls_json = switch_dictset_to_dictlist(atc2umls)
    
    json.dump(db2umls_json, open(direct+'db2umls.json','w'))
    json.dump(umls2db_json, open(direct+'umls2db.json','w'))
    json.dump(umls2atc_json, open(direct+'umls2atc.json','w'))
    json.dump(atc2umls_json, open(direct+'atc2umls.json','w'))
    json.dump(umls2mesh_json, open('output/otherMappings/umls2mesh.json','w'))
    json.dump(mesh2umls_json, open('output/otherMappings/mesh2umls.json','w'))
    json.dump(db2mesh_json, open('output/compound2compound/db2mesh.json','w'))
    json.dump(mesh2db_json, open('output/compound2compound/mesh2db.json','w'))
    if verbose:
        print(f'UMLS to MeSH mappings: {len(umls2mesh)}')
        print(f'MeSH to UMLS mappings: {len(mesh2umls)}')
        print(f'DrugBank to UMLS mappings: {len(db2umls)}')
        print(f'UMLS to DrugBank mappings: {len(umls2db)}')
        print(f'DrugBank to MeSH mappings: {len(db2mesh)}')
        print(f'MeSH to DrugBank mappings: {len(mesh2db)}')
        

def extract_drugbank_xml():
    zip_path = 'input/drugbank_all_full_database.xml.zip'
    destination_folder = 'input/'

    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(destination_folder)

        
def parse_drugbank_xml():
    print('Parsing DrugBank XML')
    tree = ET.parse('input/full database.xml')
    root = tree.getroot()
    return root 


def convert_and_export(x_to_y, out_path):
    dict_val_type =  type(list(x_to_y.values())[0])
    x_to_y_copy = x_to_y.copy()
    if dict_val_type == set:
        x_to_y_copy = switch_dictset_to_dictlist(x_to_y_copy)
    with open(out_path, 'w') as fout:
        json.dump(x_to_y_copy, fout)


def map_drugbank_to_external_ids(root, direct='output/compound2compound/',
                                 verbose=False):
    db2cas, db2unii, db2onts = dict(), dict(), dict()
    cas2db, unii2db = dict(), dict()      
    
    dbs = set()


    '''DrugBank -is _______'''
    db2atc, atc2db = dict(), dict()
    db2chembl, db2chebi, db2pubchemcomp, db2chemspider, db2pubchemsub = dict(), dict(), dict(), dict(), dict()
    db2keggdrug, db2keggcompound, db2oldttd, db2uniprot, db2pharmgkb, db2wikipedia = dict(), dict(), dict(), dict(), dict(), dict()
    chembl2db, chebi2db, pubchemcomp2db, chemspider2db, pubchemsub2db = dict(), dict(), dict(), dict(), dict()
    keggdrug2db, keggcompound2db, oldttd2db, uniprot2db, pharmgkb2db, wikipedia2db = dict(), dict(), dict(), dict(), dict(), dict()


    for ele in root:
        db = ParseXML.getID(ele)
        dbs.add(db)

        # DrugBank -is- CAS
        cas = ele.find('{http://www.drugbank.ca}cas-number').text
        if cas != None:
            db2cas[db] = [cas]                   # 1 DrugBank is 1 CAS
            cas2db.setdefault(cas,set()).add(db) # 1 CAS is 1+ DrugBank (usually 1)

        # DrugBank -is- UNII
        unii = ParseXML.getUNII(ele)
        if(unii != None):
            db2unii[db] = [unii]                   # 1 to 1
            unii2db.setdefault(unii,set()).add(db) # 1 to 1

        # DrugBank -is- ATC
        atcs = ParseXML.getATCCode(ele)
        if(atcs != []):
            db2atc[db] = atcs
            for atc in atcs:
                atc2db.setdefault(atc,set()).add(db)

        # DrugBank -is- Other ontologies
        ext_id_dicts = ParseXML.getExternalIDs(ele)
        db2onts[db] = ext_id_dicts


    # Ontology to DrugBank Dictionary
    allonts_d, allonts = dict(), set()
    for drug, onts in db2onts.items():
        for ont in onts:
            allonts_d.setdefault(ont,set()).add(drug)   
            allonts.add(ont)

    for db, onts in db2onts.items():
        # DrugBank -is- ChEMBL
        if 'ChEMBL' in onts:
            db2chembl[db] = onts['ChEMBL']
            for chembl in onts['ChEMBL']:
                chembl2db.setdefault(chembl,set()).add(db)

        # DrugBank -is- ChEBI
        if 'ChEBI' in onts:
            db2chebi[db] = onts['ChEBI']
            for chebi in onts['ChEBI']:
                chebi2db.setdefault(chebi,set()).add(db)

        # DrugBank -is- PubChem Compound
        if 'PubChem Compound' in onts:
            db2pubchemcomp[db] = onts['PubChem Compound']
            for pubchemcomp in onts['PubChem Compound']:
                pubchemcomp2db.setdefault(pubchemcomp,set()).add(db)

        # DrugBank -is- PubChem Substance
        if 'PubChem Substance' in onts:
            db2pubchemsub[db] = onts['PubChem Substance']
            for pubchemsub in onts['PubChem Substance']:
                pubchemsub2db.setdefault(pubchemsub,set()).add(db)


        # DrugBank -is- ChemSpider
        if 'ChemSpider' in onts:
            db2chemspider[db] = onts['ChemSpider']
            for chemspider in onts['ChemSpider']:
                chemspider2db.setdefault(chemspider,set()).add(db)

        # DrugBank -is- KEGG Drug
        if 'KEGG Drug' in onts:
            db2keggdrug[db] = onts['KEGG Drug']
            for keggdrug in onts['KEGG Drug']:
                keggdrug2db.setdefault(keggdrug,set()).add(db)

        # DrugBank -is- KEGG Compound
        if 'KEGG Compound' in onts:
            db2keggcompound[db] = onts['KEGG Compound']
            for keggcomp in onts['KEGG Compound']:
                keggcompound2db.setdefault(keggcomp,set()).add(db)

        # DrugBank -is- Therapeutic Targets Database
        if 'Therapeutic Targets Database' in onts:
            db2oldttd[db] = onts['Therapeutic Targets Database']
            for oldttd in onts['Therapeutic Targets Database']:
                oldttd2db.setdefault(oldttd,set()).add(db)

        # DrugBank -is- UniProt
        if 'UniProtKB' in onts:
            db2uniprot[db] = onts['UniProtKB']
            for uniprot in onts['UniProtKB']:
                uniprot2db.setdefault(uniprot,set()).add(db)

        # DrugBank -is- PharmGKB
        if 'PharmGKB' in onts:
            db2pharmgkb[db] = onts['PharmGKB']  
            for pharm in onts['PharmGKB']:
                pharmgkb2db.setdefault(pharm,set()).add(db)

        # DrugBank -is- Wikipedia
        if 'Wikipedia' in onts:
            db2wikipedia[db] = onts['Wikipedia']
            for wiki in onts['Wikipedia']:
                wikipedia2db.setdefault(wiki,set()).add(db)


    '''Print numbers'''  
    # Display the number of DrugBank-provided xref alignments
    print('DrugBank-provided external mappings:')
    for ont,drugs in allonts_d.items():
        print(len(drugs),'DrugBank -is-',ont)
    print(len(db2cas), 'DrugBank -is- CAS')
    print(len(cas2db), 'CAS -is/are- DrugBank')
    print(len(db2unii), 'DrugBank -is- UNII')
    print(len(db2atc), 'DrugBank -is/are- ATC')
    print(len(atc2db), 'ATC -is/are- DrugBank')


    '''Export alginments'''    
    convert_and_export(chembl2db, os.path.join(direct, 'chembl2db.json'))
    convert_and_export(chebi2db, os.path.join(direct, 'chebi2db.json'))
    convert_and_export(db2chebi, os.path.join(direct, 'db2chebi.json'))
    convert_and_export(db2chembl, os.path.join(direct, 'db2chembl.json'))
    convert_and_export(pubchemcomp2db, os.path.join(direct, 'pubchemcomp2db.json'))
    convert_and_export(db2pubchemcomp, os.path.join(direct, 'db2pubchemcomp.json'))
    convert_and_export(db2pubchemsub, os.path.join(direct, 'db2pubchemsub.json'))
    convert_and_export(chemspider2db, os.path.join(direct, 'chemspider2db.json'))
    convert_and_export(pubchemsub2db, os.path.join(direct, 'pubchemsub2db.json'))
    convert_and_export(keggdrug2db, os.path.join(direct, 'keggdrug2db.json'))
    convert_and_export(keggcompound2db, os.path.join(direct, 'keggcompound2db.json'))
    convert_and_export(wikipedia2db, os.path.join(direct, 'wikipedia2db.json'))
    convert_and_export(oldttd2db, os.path.join(direct, 'oldttd2db.json'))
    convert_and_export(uniprot2db, os.path.join(direct, 'uniprot2db.json'))
    convert_and_export(pharmgkb2db, os.path.join(direct, 'pharmgkb2db.json'))
    convert_and_export(unii2db, os.path.join(direct, 'unii2db.json'))
    convert_and_export(db2unii, os.path.join(direct, 'db2unii.json'))
    convert_and_export(unii2db, os.path.join(direct, 'unii2db.json'))
    convert_and_export(atc2db, os.path.join(direct, 'atc2db.json'))
    convert_and_export(db2atc, os.path.join(direct, 'db2atc.json'))
    convert_and_export(cas2db, os.path.join(direct, 'cas2db.json'))
    convert_and_export(db2cas, os.path.join(direct, 'db2cas.json'))
    convert_and_export(db2oldttd, os.path.join(direct, 'db2oldttd.json'))
    
    output_db_to_atc_mapping()
    
    return dbs

def output_db_to_atc_mapping():
    db_to_atc = json.load(open('output/compound2compound/db2atc.json'))

    with open('output/compound2compound/drugbank_drug_to_atc_class.csv','w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Compound (DrugBank)','Compound (ATC)','-is->'])
        rel = '-is-'

        for drugbank, atcs in db_to_atc.items():
            for atc in atcs:
                writer.writerow([drugbank, atc, rel])
                
                
def align_db_to_mesh_via_atc():
    db2atc = json.load(open('output/compound2compound/db2atc.json'))
    db2mesh = json.load(open('output/compound2compound/db2mesh.json'))
    mesh2db = json.load(open('output/compound2compound/mesh2db.json'))
    db2umls = json.load(open('output/compound2compound/db2umls.json'))
    umls2mesh = json.load(open('output/otherMappings/umls2mesh.json'))
    atc2umls = json.load(open('output/compound2compound/atc2umls.json'))
    
    db2umls = switch_dictlist_to_dictset(db2umls)
    
    for db, atcs in db2atc.items():
        for atc in atcs:
            try:
                umlses = atc2umls[atc]
            except:
                continue
            for umls in umlses:
                db2umls.setdefault(db, set()).add(umls)

    db2mesh = switch_dictlist_to_dictset(db2mesh)
    mesh2db = switch_dictlist_to_dictset(mesh2db)
    for db, umlses in db2umls.items():
        for umls in umlses:
            try:
                meshes = umls2mesh[umls]
            except:
                continue
            for mesh in meshes:
                db2mesh.setdefault(db, set()).add(mesh)
                mesh2db.setdefault(mesh, set()).add(db)

    db2mesh = switch_dictset_to_dictlist(db2mesh)
    mesh2db = switch_dictset_to_dictlist(mesh2db)
    json.dump(db2mesh, open('output/compound2compound/db2mesh.json','w'))
    json.dump(mesh2db, open('output/compound2compound/mesh2db.json','w'))
    

                
def map_unii_to_mesh_via_mesh_api(direct='output/compound2compound/'):
    '''MeSH API/Scraper'''
    print('Mapping UNII to MeSH')
    db2unii = json.load(open(direct+'db2unii.json'))
    unii2db = json.load(open(direct+'unii2db.json'))
    
    multiprocess_a_dict_of_lists_values(db2unii, unii2mesh_meshapi)

    ### Make unii2mesh: Merge the temp files into a dictionary
    unii2mesh, mesh2unii = dict(), dict()
    procs = cpu_count()

    for b_id in range(procs):
        temp_path = direct+'temp_unii2mesh'+str(b_id)+'.txt'
        with open(temp_path) as fin:
            for line in fin:

                # UNII
                unii = line.split('|')[0].strip()
                # MeSH
                mesh = line.split('|')[1].strip()    
                # UNII -is- MeSH
                unii2mesh.setdefault(unii, set()).add(mesh)
                # MeSH -is- UNII
                mesh2unii.setdefault(mesh, set()).add(unii)

        os.remove(temp_path)

    print('UNII to MeSH:', len(unii2mesh))
    print('MeSH to UNII:', len(mesh2unii))
    
    return unii2mesh, mesh2unii            
   
     
import time # delaying request because the server doesn't like too many at once. It would
# be better to make this a batched POST request. Check if other APIs can be improved this way
def unii2mesh_mycheminfoapi(b_id, unii_batch):
    '''
    FUNCTION:
    - Uses the MyChem.info API to align UNII with MeSH IDs
    
    PARAMS:
    - b_id: batch ID (integer number ranging from 0 to number of processors)
    - unii_batch: UNII IDs used as input to be mapped to MeSH
    '''
    
    temp_path = 'output/compound2compound/temp_unii2mesh_mycheminfo'+str(b_id)+'.txt'
    open(temp_path,'w').write('')
    tot = str(len(unii_batch))
    
    for count, unii in enumerate(unii_batch):  
        unii = unii[0]
        time.sleep(1)
        
        # Print progress(only prints one batch, threadsafe)
        if b_id == 1:
            print('Progress of Batch 1: '+str(count)+'/'+tot, end='\r')     
        
        # Use API 
        try: 
            res = requests.get('http://mychem.info/v1/chem/%s'%unii).json()
            mesh = res['umls']['mesh'].strip()
        except: 
            try: mesh = res['drugcentral']['xrefs']['mesh_descriptor_ui'].strip()
            except: continue     
        
        # Save results
        with open(temp_path,'a') as fout:
            fout.write(unii+'|'+mesh+'\n')
    
def map_unii_to_mesh_via_mycheminfo(db2unii, unii2mesh, mesh2unii, 
                                    direct='output/compound2compound/'):
    '''MyChem.info API'''
    print('Mapping UNII to MeSH via MyChem.info')
    multiprocess_a_dict_of_lists_values(db2unii, unii2mesh_mycheminfoapi)

    ### Merge the temp unii2mesh files into a dictionary
    procs = cpu_count()
    for b_id in range(procs):
        temp_path = direct+'temp_unii2mesh_mycheminfo'+str(b_id)+'.txt'
        with open(temp_path) as fin:
            for line in fin:

                unii = line.split('|')[0].strip()
                mesh = line.split('|')[1].strip()            
                unii2mesh.setdefault(unii, set()).add(mesh)
                mesh2unii.setdefault(mesh, set()).add(unii)

        os.remove(temp_path)

    print(len(unii2mesh))
    print(len(mesh2unii))

    json.dump(switch_dictset_to_dictlist(unii2mesh), open(direct+'unii2mesh.json','w'))
    json.dump(switch_dictset_to_dictlist(mesh2unii), open(direct+'mesh2unii.json','w'))
    json.dump(switch_dictset_to_dictlist(db2unii), open(direct+'db2unii.json','w'))
    
    
def map_drugbank_to_mesh_drugbank_to_cas_via_inxight_drugs(dbs, 
                                                           direct='output/compound2compound/'):
    '''Inxight Drugs API/Scraper'''
    multiprocess_a_list(dbs, db2meshorcas_inxightAPI)

    '''Merge Files: DrugBank-is-CAS, DrugBank-is-MeSH'''
    db2mesh = json.load(open(os.path.join(direct, 'db2mesh.json')))
    mesh2db = json.load(open(os.path.join(direct, 'mesh2db.json')))
    db2mesh = switch_dictlist_to_dictset(db2mesh)
    mesh2db = switch_dictlist_to_dictset(mesh2db)


    mesh2cas, cas2mesh = dict(), dict()
    db2cas = json.load(open(os.path.join(direct,'db2cas.json')))
    cas2db = json.load(open(os.path.join(direct,'cas2db.json')))
    db2cas = switch_dictlist_to_dictset(db2cas)
    cas2db = switch_dictlist_to_dictset(cas2db)

    procs = cpu_count()
    for b_id in range(procs):

        ### DrugBank-is-MeSH
        temp_db2mesh_path = os.path.join(direct,f'temp_db2mesh_inxightapi{str(b_id)}.txt')
        with open(temp_db2mesh_path) as fin1:
            for line in fin1:            
                db = line.split('|')[0].strip()
                mesh = line.split('|')[1].strip() 
                db2mesh.setdefault(db, set()).add(mesh)
                mesh2db.setdefault(mesh, set()).add(db)
            os.remove(temp_db2mesh_path)

        ### DrugBank-is-CAS
        temp_db2cas_path = os.path.join(direct,f'db2cas_inxight{str(b_id)}.txt')
        with open(temp_db2cas_path) as fin2:
            for line in fin2:
                db = line.split('|')[0].strip()
                cases = line.split('|')[1:]
                for cas in cases:
                    cas = cas.strip()
                    cas2db.setdefault(cas, set()).add(db)
                    db2cas.setdefault(db, set()).add(cas)
            os.remove(temp_db2cas_path)

    # Export
    
    json.dump(switch_dictset_to_dictlist(db2cas), open(os.path.join(direct,'db2cas.json'),'w'))
    json.dump(switch_dictset_to_dictlist(cas2db), open(os.path.join(direct,'cas2db.json'),'w'))
    json.dump(switch_dictset_to_dictlist(db2mesh), open(os.path.join(direct,'db2mesh.json'),'w'))
    json.dump(switch_dictset_to_dictlist(mesh2db), open(os.path.join(direct,'mesh2db.json'),'w'))
    print('DrugBank-to-MeSH:', len(db2mesh))
    print('MeSH-to-DrugBank:', len(mesh2db))
    print('DrugBank-to-CAS:', len(db2cas))
    
    
def map_cas_to_mesh_via_mesh_api(db2cas, cas2mesh, mesh2cas):
    '''MeSH API/Scraper: CAS-is-MeSH'''
    multiprocess_a_dict_of_lists_values(db2cas, cas2mesh_meshapi)

    ### Make CAS-is-MeSH: Merge the temp cas2mesh files into a dictionary
    cas2mesh = switch_dictlist_to_dictset(cas2mesh)
    mesh2cas = switch_dictlist_to_dictset(mesh2cas)

    procs = cpu_count()
    for b_id in range(procs):
        temp_path = os.path.join(direct, f'temp_cas2mesh{str(b_id)}.txt')
        with open(temp_path) as fin:
            for line in fin:

                cas = line.split('|')[0].strip()
                mesh = line.split('|')[1].strip()            
                cas2mesh.setdefault(cas, set()).add(mesh)
                mesh2cas.setdefault(mesh, set()).add(cas)

        os.remove(temp_path)

    print('CAS-to-MeSH', len(cas2mesh))
    print('MeSH-to-CAS', len(mesh2cas))

    # Prepare to export
    cas2mesh = switch_dictset_to_dictlist(cas2mesh)
    mesh2cas = switch_dictset_to_dictlist(mesh2cas)

    # Export
    json.dump(cas2mesh, open(os.path.join(direct,'cas2mesh.json'),'w'))
    json.dump(mesh2cas, open(os.path.join(direct+'mesh2cas.json'),'w'))
    
    
def align_chebi_to_mesh(direct='output/compound2compound/'):
    chebi2db = json.load(open(direct+'chebi2db.json'))

    # Add items
    chebi2mesh, mesh2chebi = dict(), dict()
    db2mesh = json.load(open(os.path.join(direct, 'db2mesh.json')))
    key2mesh_via_db(chebi2db, db2mesh, chebi2mesh, mesh2chebi)
    print('ChEBI-is-MeSH', len(chebi2mesh))

    # Export
    json.dump(switch_dictset_to_dictlist(chebi2mesh), open(direct+'chebi2mesh.json','w'))
    json.dump(switch_dictset_to_dictlist(mesh2chebi), open(direct+'mesh2chebi.json','w'))
    

def align_cas_to_mesh(direct='output/compound2compound/'):
    try:
        db2cas = json.load(open(os.path.join(direct, 'db2cas.json')))
        cas2db = json.load(open(os.path.join(direct, 'cas2db.json')))
        cas2mesh = json.load(open(direct+'cas2mesh.json'))
    except:
        pass
    cas2mesh = switch_dictlist_to_dictset(cas2mesh)
    print('CAS-is-MeSH',len(cas2mesh), len(mesh2cas))

    # Add items
    key2mesh_via_db(cas2db, db2mesh, cas2mesh, mesh2cas)
    print('CAS-is-MeSH',len(cas2mesh), len(mesh2cas))

    # Export
    json.dump(switch_dictset_to_dictlist(cas2mesh), open(direct+'cas2mesh.json','w'))
    json.dump(switch_dictset_to_dictlist(mesh2cas), open(direct+'mesh2cas.json','w'))
    
    
def align_inchi_and_smiles_to_drugbank():
    # INSTRUCTIONS: Download from https://go.drugbank.com/releases/5-1-9/downloads/all-structure-links
    # You need a DrugBank account to download this file
    # ! unzip 'input/drugbank_all_structure_links.csv.zip'
    # moved the file to input
    structure_df = pd.read_csv('input/structure links.csv')

    '''InchiKey & SMILES -is- DrugBank'''
    db2inchikey, inchikey2db = dict(), dict()
    db2smiles, smiles2db = dict(), dict()

    for i in range(len(structure_df)):

        # DrugBank ID
        db = structure_df['DrugBank ID'].iloc[i]

        # Inchi Key
        inchikey = structure_df['InChIKey'].iloc[i]
        if type(inchikey) == str:
            db2inchikey.setdefault(db,set()).add(inchikey)
            inchikey2db.setdefault(inchikey,set()).add(db)

        # SMILES
        smiles = structure_df['SMILES'].iloc[i]
        if type(smiles) == str:
            db2smiles.setdefault(db,set()).add(smiles)
            smiles2db.setdefault(smiles, set()).add(db)

    # Export
    json.dump(switch_dictset_to_dictlist(db2inchikey), open(direct+'db2inchikey.json','w'))
    json.dump(switch_dictset_to_dictlist(inchikey2db), open(direct+'inchikey2db.json','w'))
    json.dump(switch_dictset_to_dictlist(db2smiles), open(direct+'db2smiles.json','w'))
    json.dump(switch_dictset_to_dictlist(smiles2db), open(direct+'smiles2db.json','w'))


def align_inchi_and_smiles_to_mesh(direct='output/compound2compound/'):
    '''InchiKey & SMILES -is- MeSH'''
    db2mesh = json.load(open(os.path.join(direct, 'db2mesh.json')))
    
    # InchiKey -is- MeSH
    inchis = list(json.load(open(direct+'inchikey2db.json')).keys())
    inchikey2db = json.load(open(os.path.join(direct, 'inchikey2db.json')))
    inchikey2mesh, mesh2inchikey = dict(), dict()
    inchi2mesh_mycheminfoAPI(inchis, inchikey2mesh, mesh2inchikey)
    key2mesh_via_db(inchikey2db, db2mesh, inchikey2mesh, mesh2inchikey)
    print('InchiKey-is-MeSH', len(inchikey2mesh), len(mesh2inchikey))

    # SMILES -is- MeSH
    smiles2mesh, mesh2smiles = dict(), dict()
    smiles2db = json.load(open('output/compound2compound/smiles2db.json'))
    key2mesh_via_db(smiles2db, db2mesh, smiles2mesh, mesh2smiles)
    print('SMILES-is-MeSH', len(smiles2mesh), len(mesh2smiles))

    # Export
    json.dump(switch_dictset_to_dictlist(smiles2mesh), open(direct+'smiles2mesh.json','w'))
    json.dump(switch_dictset_to_dictlist(mesh2smiles), open(direct+'mesh2smiles.json','w'))
    json.dump(switch_dictset_to_dictlist(inchikey2mesh), open(direct+'inchikey2mesh.json','w'))
    json.dump(switch_dictset_to_dictlist(mesh2inchikey), open(direct+'mesh2inchikey.json','w'))
    
    
def align_x_to_mesh_via_db(x2db, x_name, direct='output/compound2compound/'):

    # Import
    x2db = json.load(open(direct+f'{x_name.lower()}2db.json'))
    db2mesh = json.load(open(os.path.join(direct, 'db2mesh.json')))
    print(len(db2mesh))
    
    # Add items
    x2mesh, mesh2x = dict(), dict()
    key2mesh_via_db(x2db, db2mesh, x2mesh, mesh2x)
    print(f'{x_name}-is-MeSH {len(x2mesh)}')

    # Export
    with open(os.path.join(direct, f'{x_name.lower()}2mesh.json'),'w') as fout:
        json.dump(switch_dictset_to_dictlist(x2mesh), fout)
    with open(os.path.join(direct, f'mesh2{x_name.lower()}.json'),'w') as fout:
        json.dump(switch_dictset_to_dictlist(mesh2x), fout)
    
    db2mesh = switch_dictset_to_dictlist(db2mesh)
    json.dump(db2mesh, open('output/compound2compound/db2mesh.json','w'))
        
    
    
    
def align_kegg_drugs_and_compounds():
    keggcompound2mesh = json.load(open('output/compound2compound/keggcompound2mesh.json'))
    mesh2keggcompound = json.load(open('output/compound2compound/mesh2keggcompound.json'))

    keggdrug2mesh = json.load(open('output/compound2compound/keggdrug2mesh.json'))
    mesh2keggdrug = json.load(open('output/compound2compound/mesh2keggdrug.json'))

    keggdrug2db = json.load(open('output/compound2compound/keggdrug2db.json'))
    keggcompound2db = json.load(open('output/compound2compound/keggcompound2db.json'))

    ''' KEGG Compund -is- KEGG Drug '''
    keggcompound2keggdrug, keggdrug2keggcompound = dict(), dict()
    for mesh, kegg_drugs in mesh2keggdrug.items():
        try:
            kegg_cpds = mesh2keggcompound[mesh]
            for kegg_cpd in kegg_cpds:
                for kegg_drug in kegg_drugs:
                    keggcompound2keggdrug.setdefault(kegg_cpd,set()).add(kegg_drug)
                    keggdrug2keggcompound.setdefault(kegg_drug,set()).add(kegg_cpd)
        except:
            continue

    for mesh, kegg_compounds in mesh2keggcompound.items():
        try:
            kegg_drugs = mesh2keggdrug[mesh]
            for kegg_drug in kegg_drugs:
                for kegg_compound in kegg_compounds:
                    keggdrug2keggcompound.setdefault(kegg_drug,set()).add(kegg_compound)
                    keggcompound2keggdrug.setdefault(kegg_compound,set()).add(kegg_drug)
        except:
            continue

    print(len(keggcompound2keggdrug), len(keggdrug2keggcompound), 'KEGG Compounds aligned to Drug / Compounds aligned to Drugs with MeSH alignments')
    
    
    '''Export'''
    os.system('curl https://rest.kegg.jp/link/cpd/drug > input/KEGG/kegg_drug_to_cpd.tsv')
    for line in open('input/KEGG/kegg_drug_to_cpd.tsv'):
        line = line.split('\t')

        # KEGG Drug, KEGG Compund
        kegg_drug = line[0].split('dr:')[1]
        kegg_cpd = line[1].strip().split('cpd:')[1]

        # KEGG Drug -is- KEGG Compound
        keggcompound2keggdrug.setdefault(kegg_cpd,set()).add(kegg_drug)
        keggdrug2keggcompound.setdefault(kegg_drug,set()).add(kegg_cpd)

    print(len(keggcompound2keggdrug), len(keggdrug2keggcompound), 'KEGG Compounds aligned to Drug / Compounds aligned to Drugs')

    json.dump(switch_dictset_to_dictlist(keggcompound2keggdrug), open('output/compound2compound/keggcompound2keggdrug.json','w'))
    
    
            
def getOld2NewTTDdrug(b_id, old_ttd_list):
    '''
    Function: 
    - Saves JSON files of {old TTD ID : [new TTD ID], ...}
    
    Params: 
    - b_id: batch ID from the multiprocess process
    - old_ttd_list: list of deprecated TTD IDs / old naming scheme
    '''
    old_ttd2newttd = dict() 
    no_newttd = set()
    tot = str(len(old_ttd_list))
    
    # Find drug's updated TTD ID
    for i, old_ttd in enumerate(old_ttd_list):
        old_ttd = old_ttd[0]
        
        # Print progress
        if b_id == 1 or int(b_id) > 999:
            print('Progress of Batch 1'+str(i)+'/'+tot, end='\r')
        
        try:
            res = requests.get('http://db.idrblab.net/web/drug/'+old_ttd)
            new_ttd = res.url.split('/')[-1].upper()
            old_ttd2newttd.setdefault(old_ttd, set()).add(new_ttd)
        except:
            no_newttd.add(old_ttd)
            continue
     
    # Export Old TTD -is- New TTD
    old_ttd2newttd = switch_dictset_to_dictlist(old_ttd2newttd)
    fout1 = open('output/compound2compound/temp_newTTDfound_'+str(b_id)+'.json', 'w')
    json.dump(old_ttd2newttd, fout1)
    
    # Export Old TTD without New TTD alignments
    no_newttd = list(no_newttd)
    fout2 = open('output/compound2compound/temp_no_newTTDfound_'+str(b_id)+'.json','w')
    json.dump(no_newttd, fout2)


def align_drugbank_to_ttd():
    '''DrugBank provides old TTD IDs. Use DrugBank and TTD together.'''
    ''' Get DrugBank to Alt IDs: CAS, Chebi, PubChem Substance, PubChem Compound, oldTTD '''
    db2oldttd = json.load(open('output/compound2compound/db2oldttd.json'))
    db2chebi = json.load(open('output/compound2compound/db2chebi.json'))
    db2pubchemC = json.load(open('output/compound2compound/db2pubchemcomp.json'))
    db2pubchemS = json.load(open('output/compound2compound/db2pubchemsub.json'))
    db2cas = json.load(open('output/compound2compound/db2cas.json'))
    db2ChEMBL = json.load(open('output/compound2compound/db2chembl.json'))

    # Download TTD's External Mapping Drug file: http://db.idrblab.net/ttd/full-data-download
    os.system('wget -P input/ http://db.idrblab.net/ttd/sites/default/files/ttd_database/P1-03-TTD_crossmatching.txt')

    ''' Align New TTD Drug ID -is- xref Alt Drug ID '''
    ttd2altids = dict(); info = []
    with open('input/P1-03-TTD_crossmatching.txt','r') as fin:
        for count, line in enumerate(fin):
            if count < 28:
                continue
            if line == '\t\t\n':
                i = getTTD_IDs(info)
                ttd2altids[list(i.keys())[0]] = list(i.values())[0]
                info = []
                continue
            info.append(line.strip().split('\t'))

    cas2ttd, pubchemc2ttd, pubchemsubs2ttd, chebi2ttd = dict(), dict(), dict(), dict()
    align_ttd_to_xrefs(ttd2altids, cas2ttd, pubchemc2ttd, pubchemsubs2ttd, chebi2ttd)


    ''' Align New TTD -is- DrugBank '''
    ''' (NewTTD-is-AltID + AltID-is-DrugBank) '''
    ttd2db = dict()
    # Map DrugBank-PubChemC-TTD
    key2value2ttd(db2pubchemC, pubchemc2ttd, ttd2db)

    # Map DrugBank-PubChemSubs-TTD
    key2value2ttd(db2pubchemS, pubchemsubs2ttd, ttd2db)

    # Map DrugBank-CAS-TTD
    key2value2ttd(db2cas, cas2ttd, ttd2db)

    # Map DrugBank-Chebi-TTD
    key2value2ttd(db2chebi, chebi2ttd, ttd2db)


    ''' Align OldTTD -is- NewTTD via DrugBank '''
    #NOTE: Fix this part
    multiprocess_a_dict_of_lists_values(db2oldttd, getOld2NewTTDdrug)

    # Merge temp files
    oldttd2newttd = switch_dictset_to_dictlist(merge_oldttd2newttd())
    json.dump(oldttd2newttd, open('output/compound2compound/oldTTD2newTTD.json','w'))
    print(len(oldttd2newttd))
    # ------------------------------------

    ''' Align DrugBank -is- New TTD '''
    ttd2db = dict()
    for db, oldttd in db2oldttd.items():
        try:
            newttd = oldttd2newttd[oldttd[0]][0]
            ttd2db.setdefault(newttd,set()).add(db)
        except:
            continue
    ttd2db = switch_dictset_to_dictlist(ttd2db)


    print('New TTD IDs mapping to multiple DrugBank IDs')
    print_if_len_value_greater_than_1(ttd2db)

    json.dump(ttd2db, open('output/compound2compound/ttd2db.json','w'))
    
    
    
    ''' Align DrugBank -is- New TTD '''
    ttd2db = dict()
    for db, oldttd in db2oldttd.items():
        try:
            newttd = oldttd2newttd[oldttd[0]][0]
            ttd2db.setdefault(newttd,set()).add(db)
        except:
            continue
    ttd2db = switch_dictset_to_dictlist(ttd2db)


    print('New TTD IDs mapping to multiple DrugBank IDs')
    print_if_len_value_greater_than_1(ttd2db)
    
    
def output_final_db_to_mesh_edges():
    db2mesh = json.load(open(os.path.join('output/compound2compound/db2mesh.json')))
    mesh2db = json.load(open(os.path.join('output/compound2compound/mesh2db.json'))) 
    print(f'{len(db2mesh)} DrugBank-to-MeSH {len(mesh2db)} MeSH to DrugBank')

    output_edgefile_onerel_noweight(outpath = 'output/compound2compound/Compound_(DrugBank)_2_Compound_(MeSH).csv',
                                   columns = ['Compound (DrugBank)','Compound (MeSH)','Relationship'],
                                   dictionary = db2mesh,
                                   rel = '-is-',
                                   prefix_col1 = 'DrugBank_Compound:',
                                   prefix_col2 = 'MeSH_Compound:',
                                   )
    
'''
INSTRUCTIONS: 
Download the zipped MRCONSO.RRF file from UMLS (after you log in with your account)
into the input folder https://www.nlm.nih.gov/research/umls/licensedcontent/umlsknowledgesources.html


 File:  https://go.drugbank.com/releases/5-1-9/downloads/all-full-database
- Create a DrugBank account and request an academic license. Then, log in and download the above file in the latest version

'''





if __name__ == '__main__':
    # Align UMLS, DrugBank, MeSH (Source: MRCONSO, DrugBank, MeSH)
    create_empty_folders_for_data()
    extract_mr_conso_file()


    extract_mesh_umls_drugbank_from_mr_conso()
    extract_drugbank_xml()
    root = parse_drugbank_xml()
    dbs = map_drugbank_to_external_ids(root)  
    del root




    ''' Align identifiers to DrugBank via ___'''
    # via ATC
    align_db_to_mesh_via_atc()

    # via UNII
    unii2mesh, mesh2unii = map_unii_to_mesh_via_mesh_api()
    db2unii = json.load(open('output/compound2compound/db2unii.json'))
    map_unii_to_mesh_via_mycheminfo(db2unii, unii2mesh, mesh2unii)
    print(len(json.load(open('output/compound2compound/unii2mesh.json'))), 'UNII to MeSH')

    db2unii = json.load(open('output/compound2compound/db2unii.json'))
    unii2mesh = json.load(open('output/compound2compound/unii2mesh.json'))
    db2mesh = json.load(open(os.path.join(direct, 'db2mesh.json')))
    mesh2db = json.load(open(os.path.join(direct, 'mesh2db.json')))
    align_db_to_mesh_via_x_to_mesh_list(db2unii, unii2mesh, db2mesh, mesh2db)

    # via CAS
    map_drugbank_to_mesh_drugbank_to_cas_via_inxight_drugs(dbs)
    db2cas = json.load(open('output/compound2compound/db2cas.json'))
    cas2mesh, mesh2cas = {}, {}
    map_cas_to_mesh_via_mesh_api(db2cas, cas2mesh, mesh2cas)
    mesh2cas = json.load(open('output/compound2compound/mesh2cas.json'))
    cas2mesh = json.load(open('output/compound2compound/cas2mesh.json'))

    align_db_to_mesh_via_x_to_mesh_list(db2cas, cas2mesh, db2mesh, mesh2db)

    # via UMLS
    db2umls = json.load(open('output/compound2compound/db2umls.json'))
    umls2db = json.load(open('output/compound2compound/umls2db.json'))


    del cas2mesh
    del mesh2cas
    del unii2mesh
    del mesh2unii
    del db2cas

    umls2mesh = json.load(open('output/otherMappings/umls2mesh.json'))
    align_db_to_mesh_via_x_to_mesh_list(db2umls, umls2mesh, db2mesh, mesh2db)

    # Inchi & SMILES -to- DrugBank & MeSH
    os.system('unzip input/drugbank_all_structure_links.csv.zip')
    os.system("mv 'structure links.csv' 'input/structure links.csv'")
    align_inchi_and_smiles_to_drugbank()
    align_inchi_and_smiles_to_mesh()

    # via SMILES
    print(len(db2mesh), len(mesh2db))
    db2smiles = json.load(open('output/compound2compound/db2smiles.json'))
    db2mesh = json.load(open(os.path.join(direct, 'db2mesh.json')))
    mesh2db = json.load(open(os.path.join(direct, 'mesh2db.json')))
    smiles2mesh = json.load(open('output/compound2compound/smiles2mesh.json'))
    align_db_to_mesh_via_x_to_mesh_list(db2smiles, smiles2mesh, db2mesh, mesh2db)
    print(len(db2mesh), len(mesh2db))
    del smiles2mesh
    del db2smiles

    # via InChI
    db2inchi = json.load(open('output/compound2compound/db2inchikey.json'))
    inchi2mesh = json.load(open(direct+'inchikey2mesh.json'))
    align_db_to_mesh_via_x_to_mesh_list(db2inchi, inchi2mesh, db2mesh, mesh2db)
    print('MeSH-is-DrugBank:',len(mesh2db))
    print('Drugbank-is-MeSH:',len(db2mesh))


    ''' Align identifiers to MeSH via ___ to DrugBank '''
    # via ChEBI
    chebi2db = json.load(open(os.path.join(direct, 'chebi2db.json')))
    align_x_to_mesh_via_db(chebi2db, 'ChEBI')
    del chebi2db

    # via CAS
    cas2db = json.load(open(os.path.join(direct, 'cas2db.json')))
    align_x_to_mesh_via_db(cas2db, 'CAS')
    del cas2db

    # via KEGG Compound
    keggcompound2db = json.load(open(os.path.join(direct, 'keggcompound2db.json')))
    align_x_to_mesh_via_db(keggcompound2db, 'KEGGcompound')
    del keggcompound2db

    # via KEGG Drug
    keggdrug2db = json.load(open(os.path.join(direct, 'keggdrug2db.json')))
    align_x_to_mesh_via_db(keggdrug2db, 'KEGGdrug')
    del keggdrug2db

    # via SMILES
    print(len(db2mesh), len(mesh2db))
    db2smiles = json.load(open('output/compound2compound/db2smiles.json'))
    db2mesh = json.load(open(os.path.join(direct, 'db2mesh.json')))
    mesh2db = json.load(open(os.path.join(direct, 'mesh2db.json')))
    smiles2mesh = json.load(open('output/compound2compound/smiles2mesh.json'))
    align_db_to_mesh_via_x_to_mesh_list(db2smiles, smiles2mesh, db2mesh, mesh2db)
    print(len(db2mesh), len(mesh2db))
    del smiles2mesh
    del db2smiles

    # via InChI
    db2inchi = json.load(open('output/compound2compound/db2inchikey.json'))
    inchi2mesh = json.load(open(direct+'inchikey2mesh.json'))
    align_db_to_mesh_via_x_to_mesh_list(db2inchi, inchi2mesh, db2mesh, mesh2db)
    print('MeSH-is-DrugBank:',len(mesh2db))
    print('Drugbank-is-MeSH:',len(db2mesh))

    ''' Align identifiers to MeSH via ___ to DrugBank '''
    # via ChEBI
    chebi2db = json.load(open(os.path.join(direct, 'chebi2db.json')))
    align_x_to_mesh_via_db(chebi2db, 'ChEBI')
    del chebi2db

    # via CAS
    cas2db = json.load(open(os.path.join(direct, 'cas2db.json')))
    align_x_to_mesh_via_db(cas2db, 'CAS')
    del cas2db

    # via KEGG Compound
    keggcompound2db = json.load(open(os.path.join(direct, 'keggcompound2db.json')))
    align_x_to_mesh_via_db(keggcompound2db, 'KEGGcompound')
    del keggcompound2db

    # via KEGG Drug
    keggdrug2db = json.load(open(os.path.join(direct, 'keggdrug2db.json')))
    align_x_to_mesh_via_db(keggdrug2db, 'KEGGdrug')
    del keggdrug2db

    # KEGG Compound to KEGG Drug
    align_kegg_drugs_and_compounds()

    output_final_db_to_mesh_edges()

    align_drugbank_to_ttd()
    print('Finished')
