from biomedkg_utils import multiprocess_a_list, output_edgefile_onerel_noweight, switch_dictset_to_dictlist
import os
import pandas as pd
import csv
import json
import requests
from bs4 import BeautifulSoup
from multiprocessing import cpu_count
import xml.etree.ElementTree as ET
from parse_xml import *

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




def download_reactome_compound_to_pathway():
    os.system('wget -N -P input/ https://reactome.org/download/current/ChEBI2Reactome_PE_Pathway.txt')
    chebi2reactome = pd.read_table('input/ChEBI2Reactome_PE_Pathway.txt', header=None)
    chebi2reactome.columns=['ChEBI Compound ID','Reactome Compound ID','Reactome Compound Name','Reactome Pathway', 'Link', 'Pathway Name','Evidence','Species']
    chebi2reactome = chebi2reactome[chebi2reactome['Species']=='Homo sapiens']
    chebi2reactome.to_csv('output/compound2pathway/chebi2reactome_df.csv',index=False)
    print(len(chebi2reactome), 'rows')


def reactome_compound_type_scraper_api(b_id, reactome_comp_batch):
    '''
    FUNCTION:
    - Access compounds' Reactome webpage. Use a scraper to 
      find the type of compound (e.g., ChemicalDrug, Chemical Compound
      Polymer, Genes and Transcripts, OtherEntity, etc.)
    
    PARAMS:
    - b_id: batch ID (integer number ranging from 0 to number of processors)
    - reactome_comp_batch: a batch of Reactome compound IDs
    '''
    
    temp_path = 'output/compound2compound/temp_reactomecomp'+str(b_id)+'.txt'
    open(temp_path,'w').write('')
    
    tot = str(len(reactome_comp_batch)-1)

    for i, reactome_comp in enumerate(reactome_comp_batch):
        
        # Print progress
        if b_id == 1:
            print('Progress of Batch 1:'+str(i)+'/'+tot, end='\r')
        
        # API
        req = requests.get('https://reactome.org/content/detail/'+reactome_comp+'#Homo%20sapiens.html')
        soup = BeautifulSoup(req.content, 'html.parser')
        comp_type = str(soup.find_all('div', class_ ="details-field favth-col-lg-10 favth-col-md-9 favth-col-sm-9 favth-col-xs-12")[1].find_all('span')[0]).split('\n')[0].split('>')[-1]
        
        # Export
        with open(temp_path,'a') as fout:
            fout.write(reactome_comp+'|'+comp_type+'\n')
            
            
def merge_reactome_compound_type_results():
    '''Merge the files from the scraper above'''
    type2reactomecomp, reactomecomp2type = dict(),dict()

    procs = cpu_count()
    print(procs, 'processors')
    for b_id in range(procs):
        temp_path = 'output/compound2compound/temp_reactomecomp'+str(b_id)+'.txt'

        try:
            # Read in temp file
            with open(temp_path) as fin:
                for line in fin:

                    # Reactome Compound, Compound Type
                    reactome_comp = line.split('|')[0]
                    compound_type = line.split('|')[1].strip()

                    # Reactome Compound -is-> Compound Type
                    type2reactomecomp.setdefault(compound_type, set()).add(reactome_comp)
                    reactomecomp2type.setdefault(reactome_comp, set()).add(compound_type)
            os.remove(temp_path)
        except:
            print('No path called', temp_path)
        
    json.dump(switch_dictset_to_dictlist(type2reactomecomp), open('output/compound2compound/compoundtype2reactomecompound.json','w'))
    json.dump(switch_dictset_to_dictlist(reactomecomp2type), open('output/compound2compound/reactomecompound2compoundtype.json','w'))
    #type2reactomecomp = json.load(open('output/compound2compound/compoundtype2reactomecompound.json'))

    
def get_reactome_compound_types():
    chebi2reactome = pd.read_csv('output/compound2pathway/chebi2reactome_df.csv')
    # Reactome Compound IDs
    reactome_comps = list(set(chebi2reactome['Reactome Compound ID']))
    # Use scraper to find Reactome Compounds' types 
    multiprocess_a_list(reactome_comps, reactome_compound_type_scraper_api)
    
    
def display_summary_of_compound_to_reactome_pathway():
    chebi2reactome = pd.read_csv('output/compound2pathway/chebi2reactome_df.csv')
    chebis = [str(chebi) for chebi in list(set(chebi2reactome['ChEBI Compound ID']))]
    chebi2mesh = json.load(open('output/compound2compound/chebi2mesh.json'))
    chebi2db = json.load(open('output/compound2compound/chebi2db.json'))

    yes,no = 0,0
    for chebi in chebis:
        if chebi in chebi2mesh: yes += 1
        else: no +=1
    print(str(yes)+'/'+str(yes+no)+' ChEBI -is- DrugBank (in Reactome)')

    yes,no = 0,0
    for chebi in chebis:
        if chebi in chebi2db:  yes += 1
        else:  no +=1

    print(str(yes)+'/'+str(yes+no)+' ChEBI -is- MeSH (in Reactome)')


    yes,no = 0,0
    for chebi in chebis:
        if chebi in chebi2db and chebi in chebi2mesh:  yes += 1
        else:  no +=1

    print(str(yes)+'/'+str(yes+no)+' ChEBI -is- MeSH&DrugBank (in Reactome)')
    
    
def export_compound_to_reactome_pathway():
    chebi2mesh = json.load(open('output/compound2compound/chebi2mesh.json'))
    chebi2db = json.load(open('output/compound2compound/chebi2db.json'))
    chebi2reactome = pd.read_csv('output/compound2pathway/chebi2reactome_df.csv')
    reactomecomp2type = json.load(open('output/compound2compound/reactomecompound2compoundtype.json'))
    
    ### MeSH Compound
    with open('output/compound2pathway/edges_meshCompound-participates_in->reactomePathway.csv','w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Compound (MeSH)','Pathway (Reactome)','Relationship'])

        for i in range(0,len(chebi2reactome)):

            # Compound, Compound Type, Reactome Pathway
            chebi_comp = str(chebi2reactome['ChEBI Compound ID'].iloc[i])
            try:
                mesh_comps = chebi2mesh[chebi_comp]
            except:
                continue
            react_comp = str(chebi2reactome['Reactome Compound ID'].iloc[i])
            comp_type  = list(reactomecomp2type[react_comp])[0]
            react_pway = chebi2reactome['Reactome Pathway'].iloc[i]

            # Write edge for drugs
            if 'drug' in comp_type.lower():
                for mesh_comp in mesh_comps:
                    writer.writerow(['MeSH_Compound:'+mesh_comp,'Reactome_Pathway:'+react_pway,'-drug_participates_in_pathway->'])

    mesh2reactomepw_df = pd.read_csv('output/compound2pathway/edges_meshCompound-participates_in->reactomePathway.csv')
    mesh2reactomepw_df.to_csv('output/edges/edges_meshCompound-participates_in->reactomePathway.csv', index=False)
  #mesh2reactomepw_df.to_csv('output/edges_to_use/Compound_(MeSH)_to_Pathway_(Reactome).csv', index=False)
    print(len(set(mesh2reactomepw_df['Compound (MeSH)'])), 'MeSH Compounds')
    print(len(set(mesh2reactomepw_df['Pathway (Reactome)'])), 'Reactome Pathways')
    print(len(mesh2reactomepw_df), 'Edges')

    
    ### DrugBank Compound
    with open('output/compound2pathway/edges_drugbankCompound-participates_in->reactomePathway.csv','w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Compound (DrugBank)','Pathway (Reactome)','Relationship'])

        for i in range(0,len(chebi2reactome)):
            chebi_comp = str(chebi2reactome['ChEBI Compound ID'].iloc[i])

            try:
                # Compound, Compound Type, Reactome Pathway
                db_comps = chebi2db[chebi_comp]
                comp_type  = list(reactomecomp2type[react_comp])[0]
                react_pway = chebi2reactome['Reactome Pathway'].iloc[i]

                # Write edge for drugs
                if 'Drug' in comp_type:
                    for db_comp in db_comps:
                        writer.writerow(['DrugBank_Compound:'+db_comp,'Reactome_Pathway:'+react_pway,'-drug_participates_in_pathway->'])
            except:
                continue

    drugbank2reactomepw_df = pd.read_csv('output/compound2pathway/edges_drugbankCompound-participates_in->reactomePathway.csv')    
    drugbank2reactomepw_df.to_csv('output/edges/edges_drugbankCompound-participates_in->reactomePathway.csv', index=False)
    drugbank2reactomepw_df.to_csv('output/edges_to_use/Compound_(DrugBank)_to_Pathway_(Reactome).csv', index=False)
    print(len(set(drugbank2reactomepw_df['Compound (DrugBank)'])), 'DrugBank Compounds')
    print(len(set(drugbank2reactomepw_df['Pathway (Reactome)'])), 'Reactome Pathways')
    print(len(drugbank2reactomepw_df), 'Edges')
    
    
def map_compounds_to_smpdb_pathway(root):    
    # DrugBank Compound
    with open('output/compound2pathway/edges_drugbankCompound-participates_in-smpdbPathway.csv','w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Compound (DrugBank)','Pathway (SMPDB)','Relationship'])
        smpdb_pws = dict()

        for count, ele in enumerate(root):
            ID = ParseXML.getID(ele)
            pathways = ParseXML.getPathways(ele)

            if len(pathways) > 0:
                for pw in pathways:
                    smpdb_pws.setdefault(ID,set()).add(pw['smpdb_id'])
                    writer.writerow(['DrugBank_Compound:'+ID,'SMPDB_Pathway:'+pw['smpdb_id'],'-drug_participates_in_pathway->'])
                smpdb_pws[ID] = list(smpdb_pws[ID])
    #smpdb_pws
    df = pd.read_csv('output/compound2pathway/edges_drugbankCompound-participates_in-smpdbPathway.csv').drop_duplicates()
    df.to_csv('output/edges/edges_drugbankCompound-participates_in-smpdbPathway.csv', index=False)
    df.to_csv('output/edges_to_use/Compound_(DrugBank)_to_Pathway_(SMPDB).csv', index=False)

    
    # MeSH Compound
    db2mesh = json.load(open('output/compound2compound/db2mesh.json'))
    with open('output/compound2pathway/edges_MeSHCompound-participates_in-smpdbPathway.csv','w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Compound (MeSH)','Pathway (SMPDB)','Relationship'])
        smpdb_pws_mesh = dict()

        for count, ele in enumerate(root):
            ID = ParseXML.getID(ele)
            pathways = ParseXML.getPathways(ele)

            try:  mesh_compounds = db2mesh[ID]
            except: continue

            if len(pathways) > 0:
                for pw in pathways:
                    for mesh_compound in mesh_compounds:
                        smpdb_pws_mesh.setdefault(mesh_compound, set()).add(pw['smpdb_id'])
                        writer.writerow(['MeSH_Compound:'+mesh_compound, 'SMPDB_Pathway:'+pw['smpdb_id'],'-drug_participates_in_pathway->'])
                #smpdb_pws_mesh[mesh_compound] = list(smpdb_pws_mesh[mesh_compound])
    #smpdb_pws_mesh
    df = pd.read_csv('output/compound2pathway/edges_MeSHCompound-participates_in-smpdbPathway.csv').drop_duplicates()
    df.to_csv('output/edges/edges_MeSHCompound-participates_in-smpdbPathway.csv', index=False)
    #df.to_csv('output/edges_to_use/Compound_(MeSH)_to_Pathway_(SMPDB).csv', index=False)
    json.dump(smpdb_pws,open("output/compound2pathway/compound_to_pathway.json","w"))
    
    
def map_compound_to_kegg_pathway():
    keggcompound2mesh = json.load(open('output/compound2compound/keggcompound2mesh.json'))
    keggcompound2keggdrug = json.load(open('output/compound2compound/keggcompound2keggdrug.json'))
    keggcompound2db = json.load(open('output/compound2compound/keggcompound2db.json'))
    keggdrug2mesh = json.load(open('output/compound2compound/keggdrug2mesh.json'))
    os.system('curl https://rest.kegg.jp/link/pathway/cpd > input/KEGG/kegg_cpd_to_pathway.tsv')

    '''Compound-Pathway (Non-drug)'''
    keggpathway2meshcompound = dict()
    meshcompound2keggpathway = dict()

    for line in open('input/KEGG/kegg_cpd_to_pathway.tsv'):
        line = line.split('\t')

        try:
            # MeSH Compound, KEGG Pathway
            kegg_cpd = line[0].split('cpd:')[1]
            mesh_cpds = keggcompound2mesh[kegg_cpd]
            kegg_pathway = line[1].strip().replace('path:','path_')
            if '_map' in kegg_pathway:
                kegg_pathway = kegg_pathway.replace('_map','_hsa')

            # MeSH Compound -> KEGG Pathway
            if kegg_cpd not in keggcompound2keggdrug:
                for mesh_cpd in mesh_cpds:
                    meshcompound2keggpathway.setdefault(mesh_cpd, set()).add(kegg_pathway)
                    keggpathway2meshcompound.setdefault(kegg_pathway,set()).add(mesh_cpd)
        except:
            continue

    print(len(keggpathway2meshcompound), 'KEGG Pathways')
    print(len(meshcompound2keggpathway), 'MeSH Compounds')

    file = 'Compound_(MeSH)_compound_participates_in_Pathway_(KEGG).csv'
    outpath = os.path.join('output/compound2pathway',file)
    output_edgefile_onerel_noweight(
        outpath = outpath,
        columns = ['Compound (MeSH)','Pathway (KEGG)','Relationship'],
        dictionary = meshcompound2keggpathway,
        rel = '-compound_participates_in->',
        prefix_col1='MeSH_Compound:',
        prefix_col2='KEGG_Pathway:',
        edges_to_use_folder=False,
    )
    df = pd.read_csv(outpath)
    
    
    ''' DrugBank Compound-Pathway (Non-drug)'''
    keggpathway2dbcompound = dict()
    dbcompound2keggpathway = dict()

    for line in open('input/KEGG/kegg_cpd_to_pathway.tsv'):
        line = line.split('\t')

        try:
            # db Compound, KEGG Pathway
            kegg_cpd = line[0].split('cpd:')[1]
            db_cpds = keggcompound2db[kegg_cpd]
            kegg_pathway = line[1].strip().replace('path:','path_')
            if '_map' in kegg_pathway:
                kegg_pathway = kegg_pathway.replace('_map','_hsa')
            
            # db Compound -> KEGG Pathway
            if kegg_cpd not in keggcompound2keggdrug:
                for db_cpd in db_cpds:
                    dbcompound2keggpathway.setdefault(db_cpd, set()).add(kegg_pathway)
                    keggpathway2dbcompound.setdefault(kegg_pathway,set()).add(db_cpd)
        except:
            continue

    print(len(keggpathway2dbcompound), 'KEGG Pathways')
    print(len(dbcompound2keggpathway), 'DrugBank Compounds')


    file = 'Compound_(DrugBank)_compound_participates_in_Pathway_(KEGG).csv'
    outpath = os.path.join('output/compound2pathway/',file)
    output_edgefile_onerel_noweight(
        outpath = outpath,
        columns = ['Compound (DrugBank)','Pathway (KEGG)','Relationship'],
        dictionary = dbcompound2keggpathway,
        rel = '-compound_participates_in->',
        prefix_col1='DrugBank_Compound:',
        prefix_col2='KEGG_Pathway:',
    )
    
    
def map_drug_to_kegg_pathway():
    os.system('curl https://rest.kegg.jp/link/pathway/drug > input/KEGG/kegg_drug_to_pathway.tsv')
    keggdrug2mesh = json.load(open('output/compound2compound/keggdrug2mesh.json'))

    '''MeSH Compound-Pathway (Drug)'''
    keggpathway2meshdrug = dict()
    meshdrug2keggpathway = dict()

    for line in open('input/KEGG/kegg_drug_to_pathway.tsv'):
        line = line.split('\t')

        try:
            # MeSH drug, KEGG Pathway
            kegg_drug = line[0].split('dr:')[1]
            mesh_drugs = keggdrug2mesh[kegg_drug]
            kegg_pathway = line[1].strip().replace('path:','path_')
            if '_map' in kegg_pathway:
                kegg_pathway = kegg_pathway.replace('_map','_hsa') ### NOTE: Temp fix

            # MeSH drug -> KEGG Pathway
            for mesh_drug in mesh_drugs:
                meshdrug2keggpathway.setdefault(mesh_drug, set()).add(kegg_pathway)
                keggpathway2meshdrug.setdefault(kegg_drug, set()).add(mesh_drug)
        except:
            continue

    print(len(keggpathway2meshdrug), 'KEGG Pathways')
    print(len(meshdrug2keggpathway), 'MeSH Drugs')

    file = 'Compound_(MeSH)_drug_participates_in_Pathway_(KEGG).csv'
    outpath = os.path.join('output/compound2pathway',file)
    output_edgefile_onerel_noweight(
        outpath = outpath,
        columns = ['Compound (MeSH)','Pathway (KEGG)','Relationship'],
        dictionary = meshdrug2keggpathway,
        rel = '-drug_participates_in->',
        prefix_col1='MeSH_Compound:',
        prefix_col2='KEGG_Pathway:',
        edges_to_use_folder=False,
    )
    df = pd.read_csv(outpath)
    
    
    # DrugBank
    keggdrug2db = json.load(open('output/compound2compound/keggdrug2db.json'))

    '''DrugBank Compound-Pathway (Drug)'''
    keggpathway2dbdrug = dict()
    dbdrug2keggpathway = dict()

    for line in open('input/KEGG/kegg_drug_to_pathway.tsv'):
        line = line.split('\t')

        try:
            # db drug, KEGG Pathway
            kegg_drug = line[0].split('dr:')[1]
            db_drugs = keggdrug2db[kegg_drug]
            kegg_pathway = line[1].strip().replace('path:','path_')
            if '_map' in kegg_pathway:
                kegg_pathway = kegg_pathway.replace('_map','_hsa') ### NOTE: Temp fix

        # db drug -> KEGG Pathway
            for db_drug in db_drugs:
                dbdrug2keggpathway.setdefault(db_drug, set()).add(kegg_pathway)
                keggpathway2dbdrug.setdefault(kegg_drug, set()).add(db_drug)
        except:
            continue

    print(len(keggpathway2dbdrug), 'KEGG Pathways')
    print(len(dbdrug2keggpathway), 'DrugBank Drugs')

    file = 'Compound_(DrugBank)_drug_participates_in_Pathway_(KEGG).csv'
    outpath = os.path.join('output/compound2pathway', file)
    output_edgefile_onerel_noweight(
        outpath = outpath,
        columns = ['Compound (DrugBank)','Pathway (KEGG)','Relationship'],
        dictionary = dbdrug2keggpathway,
        rel = '-drug_participates_in->',
        prefix_col1='DrugBank_Compound:',
        prefix_col2='KEGG_Pathway:'
    )
    df = pd.read_csv(outpath)
    df.to_csv(os.path.join('output/edges',file),index=False)
    df.to_csv(os.path.join('output/edges_to_use/',file),index=False)
    
    
if __name__ == '__main__':
    try:
        root = parse_drugbank_xml()
    except:
        extract_drugbank_xml()
        root = parse_drugbank_xml() 
    download_reactome_compound_to_pathway()
    get_reactome_compound_types()
    merge_reactome_compound_type_results()
    export_compound_to_reactome_pathway()    
    display_summary_of_compound_to_reactome_pathway()
    map_compounds_to_smpdb_pathway(root)
    map_compound_to_kegg_pathway()
    map_drug_to_kegg_pathway()