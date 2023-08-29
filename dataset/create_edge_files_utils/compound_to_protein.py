import xml.etree.ElementTree as ET
import pandas as pd
from parse_xml import *
import csv
import json
import os
from biomedkg_utils import switch_dictset_to_dictlist
from biomed_apis import *

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


def map_compound_to_protein_df(root):
    compoundprotein_df = pd.DataFrame(columns=['ID','Name','Synonyms','ATC','Pathways','Indications','Drug-->Protein'])
    for i, ele in enumerate(root):
        print(i, end = '\r')

        name = ParseXML.getName(ele)                               # Drug name
        syn = [drug.lower() for drug in ParseXML.getSynonyms(ele)] # Drug synonyms
        ID = ParseXML.getID(ele)    # DrugBank ID
        indications = ParseXML.getIndication(ele)                  # Description of indications
        pathways = ParseXML.getPathways(ele)                       # SMPDB Pathway information (drug targets pathways)
        atc = ParseXML.getATCCode(ele)                             # ATC, drug classification code

        # Drug-[relationship]->Protein dictionary
        d2p = []
        for relation in ['carriers', 'targets', 'transporters' , 'enzymes']:
            d2p.append({relation:ParseXML.getEntities(ele, relation)})

        compoundprotein_df.loc[len(compoundprotein_df.index)] = [ID, name, syn, atc, pathways, indications,d2p]
    
    drugbankIDs = list(compoundprotein_df['ID'])
    compoundprotein_df['ID'].to_csv('output/compound2compound/drugbankIDs')

    return compoundprotein_df


def map_compounds_to_targets(root):
    # Option 2: Use Dictionaries for each relationship: {Drug: [upids], ...}
    count = 0
    drug_targets_prot, drug_carry_prot, drug_transport_prot, drug_enzyme_prot = dict(), dict(), dict(), dict()
    drugTargetsProtrelationships = 0
    with open('output/compound2protein/edges_compoundprotein_drugbank.csv', 'w') as fout,\
        open('output/compound2protein/edges_compoundprotein_mesh.csv','w') as fout1:
        writer = csv.writer(fout)
        writer.writerow(['Compound (DrugBank)', 'Protein (UniProt)','Relationship'])

        writer2 = csv.writer(fout1)
        writer2.writerow(['Compound (MeSH)', 'Protein (UniProt)', 'Relationship'])

        for count, ele in enumerate(root):
            ID = ParseXML.getID(ele)
            for relation in ['carriers', 'targets', 'transporters' , 'enzymes']:
                prots = ParseXML.getEntities(ele, relation)

                for prot in prots:
                    upid = prot['uniprot_id']
                    if(upid != 'Null'):

                        if(relation == 'carriers'):
                            drug_carry_prot.setdefault(ID,[]).append(upid)

                        elif(relation == 'targets'):
                            drug_targets_prot.setdefault(ID,[]).append(upid)
                            drugTargetsProtrelationships += 1

                        elif(relation == 'transporters'):
                            drug_transport_prot.setdefault(ID,[]).append(upid)

                        elif(relation == 'enzymes'):
                            drug_enzyme_prot.setdefault(ID,[]).append(upid)

                for prot in prots:
                    upid = prot['uniprot_id']
                    if(upid != 'Null'):
                        if relation == 'targets':
                            writer.writerow(['DrugBank_Compound:'+ID, 'UniProt:'+upid, '-drug_%s_protein->'%relation])
                            try:
                                mesh_compounds = db2mesh[ID]
                                for mesh in mesh_compounds:
                                    writer2.writerow(['MeSH_Compound:'+mesh, 'UniProt:'+upid, '-drug_%s_protein->'%relation])
                            except:
                                continue

                        else:
                            writer.writerow(['DrugBank_Compound:'+ID, 'UniProt:'+upid, '-drug_uses_protein_as_%s-'%relation])
                            try:
                                mesh_compounds = db2mesh[ID]
                                for mesh in mesh_compounds:
                                    writer2.writerow(['MeSH_Compound:'+mesh, 'UniProt:'+upid, '-drug_uses_protein_as_%s-'%relation])
                            except:
                                continue
    print('drugTargetsProtrelationships', drugTargetsProtrelationships)


    drugwithtarget = 0
    num_actions, num_prots, num_drugs = 0, 0, 0 
    uniq_acts, uniq_prot, unique_drug = dict(), set(), set()
    with open('output/compound2protein/edges_compoundsometargetsdetailed_drugbank.csv', 'w') as fout,\
        open('output/compound2protein/edges_compoundsometargetsdetailed_mesh.csv', 'w') as fout1:

        writer = csv.writer(fout)
        writer.writerow(['Compound (DrugBank)','Protein (UniProt)', 'Relationship'])

        writer2 = csv.writer(fout1)
        writer2.writerow(['Compound (MeSH)', 'Protein (UniProt)', 'Relationship'])

        for ele in root:
            ID = ParseXML.getID(ele)
            prots = ParseXML.getEntities(ele, 'targets')

            if len(prots) > 0:
                drugwithtarget += 1
                has_action = False
                for prot in prots:
                    upid = prot['uniprot_id']
                    actions = prot['actions']
                    uniq_prot.add(upid)
                    num_prots += 1

                    if(upid != 'Null'):
                        if len(actions) > 0:
                            has_action = True
                            for action in actions:

                                # DrugBank Compound
                                writer.writerow(['DrugBank_Compound:'+ID, 'UniProt:'+upid, '-%s->'%action.replace(' ','_')])
                                uniq_acts.setdefault(action,0)
                                num_actions += 1
                                uniq_acts[action] += 1

                                # MeSH Compound
                                try:
                                    mesh_compounds = db2mesh[ID]
                                    for mesh in mesh_compounds:
                                        writer2.writerow(['MeSH_Compound:'+mesh, 'UniProt:'+upid, '-%s->'%action.replace(' ','_')])
                                except:
                                    continue
                        else:
                            # DrugBank Compound
                            writer.writerow(['DrugBank_Compound:'+ID, 'UniProt:'+upid, '-drug_targets_protein->'])

                            # MeSH Compound
                            try:
                                mesh_compounds = db2mesh[ID]
                                for mesh in mesh_compounds:
                                    writer2.writerow(['MeSH_Compound:'+mesh, 'UniProt:'+upid, '-drug_targets_protein->'])
                            except:
                                continue


                if has_action == True:
                    num_drugs += 1
                    unique_drug.add(ID)

    os.system('cp output/compound2protein/edges_compoundsometargetsdetailed_drugbank.csv output/edges/edges_compoundsometargetsdetailed_drugbank.csv')
    os.system('cp output/compound2protein/edges_compoundsometargetsdetailed_mesh.csv output/edges/edges_compoundsometargetsdetailed_mesh.csv')
    os.system('cp output/compound2protein/edges_compoundprotein_drugbank.csv output/edges/edges_compoundprotein_drugbank.csv')
    os.system('cp output/compound2protein/edges_compoundprotein_mesh.csv output/edges/edges_compoundprotein_mesh.csv')

    print('Drug with Target:', drugwithtarget)
    print('Total Drug-action-Target Relationships:',num_actions, num_prots, num_drugs)
    print('Unique Actions:',len(uniq_acts), 'Unique Proteins:',len(uniq_prot), 'Unique Drugs:',len(unique_drug))
    print('Actions:', uniq_acts)

    
    
    d_proteins = set()
    for count, ele in enumerate(root):
        for relation in ['carriers', 'targets', 'transporters' , 'enzymes']:

            prots = ParseXML.getEntities(ele, relation)
            for prot in prots:
                upid = prot['uniprot_id']
                if(upid != 'Null'):
                    d_proteins.add(upid) 
                    
    # Export Drug-Protein
    json.dump(drug_targets_prot, open("output/compound2protein/compound_targets_prot.json","w"))
    json.dump(drug_carry_prot, open("output/compound2protein/compound_carry_prot.json","w"))
    json.dump(drug_transport_prot,open("output/compound2protein/compound_transport_prot.json","w"))
    json.dump(drug_enzyme_prot,open("output/compound2protein/compound_enzyme_prot.json","w"))
    json.dump(list(d_proteins),open("output/compound2protein/drugs_proteins.json","w"))    
    

def map_ttd_drug_to_ttd_target():
    '''Download TTD drug-to-protein data'''
    os.system('wget -N -P input/ http://db.idrblab.net/ttd/sites/default/files/ttd_database/P1-07-Drug-TargetMapping.xlsx')
    target2compound_df = pd.read_excel('input/P1-07-Drug-TargetMapping.xlsx')
    ttd2db = json.load(open('output/compound2compound/ttd2db.json'))

    ''' TTD Drug-targets-> TTD Protein '''
    ttd2ttdtarget, ttd2ttdtarget_allstatuses = dict(), dict()
    all_ttdtargets = set()
    for i in range(0,len(target2compound_df)):

        # TTD Drug
        drug = target2compound_df['DrugID'][i]
        if drug in ttd2db:

            # Protein target
            ttd_target = target2compound_df['TargetID'][i]

            # Approval status (e.g., Approved, Patented)
            status = target2compound_df['Highest_status'][i]

            # Save drug->target mapping
            ttd2ttdtarget_allstatuses.setdefault(drug,[]).append(ttd_target)
            all_ttdtargets.add(ttd_target)
            if status == 'Patented' or status == 'Approved' or status == 'Phase 3':
                ttd2ttdtarget.setdefault(drug,[]).append(ttd_target)

    return ttd2ttdtarget
    #len(all_ttdtargets)
    
    
# Read through TTD-provided file
def get_ttd_to_uniprot_name(info):
    '''
    FUNCTION:
    - Align TTD Target ID -is- UniProt Name
    
    PARAMS:
    - info: line in TTD input file
    '''
    up_name, targetID = '',''
    for line in info:
        if 'TARGETID' in line:
            targetID = line[2]
        if 'UNIPROID' in line:
            up_name = line[2]

    if up_name == '' or targetID == '':
        return '',''
    else:
        return targetID, up_name

    
def align_ttd_target_id_to_uniprot_target_name():
    os.system('wget -N -P input/ http://db.idrblab.net/ttd/sites/default/files/ttd_database/P1-01-TTD_target_download.txt')

    ''' Align UniProt name to TTD Target ID '''
    info = []
    up_name2ttdtarg, ttdtarg2up_name = dict(), dict()
    with open('input/P1-01-TTD_target_download.txt','r') as fin:
        for i, line in enumerate(fin):
            # Start when table starts
            if i < 40:
                continue

            ### TTD ID -is- UniProt Name
            if ((line == '\n' or line.startswith('\t') or line == '') and info != []):
                ttdid, uniprot = get_ttd_to_uniprot_name(info)
                uniprots = uniprot.split(';')
                for uniprot in uniprots:
                    #if '_HUMAN' in uniprot:
                    up_name2ttdtarg.setdefault(uniprot, set()).add(ttdid)
                    ttdtarg2up_name.setdefault(ttdid,set()).add(uniprot)

                info = []
                continue
            info.append(line.strip().split('\t'))

    # Export (checkpoint)
    json.dump(switch_dictset_to_dictlist(up_name2ttdtarg), open('output/protein2protein/up_name2ttdtarg.json','w'))
    json.dump(switch_dictset_to_dictlist(ttdtarg2up_name), open('output/protein2protein/ttdtarg2up_name.json','w'))
    
    
def align_uniprot_name_to_uniprot_id():
    up_name2ttdtarg = json.load(open('output/protein2protein/up_name2ttdtarg.json'))
    up_name2ttdtarg.pop('')
    protein_names = [name.split('_HUMAN')[0].strip().split(' ')[0] for name in list(up_name2ttdtarg.keys()) if '_HUMAN' in name]
    with open('submit_this_to_uniprot.txt','w') as fout:
        for protein_name in protein_names:
            fout.write(protein_name+', ')


    ''' Proteins to map to ID '''
    up_name2ttdtarg = json.load(open('output/protein2protein/up_name2ttdtarg.json'))
    protein_names = [name.split('_HUMAN')[0].strip().split(' ')[0] for name in list(up_name2ttdtarg.keys()) if '_HUMAN' in name]

    ''' Submit API mapping request'''
    job_id = submit_id_mapping_UniProtAPI(
                      from_db = 'Gene_Name',
                      to_db = 'UniProtKB-Swiss-Prot', 
                      ids = protein_names)

    ''' Get API call's results'''
    if check_id_mapping_results_ready_UniProtAPI(job_id):
        link = get_id_mapping_results_link_UniProtAPI(job_id)
        results = get_id_mapping_results_search_UniProtAPI(link)

    ''' Print summary '''
    try: 
        print(len(results['failedIds']), 'failed alignments')
    except: 
        results['failedIds'] = []
    try: 
        print(len(results['results']), 'successful alignments')
    except: 
        print()
    print('\nUnaligned IDs', len(set(results['failedIds'])))
 

    ''' Store mapping results '''
    name2up_id = get_to_uniprotid_from_genename_mapping_dict_UniProtAPI(results, [9606])
    name2up_id = switch_dictset_to_dictlist(name2up_id)
    json.dump(name2up_id, open('output/protein2protein/uniprotname2uniprotid_from_ttd.json','w'))
    print('Aligned Human IDs', len(name2up_id))   
    
    # Checks for an unexpected error in submitted names
    ID_not_in_submitted_IDs = set()
    for ID in results['failedIds']:
        if ID not in protein_names:
            ID_not_in_submitted_IDs.add(ID)
    if len(ID_not_in_submitted_IDs) > 0:
        print('These IDs were not submitted but were failedIds')
        print('There is probably a bad name submitted, e.g., ALT (123) instead of ALT')
        display(ID_not_in_submitted_IDs)
        
        
def align_ttd_target_id_to_uniprot_target_id():
    ttdtarg2up_id = dict()
    up_name2ttdtarg = json.load(open('output/protein2protein/up_name2ttdtarg.json'))

    for up_name, ttdtargs in up_name2ttdtarg.items():
        for ttdtarg in ttdtargs:
            try:
                up_name = up_name.split('_HUMAN')[0]
                up_ids = merged_name2up_id[up_name]
                for up_id in up_ids:
                    ttdtarg2up_id.setdefault(ttdtarg,set()).add(up_id)
            except:
                continue

    # Switch set to list
    ttdtarg2up_id = switch_dictset_to_dictlist(ttdtarg2up_id)
    
    with open('output/compound2compound/ttdtarg2up_id.json','w') as fout:
        json.dump(ttdtarg2up_id, fout)
    
    
def drugbank_id_targets_uniprot_id():
    db2up = dict()
    ttd2db = json.load(open('output/compound2compound/ttd2db.json'))
    ttdtarg2up_id = json.load(open('output/compound2compound/ttdtarg2up_id.json'))
    for ttddrug, ttdtargets in ttd2ttdtarget.items():

        # DrugBank Drugs -are- TTD Drug
        dbs = ttd2db[ttddrug]
        for ttdtarg in ttdtargets:
            try:
                # UniProt Proteins -are- TTD Target
                ups = ttdtarg2up_id[ttdtarg]

                # DrugBank Drug -targets- UniProt Protein
                for db in dbs:
                    for up in ups:
                        db2up.setdefault(db,set()).add(up)
            except:
                continue
                
                
    alignable_ttd_rows = list()
    target2compound_df = pd.read_excel('input/P1-07-Drug-TargetMapping.xlsx')
    for i in range(0,len(target2compound_df['DrugID'])):
        TTD_drug = target2compound_df['DrugID'].iloc[i]
        TTD_targ = target2compound_df['TargetID'].iloc[i]

        try:
            dbs = ttd2db[TTD_drug]
            for db in dbs:
                if db in db2up:

                    ups = ttdtarg2up_id[TTD_targ]
                    for up in ups:
                        if up in db2up[db]:
                            alignable_ttd_rows.append(i)
        except:
            continue
            
    print('target2compound_df', len(target2compound_df))
    print(len(alignable_ttd_rows), 'alignable_ttd_rows')
    alignable_targ2comp_df = target2compound_df.iloc[alignable_ttd_rows]
    
    return alignable_targ2comp_df


def export_drug_to_protein_from_ttd(alignable_targ2comp_df):
    ''' Output DrugBank-[moa]->Protein (Source: TTD)'''
    ttd2db = json.load(open('output/compound2compound/ttd2db.json'))
    ttdtarg2up_id = json.load(open('output/compound2compound/ttdtarg2up_id.json'))
    with open('output/compound2protein/edges_drugbank-moa->protein_ttd.csv','w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Compound (DrugBank)','Protein (UniProt)','Relationship'])

        for i in range(0, len(alignable_targ2comp_df)):
            dbs = ttd2db[alignable_targ2comp_df['DrugID'].iloc[i]]
            ups = ttdtarg2up_id[alignable_targ2comp_df['TargetID'].iloc[i]]
            moa = alignable_targ2comp_df['MOA'].iloc[i].lower()
            if moa == '.':
                moa = 'drug_targets_protein'

            # Write edge
            for db in dbs:
                for up in ups:
                    writer.writerow(['DrugBank_Compound:'+db, 'UniProt:'+up, '-'+moa+'->'])

    os.system('cp output/compound2protein/edges_drugbank-moa->protein_ttd.csv output/edges/edges_drugbank-moa->protein_ttd.csv')

    ''' Output MeSH-[moa]->Protein (Source: TTD)'''
    db2mesh = json.load(open('output/compound2compound/db2mesh.json'))

    with open('output/compound2protein/edges_MeSH-moa->protein_ttd.csv','w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Compound (MeSH)','Protein (UniProt)','Relationship'])

        for i in range(0, len(alignable_targ2comp_df)):
            dbs = ttd2db[alignable_targ2comp_df['DrugID'].iloc[i]]
            ups = ttdtarg2up_id[alignable_targ2comp_df['TargetID'].iloc[i]]
            moa = alignable_targ2comp_df['MOA'].iloc[i].lower()
            if moa == '.':
                moa = 'drug_targets_protein'

            # Write edge
            for db in dbs:
                try:
                    meshes = db2mesh[db]
                except:
                    continue
                for mesh in meshes:
                    for up in ups:
                        writer.writerow(['MeSH_Compound:'+mesh, 'UniProt:'+up, '-'+moa+'->'])

    os.system('output/compound2protein/edges_MeSH-moa->protein.csv output/edges/edges_MeSH-moa->protein.csv')
    
    
def display_and_save_relationship_details():
    d2p_TTD = pd.read_csv('output/compound2protein/edges_drugbank-moa->protein_ttd.csv')
    d2p_DB = pd.read_csv('output/compound2protein/edges_compoundsometargetsdetailed_drugbank.csv')
    d2p = pd.concat([d2p_DB, d2p_TTD]).drop_duplicates()

    print(len(d2p_DB), 'rows from DrugBank')
    print(len(d2p_TTD), 'rows from TTD')
    print(len(d2p), 'combined rows')
    #print(d2p['Relationship'].value_counts())

    d2p.to_csv('output/compound2protein/edges_drugbank-alltargetmoa->protein.csv', index=False)
    d2p.to_csv('output/edges/edges_drugbank-alltargetmoa->protein.csv', index=False)
    d2p.to_csv('output/edges_to_use/Compound_(DrugBank)_2_Protein_(UniProt).csv', index=False)

if __name__ == '__main__':
    try:
        root = parse_drugbank_xml()
    except:
        extract_drugbank_xml()
        root = parse_drugbank_xml()    
    map_compound_to_protein_df(root)
    map_compounds_to_targets(root)
    align_ttd_target_id_to_uniprot_target_name()
    align_ttd_target_id_to_uniprot_target_id()
    ttd2ttdtarget = map_ttd_drug_to_ttd_target()
    #align_ttd_target_to_uniprot_target()
    align_uniprot_name_to_uniprot_id()
    align_ttd_target_id_to_uniprot_target_id()
    target_to_compound = drugbank_id_targets_uniprot_id()
    print('alignable_targ2comp_df', len(target_to_compound))
    export_drug_to_protein_from_ttd(target_to_compound)
    display_and_save_relationship_details()