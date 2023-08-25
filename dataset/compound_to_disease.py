# Some of the source code: https://github.com/mmayers12/metapaths/blob/main/1_code/08_CTD.ipynb
import json
import pandas as pd
import csv
import numpy as np
import os
from biomedkg_utils import *

def download_ctd_compound_to_disease_associations():
    #os.remove('input/CTD_chemicals_diseases.csv.gz')
    try: 
        os.remove('input/CTD_chemicals_diseases.csv')
    except:
        pass
    os.system('wget -N -P input/ http://ctdbase.org/reports/CTD_chemicals_diseases.csv.gz')
    os.system('gunzip input/CTD_chemicals_diseases.csv.gz .')
    

def process_ctd_compound_to_disease_associations():
    
    # Extract as dataframe
    rows = [] 
    with open('input/CTD_chemicals_diseases.csv','rt') as fin:
        data = csv.reader(fin)

        for i, row in enumerate(data):
            if i < 27:
                continue
            elif i == 27:  #NOTE: This may change with new downloads
                headers = row
            else:
                rows.append(row)

    df = pd.DataFrame(rows, columns=headers)

    # Filter by 'treats' and 'biomarker' relationships
    df2 = df.DirectEvidence.replace(r'^\s*$', np.nan, regex=True)
    df['DirectEvidence'] = df2
    df['DirectEvidence'].value_counts()

    cd_edges = df.dropna(subset=['DirectEvidence'])[['ChemicalID', 'DiseaseID', 'DirectEvidence']]
    temp = cd_edges.DiseaseID.str.replace('MESH:', '')
    cd_edges.loc[:,'DiseaseID'] = temp
    grouped = cd_edges.groupby(['DirectEvidence'])
    print(len(cd_edges), 'rows')
    
    therapeutic = grouped.get_group('therapeutic')
    therapeutic.to_csv('output/compound2disease/chem2disease_treats.csv')
    marker_or_mechanism = grouped.get_group('marker/mechanism')
    marker_or_mechanism.to_csv('output/compound2disease/chem2disease_marker.csv')

    db2mesh = json.load(open('output/compound2compound/db2mesh.json'))
    meshDrugs = list(set([mesh[0] for mesh in list((db2mesh.values()))]))
    tChem = list(set(therapeutic['ChemicalID']))
    mChem = list(set(marker_or_mechanism['ChemicalID']))
    t,m,tot  = 0,0,0
    for chem in meshDrugs:
        if chem in tChem:
            t += 1
        if chem in mChem:
            m += 1
        if chem in mChem or chem in tChem:
            tot += 1
    print(t,'/', len(tChem), 'Drugs with a drug-treats-disease relationship that have a with DrugBank2MeSH ID')
    print(m,'/', len(mChem), 'Drugs with a drug-marker-disease relationship that have a with DrugBank2MeSH ID')
    print(tot,'in either')
    print(len(meshDrugs),'/',len(db2mesh), 'MeSH drugs of drugs with Drugbank-MeSH ID')

    drugInd_ther = 0
    for mesh in therapeutic['ChemicalID']:
        if mesh in meshDrugs:
            drugInd_ther += 1
    print(drugInd_ther, 'drug-treats-disease relationships')

    drugInd_mark = 0
    for mesh in marker_or_mechanism['ChemicalID']:
        if mesh in meshDrugs:
            drugInd_mark += 1
    print(drugInd_mark, 'drug-marker-disease relationships')

    return cd_edges
    

def export_ctd_compound_to_disease_associations(cd_edges):
    # Export Drug-Disease 
    compTREATSdis_ctd, compMARKERdis_ctd = dict(), dict()

    for i in range(0,len(cd_edges)):
        comp = cd_edges['ChemicalID'].iloc[i]
        dis = cd_edges['DiseaseID'].iloc[i]
        rel = cd_edges['DirectEvidence'].iloc[i]

        if rel == 'therapeutic':
            compTREATSdis_ctd.setdefault(comp,[]).append(dis)
        elif rel == 'marker/mechanism':
            compMARKERdis_ctd.setdefault(comp,[]).append(dis)

    with open('output/compound2disease/meshCompound-TREATS-meshDisease_ctd.json','w') as fout:
        json.dump(compTREATSdis_ctd, fout)

    with open('output/compound2disease/meshCompound-MARKER_OF-meshDisease_ctd.json','w') as fout:
        json.dump(compMARKERdis_ctd, fout)                            
        
        
# INSTRUCTIONS: download PathFX
def map_pathfx_compound_treats_disease():
    # Note: All drugbank drugs here have at least one UMLS CUI for itself
    os.system('wget -N -P input/ https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6285459/bin/pcbi.1006614.s001.xlsx')
    pathfx = pd.read_excel('input/pcbi.1006614.s001.xlsx')[['DrugBankID','Indication CUI']]

    umls2mesh = json.load(open('output/otherMappings/umls2mesh.json'))
    db2mesh = json.load(open('output/compound2compound/db2mesh.json'))
    
    # Map MeSH drug to MeSH disease indication
    compTREATSdis_pathfx = dict()
    for i in range(0,len(pathfx)):

        try:
            # Drug/Compound
            db_drug = pathfx['DrugBankID'][i]
            mesh_comps = db2mesh[db_drug]

            # Disease
            umls_dis = pathfx['Indication CUI'][i]
            meshes_dis = umls2mesh[umls_dis]

            # Drug-treats-disease
            for mesh_dis in meshes_dis:
                for mesh_comp in mesh_comps:
                    compTREATSdis_pathfx.setdefault(mesh_comp,set()).add(mesh_dis)

        except:
            continue


    for k in compTREATSdis_pathfx.copy():
        compTREATSdis_pathfx[k] = list(compTREATSdis_pathfx[k])

    with open('output/compound2disease/meshCompound-TREATS-meshDisease_pathfx.json','w') as fout:
        json.dump(compTREATSdis_pathfx, fout)
        
        
def align_doid_to_mesh_via_hetionet_disease():
    os.system('wget -N -P input/ https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/HumanDO.obo')
    umls2mesh = json.load(open('output/otherMappings/umls2mesh.json'))

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
        
        
def map_hetionets_drugbank_compound_treats_mesh_disease():
    os.system('wget -N -P input/ https://raw.githubusercontent.com/dhimmel/indications/49cf8e3858f124e552051fc1e8c1c229a43dadf6/curation/results-three-curators.tsv')
    compTREATSdis_hetionet = dict()
    doid2mesh = json.load(open('output/disease2disease/doid2mesh.json'))
    db2mesh = json.load(open('output/compound2compound/db2mesh.json'))
    
    with open('input/results-three-curators.tsv') as fin:
        
        for line in fin:
            line = line.split('\t')
            try:
                diseases = doid2mesh[line[8]]
                comps = db2mesh[line[9]]
                moa = line[10]

                if 'DM' in moa: #if the drug treats the disease
                    for comp in comps:
                        for disease in diseases:
                            compTREATSdis_hetionet.setdefault(comp,set()).add(disease)
            except:
                continue

    # Compound-treats-Disease from Hetionet
    count = 0
    for k,v in compTREATSdis_hetionet.items():
        for V in v:
            count += 1
    print(count, 'Compound-treats-disease from Hetionet')
    compTREATSdis_hetionet = switch_dictset_to_dictlist(compTREATSdis_hetionet)
    
    json.dump(compTREATSdis_hetionet, open('output/compound2disease/hetionet_comp_to_dis.json','w'))
    
        
''' Merge compound-treats-disease dictionaries '''
def merge_dict_values_of_lists(dMain, d2):
    for key,vals in d2.items():
        if key in dMain:
            new_vals = set(dMain[key])
            for val in vals:
                new_vals.add(val)
            dMain[key] = list(new_vals)
        else:
            dMain[key] = vals

def merge_compound_treats_disease_associations():
    ctd_ctd = json.load(open('output/compound2disease/meshCompound-TREATS-meshDisease_ctd.json'))
    pathfx_ctd = json.load(open('output/compound2disease/meshCompound-TREATS-meshDisease_pathfx.json'))
    hetionet_ctd = json.load(open('output/compound2disease/hetionet_comp_to_dis.json'))
    ctd = ctd_ctd.copy()
    
    merge_dict_values_of_lists(ctd, pathfx_ctd)
    merge_dict_values_of_lists(ctd, hetionet_ctd)

    ''' Print summary counts '''
    rels = 0
    for compound, diseases in ctd.items():
        rels += len(diseases)
    print(rels, 'compound-treats-disease relationships')
    
    return ctd

def export_all_compound_treats_disease_associations(compound_treats_disease):
    # MeSH Compound -TREATS-> MeSH Disease
    output_edgefile_onerel_noweight(
                outpath = 'output/compound2disease/edges_meshCompound-TREATS->meshDisease.csv',
                columns = ['Compound (MeSH)','Disease (MeSH)','Relationship'],
                dictionary = compound_treats_disease,
                rel = '-treats->',
                prefix_col1 = 'MeSH_Compound:',
                prefix_col2 = 'MeSH_Disease:')

    df = pd.read_csv('output/edges/edges_meshCompound-TREATS->meshDisease.csv')
    df.to_csv('output/edges_to_use/Compound_(MeSH)_treats_Disease_(MeSH).csv', index=False)

    dbcompTREATSdis = dict()
    for comp, diseases in compound_treats_disease.items():
        try:
            db_comps = mesh2db[comp]
            for db in db_comps:
                for disease in diseases:
                    dbcompTREATSdis.setdefault(db, list()).append(disease)
        except:
            continue


    ### DrugBank Compound -TREATS-> drugbank Disease
    output_edgefile_onerel_noweight(outpath = 'output/compound2disease/edges_drugbankCompound-TREATS->meshDisease.csv',
                                   columns = ['Compound (DrugBank)','Disease (MeSH)','Relationship'],
                                   dictionary = dbcompTREATSdis,
                                   rel = '-treats->',
                                   prefix_col1 = 'DrugBank_Compound:',
                                   prefix_col2 = 'MeSH_Disease:')

    df = pd.read_csv('output/edges/edges_drugbankCompound-TREATS->meshDisease.csv')
    df.to_csv('output/edges_to_use/Compound_(DrugBank)_treats_Disease_(MeSH).csv', index=False)
    
    
    
if __name__ == '__main__':
    # CTD
    download_ctd_compound_to_disease_associations()        
    ctd_cd_edges = process_ctd_compound_to_disease_associations()
    export_ctd_compound_to_disease_associations(ctd_cd_edges)

    # PathFX (Download this)
    map_pathfx_compound_treats_disease()

    # Hetionet
    align_doid_to_mesh_via_hetionet_disease()
    map_hetionets_drugbank_compound_treats_mesh_disease()

    # All
    compound_treats_disease = merge_compound_treats_disease_associations()
    export_all_compound_treats_disease_associations(compound_treats_disease)