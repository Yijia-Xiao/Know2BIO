import os
import json
import csv
import pandas as pd


def download_drug_to_side_effect():
    os.system('wget -N -P input/ http://sideeffects.embl.de/media/download/meddra_all_se.tsv.gz')
    try:
        os.remove('input/meddra_all_se.tsv')
    except:
        pass
    os.system('gunzip input/meddra_all_se.tsv.gz')
    
    
def map_pubchem_compound_to_side_effect():
    # 1 & 2: STITCH compound ids (flat/stereo, see above)
    # 3: UMLS concept id as it was found on the label
    # 4: MedDRA concept type (LLT = lowest level term, PT = preferred term; in a few cases the term is neither LLT nor PT)
    # 5: UMLS concept id for MedDRA term
    # 6: side effect name
    #All side effects found on the labels are given as LLT. Additionally, the PT is shown. There is at least one
    #PT for every LLT, but sometimes the PT is the same as the LLT.
    umls2mesh = json.load(open('output/otherMappings/umls2mesh.json'))
    pubchem2side = dict()
    with open('input/meddra_all_se.tsv') as fin:
        for i,line in enumerate(fin):
            line = line.split('\t')
            stitch_sterio = line[1] # Hetionet did stereo https://github.com/dhimmel/SIDER4/blob/2acca0b065e736bc99702906024efd4718e502ee/SIDER4.ipynb
            pubchem = stitch_sterio[3:].strip('0')
            side = line[2]
            meddra_type = line[3]
            side2 = line[4]
            if meddra_type == 'PT':   # Only reduces side effects from 146,518 to 145,617. 1,550 drugs either way.
                try: mesh_sides = umls2mesh[side]
                except: continue
                for mesh_side in mesh_sides:
                    pubchem2side.setdefault(pubchem,set()).add(mesh_side)

    for pc,se in pubchem2side.items():
        pubchem2side[pc] = list(se)  
        
    json.dump(pubchem2side, open('output/compound2disease/pubchem2side.json','w'))

    
def map_drug_to_side_effect():
    pubchem2side = json.load(open('output/compound2disease/pubchem2side.json'))
    pubchem2db = json.load(open('output/compound2compound/pubchemcomp2db.json'))
    pubchems = list(pubchem2db.keys())
    sider_pc = list(pubchem2side.keys())
    db2se = dict()
    i = 0
    for pc in pubchems:
        if pc in sider_pc:
            ses = pubchem2side[pc]
            dbs = pubchem2db[pc]
            for db in dbs:
                for se in ses:
                    db2se.setdefault(db,[]).append(se)
                    
                    
            
    dbs = list(db2se.keys())
    print('Drugs', len(dbs))
    se = list(db2se.values())
    se1 = []
    for s in se:
        se1+=s 

    se2 = set()
    for s in se:
        se2 = se2.union(set(s))
    print('Side Effects',len(se2))
    print('Edges',len(se1))
    return db2se


def export_drug_to_side_effect():
   
    ''' DrugBank Drug '''
    pubchem2side = json.load(open('output/compound2disease/pubchem2side.json'))
    pubchem2db = json.load(open('output/compound2compound/pubchemcomp2db.json'))
    
    pubchems = list(pubchem2db.keys())
    sider_pc = list(pubchem2side.keys())

    with open('output/compound2disease/edges_drugbank2meshsideeffect.csv','w') as fout:
        writer = csv.writer(fout) 
        writer.writerow(['Compound (DrugBank)','Disease (MeSH)','Relationship'])
        db2se = dict()
        for pc in pubchems:
            if pc in sider_pc:

                ses = pubchem2side[pc] # Side effects
                dbs = pubchem2db[pc]        # DrugBank IDs

                for db in dbs:
                    for se in ses:
                        db2se.setdefault(db,[]).append(se)
                        writer.writerow(['DrugBank_Compound:'+db,'MeSH_Disease:'+se,'-causes_side_effect->'])
    
    os.system('cp output/compound2disease/edges_drugbank2meshsideeffect.csv output/edges/edges_drugbank2meshsideeffect.csv')
    file = 'Compound_(DrugBank)_sideeffect_Disease_(MeSH).csv'
    df = pd.read_csv('output/edges/edges_drugbank2meshsideeffect.csv')
    df.to_csv(os.path.join('output/edges_to_use', file), index=False)

    db2mesh = json.load(open('output/compound2compound/db2mesh.json'))
    
    
    '''MeSH Drug'''
    pubchems = list(pubchem2db.keys())
    sider_pc = list(pubchem2side.keys())

    with open('output/compound2disease/edge_meshcompound2sideeffect.csv','w') as fout:
        writer = csv.writer(fout) 
        writer.writerow(['Compound (MeSH)','Disease (MeSH)','Relationship'])
        mesh2se = dict()
        for pc in pubchems:
            if pc in sider_pc:

                ses = pubchem2side[pc] # Side effects
                dbs = pubchem2db[pc]        # DrugBank IDs

                # DrugBank Drugs
                for db in dbs:

                    # MeSH Drugs
                    try: meshes = db2mesh[db]
                    except: continue
                    for mesh in meshes:

                        # Side Effect
                        for se in ses:
                            mesh2se.setdefault(mesh,[]).append(se)
                            writer.writerow(['MeSH_Compound:'+mesh,'MeSH_Disease:'+se,'-causes_side_effect->'])

    os.system('cp output/compound2disease/edge_meshcompound2sideeffect.csv output/edges/edge_meshcompound2sideeffect.csv')

    file = 'Compound_(DrugBank)_sideeffect_Disease_(MeSH).csv'
    df = pd.read_csv('output/edges/edge_meshcompound2sideeffect.csv')
    #df.to_csv(os.path.join('output/edges_to_use', file), index=False)

    
if __name__ == '__main__':
    download_drug_to_side_effect()
    map_pubchem_compound_to_side_effect()
    db2se = map_drug_to_side_effect()
    export_drug_to_side_effect()