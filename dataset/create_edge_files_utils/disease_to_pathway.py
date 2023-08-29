import pandas as pd
import os 
import json
import csv
import math
from biomed_apis import *
from biomedkg_utils import *

def download_kegg_data():
    
    # drug to disease edges
    os.system('curl https://rest.kegg.jp/link/drug/disease > input/KEGG/kegg_drug_to_disease.tsv')

    # gene to disease edges
    os.system('curl https://rest.kegg.jp/link/hsa/disease > input/KEGG/kegg_gene_to_disease.tsv')

    # pathway to disease edges
    os.system('curl https://rest.kegg.jp/link/pathway/disease > input/KEGG/kegg_pathway_to_disease.tsv')

    # These are all pathway nodes for human in KEGG
    os.system('curl https://rest.kegg.jp/list/pathway/hsa > input/KEGG/human_pathways.tsv')

    # These are all gene nodes for human
    os.system('curl https://rest.kegg.jp/list/hsa > input/KEGG/human_genes.tsv')

    # These are all compounds nodes in KEGG
    os.system('curl https://rest.kegg.jp/list/compound > input/KEGG/compounds.tsv')

    # These are all drugs nodes in KEGG
    os.system('curl https://rest.kegg.jp/list/drug > input/KEGG/drug.tsv')

    # disease nodes
    os.system('curl https://rest.kegg.jp/list/disease > input/KEGG/disease.tsv')

    # no protein nodes in KEGG, seems to be genes only. Convert to NCBI protein:
    os.system('curl https://rest.kegg.jp/conv/hsa/ncbi-proteinid > input/KEGG/ncbi_to_kegg_gene.tsv')

    # gene to pathway edges
    os.system('curl https://rest.kegg.jp/link/pathway/hsa > input/KEGG/kegg_gene_to_pathway.tsv')

    # compound to pathway edges
    os.system('curl https://rest.kegg.jp/link/pathway/cpd > input/KEGG/kegg_cpd_to_pathway.tsv')

    # drug to pathway edges
    os.system('curl https://rest.kegg.jp/link/pathway/drug > input/KEGG/kegg_drug_to_pathway.tsv')

    # drug to compound edges
    os.system('curl https://rest.kegg.jp/link/cpd/drug > input/KEGG/kegg_drug_to_cpd.tsv')

    # drug to gene edges
    os.system('curl https://rest.kegg.jp/link/hsa/drug > input/KEGG/kegg_drug_to_gene.tsv')

    # drug to disease edges
    os.system('curl https://rest.kegg.jp/link/drug/disease > input/KEGG/kegg_drug_to_disease.tsv')

    # gene to disease edges
    os.system('curl https://rest.kegg.jp/link/hsa/disease > input/KEGG/kegg_gene_to_disease.tsv')

    # pathway to disease edges
    os.system('curl https://rest.kegg.jp/link/pathway/disease > input/KEGG/kegg_pathway_to_disease.tsv')
    
def align_disease_kegg_to_mesh():
    kegg_dis2mesh_df = pd.read_csv('input/KEGG/kegg_disease_to_mesh_and_omim.csv')

    mesh2kegg_disease = dict()
    kegg_disease2mesh = dict()

    for row_index in range(0,len(kegg_dis2mesh_df)):

        try:
            # KEGG Disease, MeSH Disease
            kegg_disease = kegg_dis2mesh_df['KEGG Disease'].iloc[row_index]
            mesh_diseases = kegg_dis2mesh_df['MeSH'].iloc[row_index].split('; ')
        except:
            continue

        # KEGG Disease -is- MeSH Disease
        for mesh_disease in mesh_diseases:
            if type(mesh_disease) == str and type(kegg_disease) == str:
                mesh2kegg_disease.setdefault(mesh_disease,set()).add(kegg_disease)
                kegg_disease2mesh.setdefault(kegg_disease,set()).add(mesh_disease)
    print(len(mesh2kegg_disease), 'MeSH-is-KEGG diseases')
    print(len(kegg_disease2mesh), 'KEGG-is-MeSH diseases')
    
    return mesh2kegg_disease, kegg_disease2mesh


def map_disease_to_kegg_pathway(kegg_disease2mesh):
    mesh_disease2kegg_pathway = dict()
    kegg_pathway2mesh_disease = dict()

    for line in open('input/KEGG/kegg_pathway_to_disease.tsv'):
        line = line.split('\t')
        try:
            # MeSH Disease, KEGG Pathway
            kegg_disease = line[0].split(':')[1]
            mesh_diseases = kegg_disease2mesh[kegg_disease]
            kegg_pathway = line[1].strip().replace('path:','path_')

            # MeSH Disease - KEGG Pathway
            for mesh_disease in mesh_diseases:
                mesh_disease2kegg_pathway.setdefault(mesh_disease,set()).add(kegg_pathway)
                kegg_pathway2mesh_disease.setdefault(kegg_pathway,set()).add(mesh_disease)
        except:
            continue

    print(str(len(mesh_disease2kegg_pathway)), 'MeSH Diseases - '+\
         str(len(kegg_pathway2mesh_disease)), 'KEGG Pathways')

    file = 'Disease_(MeSH)_2_Pathway_(KEGG).csv'
    outpath = os.path.join('output/disease2pathway',file)
    output_edgefile_onerel_noweight(
        outpath = os.path.join('output/disease2pathway',file),
        columns = ['Disease (MeSH)','Pathway (KEGG)','Relationship'],
        dictionary = mesh_disease2kegg_pathway,
        rel = '-disease_involves_pathway->',
        prefix_col1='MeSH_Disease:',
        prefix_col2='KEGG_Pathway:'
    )
    

def map_reactome_pathway_to_disease():
    # Currently requires Reactome Graph database, but is being improved to not re
    disease_nodes = pd.read_csv('input/human_disease_with_pathways_involved.csv')[['_id', 'identifier']].drop_duplicates()

    # Reactome database ID to Disease Ontology ID
    db2doid = dict(zip(disease_nodes._id, disease_nodes.identifier))
    db2doid = {int(dbId):int(doid) for dbId, doid in db2doid.items() if not math.isnan(doid) and not math.isnan(dbId)}

    pathway_nodes = pd.read_excel('input/human_pathways_involved_in_diseases.xlsx')[['_id', 'stId','displayName', 'name']].drop_duplicates()
    edges = pd.read_excel('input/human_pathway_disease_edges.xlsx').drop_duplicates()

    db2stID = dict(zip(pathway_nodes._id, pathway_nodes.stId))
    db2reactome_id = {int(dbId):reactome_id for dbId, reactome_id in db2stID.items() if not math.isnan(dbId)}

    '''Align: Reactome Pathway to DOID'''
    # Edges
    edge_list = list()

    for i in range(0,len(edges)):
        pathway = edges['_start'].iloc[i]
        disease = edges['_end'].iloc[i]
        rel = edges['_type'].iloc[i]

        try:
            pathway = db2reactome_id[pathway]
            disease = db2doid[disease]        
            edge_list.append([pathway, disease])

        except:
            print(pathway, disease)
            continue
            
    return edge_list


def export_mesh_disease_to_reactome_pathway(reactome_pw_to_disease_edges):
    # Output Edges
    file = 'Disease_(MeSH)_2_Pathway_(Reactome).csv'
    with open(os.path.join('output/disease2pathway', file),'w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Disease (MeSH)','Pathway (Reactome)','Relationship'])
        relationship = '-disease_involves_pathway->'

        for edge in reactome_pw_to_disease_edges:
            pathway = edge[0]
            doid_disease = edge[1]

            try:
                mesh_diseases = doid2mesh['DOID:'+str(doid_disease)]
                for mesh_disease in mesh_diseases:
                    writer.writerow(['MeSH_Disease:'+mesh_disease, 'Reactome_Pathway:'+pathway, relationship])
            except:
                continue

    df = pd.read_csv(os.path.join('output/disease2pathway', file)).drop_duplicates()
    df.to_csv(os.path.join('output/edges/', file), index=False)
    df.to_csv(os.path.join('output/edges_to_use/', file), index=False)
    
    print(str(len(set(df['Pathway (Reactome)'])))+' Reactome Pathways mapped to '+\
          str(len(set(df['Disease (MeSH)'])))+' MeSH Diseases '+\
          str(len(df))+' total relationships')
    
    return df


if __name__ == '__main__':

    '''KEGG PATHWAYS'''
    download_kegg_data()
    _, kegg_disease2mesh = align_disease_kegg_to_mesh()
    map_disease_to_kegg_pathway(kegg_disease2mesh)

    
    
    '''
    Note: Relations (only ~600) from Reactome are provided in the output folder. To
    redo the curation of these edges, more manual steps are needed. Efforts are
    underway to automate this process
    
    REACTOME PATHWAYS
    Download the Reactome database as a Neo4j graph: https://reactome.org/download-data/
    https://reactome.org/download/current/reactome.graphdb.tgz or https://reactome.org/download/current/reactome.graphdb.dump

    Neo4j Query:
        MATCH p=(pw:Pathway)-[]-(dis:Disease)
        WHERE toLower(pw.speciesName) = 'homo sapiens'
        RETURN dis

    OR
    Use the older version we provide. (Reactome has CC0 license for identifier mapping and 
    interaction data files: https://reactome.org/license
    '''
    # NOTE: REPLACE WITH DOWNLOAD LINKS
    #os.system('cp ../../Drug\ Repurposing\ Knowledge\ Graph/input/human_disease_with_pathways_involved.csv input/human_disease_with_pathways_involved.csv')
    #os.system('cp ../../Drug\ Repurposing\ Knowledge\ Graph/input/human_pathways_involved_in_diseases.xlsx input/human_pathways_involved_in_diseases.xlsx')
    #os.system('cp ../../Drug\ Repurposing\ Knowledge\ Graph/input/human_pathway_disease_edges.xlsx input/human_pathway_disease_edges.xlsx')

    try:
        reactome_pw_to_disease_edges = map_reactome_pathway_to_disease()
        doid2mesh = json.load(open('output/disease2disease/doid2mesh.json'))
        df = export_mesh_disease_to_reactome_pathway(reactome_pw_to_disease_edges)
    except:
        msg = 'Could not map Reactome pathways to disease'
        file = 'Disease_(MeSH)_2_Pathway_(Reactome).csv'
        reactome_path_dis_file = os.path.join('output/disease2pathway', file)
        reactome_pathways_to_disease_exists = os.path.exists(reactome_path_dis_file)
        if reactome_pathways_to_disease_exists:
            msg += ', but edges already exist (provided previously)'
            os.system(f'cp {reactome_path_dis_file} output/edges_to_use/{file}')
        print(msg)