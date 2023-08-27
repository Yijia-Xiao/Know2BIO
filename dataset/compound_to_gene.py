import os
import pandas as pd 
import numpy as np
import csv
import json
from biomedkg_utils import *

def download_ctd_compound_to_gene():
    os.system('rm input/CTD_chem_gene_ixns.tsv')
    os.system('wget -N -P input/ http://ctdbase.org/reports/CTD_chem_gene_ixns.tsv.gz')
    os.system('gunzip input/CTD_chem_gene_ixns.tsv.gz')

def process_ctd_compound_to_gene():
    with open('input/formatted_CTD_chem_gene_ixns.tsv','w') as fout:
        with open('input/CTD_chem_gene_ixns.tsv') as fin:
            for i, line in enumerate(fin):
                if i >= 27:
                    fout.write(line)

    df = pd.read_table('input/formatted_CTD_chem_gene_ixns.tsv')
    df = df.loc[np.where(df['Organism']=='Homo sapiens')]
    
    
    gene_name2id = dict()
    for i in range(0,len(df)):
        gene_name = df['GeneSymbol'].iloc[i]
        gene_ID = df['GeneID'].iloc[i]
        gene_name2id[gene_name] = gene_ID # 1 to 1
        
        
    intact1 = set(df['InteractionActions'])
    indiv_intact = []
    for intact in intact1:
        row = intact.split('|')
        ia = [r for r in row]
        indiv_intact += ia
        
        
    acts = list({act.split('^')[1] for act in indiv_intact})
    act_dict = dict()
    for act in acts:
        act_dict[act] = 0

    for i in range(0,len(df)):
        ixs = df['InteractionActions'].iloc[i].split('|')
        act = list({act.split('^')[1] for act in ixs})
        for a in act:
            act_dict[a] += 1

    #l = sorted(act_dict.items(), key = lambda x: x[1], reverse=True)
    #d = dict(l)
    
    
    
    acts = list(indiv_intact)
    intact_dict = dict()
    for act in acts:
        intact_dict[act] = 0

    for i in range(0,len(df)):
        ixs = df['InteractionActions'].iloc[i].split('|')
        for ix in ixs:
            intact_dict[ix] += 1

    #l = sorted(intact_dict.items(), key = lambda x: x[1], reverse=True)
    #d = dict(l)

        
    
    obv = ['increases^expression','decreases^expression',
       'decreases^reaction','increases^reaction',
       'increases^activity','increases^abundance',
       'increases^abundance','decreases^activity']

    goodrow = 0

    for i in range(0,len(df)):
        ixs = df['InteractionActions'].iloc[i].split('|')
        for ix in ixs:
            if ix in obv:
                goodrow += 1
                break

    print(goodrow, len(df))
    return df
    
    
def export_ctd_compound_to_gene(df):
    # Open output file
    with open('output/compound2gene/edges_compound2gene.csv','w') as fout1:
        writer = csv.writer(fout1)
        writer.writerow(['Compound (MeSH)','Gene (Entrez)','Relationship'])

        # Iterate through table
        for i in range(0,len(df)):

            # Find rows with compound-expression-gene relationships
            intActs1 = df['InteractionActions'].iloc[i]
            if 'expression' in intActs1:
                intActs = intActs1.split('|')
                chemID = df['ChemicalID'].iloc[i]
                chemNm = df['# ChemicalName'].iloc[i]
                geneID = df['GeneID'].iloc[i]
                geneSy = df['GeneSymbol'].iloc[i]
                desc   = df['Interaction'].iloc[i]

                # Case: 'Expression' is the only relationship
                if len(intActs) == 1:
                    rel = intActs[0].split('^')[0]
                    if 'affects' in rel:
                        continue
                    writer.writerow(['MeSH_Compound:'+chemID, 'Entrez:'+str(int(geneID)), '-'+rel+'->'])

                # Case: 'Expression' + 'Reaction' relationships
                elif len(intActs) == 2 and 'creases^reaction' in intActs1 and 'mutant' not in desc and 'affects' not in intActs1:                
                    if 'increases^reaction' in intActs and 'increases^expression' in intActs:
                        rel = 'increases'#^expression'
                    elif 'increases^reaction' in intActs and 'decreases^expression' in intActs:
                        rel = 'decreases'#^expression'
                    elif 'decreases^reaction' in intActs and 'increases^expression' in intActs:
                        rel = 'decreases'#^expression'
                    elif 'decreases^reaction' in intActs and 'decreases^expression' in intActs:
                        rel = 'increases'#^expression'
                    writer.writerow(['MeSH_Compound:'+chemID, 'Entrez:'+str(int(geneID)), rel])

            # Case: 'Activity' is the only relationship
            intActs = intActs1.split('|')
            if 'activity' in intActs1 and len(intActs) == 1:
                chemID = df['ChemicalID'].iloc[i]
                geneID = df['GeneID'].iloc[i]
                desc   = df['Interaction'].iloc[i]
                rel    = intActs[0].split('^')[0]
                writer.writerow(['MeSH_Compound:'+chemID, 'Entrez:'+str(int(geneID)), '-'+rel+'->'])

    meshcompound2gene = pd.read_csv('output/compound2gene/edges_compound2gene.csv')[['Compound (MeSH)','Gene (Entrez)', 'Relationship']].drop_duplicates()
    meshcompound2gene.to_csv('output/compound2gene/edges_compound2gene.csv', index=False)
    meshcompound2gene.to_csv('output/edges/edges_compound2gene.csv', index=False)
    

def print_summary_stats_compound_to_gene():
    meshcompound2gene = pd.read_csv('output/compound2gene/edges_compound2gene.csv')[['Compound (MeSH)','Gene (Entrez)', 'Relationship']]
    mesh2db = json.load(open('output/compound2compound/mesh2db.json'))
    print('Rows', len(meshcompound2gene))
    meshdrugs = [meshdrug.split(':')[1] for meshdrug in set(meshcompound2gene['Compound (MeSH)'])]
    print('MeSH Compounds', len(meshdrugs))

    gooddbs = list()
    for m in meshdrugs:
        if m in mesh2db:
            gooddbs.append(mesh2db[m])
    print('MeSH Compounds with DrugBank Drug ID', len(gooddbs))

    dbspresent = set()
    for dbl in gooddbs:
        for db in dbl:
            dbspresent.add(db)
    print('DrugBank Drug IDs mapped from MeSH drugs/compounds', len(dbspresent))
    

def map_mesh_gene_to_compound(df):
    with open('output/compound2gene/edges_gene2compound.csv','w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Gene (Entrez)','Compound (MeSH)','Relationship'])
        for i in range(0,len(df)):
            intrAxs = df['InteractionActions'].iloc[i]

            if 'abundance' in intrAxs:
                intrAx = intrAxs.split('|')
                chemID = df['ChemicalID'].iloc[i]
                geneID = df['GeneID'].iloc[i]
                desc = df['Interaction'].iloc[i]

                # If 'inc/dec^abundance' is the only relationship
                if len(intrAx) == 1:
                    if True not in [phrase in desc for phrase in ['affects','modified form','mutant']]:
                        desc = desc.split(' ')
                        rel = intrAx[0].split('^')[0]
                        if desc[1] in ('protein','gene','mRNA'):
                            writer.writerow(['Entrez:'+str(int(geneID)), 'MeSH_Compound:'+chemID, '-'+rel+'->'])   

    meshgene2compound = pd.read_csv('output/compound2gene/edges_gene2compound.csv').drop_duplicates()[['Gene (Entrez)', 'Compound (MeSH)','Relationship']]
    meshgene2compound.to_csv('output/compound2gene/edges_gene2compound.csv',index=False)
    meshgene2compound.to_csv('output/edges/edges_gene2compound.csv',index=False)
    

def print_summary_stats_gene_to_compound():
    meshgene2compound = pd.read_csv('output/compound2gene/edges_gene2compound.csv')[['Gene (Entrez)', 'Compound (MeSH)','Relationship']]
    print('Rows', len(meshgene2compound))
    meshdrugs = [meshdrug.split(':')[1] for meshdrug in set(meshgene2compound['Compound (MeSH)'])]
    print('MeSH Compounds', len(meshdrugs))
    mesh2db = json.load(open('output/compound2compound/mesh2db.json'))
    gooddbs = list()
    for m in meshdrugs:
        if m in mesh2db:
            gooddbs.append(mesh2db[m])
    print('MeSH Compounds with DrugBank Drug ID', len(gooddbs))

    dbspresent = set()
    for dbl in gooddbs:
        for db in dbl:
            dbspresent.add(db)
    print('DrugBank Drug IDs mapped from MeSH drugs/compounds', len(dbspresent))


def map_gene_to_compound_via_kegg():
    '''Align KEGG Drug to MeSH Compound'''
    keggdrug2db = json.load(open('output/compound2compound/keggdrug2db.json'))
    db2mesh = json.load(open('output/compound2compound/db2mesh.json'))
    keggdrug2mesh = dict()
    for keggdrug, dbs in keggdrug2db.items():

        for db in dbs:
            try:
                meshes = db2mesh[db]
                for mesh in meshes:
                    keggdrug2mesh.setdefault(keggdrug, set()).add(mesh)
            except:
                continue
                

    '''Create Compound to Gene via KEGG'''
    meshcompound2gene_kegg, gene2meshcompound_kegg, drugbank2gene_kegg, gene2drugbank_kegg = dict(), dict(), dict(), dict()
    if not os.path.exists('input/KEGG/kegg_drug_to_gene.tsv'):
        os.system('curl https://rest.kegg.jp/link/hsa/drug > input/KEGG/kegg_drug_to_gene.tsv')
    
    for line in open('input/KEGG/kegg_drug_to_gene.tsv'):
        line = line.split('\t')

        # KEGG Drug
        kegg_drug = line[0].split(':')[1]

        # MeSH Drugs
        try: mesh_drugs = keggdrug2mesh[kegg_drug]
        except: mesh_drugs = ''

        # DrugBank Drugs
        try: db_drugs = keggdrug2db[kegg_drug]
        except: db_drugs = ''

        # Skip if no DrugBank or MeSH drug mapping/alignment
        if db_drugs == ''  and mesh_drugs == '':
            continue

        # Gene
        entrez_gene = str(int(line[1].strip().split(':')[1]))

        # MeSH Drug - Gene
        if mesh_drugs != '':
            for mesh_drug in mesh_drugs:
                meshcompound2gene_kegg.setdefault(mesh_drug, set()).add(entrez_gene)
                gene2meshcompound_kegg.setdefault(entrez_gene, set()).add(mesh_drug)

        # DrugBank Drug - Gene
        if db_drugs != '':
            for db_drug in db_drugs:
                drugbank2gene_kegg.setdefault(db_drug, set()).add(entrez_gene)
                gene2drugbank_kegg.setdefault(entrez_gene, set()).add(db_drug)

    print(len(meshcompound2gene_kegg), len(gene2meshcompound_kegg), len(drugbank2gene_kegg), len(gene2drugbank_kegg))
    
    
    output_edgefile_onerel_noweight('output/compound2gene/edges_drugbankcompound2gene_kegg.csv',
                               ['Compound (DrugBank)','Gene (Entrez)','Relationship'],
                               drugbank2gene_kegg,
                               '-associated_with->',
                                'DrugBank_Compound:',
                                'Entrez:',
                                edges_to_use_folder=False)

    output_edgefile_onerel_noweight('output/compound2gene/edges_meshcompound2gene_kegg.csv',
                               ['Compound (MeSH)','Gene (Entrez)','Relationship'],
                               meshcompound2gene_kegg,
                               '-associated_with->',
                                'MeSH_Compound:',
                                'Entrez:',
                                edges_to_use_folder=False)
    
def export_all_merged_gene_to_compound_relationships():
    gene2comp_df1 = pd.read_csv('output/edges/edges_gene2compound.csv')
    gene2comp_df1.to_csv('output/edges_to_use/Gene_(Entrez)_2_Compound_(MeSH).csv',
                         index=False)
    gene2comp_df2 = pd.read_csv('output/edges/edges_gene2compound.csv')
    gene2comp_df2.to_csv('output/edges_to_use/Gene_(Entrez)_2_Compound_(DrugBank).csv', 
                         index=False)

    comp2gene_df1 = pd.read_csv('output/edges/edges_compound2gene.csv')
    comp2gene_df2 = pd.read_csv('output/compound2gene/edges_meshcompound2gene_kegg.csv')
    comp2gene_df = pd.concat([comp2gene_df2, comp2gene_df1]).drop_duplicates()
    comp2gene_df.to_csv('output/edges_to_use/Compound_(MeSH)_2_Gene_(Entrez).csv', index=False)
    
if __name__ == '__main__':    
    download_ctd_compound_to_gene()
    df = process_ctd_compound_to_gene()
    export_ctd_compound_to_gene(df)
    print_summary_stats_compound_to_gene()
    map_mesh_gene_to_compound(df)
    print_summary_stats_gene_to_compound()
    map_gene_to_compound_via_kegg()
    export_all_merged_gene_to_compound_relationships()