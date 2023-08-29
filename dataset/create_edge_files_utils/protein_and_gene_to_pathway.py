import os
import csv
import pandas as pd
import json
from biomedkg_utils import output_edgefile_onerel_noweight


def map_protein_to_reactome_pathway():
    os.system('wget -N -P input/ https://reactome.org/download/current/UniProt2Reactome_PE_Pathway.txt')

    protein2pathway, pathway2protein = dict(), dict()
    file = 'Protein_(UniProt)_2_Pathway_(Reactome).csv'

    with open(os.path.join('output/pathway2protein/', file),'w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Protein (UniProt)','Pathway (Reactome)', 'Relationship'])

        issue_lines = []
        with open('input/UniProt2Reactome_PE_Pathway.txt') as fin:
            for line in fin:
                nline = line.strip().split('\t')

                # Protein -involved in-> Pathway
                try:  
                    species, protein, pathway = nline[7], nline[0], nline[3]

                    if species.lower() == 'homo sapiens':
                        writer.writerow(['UniProt:'+protein, 'Reactome_Pathway:'+pathway, '-may_participate_in->'])
                        protein2pathway.setdefault(protein, set()).add(pathway)
                        pathway2protein.setdefault(pathway, set()).add(protein)

                # Save other lines with problems
                except: 
                    if 'homo sapiens' in line.lower():
                        issue_lines.append(line)          

    # Drop duplicates
    df = pd.read_csv(open(os.path.join('output/pathway2protein/', file))).drop_duplicates().reset_index(drop=True)
    df.to_csv(os.path.join('output/edges_to_use/', file), index=False)
    print('Human Reactome Protein-Pathway:', len(df))

    # Save proteins from protein-pathway relationships
    prot2pw = df.copy(deep=True)
    prots = [prot.split(':')[1] for prot in set(prot2pw['Protein (UniProt)'])]
    json.dump(prots, open('output/pathway2protein/protein_to_reactome_pathway.json','w'))

    
def map_gene_to_reactome_pathway():
    '''Load Gene-Pathway Data'''
    os.system('wget -N -P input/ https://reactome.org/download/current/NCBI2Reactome_PE_Pathway.txt')
    gene2pathway_df = pd.read_table('input/NCBI2Reactome_PE_Pathway.txt', low_memory=False)
    gene2pathway_df.columns = ['Entrez Gene ID','Reactome Gene ID','Reactome Pathway Name','Reactome Pathway', 'Link', 'Pathway Name','Evidence','Species']
    gene2pathway_df = gene2pathway_df[gene2pathway_df['Species']=='Homo sapiens'][['Entrez Gene ID','Reactome Gene ID', 'Reactome Pathway Name','Reactome Pathway']].drop_duplicates()
    print(len(gene2pathway_df),'rows')

    '''Get Gene-Pathway Data'''
    gene2pathway, pathway2gene = dict(), dict()

    for i in range(0,len(gene2pathway_df)):
        entrez = gene2pathway_df['Entrez Gene ID'].iloc[i]
        pathway = gene2pathway_df['Reactome Pathway'].iloc[i]

        gene2pathway.setdefault(entrez, set()).add(pathway)
        pathway2gene.setdefault(pathway, set()).add(entrez)

    '''Export Edges'''
    file = 'Gene_(Entrez)_2_Pathway_(Reactome).csv'
    outpath = os.path.join('output/pathway2gene',file)
    output_edgefile_onerel_noweight(
        outpath = outpath,           
        columns = ['Gene (Entrez)','Pathway (Reactome)','Relationship'], 
        dictionary = gene2pathway, 
        rel = '-may_participate_in->', 
        prefix_col1='Entrez:', 
        prefix_col2='Reactome_Pathway:')

    df = pd.read_csv(os.path.join('output/pathway2gene/', file))
    df.to_csv(os.path.join('output/edges',file),index=False)
    df.to_csv(os.path.join('output/edges_to_use/',file),index=False)

    print('Human Reactome Gene-Pathway:', len(df))

    
def map_gene_to_kegg_pathway():
    os.system('curl https://rest.kegg.jp/link/pathway/hsa > input/kegg_gene_to_pathway.tsv')

    ''' Entrez Gene -in- KEGG Pathway '''
    kegg_pathway_HAS_entrez_gene = dict()
    entrez_gene_IN_kegg_pathway = dict()

    with open('input/kegg_gene_to_pathway.tsv') as fin:
        for line in fin:
            line = line.split('\t')
            kegg_gene = line[0].strip()
            kegg_pathway = line[1].strip().replace('path:','path_')

            try:
                entrez_gene = kegg_gene.split(':')[1]
                kegg_pathway_HAS_entrez_gene.setdefault(kegg_pathway,set()).add(entrez_gene)
                entrez_gene_IN_kegg_pathway.setdefault(entrez_gene,set()).add(kegg_pathway)
            except:
                continue        

    print(len(kegg_pathway_HAS_entrez_gene), 'KEGG Pathways with Entrez Genes')
    print(len(entrez_gene_IN_kegg_pathway), 'Entrez Genes with KEGG Pathways')

    file = 'Gene_(Entrez)_2_Pathway_(KEGG).csv'
    outpath = os.path.join('output/pathway2gene/',file)
    output_edgefile_onerel_noweight(
        outpath = outpath,
        columns = ['Gene (Entrez)','Pathway (KEGG)','Relationship'],
        dictionary = entrez_gene_IN_kegg_pathway,
        rel = '-may_participate_in->',
        prefix_col1='Entrez:',
        prefix_col2='KEGG_Pathway:')

    df = pd.read_csv(outpath)
    df.to_csv(os.path.join('output/edges',file),index=False)
    df.to_csv(os.path.join('output/edges_to_use',file), index=False)
    

def map_protein_to_smpdb_pathway():
    os.system('wget -N -P input/smpdb_files https://smpdb.ca/downloads/smpdb_proteins.csv.zip')
    os.system('unzip input/smpdb_files/smpdb_proteins.csv.zip -d input/smpdb_files')

    # Get all SMPDB pathway names
    files = 0
    all_smpdb_pws = list()
    directory = 'input/smpdb_files'
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        if os.path.isfile(f):
            files += 1
            all_smpdb_pws.append(filename.split('_proteins.csv')[0])
    pw_notfound = []

    # Get all SMPDBPathway-Proteins for Pathways targeted by DrugBank
    with open('output/pathway2protein/edges_DrugBankTargetprotein2SMPDBpathway.csv','w',newline='') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Protein (UniProt)','Pathway (SMPDB)','Relationship'])
        pw2prot = dict()
        for pw in all_smpdb_pws:
            try: df = pd.read_csv('input/smpdb_files/'+pw+'_proteins.csv')
            except: continue
            prots = []
            for i in range(0,len(df)):
                prot = df['Uniprot ID'][i]
                if prot == 'Unknown' or type(prot) == float:
                    continue
                prots.append(prot)
                writer.writerow([prot,pw,'-participates_in->'])
            pw2prot[pw] = prots

    # Get all Protein-SMPDBPathway Dictionary
    count = 0
    all_prot2pws = dict()
    directory = 'input/smpdb_files'
    for filename in os.listdir(directory):
        print(count, end='\r'); count += 1
        f = os.path.join(directory, filename)
        if os.path.isfile(f) and f.endswith('_proteins.csv'):
            df = pd.read_csv(f)
            pw = filename.split('_proteins.csv')[0]
            for i in range(0,len(df)):
                prot = df['Uniprot ID'][i]
                if prot == 'Unknown' or type(prot) == float: 
                    continue
                all_prot2pws.setdefault(prot,[]).append(pw)        

    # Get all SMPDBPathway-Protein Dictionary
    with open('output/pathway2protein/edges_protein-INVOLVED_IN->smpdbpathway.csv','w',newline='') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Protein (UniProt)','Pathway (SMPDB)','Relationship'])
        all_pw2prots = dict()
        for prot,pws in all_prot2pws.items():
            for pw in pws:
                all_pw2prots.setdefault(pw,[]).append(prot)
                writer.writerow(['UniProt:'+prot,'SMPDB_Pathway:'+pw,'-participates_in->'])

    df = pd.read_csv('output/pathway2protein/edges_protein-INVOLVED_IN->smpdbpathway.csv').drop_duplicates()
    df.to_csv('output/edges/edges_protein-INVOLVED_IN->smpdbpathway.csv', index=False)
    file = 'Protein_(UniProt)_2_Pathway_(SMPDB).csv'
    df.to_csv(os.path.join('output/edges_to_use', file), index=False)


    json.dump(all_prot2pws, open('output/pathway2protein/proteins_to_pathways.json','w'))
    json.dump(all_pw2prots, open('output/pathway2protein/pathway_to_proteins.json', 'w'))

    # Deletes the ~48,000 unzipped SMPDB files
    directory = 'input/smpdb_files'
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        if os.path.isfile(f):
            if(f.endswith('_proteins.csv')):
                os.remove(f)

    prot2smpdb = df.copy(deep=True)
    prots = [prot.split('UniProt:')[1] for prot in set(prot2smpdb['Protein (UniProt)'])][1:]
    json.dump(prots, open('output/protein2protein/protein_to_smpdb_pathway.json','w'))
   
    
if __name__ == '__main__':
    map_protein_to_reactome_pathway()
    map_gene_to_reactome_pathway()
    map_gene_to_kegg_pathway()
    map_protein_to_smpdb_pathway()    
    print('Finished mapping proteins and gees to pathways')