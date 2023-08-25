import os
import pandas as pd
import csv
import requests
import json
from biomedkg_utils import output_edgefile_onerel_noweight


def download_human_proteins_no_prefix():
    # Reviewed human proteome
    url = 'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names&format=json&query=%28reviewed%3Atrue%20AND%20proteome%3Aup000005640%29'
    res = requests.get(url)
    r = res.json()
    protein_set = set([entry['primaryAccession'] for entry in r['results']])
    print(len(protein_set))
    return protein_set


def map_protein_to_reaction():
    protein2reaction, reaction2protein = dict(), dict()
    proteins = set()

    os.system('wget -N -P input/ https://reactome.org/download/current/ProteinRoleReaction.txt')          
    reviewed_proteins = download_human_proteins_no_prefix()        

    file = 'Protein_(UniProt)_2_Reaction_(Reactome).csv'
    with open(os.path.join('output/protein2reaction',file),'w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Protein (UniProt)','Reaction (Reactome)', 'Relationship'])

        for line in open('input/ProteinRoleReaction.txt'):
            line = line.split('\t')

            # Protein, reaction, reaction type
            protein = line[0]
            react_type = line[1]
            reaction = line[2].strip()

            # Human reaction?
            human_reaction = reaction.startswith('R-HSA')
            if human_reaction:

                # Human protein?
                if protein in reviewed_proteins:

                    # Protein - Reaction
                    protein2reaction.setdefault(protein, set()).add(reaction)
                    reaction2protein.setdefault(reaction, set()).add(protein)
                    writer.writerow(['UniProt:'+protein, 'Reactome_Reaction:'+reaction, '-'+react_type+'->'])

        df = pd.read_csv(os.path.join('output/protein2reaction',file)).drop_duplicates()
        df.to_csv(os.path.join('output/edges',file),index=False)
        df.to_csv(os.path.join('output/edges_to_use',file),index=False)

        print(len(protein2reaction), 'Proteins')
        print(len(reaction2protein), 'Reactions')

        
def map_compound_to_reaction():
    os.system('wget -N -P input/ https://reactome.org/download/current/ChEBI2Reactome_PE_Reactions.txt')

    '''Get compounds'''
    chebi2db = json.load(open('output/compound2compound/chebi2db.json'))
    chebi2mesh = json.load(open('output/compound2compound/chebi2mesh.json'))

    mesh_compound2reactome_reaction, reactome_reaction2mesh_compound = dict(), dict()
    db_compound2reactome_reaction, reactome_reaction2db_compound = dict(), dict()
    chebi2reaction, reaction2chebi = dict(), dict()

    compounds = set()

    for line in open('input/ChEBI2Reactome_PE_Reactions.txt'):
        line = line.split('\t')

        # Compound, Reaction
        chebi = line[0]
        reaction = line[3]

        # Human reaction?
        human_reaction = reaction.startswith('R-HSA')
        if human_reaction:

            # ChEBI Compound
            chebi2reaction.setdefault(chebi, set()).add(reaction)
            reaction2chebi.setdefault(reaction, set()).add(chebi)


            # MeSH Compound
            try:
                mesh_compounds = chebi2mesh[chebi]
                for mesh_compound in mesh_compounds:
                    mesh_compound2reactome_reaction.setdefault(mesh_compound, set()).add(reaction)
                    reactome_reaction2mesh_compound.setdefault(reaction, set()).add(mesh_compound)

            except:
                pass


            # DrugBank Compound
            try:
                db_compounds = chebi2db[chebi]
                for db_compound in db_compounds:
                    db_compound2reactome_reaction.setdefault(db_compound, set()).add(reaction)
                    reactome_reaction2db_compound.setdefault(reaction, set()).add(db_compound)

            except:
                pass

    print('DrugBank Compounds:', len(db_compound2reactome_reaction), 
          '\nMeSH Compounds:', len(mesh_compound2reactome_reaction), 
          '\nChEBI Compounds:', len(chebi2reaction))
    print('Reactions to MeSH Compounds:', len(reactome_reaction2mesh_compound), 
          '\nReactions to DrugBank Drugs',len(reactome_reaction2db_compound), 
          '\nReactions to ChEBI:',len(reaction2chebi))
    

    file = 'Compound_(DrugBank)_2_Reaction_(Reactome).csv'
    outpath = os.path.join('output/compound2reaction',file)
    output_edgefile_onerel_noweight(outpath,
                                   ['Compound (DrugBank)', 'Reaction (Reactome)', 'Relationship'],
                                   db_compound2reactome_reaction,
                                    '-participates_in->',
                                   'DrugBank_Compound:',
                                   'Reactome_Reaction:')
    df = pd.read_csv(os.path.join('output/compound2reaction',file)).drop_duplicates()
    df.to_csv(os.path.join('output/edges',file), index=False)
    df.to_csv(os.path.join('output/edges_to_use',file), index=False)
    
    file = 'Compound_(MeSH)_2_Reaction_(Reactome).csv'
    outpath = os.path.join('output/compound2reaction',file)
    output_edgefile_onerel_noweight(outpath,
                                   ['Compound (MeSH)', 'Reaction (Reactome)', 'Relationship'],
                                   mesh_compound2reactome_reaction,
                                    '-participates_in->',
                                   'MeSH_Compound:',
                                   'Reactome_Reaction:')
    df = pd.read_csv(os.path.join('output/compound2reaction',file)).drop_duplicates()
    df.to_csv(os.path.join('output/edges',file), index=False)
    df.to_csv(os.path.join('output/edges_to_use',file), index=False)
    
    
if __name__ == '__main__':
    map_protein_to_reaction()
    map_compound_to_reaction()
    print('Finished mapping protein and compound to reaction')