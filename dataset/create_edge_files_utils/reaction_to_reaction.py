import json
import csv
from joblib import Parallel, delayed
from multiprocessing import cpu_count
from gene_to_anatomy import batch_iterator
import requests
import os
import pandas as pd

def map_reaction_to_preceding_reaction(reactions):
    reaction_to_reaction = {}
    
    url = 'https://reactome.org/ContentService/data/query/ids/map'
    headers = {'accept': 'application/json', 
               'Content-Type': 'text/plain'}
    BATCH_SIZE = 20

    for idx, id_batch in enumerate(batch_iterator(reactions, BATCH_SIZE)):

        # API
        reaction_ids = ','.join(id_batch)
        response = requests.post(url, 
                                 headers=headers, 
                                 data=reaction_ids)
        reaction_dict = response.json()

        # Map API results
        reaction_ids = reaction_ids.split(',')
        for reaction_id in reaction_ids:
            try:
                preceding_events = reaction_dict[reaction_id]['precedingEvent']
            except:
                continue
            for preceding_event in preceding_events:
                if 'reaction' in preceding_event['className'].lower():
                    preceding_reaction = preceding_event['dbId']
                    reaction_to_reaction.setdefault(reaction_id, 
                                                    set()).add(preceding_reaction)

    return reaction_to_reaction
        
        
def output_reaction_to_reaction(reaction_to_reaction):
    file = 'Reaction_(Reactome)_to_Reaction_(Reactome).csv'
    with open(os.path.join('output/reaction2reaction/',file), 'w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Reaction (Reactome)', 'Reaction (Reactome)', 'Relationship'])
    
        for reaction, preceding_reactions in reaction_to_reaction.items():
            for preceding_reaction in preceding_reactions:
                row = ['Reactome_Reaction:'+reaction, 
                      'Reactome_Reaction:R-HSA-'+str(preceding_reaction),
                      '-head_reaction_precedes_tail_reaction->']
                writer.writerow(row)

    reaction_to_reaction_df = pd.read_csv(os.path.join('output/reaction2reaction/',file))
    reaction_to_reaction_df.to_csv(os.path.join('output/edges',file), index=False)
    reaction_to_reaction_df.to_csv(os.path.join('output/edges_to_use',file),
                                   index=False)
        
        
def get_reaction_to_reaction_dict():
    BATCH_SIZE = cpu_count()
    reactions = json.load(open('input/reactome_reactions.json'))
    reaction_batches = [reactions[i:i+BATCH_SIZE] for i in range(0, len(reactions), BATCH_SIZE)]
    reaction_to_reaction_dicts = Parallel(n_jobs=-1)(delayed(map_reaction_to_preceding_reaction)(reaction_batch) 
                                                     for reaction_batch in reaction_batches)        
    reaction_to_reaction = {}
    for reaction_to_reaction_batch in reaction_to_reaction_dicts:
        reaction_to_reaction.update(reaction_to_reaction_batch)
        
    return reaction_to_reaction


if __name__ == '__main__':
    reaction_to_reaction = get_reaction_to_reaction_dict()
    output_reaction_to_reaction(reaction_to_reaction)