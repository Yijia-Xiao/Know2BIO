'''Note: Currently, this requires you to use Neo4j. Efforts to automate
the process of creating edges (via Reactome API) are underway.'''

import pandas as pd
import csv
import os
import math

def map_reaction_like_events_to_reaction_like_events():
    # Note that Reactions are a proper subset of Reaction-Like-Events
    nodes = pd.read_csv('input/human_reaction-like_events.csv', low_memory=False)[['_id', 'stId','displayName', 'name', 'category']].drop_duplicates()
    edges = pd.read_excel('input/human_reaction-like_event_edges.xlsx').drop_duplicates()

    db2stID = dict(zip(nodes._id, nodes.stId))
    db2reactome_id = {int(dbId):reactome_id for dbId,reactome_id in db2stID.items() if not math.isnan(dbId)}

    edge_list = list()

    for i in range(0,len(edges)):
        start = edges['_start'].iloc[i]
        end = edges['_end'].iloc[i]
        rel = edges['_type'].iloc[i]

        try:
            start = db2reactome_id[start]
            end = db2reactome_id[end]        
            edge_list.append([start, end])

        except:
            continue

    file = 'Reaction_(Reactome)_precedes_Reaction_(Reactome).csv'
    with open(os.path.join('output/reaction2reaction/',file), 'w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Reaction (Reactome)','Reaction (Reactome)','Relationship'])
        relationship = '-head_happens_before_that_tail_reaction->'

        for edge in edge_list:
            first_event = edge[1]
            second_event = edge[0]
            writer.writerow(['Reactome_Reaction:'+first_event, 'Reactome_Reaction:'+second_event, relationship])

    df = pd.read_csv(os.path.join('output/reaction2reaction/',file))
    df.to_csv(os.path.join('output/edges/',file), index=False)
    df.to_csv(os.path.join('output/edges_to_use/',file), index=False)
    
if __name__ == '__main__':
    '''
    Neo4j Query:
    MATCH (r1:ReactionLikeEvent)-[edge:precedingEvent]-(r2:ReactionLikeEvent)
    WHERE toLower(r1.speciesName) = 'homo sapiens' and toLower(r2.speciesName) = 'homo sapiens'
    WITH collect(r1) as reactions1, collect(r2) as reactions2, collect(edge) as relation
    CALL apoc.export.csv.data(reactions1 + reactions2, relation, 'human_reaction-like_events.csv', {})
    YIELD file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data
    RETURN file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data
    '''
    map_reaction_like_events_to_reaction_like_events()