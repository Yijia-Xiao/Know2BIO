import pandas as pd
import math
import os
import csv

# INSTRUCTIONS (Download the Reactome Neo4j graph):
# Download the Reactome database as a Neo4j graph: https://reactome.org/download-data/ is the main page
# These are the possible files: https://reactome.org/download/current/reactome.graphdb.tgz 
# or https://reactome.org/download/current/reactome.graphdb.dump\n",
'''
Neo4j Command:
    MATCH p=(pw:Pathway)-[edge:hasEvent]-(r:ReactionLikeEvent)
    WHERE toLower(pw.speciesName) = 'homo sapiens'
    WITH collect(pw) as pathways, collect(r) as RLEs , collect(edge) as relation
    CALL apoc.export.csv.data(pathways + RLEs, relation, 'human_pathway_reaction.csv', {})
    YIELD file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data
    RETURN file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data
'''


def map_reaction_to_pathway():
    nodes = pd.read_csv('input/human_pathway_reaction.csv', low_memory=False)[['_id', 'stId','displayName', 'name']]
    edges = pd.read_excel('input/human_pathways_reaction-like-event_edges.xlsx')
    db2reactome_id = dict(zip(nodes._id, nodes.stId))
    db2reactome_id = {int(db):reactome_id for db, reactome_id in db2reactome_id.items() if not math.isnan(db)}
    #db2reactome_id[931828]

    
    # Edges
    edge_list = list()

    for i in range(0,len(edges)):
        pathway = edges['_start'].iloc[i]
        reaction = edges['_end'].iloc[i]
        rel = edges['_type'].iloc[i]

        try:
            pathway = db2reactome_id[pathway]
            reaction = db2reactome_id[reaction]        
            edge_list.append([pathway, reaction])

        except:
            continue
            
    # Output Edges
    file = 'Reaction_(Reactome)_2_Pathway_(Reactome).csv'
    with open(os.path.join('output/pathway2reaction',file),'w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Reaction (Reactome)','Pathway (Reactome)','Relationship'])
        relationship = '-involved_in->'


        for edge in edge_list:
            pathway = edge[0]
            reaction = edge[1]
            writer.writerow(['Reactome_Reaction:'+reaction, 'Reactome_Pathway:'+pathway, relationship])

    df = pd.read_csv(os.path.join('output/pathway2reaction',file))
    df.to_csv(os.path.join('output/edges/', file), index=False)
    df.to_csv(os.path.join('output/edges_to_use/',file), index=False)
    

if __name__ == '__main__':
    map_reaction_to_pathway()
    print('Finished mapping reaction to pathway')