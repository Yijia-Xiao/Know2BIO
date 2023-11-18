import pandas as pd
from neo4j import GraphDatabase


def clear_graph(session, debug=True):
    """
    This function deletes all nodes and edges in the graph.
    """
    
    query = """
    MATCH (n)
    DETACH DELETE n
    """
    session.run(query)
    
    if debug:
        count_nodes_and_relationships(session, debug=True)
    

def count_nodes_and_relationships(session, debug=True):
    """
    This function counts the number of nodes and edges in the graph.
    """
    # Count nodes
    node_query = "MATCH (n) RETURN COUNT(n) AS nodeCount"
    nodes = session.run(node_query).single()["nodeCount"]

    # Count relationships
    relationship_query = "MATCH ()-->() RETURN COUNT(*) AS relationshipCount"
    rels = session.run(relationship_query).single()["relationshipCount"]
    
    if debug:
        print(f"Number of nodes: {nodes}")
        print(f"Number of relationships: {rels}")

    return nodes, rels


def create_node(session, node_name, node_type):
    """
    This function creates a new node with node_name and type node_type.
    """
    query = """
    CREATE (n:%s {name:"%s"})
    RETURN n
    """ % (node_type, ":".join([node_type, node_name]))
    
    session.run(query)
    
    
def create_relationship(session, entity_1, entity_2, relation):
    """
    This function creates a new relationship between entity_1 and entity_2.
    """
    e1_type = entity_1.split(":")[0]
    e2_type = entity_2.split(":")[0]
    
    query = """
    MATCH (e1:%s {name: "%s"})
    MATCH (e2:%s {name: "%s"})
    WHERE e1 IS NOT NULL AND e2 IS NOT NULL
    MERGE (e1)-[:`%s`]->(e2)
    """ % (e1_type, entity_1, e2_type, entity_2, relation)
    result = session.run(query)
    
    # Consume the result to get the query summary
    summary = result.consume()
    
    # Check if the relationship was created
    if summary.counters.relationships_created == 0:
        print("Warning: One or both persons do not exist in the graph.")


def load_kg(kg_file):
    df = pd.read_csv(kg_file, sep="\t", header=None)
    df.columns = ['h', 'r', 't']

    nodes = set(df['h']).union(set(df['t']))
    uniq_edges = set(df['r'])
    relationships = list(zip(df['h'], df['r'], df['t']))

    print("%d nodes and %d edges (%d unique edges)" % (len(nodes), len(relationships), len(uniq_edges)))
    return df, nodes, relationships


# Class handling the database connection
class Neo4jConnector:
    def __init__(self, uri, user, password):
        self._driver = GraphDatabase.driver(uri, auth=(user, password))

    def close(self):
        self._driver.close()        

# Replace with your own credentials
NEO4J_URI = "bolt://localhost:7687"
NEO4J_USER = "neo4j"
NEO4J_PASSWORD = "password"

# Connect to Neo4j Server
print("Connecting to Neo4j Server")
connector = Neo4jConnector(NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD)

# Load Know2BIO Data
print("Loading Know2BIO edges")
input_file = '../know2bio/whole_kg.txt'
# input_file = '../know2bio/whole_kg_sampled.txt' # Smaller dataset for testing purposes
_, nodes, relationships = load_kg(input_file)

with connector._driver.session() as session:
    # Clear graph
    print("Clearing Graph")
    clear_graph(session, debug=False)
    
    # Create nodes
    print("Creating nodes")
    count = 0
    for node in nodes:
        node_type = node.split(":")[0]
        node_name = ":".join(node.split(":")[1:])
        create_node(session, node_name, node_type)
        count += 1
        print(f'Nodes completed: %d out of %d' % (count,len(nodes)), end='\r', flush=True)
    print("%d nodes created.                        " % count)
    
    # Create relationships
    print("Creating relationships")
    count = 0
    for h, r, t in relationships:
        create_relationship(session, h, t, r)
        count += 1
        print(f'Relationships completed: %d out of %d' % (count, len(relationships)), end='\r', flush=True)
    print("%d relationships created.                " % count)
    
    # Print summary
    node_count, relationship_count = count_nodes_and_relationships(session)
