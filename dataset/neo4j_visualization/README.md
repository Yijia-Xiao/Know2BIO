This folder contains the software code to query the constructed knowledge graph.

Steps:
1. Download Neo4j Desktop client from https://neo4j.com/
2. Create a local DBMS and start the service.
3. Edit create_know2bio_neo4j.py with login credentials (default values are included)
4. Run `python create_know2bio_neo4j.py`
5. Start up Neo4j Browser from the Desktop client and enter Cypher query. Refer to Neo4j tutorial for more help: https://neo4j.com/developer/get-started/

![An example of a Neo4j query, once the knowledge graph is loaded.](./query_example.png)