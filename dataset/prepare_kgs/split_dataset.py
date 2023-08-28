import sys
from grape import Graph
import os
import pandas as pd
import numpy as np

def read_df(file_name, overwrite=True):
    df = pd.read_csv(file_name,sep='\t',header=None)
    df.columns = ['h','r','t']
    df = df.drop_duplicates()
    
    if overwrite:
        df.to_csv(file_name,sep='\t',header=False,index=False)
		
    return df


def export_grape_files(file_name):
    grape_edges_name = file_name[:-4]+"_grape.txt"
    grape_nodes_name = file_name[:-4]+"_nodes.txt"

    df = read_df(file_name)

    # output grape edges
    df['weight'] = 1.0
    df.to_csv(grape_edges_name,index=False, sep='\t')

    # output grape nodes
    nodes = list(set(df['h']).union(set(df['t'])))

    nodes_df = pd.DataFrame({'node': nodes,
        'node_type': [n.split(":")[0] for n in nodes]})
    nodes_df.to_csv(grape_nodes_name,index=False, sep='\t')

    return grape_edges_name, grape_nodes_name


def get_subgraphs(g, thresh=10, debug=True):
    '''
    Returns subgraphs with over thresh nodes in the subgraph
    '''

    # Calculate connected components
    ccs, n_comps, comp_min, comp_max = g.get_connected_components()
    
    print("%d connected components"%n_comps)

    # Gather nodes in each subgraph
    comp_to_nodes = {}
    for idx, group in enumerate(list(ccs)):
        if group in comp_to_nodes:
            comp_to_nodes[group].append(idx)
        else:
            comp_to_nodes[group] = [idx]

    # Extract subgraph larger than thresh size
    sub_gs = []
    for idx, nodes in comp_to_nodes.items():
        if len(nodes) > thresh:
            sub_g = g.filter_from_ids(nodes)
            sub_gs += [sub_g]
    
    if debug:
        print("Original graph: %d nodes %d edges"%(g.get_number_of_nodes(), g.get_number_of_edges()))
        print("%d graphs with over %d nodes"%(len(sub_gs),thresh))
        if len(sub_gs) < 10:
            for idx, sub_g in enumerate(sub_gs):
                print("SCC %d: %d nodes %d edges"%(idx, sub_g.get_number_of_nodes(), sub_g.get_number_of_edges()))
    return sub_gs    

def train_test_split(edges_file,nodes_file,training_size=0.8, include_headers=False, include_weights=False, debug=True):
    # load graph as grape graph
    g = Graph.from_csv(
      directed=False, 
      node_path=nodes_file,
      edge_path=edges_file,
      verbose=True,
      nodes_column='node',
      node_list_node_types_column='node_type',
      default_node_type='None',
      sources_column='h',
      destinations_column='t',
      edge_list_edge_types_column='r',
      weights_column='weight',
      node_list_separator='\t',
      edge_list_separator='\t',
      name="mykg"
    )

    components = get_subgraphs(g)
    train_files = []
    test_files = []
    for idx, sub_g in enumerate(components):
        assert sub_g.is_connected()
    
        try:
            train, test = sub_g.connected_holdout(train_size=training_size)    
            train_edges = (edges_file.replace(".txt","_train_scc_%d.txt"%(idx)))
            test_edges = (edges_file.replace(".txt","_test_scc_%d.txt"%(idx)))
            train.dump_edges(train_edges, separator="\t",
                        sources_column='h',sources_column_number=0,
                        edge_type_column='r',edge_types_column_number=1,
                        destinations_column='t',destinations_column_number=2)
            test.dump_edges(test_edges, separator="\t",
                        sources_column='h',sources_column_number=0,
                        edge_type_column='r',edge_types_column_number=1,
                        destinations_column='t',destinations_column_number=2)
            train_files += [train_edges]
            test_files += [test_edges]
            if debug:
                print("SCC %d: %d nodes %d edges. Train: %d nodes, %d edges; Test: %d nodes, %d edges"%(idx,
                             sub_g.get_number_of_nodes(), sub_g.get_number_of_edges(),
                             train.get_number_of_nodes(), train.get_number_of_edges(),
                             test.get_number_of_nodes(), test.get_number_of_edges()))
        except Exception as error:
            #print(error)
            rate = float(str(error).split("train rate of ")[1][:-1])
            if debug:
                print("SCC %d was unsuccessful. Retrying with new train rate %d"%(idx,rate))
            train_edges = (edges_file.replace(".txt","_train_scc_%d.txt"%(idx)))
            test_edges = (edges_file.replace(".txt","_test_scc_%d.txt"%(idx)))
            if rate == 1.: # Don't have to split, it's all training set
                sub_g.dump_edges(train_edges, separator="\t", 
                        sources_column='h',sources_column_number=0,
                        edge_type_column='r',edge_types_column_number=1,
                        destinations_column='t',destinations_column_number=2)
                train_files += [train_edges]
                if debug:
                    print("SCC %d: %d nodes %d edges. Added to train set"%(idx,sub_g.get_number_of_nodes(), sub_g.get_number_of_edges()))
            else:
                train, test = sub_g.connected_holdout(train_size=rate)    
                train.dump_edges(train_edges, separator="\t",
                            sources_column='h',sources_column_number=0,
                            edge_type_column='r',edge_types_column_number=1,
                            destinations_column='t',destinations_column_number=2)
                test.dump_edges(test_edges, separator="\t",
                            sources_column='h',sources_column_number=0,
                            edge_type_column='r',edge_types_column_number=1,
                            destinations_column='t',destinations_column_number=2)
                train_files += [train_edges]
                test_files += [test_edges]
                if debug:
                    print("SCC %d: %d nodes %d edges. Train: %d nodes, %d edges; Test: %d nodes, %d edges"%(idx,
                                 sub_g.get_number_of_nodes(), sub_g.get_number_of_edges(),
                                 train.get_number_of_nodes(), train.get_number_of_edges(),
                                 test.get_number_of_nodes(), test.get_number_of_edges()))

        
    # Combine files and split validation from test set
    train_dfs = [pd.read_csv(f,sep='\t') for f in train_files]
    train_df = pd.concat(train_dfs, ignore_index=True)

    test_dfs = [pd.read_csv(f,sep='\t') for f in test_files]
    test_valid_df = pd.concat(test_dfs, ignore_index=True)
    test_df, valid_df = split_df(test_valid_df)

    # Save files
    train_out_file = edges_file.replace("_grape.txt","_train.txt")
    test_out_file = edges_file.replace("_grape.txt","_test.txt")
    valid_out_file = edges_file.replace("_grape.txt","_valid.txt")
    train_df.to_csv(train_out_file,index=False,header=False)
    test_df.to_csv(test_out_file,index=False,header=False)
    valid_df.to_csv(valid_out_file,index=False,header=False)
 
    print("Original graph: %d nodes %d edges"%(g.get_number_of_nodes(),g.get_number_of_edges()))
    train_nodes,train_edges = get_nodes_edges_from_df(train_df)
    test_nodes,test_edges = get_nodes_edges_from_df(test_df)
    valid_nodes,valid_edges = get_nodes_edges_from_df(valid_df)
    print("Train graph: %d nodes %d edges"%(len(train_nodes),len(train_edges)))
    print("Test graph: %d nodes %d edges"%(len(test_nodes),len(test_edges)))
    print("Valid graph: %d nodes %d edges"%(len(valid_nodes),len(valid_edges)))

    print("Checking if all nodes appear in train")
    assert len(test_nodes.difference(train_nodes)) == 0
    assert len(valid_nodes.difference(train_nodes)) == 0
    print(True)

def get_nodes_edges_from_df(df):
    nodes = set(df['h']).union(set(df['t']))
    edges = df['r']
    return nodes,edges

def split_df(df):
    shuffled_indices = np.random.permutation(df.index)
    split_point = len(shuffled_indices) // 2
    first_half_indices = shuffled_indices[:split_point]
    second_half_indices = shuffled_indices[split_point:]

    first_half_df = df.loc[first_half_indices]
    second_half_df = df.loc[second_half_indices]

    return first_half_df, second_half_df

in_files = ['kg1f_instances.txt','kg2f_ontologies.txt','alignf_bridges.txt']
#in_file = 'kg1f_instances.txt'
#in_file = 'kg2f_ontologies.txt'
#in_file = 'alignf_bridges.txt'
for in_file in in_files:
    print(in_file)
    edges_file, nodes_file = export_grape_files(in_file)
    #edges_file = in_file[:-4]+"_grape.txt"
    #nodes_file = in_file[:-4]+"_nodes.txt"
    train_test_split(edges_file, nodes_file,debug=False)





