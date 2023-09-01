import sys
import os
import pandas as pd


def read_df(input_file):
    # read as pandas df
    df = pd.read_csv(input_file)
    # keep track of node type, indicated by header
    node_to_type = {}
    node_types = list(df.columns)[0:2]
    for node in df[node_types[0]]:
        node_to_type[node] = node_types[0]
    for node in df[node_types[1]]:
        node_to_type[node] = node_types[1].replace(".1","") # tag added to prevent duplicate names
    # rename header to hrt weight convetion
    if df.shape[1] == 3:
        df.columns = ['h','t','r']
    elif df.shape[1] == 4:
        df.columns = ['h','t','r','weight']
    else:
        print("Error!",input_file)
        print(df.shape)
        print(df.head())    

    return df, node_to_type


def read_and_merge_files(input_directory, files_to_merge, default_weight=1.0, rename_ontology=False):
    merged_df = pd.DataFrame({'h': pd.Series(dtype='str'),
                   'r': pd.Series(dtype='str'),
                   't': pd.Series(dtype='str'),
                   'weight': pd.Series(dtype='float')})
    node_to_types = {}    

    # read and merge files
    for f in files_to_merge:
        file_path = os.path.join(input_directory,f)
        df,node_to_type = read_df(file_path)
        merged_df = pd.concat([merged_df,df],sort=False)
        node_to_types = {**node_to_types, **node_to_type}

    # fill NA
    edges_df =  merged_df.fillna(value={'weight':default_weight})
    # filter out 0 weights
    edges_df = edges_df[edges_df['weight'] != 0]

    # convert nodes_to_types to df
    nodes_df = pd.DataFrame({'node':list(node_to_types.keys()),
                             'node_type':list(node_to_types.values())})    

    # drop duplicate nodes
    edges_df = edges_df.drop_duplicates()

    return edges_df, nodes_df


def parse_input_list(file_name):
    lines = [l.strip("\n") for l in open(file_name,"r").readlines()]
    name_to_files = {}
    name = ""
    for l in lines:
        if "#" in l: # header
            s = l.split(" ")
            name = s[1]
            name_to_files[name] = []
        elif len(l) > 0: # file
            name_to_files[name] += [l]
    return name_to_files['ontologies'], name_to_files['bridges'], name_to_files['instances']


### Parse arguments ###
if len(sys.argv) < 4:
    print("Incorrect usage.")
    print("python prepare_kgs.py <input_directory> <output_directory> <input_lists>")
    print("Example: python prepare_kgs.py Know2BIO_files Know2BIO_kgs ont_bridge_inst_lists.txt")
    sys.exit(0)

# Parse arguments
input_directory = sys.argv[1]
output_directory = sys.argv[2]
input_lists = sys.argv[3]

# Make directory if not exist
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Output file paths
output_file_names = ['kg2f_ontologies.txt','alignf_bridges.txt','kg1f_instances.txt']
output_file_names = [os.path.join(output_directory,f) for f in output_file_names]

# Parse file lists
ontologies, bridges, instances = parse_input_list(input_lists)

# Change flags for different formatting.
include_headers=False
include_weights=False

# Read files for each kg view
ontologies_edges_nodes = read_and_merge_files(input_directory, ontologies)
bridges_edges_nodes = read_and_merge_files(input_directory, bridges)
instances_edges_nodes = read_and_merge_files(input_directory, instances)

dfs = []
# Save files
for out_file_name, edges_nodes in zip(output_file_names,[ontologies_edges_nodes,bridges_edges_nodes,instances_edges_nodes]):
    edges_df, _ = edges_nodes
    
    # write edges file
    if include_weights:
        out_edges_df = edges_df
    else:
        out_edges_df = edges_df[['h','r','t']]
    out_edges_df.to_csv(out_file_name, index=False, sep="\t", header=include_headers)

    print(out_file_name, out_edges_df.shape)

    dfs += [out_edges_df]

whole_df_outfile = os.path.join(output_directory,'whole_kg.txt')
whole_df = pd.concat(dfs)  
whole_df = whole_df.drop_duplicates()
whole_df.to_csv(whole_df_outfile, index=False, sep="\t", header=include_headers)
print(whole_df_outfile, whole_df.shape)
