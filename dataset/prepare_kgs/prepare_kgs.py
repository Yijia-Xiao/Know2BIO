import sys
from grape import Graph
import os
import pandas as pd

def read_df(input_file):
	# read as pandas df
	df = pd.read_csv(input_file,index_col=0)
	if df.shape[1] == 2:
		# some of the files didn't have index column
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


def rename_ontology_df(df):
	print("Renaming ontology")
	# rename the key relation to (new_relation, flip)
	# if flip is True, will flip the relation
	remap_relations_flip = {'-parent_of->':('isa',True), 
							'-part_of->':('isa',False), 
							'is_subclass_of->':('isa',False), 
							'-subclass_of->':('isa',False),
							'-is_a->':('isa',False)}
	# keep track of unused relations to append at end
	unused_rels = set(df['r'])
	# keep track of relabeled rels
	ret_df = None
	# loop through each rel in remap function
	for rel, rf in remap_relations_flip.items():
		new_rel, flip = rf
		
		# get just the edges w/ rel
		sub_df = df[df['r'] == rel]
		
		if flip: # flip head and tail position, then rename
			dat = [list(sub_df['t']), 
				   [new_rel for i in range(sub_df.shape[0])],
				   list(sub_df['h'])]
			new_df = pd.DataFrame(dat).T
			new_df.columns=['h','r','t']
			if 'weight' in sub_df.columns:
				new_df['weight'] = list(sub_df['weight'])
		else: # only rename
			new_df = sub_df.copy(deep=True)
			new_df['r'] = new_rel
		# append to ret_Df
		if (ret_df is None):
			ret_df = new_df
		else:
			ret_df = pd.concat([ret_df,new_df])
		#	ret_df = ret_df.append(new_df,ignore_index=True)
		unused_rels.remove(rel)
	# append unused rels to the dataframe
	for unused_rel in unused_rels:
		sub_df = df[df['r'] == unused_rel]
		#ret_df = ret_df.append(sub_df,ignore_index=True)
		ret_df = pd.concat([ret_df,sub_df])
	# error checking
	if df.shape != ret_df.shape:
		print("New dataframe is not the same shape!!!")
		print(df.shape)
		print(ret_df.shape)
	return ret_df


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

	# rename ontology
	if rename_ontology:
		edges_df = rename_ontology_df(edges_df)

	# convert nodes_to_types to df
	nodes_df = pd.DataFrame({'node':list(node_to_types.keys()),
							 'node_type':list(node_to_types.values())})	

	return edges_df, nodes_df


def train_test_split(edges_file,nodes_file,training_size=0.7, include_headers=False, include_weights=False):
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
	g = g.remove_disconnected_nodes()

	# split train and test, write to file
	train, test = g.connected_holdout(train_size=training_size)
	train_edges = (edges_file.replace(".txt","")+"_train.txt").replace("grape_","")
	test_edges = (edges_file.replace(".txt","")+"_test.txt").replace("grape_","")
	train.dump_edges(train_edges, separator="\t",
					sources_column='h',sources_column_number=0,
					edge_type_column='r',edge_types_column_number=1,
					destinations_column='t',destinations_column_number=2,
					weights_column='weight',weights_column_number=3)
	test.dump_edges(test_edges, separator="\t",
					sources_column='h',sources_column_number=0,
					edge_type_column='r',edge_types_column_number=1,
					destinations_column='t',destinations_column_number=2,
					weights_column='weight',weights_column_number=3)

	# report results
	full_size = (g.get_number_of_nodes(),g.get_number_of_edges())
	train_size = (train.get_number_of_nodes(),train.get_number_of_edges())
	test_size = (test.get_number_of_nodes(),test.get_number_of_edges())
	
	# formatting
	for outname in [train_edges,test_edges]:
		df = pd.read_csv(outname,header=0, sep="\t")
		if not include_weights:
			df = df[['h','r','t']]
		print(df.shape,outname)
		df.to_csv(outname, index=False, sep="\t", header=include_headers)
	return (full_size,train_size,test_size)	


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

input_directory = sys.argv[1]
output_directory = sys.argv[2]
input_lists = sys.argv[3]

output_file_names = ['kg2f_ontologies.txt','alignf_bridges.txt','kg1f_instances.txt']
output_file_names = [os.path.join(output_directory,f) for f in output_file_names]
ontologies, bridges, instances = parse_input_list(input_lists)

include_headers=False
include_weights=False

ontologies_edges_nodes = read_and_merge_files(input_directory, ontologies, rename_ontology=True)
bridges_edges_nodes = read_and_merge_files(input_directory, bridges)
instances_edges_nodes = read_and_merge_files(input_directory, instances)

summary_out = []
for out_file_name, edges_nodes in zip(output_file_names,[ontologies_edges_nodes,bridges_edges_nodes,instances_edges_nodes]):
	edges_df, nodes_df = edges_nodes
	
	# write edges file
	if include_weights:
		out_edges_df = edges_df
	else:
		out_edges_df = edges_df[['h','r','t']]
	out_edges_df.to_csv(out_file_name, index=False, sep="\t", header=include_headers)
	print(out_file_name, edges_df.shape)

	# write edges file for grape
	grape_edges_name = out_file_name[:-4]+"_grape.txt"
	edges_df.to_csv(grape_edges_name, index=False, sep="\t")	
	# write nodes file for grape
	nodes_outfile_name = out_file_name.replace(".txt","")+"_nodes.txt"
	nodes_df.to_csv(nodes_outfile_name,index=False,sep="\t")
	print(nodes_outfile_name,nodes_df.shape)	
	
	# train-test split
	split_sizes = train_test_split(grape_edges_name,nodes_outfile_name)
	summary_out += [(out_file_name,split_sizes)]

print("### Summary ###")
for out_file_name, split_sizes in summary_out:
	print(out_file_name,split_sizes)
