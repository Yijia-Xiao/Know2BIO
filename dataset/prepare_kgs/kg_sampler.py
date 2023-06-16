import sys
from grape import Graph
import os
import pandas as pd

def read_df(input_file):
	# read as pandas df
	df = pd.read_csv(input_file,sep="\t")
	if df.shape[1] == 2:
		# some of the files didn't have index column
		df = pd.read_csv(input_file)
	# rename header to hrt weight convetion
	if df.shape[1] == 3:
		df.columns = ['h','r','t']
	elif df.shape[1] == 4:
		df.columns = ['h','r','t','weight']
	else:
		print("Error!",input_file)
		print(df.shape)
		print(df.head())	
	# keep track of node type, indicated by header
	node_to_type = {}
	nodes = list(set(df['h']).union(set(df['t'])))
	for n in nodes:
		if ":" in n:
			t = n.split(":")[0]
		else:
			t = "Undefined"
		node_to_type[n] = t
	nodes_df = pd.DataFrame({'node':list(node_to_type.keys()),
							 'node_type':list(node_to_type.values())})	
	return df, nodes_df


# repurposed train_test_split to sample 1% of the kg
def train_test_split(edges_file,nodes_file,training_size=0.99, include_headers=False, include_weights=False):
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
	  node_list_separator='\t',
	  edge_list_separator='\t',
	  name="mykg"
	)
	g = g.remove_disconnected_nodes()

	# split train and test, write to file
	train, test = g.connected_holdout(train_size=training_size)
	#train_edges = (edges_file.replace(".txt","")+"_train.txt").replace("grape_","")
	sampled_edges = (edges_file.replace(".txt","")+"_sampled.txt").replace("grape_","")
#	train.dump_edges(train_edges, separator="\t",
#					sources_column='h',sources_column_number=0,
#					edge_type_column='r',edge_types_column_number=1,
#					destinations_column='t',destinations_column_number=2,
#					weights_column='weight',weights_column_number=3)
	test.dump_edges(sampled_edges, separator="\t",
					sources_column='h',sources_column_number=0,
					edge_type_column='r',edge_types_column_number=1,
					destinations_column='t',destinations_column_number=2)
					#weights_column='weight',weights_column_number=3)

	# report results
	full_size = (g.get_number_of_nodes(),g.get_number_of_edges())
	train_size = (train.get_number_of_nodes(),train.get_number_of_edges())
	test_size = (test.get_number_of_nodes(),test.get_number_of_edges())
	
	# formatting
	for outname in [sampled_edges]:
		df = pd.read_csv(outname,header=0, sep="\t")
		print(df.shape,outname)
		df.to_csv(outname, index=False, sep="\t", header=include_headers)
	return (full_size,train_size,test_size)	


### Parse arguments ###
if len(sys.argv) < 3:
	print("Incorrect usage.")
	print("python kg_sampler.py <kg_file> <output_name>")
	print("Example: python kg_sampler.py ../dgs_input/kg1f_instances_test.txt ../dgs_input/sampled_kg1f_instances_test.txt")
	sys.exit(0)

input_kg_name = sys.argv[1]
output_kg_name = sys.argv[2]

print(input_kg_name,output_kg_name)

# read kg
edges_df, nodes_df = read_df(input_kg_name)

# export for grape format
# write edges file for grape
grape_edges_name = input_kg_name[:-4]+"_grape.txt"
edges_df.to_csv(grape_edges_name, index=False, sep="\t")	
# write nodes file for grape
nodes_outfile_name = input_kg_name.replace(".txt","")+"_nodes.txt"
nodes_df.to_csv(nodes_outfile_name,index=False,sep="\t")

# sample from KG, ensuring connected components 
split_sizes = train_test_split(grape_edges_name,nodes_outfile_name)
original_size, _, sampled_size = split_sizes
print("Original KG size: ",original_size)
print("Sampled KG size: ",sampled_size)
