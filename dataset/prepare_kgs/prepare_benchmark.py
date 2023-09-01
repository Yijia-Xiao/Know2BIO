import pandas as pd
import os
import sys

def is_valid_file(file_path):
    return os.path.exists(file_path) and os.path.isfile(file_path)


def get_input_files(directory, debug=False):
    '''
    This function returns a dictionary mapping of each base_filename to the dataset_name files.
    '''

    dataset_names = ["all", "train", "test", "valid", "all_sampled", "train_sampled", "test_sampled", "valid_sampled"]
    base_filenames = ["kg1f_instances","kg2f_ontologies","alignf_bridges","whole_kg"]    

    input_files = {}
    
    # Get all possible files for each dataset
    for base_filename in base_filenames:
        dataset_files = {}

        for dataset_name in dataset_names:
            if dataset_name == "all":
                filename = f"{base_filename}.txt"
            elif dataset_name == "all_sampled":
                filename = f"{base_filename}_sampled.txt"
            else:
                filename = f"{base_filename}_{dataset_name}.txt"
            
            file_path = os.path.join(directory, filename)
            if debug:
                print(file_path)    
            if is_valid_file(file_path):
                dataset_files[dataset_name] = file_path
                if debug:
                    print(f"{filename} is a valid file for {dataset_name} dataset.")
            else:
                if debug:
                    print(f"{filename} is not a valid file for {dataset_name} dataset.")

        # Make sure all datasets have base, train, test, valid. Must have 4 for all, sampled, or both.
        # Check if all required datasets are present for this base_filename
        required_datasets = ["all", "train", "test", "valid"]
        required_sampled_datasets = ["all_sampled", "train_sampled", "test_sampled", "valid_sampled"]
        if all(dataset in dataset_files for dataset in required_datasets) or \
           all(dataset in dataset_files or dataset.replace("_sampled", "") in dataset_files for dataset in required_sampled_datasets):
            input_files[base_filename] = dataset_files
        else:
            missing_datasets = required_datasets + [ds for ds in required_sampled_datasets if ds not in dataset_files]
            print(f"Missing required datasets for {base_filename}: {missing_datasets}")
    
    return input_files


def read_kg(kg_file):
    '''
    Reads and returns kg_file as Pandas DataFrame
    '''
    df = pd.read_csv(kg_file,sep='\t',header=None)
    df.columns = ['h','r','t']
    return df


def create_entity_relation_mappings(kg_file, directory=None, debug=False):
    '''
    This function prepares an entity2id and relation2id mapping.
    Outputs to provided directory provided in the below format.
        - # lines in first line
        - # name\tid
    '''

    kg_df = read_kg(kg_file)    
    nodes = list(set(kg_df['h']).union(set(kg_df['t'])))
    relations = list(set(kg_df['r']))

    entity2id = {n:idx for idx,n in enumerate(nodes)}
    relation2id = {r:idx for idx,r in enumerate(relations)}

    if debug:
        print("%d entities and %d relations parsed."%(len(entity2id),len(relation2id)))

    if directory:
        entity_out_file = os.path.join(directory,'entity2id.txt')
        with open(entity_out_file,'w') as outfile:
            outfile.write(str(len(entity2id))+"\n")
            outfile.write("\n".join(["\t".join([ent_name,str(ent_id)]) for ent_name, ent_id in entity2id.items()]))
        if debug:
            print("%d IDs written to %s"%(len(entity2id),entity_out_file))

        rel_out_file = os.path.join(directory,'relation2id.txt')
        with open(rel_out_file,'w') as outfile:
            outfile.write(str(len(relation2id))+"\n")
            outfile.write("\n".join(["\t".join([rel_name,str(rel_id)]) for rel_name, rel_id in relation2id.items()]))
        if debug:
            print("%d IDs written to %s"%(len(relation2id),rel_out_file))
    return entity2id, relation2id


def reformat_kg(kg_file, entity2id, relation2id, out_file=None, debug=False):

    # read kg
    kg_df = read_kg(kg_file)    
    if debug:
        print("%s: %s"%(kg_file,str(kg_df.shape)))

    # remap to IDs
    remapped_df = pd.DataFrame({'h':[entity2id[e] for e in kg_df['h']],
                                'r':[relation2id[r] for r in kg_df['r']],
                                't':[entity2id[e] for e in kg_df['t']]
                                })
    # correct format for benchmark data is h,t,r
    remapped_df = remapped_df[['h','t','r']]

    if out_file:
        with open(out_file,'w') as outfile:
            outfile.write(str(remapped_df.shape[0])+"\n")
            data = "\n".join(["\t".join([str(x) for x in [row['h'],row['t'],row['r']]]) for _, row in remapped_df.iterrows()])
            outfile.write(data)
        if debug:
            print("Data file saved to %s"%out_file)
    return remapped_df 


def write_benchmark_files(directory, basename, input_files, debug=False):

    # Mapping of input file to output name
    kg_out_dict = {'all':'whole2id.txt',
                    'train':'train2id.txt',
                    'test':'test2id.txt',
                    'valid':'valid2id.txt'}
    folder_path = os.path.join(directory, basename)    

    # Make output directory not exist
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        print(f"Folder '{folder_path}' created.")
    else:
        print(f"Folder '{folder_path}' already exists. Files will be overwritten.")
    
    # Create entity-relation id mappings from whole_kg.txt
    if 'all' in input_files:
        whole_kg_file = input_files['all']
    elif 'all_sampled' in input_files:
        whole_kg_file = input_files['all_sampled']
    entity2id, relation2id = create_entity_relation_mappings(whole_kg_file,directory=folder_path)

    # Create train2id.txt, test2id.txt, valid2id.txt, whole2id.txt
    for filetype, filename in input_files.items():
        filetype = filetype.replace("_sampled","")
        out_file_name = os.path.join(folder_path,kg_out_dict[filetype])
        out_df = reformat_kg(filename, entity2id, relation2id, out_file=out_file_name, debug=debug) 


def prepare_benchmarks(input_directory, output_directory):
    '''
    This function prepares the input files for benchmarking. 
    '''

    tag_dict = {'kg1f_instances':'K2BIO_INST',
                'kg2f_ontologies':'K2BIO_ONTO',
                'whole_kg':'K2BIO'} # skip bridges

    # Read input files from input_directory as dict
    input_files = get_input_files(input_directory)

    # Create a benchmark folder for each view
    for base_name, files in input_files.items():
        
        if base_name in tag_dict:
            tag = tag_dict[base_name]
            write_benchmark_files(output_directory, tag, files, debug=True)
        
        
### Parse arguments ###
if len(sys.argv) < 3:
	print("Incorrect usage.")
	print("python prepare_benchmark.py <input_directory> <output_directory>")
	print("Example: python prepare_benchmark.py ../sampled_know2bio_safe_release/ ../../benchmark/data")
	sys.exit(0)
    

input_directory = sys.argv[1]
output_directory = sys.argv[2]

prepare_benchmarks(input_directory,output_directory)

