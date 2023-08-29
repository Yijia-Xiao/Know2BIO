import json
import requests
import subprocess
import os


def runcmd(cmd, verbose = False, *args, **kwargs):

	process = subprocess.Popen(
		cmd,
		stdout = subprocess.PIPE,
		stderr = subprocess.PIPE,
		text = True,
		shell = True
	)
	std_out, std_err = process.communicate()
	if verbose:
		print(std_out.strip(), std_err)
	pass


def check_website_exists(url):
	response = requests.get(url)
	website_exists = response.status_code == 200
	return website_exists


def get_3d_structure(uid, cif_out_dir, pae_out_dir, version=4):
	
	pae_url = "https://alphafold.ebi.ac.uk/files/AF-%s-F1-predicted_aligned_error_v%d.json"%(uid,version)
	cif_url = "https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_v%d.cif"%(uid,version)
	
	cif_valid = check_website_exists(cif_url)
	pae_valid = check_website_exists(pae_url)
	
	# check if url is valid
	if cif_valid:
		# download file
		runcmd('wget %s'%(cif_url))
		# move file into directory
		base_name = cif_url.split("/")[-1]
		runcmd("mv %s %s"%(base_name,cif_out_dir))	  
	
	# check if url is valid
	if pae_valid:
		# download file
		runcmd('wget %s'%(pae_url))
		# move file into directory
		base_name = pae_url.split("/")[-1]
		runcmd("mv %s %s"%(base_name,pae_out_dir))	  
   	
	# keep track of proteins with valid structures
	return cif_valid
	
def read_prot2seq(file_name):
	pid2seq = json.loads(open(file_name,'r').read())
	# it's uniprot id to list, but the list has 100% duplicate entries
	# take the first entry only
	protein_id2sequence = {u:sl[0] for u,sl in pid2seq.items()}
	print(len(protein_id2sequence),"protein ids with sequences")
	return protein_id2sequence

if __name__ == "__main__":
	print("Downloading 3d structures")

	# make output directories
	output_directory = '../output/node_features/structures/'
	cif_directory = os.path.join(output_directory,'protein_3d_structure')
	pae_directory = os.path.join(output_directory,'protein_3d_structure_prediction_errors')
	if not os.path.exists(cif_directory):
		os.makedirs(cif_directory)
	if not os.path.exists(pae_directory):
		os.makedirs(pae_directory)
	
	# load protein ids and sequences to download
	prot2seq_file = '../output/node_features/sequences/protein_id_to_sequences.json'
	protein_id2sequence = read_prot2seq(prot2seq_file)
	
	count = 0
	missing_ids = []
	for uid in protein_id2sequence.keys():
		# download structure
		valid_structure = get_3d_structure(uid, cif_directory, pae_directory, version=4)
	
		if not valid_structure:
			missing_ids += [uid]
	
		# print progress
		count += 1
		if count % 10 == 0:
			print("%d out of %d processed. %d missing."%(count,
											len(protein_id2sequence.keys()),
											len(missing_ids)),end='\r')
	print("%d out of %d successfully downloaded."%(len(protein_id2sequence.keys())-len(missing_ids),
								len(protein_id2sequence.keys())))


