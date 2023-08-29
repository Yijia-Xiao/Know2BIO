import json
import requests
import subprocess
import os
from downloader import *


def check_length(missing_ids, protein_id2sequence, lower=16, upper=1280):
    too_short = []
    too_long = []
    for uid in missing_ids:
        seq = protein_id2sequence[uid]
        if len(seq) < lower:
            too_short += [uid]
        elif len(seq) > upper:
            too_long += [uid]
    print("%d sequences were shorter than %d AA"%(len(too_short),lower))
    print("%d sequences were longer than %d AA"%(len(too_long),upper))
    
    return too_short, too_long


def check_nonstandard(missing_ids, protein_id2sequence, nonstandard="BZXJ"):
    nonstandard_aa = []
    for uid in missing_ids:
        seq = protein_id2sequence[uid]
        found_nonstandard = False        
        for aa in nonstandard:
            if aa in seq:
                found_nonstandard = True
                break
        if found_nonstandard:
            nonstandard_aa += [uid]
    print("%d sequences contained nonstandard AA's from %s"%(len(nonstandard_aa),nonstandard))
    
    return nonstandard_aa


def check_against_uniprot_ospg(missing_ids, uniprot_ospg_file):
    uniprot_ref_oppg = set(l.strip("\n").split("|")[1] for l in open(fasta_file,"r").readlines() if ">" in l)
    uniprot_missing_ospg = set(missing_ids).difference(uniprot_ref_oppg)
    print("%d sequences were not found in Uniprot One protein Sequence per Gene"%len(uniprot_missing_ospg))
    return list(uniprot_missing_ospg)


if __name__ == "__main__":
    print("Checking downloaded files")

    redownload = False
    display_examples = True

    # directories
    output_directory = '../output/node_features/structures'
    cif_directory = os.path.join(output_directory,'protein_3d_structure')
    pae_directory = os.path.join(output_directory,'protein_3d_structure_prediction_errors')
    cif_files = os.listdir(cif_directory)
    pae_files = os.listdir(pae_directory)
    print("%d files in %s"%(len(cif_files),cif_directory))
    print("%d files in %s"%(len(pae_files),pae_directory))
    
    # load protein ids and sequences to download
    prot2seq_file = '../output/node_features/sequences/protein_id_to_sequences.json'
    protein_id2sequence = read_prot2seq(prot2seq_file)
    
    
    count = 0
    missing_ids = []
    for uid in protein_id2sequence.keys():
        uid_found = any([uid in u for u in cif_files])
        if not uid_found:
            missing_ids += [uid]
    print("%d missing id's"%len(missing_ids))
    
    if redownload:
        print("Attempting re-download")
        count = 0
        for uid in missing_ids:
            # download structure
            valid_structure = get_3d_structure(uid, cif_directory, pae_directory, version=4)
        
            if valid_structure:
                missing_ids.pop(uid)
                count += 1
        print("%d out of %d successfully downloaded."%(count,len(missing_ids)+count))

    ### Examine why these are missing ###
    ## FAQ 1: ##
    #     - Outside length range. Minimum length is 16 AA's, maximum is 2700 or 1280 AA.
    too_short, too_long = check_length(missing_ids, protein_id2sequence)

    ## FAQ 2: ##
    #    - Non-standard AA's (e.g. X)
    nonstandard_aa = check_nonstandard(missing_ids, protein_id2sequence)
        
    ## FAQ 3: ##
    #     - Not in UniProt's One protein Sequence Per Gene list
    #fasta_file = "2023-02-01_UP000005640_9606.fasta"
    #uniprot_ospg = check_against_uniprot_ospg(missing_ids, fasta_file)

    # Summary
    explained_missing = set(too_short).union(set(too_long)).union(set(nonstandard_aa))
    #explained_missing = set(too_short).union(set(too_long)).union(set(nonstandard_aa)).union(set(uniprot_ospg))
    unexplained_missing = set(missing_ids).difference(explained_missing)
    print("%d out of %d were explained by FAQ's. %d unexplained." %(len(explained_missing),len(missing_ids),len(unexplained_missing)))

    if display_examples:
        n=5
        print(too_short[:n])
        print(too_long[:n])
        print(nonstandard_aa[:n])
        #print(uniprot_ospg[:n])
        print(list(unexplained_missing)[:n])

