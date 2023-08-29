import os 
from Bio.PDB import *
import numpy as np
from tqdm import tqdm
import pdb
from os.path import exists
from absl import app, flags
import time
from multiprocessing import Pool

directory = '../output/node_features/structures/protein_3d_structure'
out_dir = '../output/node_features/structures/protein_contact_map'
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
# parser = PDBParser()
parser = MMCIFParser()

thres = 10

a = os.popen('ls '+ directory)

filenames = [line.strip() for line in a]

def main(argv):
    t0 = time.time()
    # parser = PDBParser()
    a = os.popen('ls '+ directory)
    filenames = [ line.strip() for i, line in enumerate(a)]
    print("%d files to process"%len(filenames))
    with Pool(96) as p:
        p.map(getadj, filenames)
        #print(p.map(getadj, filenames))
    print(time.time()-t0, "seconds.")
        

def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    return np.sqrt(np.sum(diff_vector * diff_vector))

def calc_dist_matrix(chain_one, chain_two) :
    """Returns a matrix of C-alpha distances between two chains"""
    answer = np.zeros((len(chain_one), len(chain_two)), float)
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            answer[row, col] = calc_residue_dist(residue_one, residue_two)
    return answer

def find_res(chain):
    res_ls = []
    for r in chain:
        if is_aa(r):
            res_ls.append(r) 
    return res_ls


def getadj(filename):
    global out_dir
    uid = filename.split('-')[1]
    if exists(os.path.join(out_dir,uid+'.cont_map.npy')):
        return None
    structure = parser.get_structure(uid, os.path.join(directory, filename))
    # # retain the longest chain
    if len(structure[0])>1:
        print(uid, len(structure[0]))
        return None
    res = structure[0]['A'].get_residues()
    res = [r for r in res if is_aa(r)]

    res_ls = find_res(res)
    dist_map = calc_dist_matrix(res_ls, res_ls)
    cont_map = dist_map < thres
    cont_map = cont_map.astype(int)
    fw = open(os.path.join(out_dir, uid+ '.cont_map.npy'),'bw')
    np.save(fw, cont_map)

if __name__ == '__main__':
    app.run(main)
