import xml.etree.ElementTree as ET
from parse_xml import *
from biomedkg_utils import output_edgefile_onerel_noweight
import json
import zipfile
import os

'''
INSTRUCTIONS:
- File:  https://go.drugbank.com/releases/5-1-9/downloads/all-full-database
- Create a DrugBank account and request an academic license. Then, log in and download the above file in the latest version
'''


def extract_drugbank_xml():
    zip_path = 'input/drugbank_all_full_database.xml.zip'
    destination_folder = 'input/'

    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(destination_folder)

        
def parse_drugbank_xml():
    print('Parsing DrugBank XML')
    tree = ET.parse('input/full database.xml')
    root = tree.getroot()
    return root 


def map_interacting_drugbank_compounds(root):    
    drug2interacting_drugs = dict()
    for i, ele in enumerate(root):
        # Main drug's DrugBank ID
        db_id = ParseXML.getID(ele)

        # Interacting drugs' DrugBank IDs
        interacting_drugs = ParseXML.getInteractingDrugs(ele)

        # Drug - Drug interaction
        drug2interacting_drugs[db_id] = interacting_drugs

    outpath_compound2compound = 'output/compound2compound/edges_MeSH-interacts_with->MeSH.csv'
    output_edgefile_onerel_noweight('output/compound2compound/edges_drugbank-interacts_with->drugbank.csv',
                                   ['Compound (DrugBank)', 'Compound (DrugBank)', 'Relationship'],
                                   drug2interacting_drugs,
                                    '-interacts_with->',
                                   'DrugBank_Compound:',
                                   'DrugBank_Compound:')
    drugbank_outfile = 'output/edges_to_use/Compound_(DrugBank)_2_Compound_(DrugBank).csv'
    output_edgefile_onerel_noweight(drugbank_outfile,
                                   ['Compound (DrugBank)', 'Compound (DrugBank)', 'Relationship'],
                                   drug2interacting_drugs,
                                    '-interacts_with->',
                                   'DrugBank_Compound:',
                                   'DrugBank_Compound:')
    return drug2interacting_drugs


def map_interacting_mesh_compounds_via_drugbank(root, drug2interacting_drugs):    
    db2mesh = json.load(open('output/compound2compound/db2mesh.json'))
    meshdrug2interacting_drugs = dict()

    for db, intdbs in drug2interacting_drugs.items():

        try:
            # Main MeSH
            mesh_mains = db2mesh[db]
            for mesh_main in mesh_mains:

                # Interacting MeSH
                for intdb in intdbs:
                    mesh_ints = db2mesh[intdb]
                    for mesh_int in mesh_ints:

                        # MeSH -interacts- MeSH Drug/Compound
                        meshdrug2interacting_drugs.setdefault(mesh_main, set()).add(mesh_int)         
        except:
            continue


    outpath_compound2compound = 'output/compound2compound/edges_MeSH-interacts_with->MeSH.csv'
    output_edgefile_onerel_noweight(outpath_compound2compound,
                                   ['Compound (MeSH)', 'Compound (MeSH)', 'Relationship'],
                                   meshdrug2interacting_drugs,
                                    '-interacts_with->',
                                   'MeSH_Compound:',
                                   'MeSH_Compound:')
    mesh_outfile = 'output/edges_to_use/Compound_(MeSH)_2_Compound_(MeSH).csv'
    output_edgefile_onerel_noweight(mesh_outfile,
                                   ['Compound (MeSH)', 'Compound (MeSH)', 'Relationship'],
                                   meshdrug2interacting_drugs,
                                    '-interacts_with->',
                                   'MeSH_Compound:',
                                   'MeSH_Compound:')

if __name__ == '__main__':
    extract_drugbank_xml()
    root = parse_drugbank_xml()
    drug2interacting_drugs = map_interacting_drugbank_compounds(root)
    map_interacting_mesh_compounds_via_drugbank(root, drug2interacting_drugs)
    print('Finished mapping compound interactions')