import os
import pandas as pd
import json
from biomedkg_utils import switch_dictset_to_dictlist, switch_dictlist_to_dictset
import requests
import zipfile
from disease_to_disease import download_mesh_xml, parse_mesh_xml
from pathway_to_pathway import remove_letter_before_kegg_pathway
from multiprocessing import cpu_count
from joblib import Parallel, delayed
from compound_to_compound_alignment import parse_drugbank_xml
from parse_xml import *
from gene_to_anatomy import batch_iterator


class NodeNames():
    
    def __init__(self):
        self.all_ids = []
        self.root = ''
    
    def get_all_node_ids(self):
        '''Get all IDs of the nodes'''
        all_ids = set()

        root = 'output/edges_to_use/'
        dfs = []
        for file in os.listdir(root):
            if file.endswith('.csv'):
                df = pd.read_csv(os.path.join(root, file))
                dfs.append(df)
                # Columns and Entries
                try:
                    cols = list(df.columns)
                    first_col_ids = set(df[cols[0]])
                    sec_col_ids = set(df[cols[1]])
                    all_ids = all_ids.union(first_col_ids, sec_col_ids)
                    #print('got IDs from', file)
                except:
                    print(file, 'has an issue')

        json.dump(list(all_ids), open('output/otherMappings/all_ids.json','w'))
        self.all_ids = list(set(all_ids))
        
        all_ids_no_prefix = list([ID.split(':')[1] for ID in all_ids])
        #duplicate_items = [item for item in all_ids_no_prefix if all_ids_no_prefix.count(item) > 1]
        assert len(all_ids) == len(all_ids_no_prefix)
        print(len(all_ids), 'nodes/IDs')
        ids_without_prefixes = [ID for ID in all_ids if ':' not in ID]
        print(len(ids_without_prefixes), 'IDs without prefixes')
        if len(ids_without_prefixes) > 0:
            print(ids_without_prefixes)
        
        self.mesh_compound_ids = list({ID.split(':')[1] for ID in self.all_ids if ID.startswith('MeSH_Compound:')})
        self.mesh_disease_ids = list({ID.split(':')[1] for ID in self.all_ids if ID.startswith('MeSH_Disease:')})
        self.mesh_anatomy_ids = list({ID.split(':')[1] for ID in self.all_ids if ID.startswith('MeSH_Anatomy:')})

        
    def map_protein_id_to_name(self):
        all_protein_ids = list({id_.split(':')[1] for id_ in self.all_ids if id_.startswith('UniProt')})
        url = 'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cprotein_name&format=tsv&query=%28*%29+AND+%28reviewed%3Atrue%29+AND+%28model_organism%3A9606%29'
        req = requests.get(url)

        rows = [row.split('\t') for row in req.text.split('\n')]
        rows.remove([''])
        protein_ids = []
        protein_names = []
        for idx in range(len(rows[1:])):
            protein_id = 'UniProt:'+rows[idx+1][0]
            protein_ids.append(protein_id)
            protein_name = rows[idx+1][1].replace(')','').split(' (')
            protein_names.append(protein_name)

        protein_id_to_name = dict(zip(protein_ids, protein_names))

        self.protein_id_to_name = remove_ids_not_in_all_ids(protein_id_to_name, self.all_ids)
        print(len(self.protein_id_to_name),'/',len(all_protein_ids),'protein IDs mapped to names')
        with open('output/node_features/natural_language_names/protein_id_to_names.json','w') as fout:
            json.dump(protein_id_to_name, fout)

            
    def map_gene_id_to_name(self):
        all_gene_ids = list({ID.split(':')[1] for ID in self.all_ids if ID.startswith('Entrez:')})
        gene_id_to_name = dict()

        BATCH_SIZE = 1000
        for i in range(int(len(all_gene_ids)/BATCH_SIZE)+1):
            headers = {'Accept': 'application/json', 
                      'species':'human',
                      'fields':'symbol,name,taxid'}
            data ={'ids': ','.join(all_gene_ids[i*BATCH_SIZE:(i+1)*BATCH_SIZE])}
            url = 'https://mygene.info/v3/gene?fields=symbol%2Cname%2Ctaxid&species=human'
            r = requests.post(url=url, headers=headers, json=data)
            gene_id_to_name_batch = process_mygeneinfo_entrez_to_name(r.json())
            for gene_id, gene_name in gene_id_to_name_batch.items():
                gene_id_to_name['Entrez:'+gene_id] = gene_name

        self.gene_id_to_name = remove_ids_not_in_all_ids(gene_id_to_name, self.all_ids)
        json.dump(gene_id_to_name, open('output/node_features/natural_language_names/gene_id_to_name.json','w'))
        print(f'{len(gene_id_to_name)}/{len(all_gene_ids)} genes mapped to names')

        
    def map_pathway_id_to_name(self):
        all_pathway_ids = list({ID.split(':')[1] for ID in self.all_ids if 'Pathway' in ID})
        self.pathway_id_to_name, self.pathway_name_to_id = {}, {}
        try:
            with open('input/ReactomePathways.txt') as fin:
                pass
        except:
            os.system('wget -N -P input/ https://reactome.org/download/current/ReactomePathways.txt')
        
        '''Reactome'''
        with open('input/ReactomePathways.txt') as fin:
            for line in fin:
                line = line.strip().split('\t')
                pathway_id = 'Reactome_Pathway:'+line[0]
                pathway_name = line[1]
                species = line[2] 

                if species.lower() == 'homo sapiens':
                    self.pathway_id_to_name.setdefault(pathway_id, set()).add(pathway_name)
                    self.pathway_name_to_id.setdefault(pathway_name, set()).add(pathway_id)

        '''KEGG'''
        with open('input/KEGG/kegg_pathway_hierarchy.csv') as fin:
            for line in fin:
                line = line.split(' ')
                first_number = line[0]
                if '.' not in first_number:
                    id_, name = remove_letter_before_kegg_pathway(line)
                    self.pathway_id_to_name.setdefault('KEGG_Pathway:path_hsa'+id_, set()).add(name)
                
        '''SMPDB'''
        try:
            smpdb_dfs = pd.read_csv('input/smpdb_files/smpdb_pathways.csv')
        except:
            os.system('wget -N -P input/smpdb_files/ http://www.smpdb.ca/downloads/smpdb_pathways.csv.zip')
            with zipfile.ZipFile('input/smpdb_files/smpdb_pathways.csv.zip') as zip_ref:
                zip_ref.extractall('input/smpdb_files')
            smpdb_dfs = pd.read_csv('input/smpdb_files/smpdb_pathways.csv')
            
        for i in range(len(smpdb_dfs)):

            # SMPDB Pathway ID, Description
            smpdb_id = 'SMPDB_Pathway:'+smpdb_dfs['SMPDB ID'].iloc[i]
            smpdb_name = smpdb_dfs['Name'].iloc[i]
            descript = smpdb_dfs['Description'].iloc[i]

            # SMPDB Pathway ID to Description
            self.pathway_id_to_name.setdefault(smpdb_id, set()).add(smpdb_name)
            self.pathway_id_to_name.setdefault(smpdb_id, set()).add(descript)
            self.pathway_name_to_id.setdefault(smpdb_name, set()).add(smpdb_id)
            self.pathway_name_to_id.setdefault(descript, set()).add(smpdb_id)
        
        print(f'{len(self.pathway_id_to_name)}/{len(all_pathway_ids)} pathways mapped to names')
        self.pathway_name_to_id = switch_dictset_to_dictlist(self.pathway_name_to_id)
        self.pathway_id_to_name = remove_ids_not_in_all_ids(switch_dictset_to_dictlist(self.pathway_id_to_name), self.all_ids)
        json.dump(self.pathway_id_to_name, open('output/node_features/natural_language_names/pathway_id_to_name.json','w'))
        json.dump(self.pathway_name_to_id, open('output/pathway2pathway/pathway_name_to_id.json','w'))

               
    def map_anatomy_id_to_name(self):
        all_anatomy_ids = list({ID.split(':')[1] for ID in self.all_ids if ID.startswith('MeSH_Anatomy:')})
        all_anatomy_trees = list({ID.split(':')[1] for ID in self.all_ids if ID.startswith('MeSH_Tree_Anatomy:')})
        if self.root == '':
            download_mesh_xml()
            self.root = parse_mesh_xml()
        self.anatomy_id_to_name, self.anatomy_name_to_id = dict(), dict()
        self.anatomy_tree_to_name, self.anatomy_name_to_tree = dict(), dict()

        for ele in self.root:
            try:
                # MeSH Tree Number
                tree_numbers = ele.find('TreeNumberList').findall('TreeNumber')
                for tree_number in tree_numbers:
                    tree_number = tree_number.text
                    if tree_number.startswith(('A01','A02','A03','A04',
                                               'A05','A06','A07','A08','A09','A10',
                                               'A11','A12','A14','A15','A16',
                                               'A17.815','A17.360.710','A17.360.421',
                                               'A17.360.296')):                    

                        # ID to Name
                        try:
                            ID = 'MeSH_Anatomy:'+ele.find('DescriptorUI').text
                            name = ele.find('DescriptorName').find('String').text
                            self.anatomy_name_to_id.setdefault(name,set()).add(ID)
                            self.anatomy_id_to_name.setdefault(ID,set()).add(name)
                        except:
                            pass

                        # Tree to Name
                        try:
                            name = ele.find('DescriptorName').find('String').text
                            tree_number = 'MeSH_Tree_Anatomy:'+tree_number
                            self.anatomy_name_to_tree.setdefault(name, set()).add(tree_number)
                            self.anatomy_tree_to_name.setdefault(tree_number, set()).add(name)
                        except:
                            pass
            except:
                continue        

        self.anatomy_name_to_id = switch_dictset_to_dictlist(self.anatomy_name_to_id)
        self.anatomy_id_to_name = remove_ids_not_in_all_ids(switch_dictset_to_dictlist(self.anatomy_id_to_name), self.all_ids)
        self.anatomy_name_to_tree = switch_dictset_to_dictlist(self.anatomy_name_to_tree)
        self.anatomy_tree_to_name = remove_ids_not_in_all_ids(switch_dictset_to_dictlist(self.anatomy_tree_to_name), self.all_ids)

                
        print(f'{len(self.anatomy_tree_to_name)}/{len(all_anatomy_trees)} anatomy tree numbers mapped to names')
        #print(f'{len(self.anatomy_id_to_name)}/{len(all_anatomy_ids)} anatomy IDs mapped to names')
        json.dump(self.anatomy_name_to_id, open('output/node_features/natural_language_names/anatomy_name_to_id.json','w'))
        json.dump(self.anatomy_id_to_name, open('output/node_features/natural_language_names/anatomy_id_to_name.json','w'))
        json.dump(self.anatomy_name_to_tree, open('output/node_features/natural_language_names/anatomy_name_to_tree.json','w'))
        json.dump(self.anatomy_tree_to_name, open('output/node_features/natural_language_names/anatomy_tree_to_name.json','w'))

    
    def map_disease_id_to_name(self):
        all_disease_ids = list({ID.split(':')[1] for ID in self.all_ids if ID.startswith('MeSH_Disease:')})
        all_disease_trees = list({ID.split(':')[1] for ID in self.all_ids if ID.startswith('MeSH_Tree_Disease:')})
   
        if self.root == '':
            download_mesh_xml()
            self.root = parse_mesh_xml()
        self.disease_id_to_name, self.disease_name_to_id = dict(), dict()
        self.disease_tree_to_name, self.disease_name_to_tree = dict(), dict()

        for ele in self.root:
            try:
                # MeSH Tree Number
                tree_numbers = ele.find('TreeNumberList').findall('TreeNumber')
                for tree_number in tree_numbers:
                    tree_number = tree_number.text
                    if tree_number.startswith(('C','F03')):                    

                        # ID to Name
                        try:
                            ID = 'MeSH_Disease:'+ele.find('DescriptorUI').text
                            name = ele.find('DescriptorName').find('String').text
                            self.disease_name_to_id.setdefault(name,set()).add(ID)
                            self.disease_id_to_name.setdefault(ID,set()).add(name)
                        except:
                            pass

                        # Tree to Name
                        try:
                            name = ele.find('DescriptorName').find('String').text
                            tree_number = 'MeSH_Tree_Disease:'+tree_number
                            self.disease_name_to_tree.setdefault(name, set()).add(tree_number)
                            self.disease_tree_to_name.setdefault(tree_number, set()).add(name)
                        except:
                            pass
            except:
                continue        

        self.disease_name_to_id = switch_dictset_to_dictlist(self.disease_name_to_id)
        self.disease_id_to_name = remove_ids_not_in_all_ids(switch_dictset_to_dictlist(self.disease_id_to_name), self.all_ids)
        self.disease_name_to_tree = switch_dictset_to_dictlist(self.disease_name_to_tree)
        self.disease_tree_to_name = remove_ids_not_in_all_ids(switch_dictset_to_dictlist(self.disease_tree_to_name), self.all_ids)
        
        print(f'{len(self.disease_tree_to_name)}/{len(all_disease_trees)} disease tree numbers mapped to names')
        #print(f'{len(self.disease_id_to_name)}/{len(all_disease_ids)} disease IDs mapped to names')
        json.dump(self.disease_name_to_id, open('output/node_features/natural_language_names/disease_name_to_id.json','w'))
        json.dump(self.disease_id_to_name, open('output/node_features/natural_language_names/disease_id_to_name.json','w'))
        json.dump(self.disease_name_to_tree, open('output/node_features/natural_language_names/disease_name_to_tree.json','w'))
        json.dump(self.disease_tree_to_name, open('output/node_features/natural_language_names/disease_tree_to_name.json','w'))

    
    def map_mesh_id_to_names(self):
        ''' Includes drugs, diseases, and anatomies'''
        mesh_ids = list({id_.split(':')[1] for id_ in self.all_ids if id_.startswith('MeSH_') and '.' not in id_})
        mesh_compound_ids = list({id_.split(':')[1] for id_ in self.all_ids if id_.startswith('MeSH_Compound:') and '.' not in id_})
        all_anatomy_ids = list({ID.split(':')[1] for ID in self.all_ids if ID.startswith('MeSH_Anatomy:')})
        all_disease_ids = list({ID.split(':')[1] for ID in self.all_ids if ID.startswith('MeSH_Disease:')})
        self.disease_id_to_name = switch_dictlist_to_dictset(self.disease_id_to_name)
        self.anatomy_id_to_name = switch_dictlist_to_dictset(self.anatomy_id_to_name)
       
        # Parse the MRCONSO file across CPUs
        with open('input/MRCONSO.RRF') as fin:
            lines = fin.readlines()
        procs = cpu_count()
        line_chunks = [lines[idx::procs] for idx in range(procs)]
        mesh_dicts = Parallel(n_jobs=-1)(delayed(process_mrconso_line)(idx, some_lines, mesh_ids) 
                                         for idx, some_lines in enumerate(line_chunks))

        # Merge the parallelized results
        mesh_id_to_name = {}
        for dictionary in mesh_dicts:
            for id_,names in dictionary.items():
                if id_ in self.mesh_compound_ids:
                    prefix = 'MeSH_Compound:'
                elif id_ in self.mesh_disease_ids:
                    prefix = 'MeSH_Disease:'
                elif id_ in self.mesh_anatomy_ids:
                    prefix = 'MeSH_Anatomy:'
                else:
                    print(id_, 'not found in your MeSH IDs')
                for name in names:
                    if prefix == 'MeSH_Disease:':
                        self.disease_id_to_name.setdefault(prefix+id_, set()).add(name)
                    if prefix == 'MeSH_Anatomy:':
                        self.anatomy_id_to_name.setdefault(prefix+id_, set()).add(name)
                    if prefix == 'MeSH_Compound:':
                        mesh_id_to_name.setdefault(prefix+id_, set()).add(name)
        # Export
        self.disease_id_to_name = switch_dictset_to_dictlist(self.disease_id_to_name)
        self.anatomy_id_to_name = switch_dictset_to_dictlist(self.anatomy_id_to_name)
       
        json.dump(self.disease_id_to_name, open('output/node_features/natural_language_names/disease_id_to_name.json','w'))
        json.dump(self.anatomy_id_to_name, open('output/node_features/natural_language_names/anatomy_id_to_name.json','w'))
        print(f'{len(self.anatomy_id_to_name)}/{len(all_anatomy_ids)} anatomy IDs mapped to names')
        print(f'{len(self.disease_id_to_name)}/{len(all_disease_ids)} disease IDs mapped to names')
        print(f'{len(mesh_id_to_name)}/{len(mesh_compound_ids)} MeSH compounds mapped to names')
        self.mesh_id_to_name = remove_ids_not_in_all_ids(switch_dictset_to_dictlist(mesh_id_to_name), self.all_ids)
        with open('output/node_features/natural_language_names/compound_mesh_id_to_name.json','w') as fout:
            json.dump(self.mesh_id_to_name, fout)

    
    def map_drugbank_compound_id_to_name(self):
        '''DrugBank'''
        all_drugbank_ids = list({ID for ID in self.all_ids if ID.startswith('DrugBank')})
        root = parse_drugbank_xml()
        self.drugbank_compound_id_to_name = {}
        self.drugbank_compound_name_to_id = {}
        for i, ele in enumerate(root):
            drug_id = 'DrugBank_Compound:'+ParseXML.getID(ele)    # DrugBank drug_id
            name = ParseXML.getName(ele)  
            synonyms = ParseXML.getSynonyms(ele)
            names = ([name] + synonyms)
            for name in names:
                self.drugbank_compound_id_to_name.setdefault(drug_id, set()).add(name)
                self.drugbank_compound_name_to_id.setdefault(name, set()).add(drug_id)

        # Export
        self.drugbank_compound_id_to_name = switch_dictset_to_dictlist(self.drugbank_compound_id_to_name)
        self.drugbank_compound_id_to_name = remove_ids_not_in_all_ids(switch_dictset_to_dictlist(self.drugbank_compound_id_to_name), self.all_ids)
        self.drugbank_compound_name_to_id = switch_dictset_to_dictlist(self.drugbank_compound_name_to_id)
        json.dump(self.drugbank_compound_id_to_name, open('output/node_features/natural_language_names/drugbank_drug_id_to_name.json','w'))
        json.dump(self.drugbank_compound_name_to_id, open('output/node_features/natural_language_names/drugbank_drug_id_to_name.json','w'))
        print(f'{len(self.drugbank_compound_id_to_name)}/{len(all_drugbank_ids)} DrugBank compounds mapped to names')

        
    def map_go_id_to_name(self):
        go_prefixes = ['biological_process', 'cellular_component', 'molecular_function']
        all_go_ids = list({ID.split(':')[1] for ID in self.all_ids if ID.startswith(tuple(go_prefixes))})
              
        ''' Convert GO obo file to dict '''
        ID = ''
        go_dict = dict()
        with open('input/go-basic.obo') as fin:
            for line in fin:
                if line.startswith('id: '):
                    ID = line.split('id: ')[1].strip('\n')
                    continue
                if ': ' in line and ID != '':
                    k = line.split(': ')[0]
                    v = line.split(': ')[1].strip('\n')
                    go_dict.setdefault(ID,dict()).setdefault(k,[]).append(v)


        ''' Get GO ID->Name mapping'''
        go_id_to_name, go_name_to_id = dict(), dict()

        for ID, values in go_dict.items():

            # Name
            if 'name' in values:
                name = values['name'][0]

            # Ontology
            ontology = values['namespace'][0]
            if ontology not in go_prefixes:
                continue        

            # Description
            try: 
                text_description = go_dict[ID]['def'][0].split('\"')[1]
                ID = ontology+':'+ID.split('GO:')[1]
                go_id_to_name.setdefault(ID, set()).add(text_description)
                go_name_to_id.setdefault(name, set()).add(ID)
            except:
                pass

            go_id_to_name.setdefault(ID, set()).add(name)
            go_name_to_id.setdefault(text_description, set()).add(ID)

        self.go_name_to_id = switch_dictset_to_dictlist(go_name_to_id)
        self.go_id_to_name = remove_ids_not_in_all_ids(switch_dictset_to_dictlist(go_id_to_name), self.all_ids)
        print(f'{len(self.go_id_to_name)}/{len(all_go_ids)} GO ids mapped to names')
        json.dump(self.go_id_to_name, open('output/node_features/natural_language_names/go_id_to_name.json','w'))
        json.dump(self.go_name_to_id, open('output/node_features/natural_language_names/go_name_to_id.json','w'))

        
    def map_reaction_id_to_name(self):
        all_reaction_ids = list({ID.split(':')[1] for ID in self.all_ids if ID.startswith('Reactome_Reaction')})
        BATCH_SIZE = cpu_count()
        reactions = json.load(open('input/reactome_reactions.json'))
        reaction_batches = [reactions[i:i+BATCH_SIZE] for i in range(0, len(reactions), BATCH_SIZE)]
        reaction_to_name_dicts = Parallel(n_jobs=-1)(delayed(api_to_map_reaction_id_to_name)(reaction_batch) 
                                                         for reaction_batch in reaction_batches)        
        reaction_to_name = {}
        for reaction_to_name_batch in reaction_to_name_dicts:
            reaction_to_name.update(reaction_to_name_batch)

        self.reaction_id_to_name = remove_ids_not_in_all_ids(reaction_to_name, self.all_ids)
        with open('output/node_features/natural_language_names/reaction_id_to_name.json','w') as fout:
            json.dump(self.reaction_id_to_name, fout)
        print(f'{len(self.reaction_id_to_name)}/{len(all_reaction_ids)} reactions mapped to names')

        
def api_to_map_reaction_id_to_name(reactions):
    reaction_id_to_name = {}
              
    url = 'https://reactome.org/ContentService/data/query/ids/map'
    headers = {'accept': 'application/json', 
               'Content-Type': 'text/plain'}
    BATCH_SIZE = 20

    for idx, id_batch in enumerate(batch_iterator(reactions, BATCH_SIZE)):
        print(f'{idx}/{round(len(reactions)/BATCH_SIZE)}', end='\r')

        # API
        reaction_ids = ','.join(id_batch)
        response = requests.post(url, 
                                 headers=headers, 
                                 data=reaction_ids)
        reaction_dict = response.json()

        # Map API results
        reaction_ids = reaction_ids.split(',')
        for reaction_id in reaction_ids:
            try:
                name = [reaction_dict[reaction_id]['displayName']]
            except:
                continue
            reaction_id_to_name['Reactome_Reaction:'+reaction_id] = name

    return reaction_id_to_name
        
    
def remove_ids_not_in_all_ids(my_dict, all_ids):
    return {key:value for key,value in my_dict.items() if key in all_ids}

        
def process_mygeneinfo_entrez_to_name(r):
    '''Process the MyGene.Info POST request of Entrez genes mapping to gene name '''
    gene_id_to_name_batch = dict()
    for entry in r:

        # Check if it's human
        try:
            if entry['taxid'] != 9606:
                continue
        except:
            pass

        # Check if it's found
        try:
            if entry['notfound']:
                continue
        except:
            pass

        # Gene ID
        gene_id = entry['query']

        # Gene name
        try: 
            gene_name = entry['name']
        except: 
            gene_name = entry['symbol']

        # Gene ID to name
        gene_id_to_name_batch[gene_id] = gene_name

    return gene_id_to_name_batch


def process_mrconso_line(idx, some_lines, mesh_ids_no_prefix):
    mesh_id_to_name = {}
    for idx, line in enumerate(some_lines):
        line = line.strip().split('|')
        ontology = line[11]
        if ontology == 'MSH':
            mesh = line[10]
            name = line[14]
            if mesh in mesh_ids_no_prefix:
                mesh_id_to_name.setdefault(mesh, set()).add(name)
    return mesh_id_to_name


def export_all_ids_to_names():
    id_to_name = {}
    path = 'output/node_features/natural_language_names/'
    for file in os.listdir(path):
        if 'id_to_' in file:
            id_to_name_dict = json.load(open(os.path.join(path,file)))
            id_to_name.update(id_to_name_dict)

    with open(os.path.join(path, 'all_ids_to_names.json'),'w') as fout:
        json.dump(id_to_name, fout)


if __name__ == '__main__':
    NN = NodeNames()
    NN.get_all_node_ids()
    NN.map_protein_id_to_name()
    NN.map_gene_id_to_name()
    NN.map_pathway_id_to_name()
    NN.map_disease_id_to_name()
    NN.map_anatomy_id_to_name()
    NN.map_go_id_to_name()
    NN.map_mesh_id_to_names()
    NN.map_drugbank_compound_id_to_name()
    NN.map_reaction_id_to_name()
    export_all_ids_to_names()