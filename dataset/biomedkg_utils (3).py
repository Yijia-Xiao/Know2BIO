import csv, multiprocessing, subprocess, os, pandas as pd, numpy as np
import requests, lxml.html as lh, json, pandas as pd, numpy as np, os
from lxml.html import fromstring
from multiprocessing import Pool, cpu_count, Process




def output_edgefile_onerel_noweight(outpath, columns, dictionary, rel, prefix_col1='', prefix_col2='', edges_folder = True):
    '''
    FUNCTION:
    - Reformat a dictionary {start node: end nodes, ...} to an edge file 
    (start node, end node, relationship)
    
    PARAMS:
    - outpath: the path to output the edge files (string)
    - columns: the row of column headers for the node types (list)
    - dictionary: the dictionary with {start node: end node, ...}
    - rel: the relationship (e.g., '-targets->')
    - prefix_col1: the prefix for the start node ID (e.g., 'UniProt:')
    - prefix_col2: the prefix for the end node ID (e.g., 'Entrez:')
    '''
    # Open output file
    with open(outpath,'w') as fout:
        writer = csv.writer(fout)
        
        # Write column headers
        writer.writerow(columns)
        
        # Write edges (start, end, relationship)
        for k,vs in dictionary.items():
            for v in vs:
                writer.writerow([prefix_col1+str(k),prefix_col2+str(v),rel])
           
    # Drop duplicates
    df = pd.read_csv(open(outpath)).drop_duplicates()
    df.to_csv(outpath, index=False)
    
    # Export to edges folder if this isn't an intermediary file but is a final file
    if edges_folder == True:
        edges_folder_outpath = 'output/edges/'+outpath.split('/')[2]
        df.to_csv(edges_folder_outpath, index=False)
             
                
                
                
    

def multiprocess_a_dict_of_lists_values(thedict, the_function):
    '''
    FUNCTION:
    - This takes a dictionary whose keys are strings and 
      values are lists of strings and splits it into input
      for separate processes. The processes then output 
      their results to temp files which are then merged. 
    
    PARARMS:
    - thedict: the dictionary to be split into input for 
      a multiprocessing function
    - the_function: the function that will use the dictionary
      as input for multiprocessing
    '''
    
    # How many processors can be used
    procs = cpu_count()

    # List of batches for multiprocessing
    batches = [[] for i in range(procs)]   

    # Create batches and send to multiprocessing
    for i,(key, values_list) in enumerate(thedict.items()):

        # Add synonym to a batch
        b_id = i%procs
        batches[b_id].append(values_list)

    # Create a list of jobs
    print("Running jobs...")
    jobs = []
    for b_id, batch in enumerate(batches):
        jobs.append(Process(target = the_function, \
                            args = [b_id, batch]))

    # Run the jobs
    for j in jobs: j.start()
    for j in jobs: j.join()
    print('Done!')
            
            
def multiprocess_a_list(thelist, the_function):
    '''
    FUNCTION:
    - This takes a list of strings and splits it into input
      for separate processes. The processes then output 
      their results to temp files which are then merged. 
    
    PARARMS:
    - thelist: the list to be split into input for 
      a multiprocessing function
    - the_function: the function that will use the list
      as input for multiprocessing
    '''
    
    # How many processors can be used
    procs = cpu_count()

    # List of batches for multiprocessing
    batches = [[] for i in range(procs)]   

    # Length of input dictionary 
    tot = len(thelist)

    # Create batches and send to multiprocessing
    for i, item in enumerate(thelist):

        # Add synonym to a batch
        b_id = i%procs
        batches[b_id].append(item)

    # Create a list of jobs 
    print("Running jobs...")
    jobs = []
    for b_id, batch in enumerate(batches):
        jobs.append(Process(target = the_function, \
                            args = [b_id, batch]))

    # Run the jobs
    for j in jobs: j.start()
    for j in jobs: j.join()
    print('Done!')
    
    
    
def switch_dictset_to_dictlist(the_dict):
    '''
    FUNCTION:
    - Make a new dictionary with values as lists 
      instead of values as sets
      
    PARAMS:
    - the_dict: The initial dict with values of sets
    '''
    
    dictlist = dict()
    
    for k in the_dict.copy():
        dictlist[k] = list(the_dict[k])
        
    return dictlist



def switch_dictlist_to_dictset(the_dict):
    '''
    FUNCTION:
    - Make a new dictionary with values as sets 
      instead of values as list
      
    PARAMS:
    - the_dict: The initial dict with values of lists
    '''
    
    dictset = dict()
    
    for k in the_dict.copy():
        dictset[k] = set(the_dict[k])
        
    return dictset


def get_all_entries_in_listvaluedict(dct, outputset):
    '''
    FUNCTION:
    - Get all dictionary entries in a dictionary with
      values of lists or sets
      
    PARAMS:
    - dct: dictionary with values of lists or sets
    - outputset: the set of values to output
    
    '''
    for k,vs in dct.items():
        for v in vs:
            outputset.add(v)

    
    
def key2mesh_via_db(key2db, db2mesh, key2mesh, mesh2key):
    '''
    FUNCTION:
    - key-is-MeSH = key-is-DrugBank + DrugBank-is-MeSH
    
    PARAMS:
    - key2db: original dictionary of some key, key-is-DrugBank
    - db2mesh: DrugBank-is-MeSH dictionary alignment/mapping
    - key2mesh: some key-is-MeSH
    - mesh2key: some MeSH-is-key
    '''

    for key, dbs in key2db.items():
        for db in dbs:
            try:
                meshes = db2mesh[db]
                for mesh in meshes:
                    key2mesh.setdefault(key,set()).add(mesh)
                    mesh2key.setdefault(mesh,set()).add(key)
            except:
                continue    
    
    
    
def gene_expressed_in_anatomy(df, gene2anat, anat2gene, uberon2mesh, ensembl_is_entrez):
    '''
    FUNCTION:
    - Gene -[over/under]expressed_in- Anatomy relationships
    - Uses Bgee data to make a dictionary of those relationships^
    
    PARAMS:
    - df: dataframe of under or overexpressed genes to anatomy from Bgee
    - gene2anat: dictionary mappings expressed genes to anatomy 
    - anat2gene: dictionary mappings anatomy to expressed genes
    - uberon2mesh: dictionary mapping UBERON -is- MeSH
    - ensembl_is_entrez: dictionary mapping Ensembl -is- Entrez ID
    '''
    for i in range(len(df)):
        
        try:
            # Gene
            gene = df['Gene ID'].iloc[i]
            entrez_gene = ensembl_is_entrez[gene]
            
            # Anatomy
            anat = df['Anatomical entity ID'].iloc[i]
            mesh_anat = uberon2mesh[anat]
            
            # Gene -[over/under]expressed in- Anatomy
            gene2anat.setdefault(entrez_gene, set()).add(mesh_anat)
            anat2gene.setdefault(mesh_anat, set()).add(entrez_gene)

        except:
            continue
            
    return gene2anat, anat2gene
    
    
    
    
class ParseXML():
    '''
    This is the class for parsing the DrugBank XML file to extract relevant information (ID, name, entities, etc.) for cardiovascular drugs. 
    '''
    def getCAS(ele):
        '''
        @param ele is the element in XML root
        @return CAS identifier
        '''
        try:
            CAS = ele.find('{http://www.drugbank.ca}cas-number').text
        except:
            CAS = None
        return CAS
    
    def getUNII(ele):
        '''
        @param ele is the element in XML root
        @return UNII identifier
        '''
        try:
            UNII = ele.find('{http://www.drugbank.ca}unii').text
        except:
            UNII = 'Null'
        return UNII
    
    def getID(ele):
        '''
        @param ele is the element in XML root
        @return ID is the DrugBank accession number/ID of the element (drug)
        '''
        try:
            ID = ele.find("{http://www.drugbank.ca}drugbank-id").text
        except:
            ID = "Null"
        return ID

    def getName(ele):
        '''
        @param ele is the element in XML root
        @return name is the name of the element (drug)
        '''
        try:
            name = ele.find("{http://www.drugbank.ca}name").text
        except:
            name = "Null"
        return name

    def getSynonyms(ele):
        '''
        @param ele is the element in XML root
        @return allsyns is the list of all the synonyms of the element (drug)
        '''
        try:
            synonyms = ele.find("{http://www.drugbank.ca}synonyms")
            allsyns = []
            for syn in synonyms:
                allsyns.append(syn.text)     
        except:
            allsyns = []
            
        return allsyns

    def getDescription(ele):
        '''
        @param ele is the element in XML root
        @return description is the description of the element (drug)
        '''
        try:
            description = ele.find("{http://www.drugbank.ca}description").text
        except:
            description = "Null"
        return description

    def getCategory(ele):
        '''
        @param ele is the element in XML root
        @return allcats is the list of all DrugBank categories of the element (drug)
        '''
        allcats = []
        try:
            category = ele.find("{http://www.drugbank.ca}categories").\
                    findall("{http://www.drugbank.ca}category")
            for cat in category:
                allcats.append(cat.find('{http://www.drugbank.ca}category').text)
        except:
            allcats = []
        return allcats

    def getATCCode(ele):
        '''
        @param ele is the element in XML root
        @return allcodes is the list of all ATC codes of the element (drug)
        '''
        allcodes = []
        try:
            atccodes = ele.find("{http://www.drugbank.ca}atc-codes").\
                    findall("{http://www.drugbank.ca}atc-code")
            for code in atccodes:
                allcodes.append(code.get('code'))
        except:
            allcodes = []
        return allcodes

    def getPathways(ele):
        '''
        @param ele is the element in XML root
        @return pathways is the list of all pathways of the element (drug). Each
        element in the list is a dictionary containing the name, smpdb_id, and category
        of the pathway.
        '''
        pathways = []
        try:
            path = list(ele.find("{http://www.drugbank.ca}pathways"))
            for term in path:
                pdict = {}
                name = term.find("{http://www.drugbank.ca}name").text
                smpdb_id = term.find("{http://www.drugbank.ca}smpdb-id").text
                category = term.find("{http://www.drugbank.ca}category").text 
                
                pdict.update({"name": name,\
                            "smpdb_id": smpdb_id,\
                            "category":category})
                
                pathways.append(pdict)
            
        except:
            pathways = "Null"
            
        return pathways

    def getIndication(ele):
        '''
        @param ele is the element in XML root
        @return indication is the indication of the element (drug)
        '''
        try:
            indication =   ele.find("{http://www.drugbank.ca}indication").text
        except:
            indication = "Null"
            
        return indication

    
    def getEntities(ele, entity_type):
        
        '''
        @param ele is the element in XML root
        @param entity_type is one of the following: 'carriers', 'targets', 'transporters' , 'enzymes'
        @return allEntities is the list of all entities (carriers, targets, transporters, or enzymes) for the element (drug). 
        Each element in the list contains the name, DrugBank accession number/ID, list of actions, and UniProt ID of the element. 
        '''
        
        allEntities =  []
        try:
            entity = ele.find('{http://www.drugbank.ca}' + entity_type).\
                            findall('{http://www.drugbank.ca}' +\
                            entity_type[:len(entity_type)-1])
            for child in entity:
                if child.find('{http://www.drugbank.ca}organism').text=='Humans':
                    entity_dict = {}
                    #find entity name-------------------------
                    try:
                        e_name = child.find('{http://www.drugbank.ca}name').text
                    except:
                        e_name = "Null"

                    # find actions------------------------------------    
                    try:
                        e_actions = []
                        for term in list(child.find('{http://www.drugbank.ca}actions')):
                            e_actions.append(term.text)
                    except:
                        e_actions = "Null"

                    # find entity drugbank ids------------------------------
                    try:
                        e_id = child.find('{http://www.drugbank.ca}id').text
                    except:
                        e_name = "Null"

                    # find entity Uniprot ids------------------------------    
                    try:
                        e_uid = child.find('{http://www.drugbank.ca}polypeptide').\
                                        get('id')
                    except:
                        e_uid = "Null"

                    entity_dict.update({"name": e_name,\
                                    "drugbank_id": e_id,\
                                    "actions" : e_actions,\
                                    "uniprot_id" : e_uid}) 
                    allEntities.append(entity_dict)
        except:
            allEntities = []
                
        return allEntities
    
    
    def getExternalIDs(ele):
        '''
        @param ele is the element in the XML root
        @return database2ids is a dictionary with the ontology as keys and 
        [ids] as the value for each ontology
        '''
        ontologies2ids = dict()
        r = ele.find('{http://www.drugbank.ca}external-identifiers')
        for e in r:
            for el in e.iter():
                if e != el:
                    if 'resource' in str(el):
                        resource = el.text
                    elif 'identifier' in str(el):
                        identifier = el.text
                        ontologies2ids.setdefault(resource,[]).append(identifier)
        return ontologies2ids
    
    
    
    
    
def get_num_entries_in_dict(dct):
    '''
    FUNCTION:
    Counts the number of entries in a dictionary

    PARAMS:
    - dct: a dictionary
    '''
    count = 0

    for k,vs in dct.items():
        for v in vs:
            count += 1

    return count




def merge_two_dictionaries_setvalues(d1, d2):
    '''
    FUNCTION:
    - Makes a new dictionary that merges to input dictionaries
    - Returns a new dictionary

    PARAMS:
    - d1: input dictionary 1
    - d2: input dictionary 2
    '''

    merged_d = dict()

    for k,vs in d1.items():
        for v in vs:
            merged_d.setdefault(k,set()).add(v)

    for k,vs in d2.items():
        for v in vs:
            merged_d.setdefault(k,set()).add(v)

    return merged_d



def key2value2ttd(key2value, value2ttd, ttd2db):
    for key, values in key2value.items():
        for value in values:
            try:
                ttds = value2ttd[value]
                for ttd in ttds:
                    ttd2db.setdefault(ttd,set()).add(key)
            except:
                continue  
                    
       
    
def merge_no_newttd():
    '''
    FUNCTION:
    - Merge the output files for old TTDs without 
      old TTD -is- new TTD alignments
    '''
    no_newttd = list()

    for b_id in range(procs):
        batch_no_newttd = json.load(open('output/compound2compound/temp_no_newTTDfound_'+str(b_id)+'.json'))
        no_newttd += batch_no_newttd
        
    return no_newttd



def merge_oldttd2newttd():
    '''
    FUNCTION:
    - Merge the output files for old TTD -is- new TTD
      alignments
    '''
    oldttd2newttd = dict()
    procs = cpu_count()
    
    for b_id in range(procs):
        outpath = 'output/compound2compound/temp_newTTDfound_'+str(b_id)+'.json'
        temp = json.load(open(outpath))
        oldttd2newttd = merge_two_dictionaries_setvalues(oldttd2newttd, temp)
        #os.remove(outpath)
    
    return oldttd2newttd



def getTTD_IDs(info):
    '''
    FUNCTION:
    - Gets the TTD ID to Alternate External ID mappings for the drugs
    
    PARAMS:
    - info: line in an input file matching TTD IDs
    '''
    entry = dict()
    for alist in info:
        if 'TTDDRUID' in alist:
            ttd_id = alist[2]
            entry[ttd_id] = dict()
        if 'ChEBI_ID' in alist:
            chebi = alist[2].split(':')[1]
            entry[ttd_id]['Chebi'] = chebi
        if 'PUBCHCID' in alist:
            pubchemc = alist[2]
            entry[ttd_id]['PubChemCompound'] = pubchemc
        if 'PUBCHSID' in alist:
            pubchemsubs = alist[2].split(';')
            entry[ttd_id]['PubChemSubstances'] = pubchemsubs
        if 'CASNUMBE' in alist:
            try: cas = alist[2].split('CAS ')[1]
            except: cas = alist[2]
            entry[ttd_id]['CAS'] = cas
    return entry



def align_ttd_to_xrefs(ttd2altids, cas2ttd, pubchemc2ttd, pubchemsubs2ttd, chebi2ttd):
    '''
    FUNCTION:
    - Align TTD drug IDs with other external reference 
      ontologies/IDs
    
    PARAMS:
    - ttd2altids: TTD aligned with all the alt IDs 
    - cas2ttd: CAS IDs aligned with TTD IDs
    - pubchemc2ttd: PubChem Compound IDs aligned with TTD IDs
    - pubchemsubs2ttd: PubChem Substance IDs aligned with TTD IDs
    - chebi2ttd: ChEBI IDs aligned with TTD IDs
    '''
    for ttd, altids in ttd2altids.items():
        if 'CAS' in altids:
            cas = altids['CAS']
            cas2ttd.setdefault(cas,[]).append(ttd)

        if 'Chebi' in altids:
            chebi = altids['Chebi']
            chebi2ttd.setdefault(chebi,[]).append(ttd)

        if 'PubChemCompound' in altids:
            pcc = altids['PubChemCompound']
            pubchemc2ttd.setdefault(pcc,[]).append(ttd)

        if 'PubChemSubstances' in altids:
            for pcs in altids['PubChemSubstances']:
                pubchemsubs2ttd.setdefault(pcs.strip(),[]).append(ttd)

                
# Read through TTD-provided file
def getTTD_to_UniProtname(info):
    '''
    FUNCTION:
    - Align TTD Target ID -is- UniProt Name
    
    PARAMS:
    - info: line in TTD input file
    '''
    up_name, targetID = '',''
    for line in info:
        if 'TARGETID' in line:
            targetID = line[2]
        if 'UNIPROID' in line:
            up_name = line[2]

    if up_name == '' or targetID == '':
        return '',''
    else:
        return targetID, up_name
    
                
def print_if_len_value_greater_than_1(dct):
    '''
    FUNCTION
    - Print key-value pairs if the value list has 
      more than one entry in it. Works for value sets too.
    
    PARAMS:
    - dct: dictionary input
    '''
    for k,v in dct.items():
        if len(v) != 1:
            print(k, v)
            
            

def get_short_full_names(names, section_key, section_values):
    '''
    FUNCTION:
    - Get all names from an entry
    
    PARAMS:
    - names (list)
    '''
    # Full name
    if section_key == 'fullName':
        name = section_values['value']
        names.append(name)

    # Short names
    if section_key == 'shortNames':
        for entry in section_values:
            name = entry['value']
            names.append(name)  
            
    # EC Numbers
    if section_key == 'ecNumbers':
         for entry in section_values:
            name = entry['value']
            names.append(name) 
            
    return names

    
def get_protein_names_from_id(IDs, gene_names=False):
    '''
    Params: 
    - IDs (list): list of UniProt IDs, 
    - gene_names (bool): Indicates if you want the gene names
   
    Function: 
    - Uses REST API to get a proteins' names
    
    Output: 
    - dictionary of IDs to names
    '''
    
    protein_id2names = dict()
    bad_ids = list()


    SIZE = str(500)

    # For each ID, get the protein names (one protein at a time)
    for index, ID in enumerate(IDs):
        print(index+1,'/',len(IDs),end='\r')
        names = list()

        # Get data via UniProt API 
        if type(ID) != str:
            continue

        # UniProt API to get gene names
        fields = 'accession,protein_name,gene_names,ec'
        url = 'https://rest.uniprot.org/uniprotkb/search?fields='+fields\
              +'&format=json&query='+ID+'&size='+SIZE
        res = requests.get(url).json()

        '''Protein Names'''
        try:
            r = res['results'][0]['proteinDescription']
        except:
            bad_ids.append(ID)
            continue

        ### Recommended name ###
        try:
            # Find sections
            rec_name_section = r['recommendedName']
            for section_key, section_values in rec_name_section.items():

                # Get names
                names = get_short_full_names(names, section_key, section_values)

            #print(ID, rec_name_section)
            #print(names)
        except: pass

        ### Alternative Names ###
        try:
            # Find sections
            alt_name_section = r['alternativeNames']
            for entry in alt_name_section:
                for section_key, section_values in entry.items():

                    # Get names
                    names = get_short_full_names(names, section_key, section_values)
        except: pass

        ### Submitted Names ###
        try:
            # Find sections
            sub_name_section = r['submissionNames']
            for entry in sub_name_section:
                for section_key, section_values in entry.items():

                    # Get names
                    names = get_short_full_names(names, section_key, section_values)     
        except: pass


        '''Gene Names'''
        if gene_names == True:
            # Data for gene names
            try: rgenes = res['results'][0]['genes'][0]
            except: continue

            # Get gene names
            for nametype, main_entry in rgenes.items():
                if type(main_entry) == dict:
                    name = main_entry['value']
                elif type(main_entry) == list:
                    for entry in main_entry:
                        name = entry['value']
                names.append(name)


        '''UniProt ID -> Protein & Gene Names'''
        for name in names:
            protein_id2names.setdefault(ID, list()).append(name)
            
    return protein_id2names



def map_pmc_pmids(batch_id, pmids):
    '''
    FUNCTION:
    - Map a batch of PMIDs to PMCs
    
    PARAMS:
    - batch_id: The batch ID
    - pmids (list): the batch of PMIDs
    '''
    url = 'https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids='
    tot_reqs = int(len(pmids)/200)+1
    temp_outfile = 'pmc_pmid_temp_'+str(batch_id)+'.txt'
    
    with open(temp_outfile,'w') as fout:
        for i in range(0, int(len(pmids)/200)+1):

            # Print progress
            if batch_id == 1:
                print(i, '/', tot_reqs, end='\r')  

            # Submit batch of up to 200 PMIDs
            if i < int(len(pmids)/200)+1:
                r = req.get(url+','.join(pmids[i*200:(i+1)*200])+'&format=json')
            else:
                r = req.get(url+','.join(pmids[i*200:])+'&format=json')

            # Map all PMID-PMCs
            records = r.json()['records']
            for record in records:
                try:
                    pmc = record['pmcid']
                    pmid = record['pmid']
                    fout.write(pmid+'|'+pmc)
                except:
                    break

                    
def merge_pmid_pmc_mapping(pmid2pmc, pmc2pmid):
    '''Merge the files output by the parallel processes
    
    PARAMS:
    - pmid2pmc: Dictionary to map the PMIDs to the PMCs
    - pmc2pmid: Dictionary to map the PMCs to the PMIDs
    '''
    procs = cpu_count()

    for batch_id in range(procs):
        temp_outfile = 'pmc_pmid_temp_'+str(batch_id)+'.txt'
        with open(temp_outfile) as fin:
            for line in fin:
                line = line.split('|')
                pmid = line[0]
                pmc = line[1].replace('PMC','')
                pmid2pmc.setdefault(pmid, set()).add(pmc)
                pmc2pmid.setdefault(pmc, set()).add(pmid)
            
    return pmid2pmc, pmc2pmid