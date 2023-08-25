import json, csv, pandas as pd, numpy as np, parse_xml, requests, subprocess, os,  urllib.parse, urllib.request, html
from bs4 import BeautifulSoup
from biomedkg_utils import switch_dictset_to_dictlist, switch_dictlist_to_dictset

    
# Source: https://github.com/HHS/uts-rest-api/blob/master/samples/python/Authentication.py
class Authentication:
    '''
    FUNCTION:
    - UMLS API authentication
    '''
    
    def __init__(self, apikey):
        self.apikey=apikey
        self.service="http://umlsks.nlm.nih.gov"

    def gettgt(self):
        params = {'apikey': self.apikey}
        h = {"Content-type": "application/x-www-form-urlencoded", "Accept": "text/plain", "User-Agent":"python" }
        r = requests.post("https://utslogin.nlm.nih.gov/cas/v1/api-key",data=params,headers=h)
        response = fromstring(r.text)
        ## extract the entire URL needed from the HTML form (action attribute) returned - looks similar to https://utslogin.nlm.nih.gov/cas/v1/tickets/TGT-36471-aYqNLN2rFIJPXKzxwdTNC5ZT7z3B3cTAKfSc5ndHQcUxeaDOLN-cas
        ## we make a POST call to this URL in the getst method
        tgt = response.xpath('//form/@action')[0]
        return tgt

    def getst(self,tgt):
        params = {'service': self.service}
        h = {"Content-type": "application/x-www-form-urlencoded", "Accept": "text/plain", "User-Agent":"python" }
        r = requests.post(tgt,data=params,headers=h)
        st = r.text
        return st
    
 
 

def mesh2umls_query(b_id, umls_lol):
    '''
    FUNCTION:
    - Queries the UMLS API to get the corresponding MeSH terms
    
    PARAMS:
    - b_id: batch ID (integer number ranging from 0 to number of processors)
    - umls_lol: list of list of UMLS strings
    '''
    # Personal authentication code for my account
    a = Authentication('5953671b-d1aa-4588-9d54-75472feb476d')  
    
    # Batch size
    tot = str(len(umls_lol))
    
    # Open temp output file
    temp_path = 'output/compound2compound/temp'+str(b_id)+'.txt'
    os.remove(temp_path)
    with open(temp_path,'a') as fout:
    
        # UMLS List for a concept (might have multiple UMLS)
        for count, umls_l in enumerate(umls_lol):
            
            # Print progress(only prints one batch, threadsafe)
            if b_id == 1:
                print('Progress of Batch 1: '+str(count)+'/'+tot, end='\r')

            # Each UMLS from that list above
            for umls in umls_l:

                # UMLS API
                tgt = Authentication.gettgt(a);
                st = Authentication.getst(a,tgt)
                res = requests.get('https://uts-ws.nlm.nih.gov/rest/content/current/CUI/%s/atoms?ticket=%s'%(umls,st)).json()

                # Find the UMLS ID's MeSH ID 
                for i in range(len(res['result'])):
                    try: link = res['result'][i]['sourceDescriptor']
                    except: link = ''
                    if 'MSH' in link:
                        mesh = link.split('/')[-1]
                        fout.write(umls+'|'+mesh+'\n')
                        break

                        
                        
def unii2mesh_meshapi(b_id, unii_batch):
    '''
    FUNCTION:
    - Uses the MeSH API to align UNII with MeSH IDs
    
    PARAMS:
    - b_id: batch ID (integer number ranging from 0 to number of processors)
    - unii_batch: UNII IDs used as input to be mapped to MeSH
    '''
    tot = str(len(unii_batch))
    temp_path = 'output/compound2compound/temp_unii2mesh'+str(b_id)+'.txt'
    open(temp_path,'w').write('')
    with open(temp_path,'a') as fout:
        
        # For each UNII ID
        for count, uniis in enumerate(unii_batch): 
            for unii in uniis:
            
                # Print progress(only prints one batch, threadsafe)
                if b_id == 1:
                    print('Progress of Batch 1: '+str(count)+'/'+tot, end='\r')

                # API, get data and save it
                try: 
                    mesh = requests.get('https://meshb.nlm.nih.gov/search?searchMethod=FullWord&searchInField=allRegistry&sort=&size=1&searchType=anyWord&from=0&q=%s'%unii).text.split('href="/record/ui?ui=')[1].split('\"')[0]
                except:
                    continue

                msg = str(str(unii)+'|'+mesh+'\n')
                fout.write(msg)
            
            

            
            
            
def db2meshorcas_inxightAPI(b_id, db_batch):
    '''
    FUNCTION:
    - Get DrugBank-is-CAS  (via Inxight Drugs)
    - Get DrugBank-is-MeSH (via Inxight Drugs)
    
    PARAMS:
    - b_id: batch ID (integer number ranging from 0 to number of processors)
    - db_batch: DrugBank IDs used as input to be mapped to MeSH
    ''' 
    
    # Prepare files
    temp_db2mesh_path = 'output/compound2compound/temp_db2mesh_inxightapi'+str(b_id)+'.txt'
    temp_db2cas_path = 'output/compound2compound/db2cas_inxight'+str(b_id)+'.txt'
    open(temp_db2mesh_path,'w').write('')
    open(temp_db2cas_path, 'w').write('')
    tot = len(db_batch)
    
    for i, db in enumerate(db_batch):

        ### Print progress
        if b_id == 1:
            print(str(i)+'/'+str(tot), end='\r')
        
        ### DrugBank-is-MeSH attempt (via Inxight Drugs)
        try:
            # API
            inx = ''
            inx = requests.get('https://drugs.ncats.io/substances?q=root_codes_DRUG%5C%20BANK:%22%5E'+db+'$%22').text.split('<a href="/substance/')[1].split('\"')[0]
            raw_data = urllib.request.urlopen('https://drugs.ncats.io/drug/'+inx).read()
            mesh = str(raw_data).split('http://www.ncbi.nlm.nih.gov/mesh/')[2].split('>')[1].split('<')[0].strip()
            
            # Export
            with open(temp_db2mesh_path,'a') as fout1:
                msg = db+'|'+mesh+'\n'
                fout1.write(msg)
                
            
        ### DrugBank-is-CAS attempt (via Inxight Drugs)
        except:
            try:
                if inx == '':
                    continue
                # API
                cas_set = set()
                raw_data = str(urllib.request.urlopen('https://drugs.ncats.io/drug/'+inx).read())
                data_list = raw_data.split('https://chem.nlm.nih.gov/chemidplus/rn/')
                for dat in data_list:
                    c = dat.split('\"')[0].strip()
                    if '\\' not in c:
                        cas_set.add(str(c))
                        
                # Export
                with open(temp_db2cas_path, 'a') as fout2:
                    if len(cas_set) > 0:
                        fout2.write(db)
                        for cas_num in cas_set:
                            fout2.write('|'+str(cas_num))
                        fout2.write('\n')

            except:
                there_are = 'no db to mesh or to cas'
            
            
            
            
            
            
            
            
def unii2mesh_mycheminfoapi(b_id, unii_batch):
    '''
    FUNCTION:
    - Uses the MyChem.info API to align UNII with MeSH IDs
    
    PARAMS:
    - b_id: batch ID (integer number ranging from 0 to number of processors)
    - unii_batch: UNII IDs used as input to be mapped to MeSH
    '''
    
    temp_path = 'output/compound2compound/temp_unii2mesh_mycheminfo'+str(b_id)+'.txt'
    open(temp_path,'w').write('')
    tot = str(len(unii_batch))
    
    for count, unii in enumerate(unii_batch):  
        unii = unii[0]
        
        # Print progress(only prints one batch, threadsafe)
        if b_id == 1:
            print('Progress of Batch 1: '+str(count)+'/'+tot, end='\r')     
        
        # Use API 
        try: 
            mesh = requests.get('http://mychem.info/v1/chem/%s'%unii).json()['umls']['mesh'].strip()
        except: 
            try: mesh = res['drugcentral']['xrefs']['mesh_descriptor_ui'].strip()
            except: continue     
        
        # Save results
        with open(temp_path,'a') as fout:
            fout.write(unii+'|'+mesh+'\n')
            
            
            
            
def cas2mesh_meshapi(b_id, cas_batch):
    '''
    FUNCTION:
    - Get CAS-is-MeSH (via MeSH API)
    
    PARAMS:
    - b_id: batch ID (integer number ranging from 0 to number of processors)
    - cas_batch: cas IDs used as input to be mapped to MeSH
    '''
    tot = str(len(cas_batch))
    temp_path = 'output/compound2compound/temp_cas2mesh'+str(b_id)+'.txt'
    open(temp_path,'w').write('')
    with open(temp_path,'a') as fout:
        
        # For each cas ID
        for count, cases in enumerate(cas_batch): 
            for cas in cases:
                
                # Print progress(only prints one batch, threadsafe)
                if b_id == 1:
                    print('Progress of Batch 1: '+str(count)+'/'+tot, end='\r')

                # API, get data and save it
                try: 
                    res = requests.get('https://meshb.nlm.nih.gov/search?searchMethod=FullWord&searchInField=allRegistry&sort=&size=1&searchType=exactMatch&from=0&q='+cas)
                    mesh = str(res.text).split('href="/record/ui?ui=')[1].split('\"')[0]
                except:
                    continue

                # Export
                msg = str(cas+'|'+mesh+'\n')
                fout.write(msg)
                
                
                
def inchi2mesh_mycheminfoAPI(inchis, inchikey2mesh, mesh2inchikey):
    '''
    FUNCTION:
    - Inchi -is- MeSH (via MyChem.info POST request API)
    
    PARAMS:
    - inchis: list of inchikey IDs
    - inchikey2mesh: Inchi -is- MeSH dictionary
    - mesh2inchikey: MeSH -is- Inchi dictionary
    '''
    for i in range(int(1+len(inchis)/900)):
        # Batch of size 900 (limit was 1,000)
        start_ind = i*900
        end_ind = (i+1)*900
        if len(inchis) < end_ind:
            end_ind = len(inchis)
        print('Trying '+str(start_ind)+'-'+str(end_ind), end='\r')

        # POST Request API
        inchi_batch = inchis[start_ind:end_ind]
        params = {'ids': str(",".join(inchi_batch)), 'fields': 'drugcentral.xrefs.mesh_descriptor_ui,umls.mesh,ginas'}
        time.sleep(3)
        res = requests.post('http://mychem.info/v1/chem', params).json()

        for r in res:
            inchi = r['query']

            # MeSH via DrugCentral
            try: 
                meshes = r['drugcentral']['xrefs']['mesh_descriptor_ui']
            except:

            # MeSH via Ginas
                try: 
                    gcodes = res['ginas']['codes']
                    for section in gcodes:
                        if section['codeSystem'] == 'MESH':
                            meshes = section['code'].strip()

            # MeSH via UMLS
                except: 
                    try: 
                        meshes = r['umls']['mesh']
                    except: 
                        continue

            if type(meshes) == list:
                for mesh in meshes:
                    inchikey2mesh.setdefault(inchi, set()).add(mesh)
                    mesh2inchikey.setdefault(mesh, set()).add(inchi)
            else:
                inchikey2mesh.setdefault(inchi, set()).add(meshes)
                mesh2inchikey.setdefault(meshes, set()).add(inchi)
                
                
def reactome_compound_type_scraper_api(b_id, reactome_comp_batch):
    '''
    FUNCTION:
    - Access compounds' Reactome webpage. Use a scraper to 
      find the type of compound (e.g., ChemicalDrug, Chemical Compound
      Polymer, Genes and Transcripts, OtherEntity, etc.)
    
    PARAMS:
    - b_id: batch ID (integer number ranging from 0 to number of processors)
    - reactome_comp_batch: a batch of Reactome compound IDs
    '''
    
    temp_path = 'output/compound2compound/temp_reactomecomp'+str(b_id)+'.txt'
    open(temp_path,'w').write('')
    
    tot = str(len(reactome_comp_batch)-1)

    for i, reactome_comp in enumerate(reactome_comp_batch):
        
        # Print progress
        if b_id == 1:
            print('Progress of Batch 1:'+str(i)+'/'+tot, end='\r')
        
        # API
        req = requests.get('https://reactome.org/content/detail/'+reactome_comp+'#Homo%20sapiens.html')
        soup = BeautifulSoup(req.content, 'html.parser')
        comp_type = str(soup.find_all('div', class_ ="details-field favth-col-lg-10 favth-col-md-9 favth-col-sm-9 favth-col-xs-12")[1].find_all('span')[0]).split('\n')[0].split('>')[-1]
        
        # Export
        with open(temp_path,'a') as fout:
            fout.write(reactome_comp+'|'+comp_type+'\n')

            
            
def getOld2NewTTDdrug(b_id, oldTTD_list):
    '''
    Function: 
    - Saves JSON files of {old TTD ID : [new TTD ID], ...}
    
    Params: 
    - b_id: batch ID from the multiprocess process
    - oldTTD_list: list of deprecated TTD IDs / old naming scheme
    '''
    oldttd2newttd = dict() 
    no_newttd = set()
    tot = str(len(oldTTD_list))
    
    # Find drug's updated TTD ID
    for i, oldTTD in enumerate(oldTTD_list):
        oldTTD = oldTTD[0]
        
        # Print progress
        if b_id == 1 or int(b_id) > 999:
            print('Progress of Batch 1'+str(i)+'/'+tot, end='\r')
        
        try:
            res = requests.get('http://db.idrblab.net/web/drug/'+oldTTD)
            split_on = '<link rel="canonical" href="http://db.idrblab.net/web/drug/'
            newTTD = res.text.split(split_on)[1].split("\"")[0].upper()
            oldttd2newttd.setdefault(oldTTD, set()).add(newTTD)
        except:
            no_newttd.add(oldTTD)
            continue
     
    # Export Old TTD -is- New TTD
    oldttd2newttd = switch_dictset_to_dictlist(oldttd2newttd)
    fout1 = open('output/compound2compound/temp_newTTDfound_'+str(b_id)+'.json', 'w')
    json.dump(oldttd2newttd, fout1)
    
    # Export Old TTD without New TTD alignments
    no_newttd = list(no_newttd)
    fout2 = open('output/compound2compound/temp_no_newTTDfound_'+str(b_id)+'.json','w')
    json.dump(no_newttd, fout2)
    
     
        
    
import re, time, json, zlib
from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
import requests
from requests.adapters import HTTPAdapter, Retry


POLLING_INTERVAL = 3

API_URL = "https://rest.uniprot.org"


retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session_UniProtAPI = requests.Session()
session_UniProtAPI.mount("https://", HTTPAdapter(max_retries=retries))


def submit_id_mapping_UniProtAPI(from_db, to_db, ids):
    '''
    FUNCITON:
    - This submits the post request to UniProt's API.
    
    PARAMS:
    - from_db (string): The database to map IDs from
    - to_db (string): The database to map IDs to
    - ids (list of strings): The IDs to map from 
    '''
    request = requests.post(
        f"{API_URL}/idmapping/run",
        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
    )
    request.raise_for_status()
    return request.json()["jobId"]


def get_next_link_UniProtAPI(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def check_id_mapping_results_ready_UniProtAPI(job_id):
    '''
    FUNCTION:
    - This checks if the submitted job is ready.
      If the job is not ready, it will tell you and try again.
      If the job is ready, it will tell you it is ready.
      
    PARAMS:
    - job_id: This is the ID of the job submitted to UniProt
    '''
    while True:
        request = session_UniProtAPI.get(f"{API_URL}/idmapping/status/{job_id}")
        request.raise_for_status()
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] == "RUNNING":
                print(f"Job still running. Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(request["jobStatus"])
        else:
            return bool(j["results"] or j["failedIds"])


def get_batch_UniProtAPI(batch_response, file_format, compressed):
    batch_url = get_next_link_UniProtAPI(batch_response.headers)
    while batch_url:
        batch_response = session_UniProtAPI.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results_UniProtAPI(batch_response, file_format, compressed)
        batch_url = get_next_link_UniProtAPI(batch_response.headers)


def combine_batches_UniProtAPI(all_results, batch_results, file_format):
    if file_format == "json":
        for key in ("results", "failedIds"):
            if key in batch_results and batch_results[key]:
                all_results[key] += batch_results[key]
    elif file_format == "tsv":
        return all_results + batch_results[1:]
    else:
        return all_results + batch_results
    return all_results


def get_id_mapping_results_link_UniProtAPI(job_id):
    '''
    FUNCTION:
    - This gets the link where the job results can
      be accessed and then later downloaded by another
      function here.
    
    PARAMS:
    - job_id: This is the ID of the job submitted to UniProt
    '''
    url = f"{API_URL}/idmapping/details/{job_id}"
    request = session_UniProtAPI.get(url)
    request.raise_for_status()
    return request.json()["redirectURL"]


def decode_results_UniProtAPI(response, file_format, compressed):
    if compressed:
        decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
        if file_format == "json":
            j = json.loads(decompressed.decode("utf-8"))
            return j
        elif file_format == "tsv":
            return [line for line in decompressed.decode("utf-8").split("\n") if line]
        elif file_format == "xlsx":
            return [decompressed]
        elif file_format == "xml":
            return [decompressed.decode("utf-8")]
        else:
            return decompressed.decode("utf-8")
    elif file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    elif file_format == "xlsx":
        return [response.content]
    elif file_format == "xml":
        return [response.text]
    return response.text


def get_xml_namespace_UniProtAPI(element):
    m = re.match(r"\{(.*)\}", element.tag)
    return m.groups()[0] if m else ""


def merge_xml_results_UniProtAPI(xml_results):
    merged_root = ElementTree.fromstring(xml_results[0])
    for result in xml_results[1:]:
        root = ElementTree.fromstring(result)
        for child in root.findall("{http://uniprot.org/uniprot}entry"):
            merged_root.insert(-1, child)
    ElementTree.register_namespace("", get_xml_namespace_UniProtAPI(merged_root[0]))
    return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)


def print_progress_batches_UniProtAPI(batch_index, size, total):
    n_fetched = min((batch_index + 1) * size, total)
    print(f"Fetched: {n_fetched} / {total}", end='\r')


def get_id_mapping_results_search_UniProtAPI(url):
    '''
    FUNCTION:
    - Download the API results from a url
    
    PARAMS:
    - url: the link where the job results can
      be accessed and downloaded here.
    '''
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    if "size" in query:
        size = int(query["size"][0])
    else:
        size = 500
        query["size"] = size
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl()
    request = session_UniProtAPI.get(url)
    request.raise_for_status()
    results = decode_results_UniProtAPI(request, file_format, compressed)
    total = int(request.headers["x-total-results"])
    print_progress_batches_UniProtAPI(0, size, total)
    for i, batch in enumerate(get_batch_UniProtAPI(request, file_format, compressed), 1):
        results = combine_batches_UniProtAPI(results, batch, file_format)
        print_progress_batches_UniProtAPI(i, size, total)
    if file_format == "xml":
        return merge_xml_results_UniProtAPI(results)
    return results


def get_id_mapping_results_stream_UniProtAPI(url):
    if "/stream/" not in url:
        url = url.replace("/results/", "/stream/")
    request = session_UniProtAPI.get(url)
    request.raise_for_status()
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    return decode_results_UniProtAPI(request, file_format, compressed)


def get_possible_fields_UniProtAPI():
    '''
    FUNCTION:
    - This prints the possible fields including the
      databases that can be mapped from and to.
    '''
    # Raw results
    all_fields = requests.get('https://rest.uniprot.org/configure/idmapping/fields').json()

    # Database Category to Database dictionary
    database_dict = dict()

    for group in all_fields['groups']:

        # Database category
        database_category = group['groupName']

        # Databases
        databases = list()
        for item in group['items']:
            databases.append(item['name'])

        # Database category to Database
        database_dict[database_category] = databases
    

    return all_fields, database_dict



def get_to_uniprotid_from_genename_mapping_dict_UniProtAPI(results, your_taxa = [''], filter_by_reviewed = True):
    '''
    FUNCTION:
    - Convert the raw API call's mapping results to a dictionary
      where the keys are the "from" IDs and the values are the
      "to" IDs
    
    PARAMS:
    - results (dict): The UniProt API's raw results from an API call
    - your_taxa (list of ints): a list of taxa you want to include
    - filter_by_reviewed (bool): indicates if you only want reviewed entries
    '''
    name2id_up = dict()
    filter_by_taxa = your_taxa != ['']
    
    for i in range(0, len(results['results'])):
        # API mapping results for this entity
        r = results['results'][i]

        # Gene ID
        from_gene_name = r['from']
        try: 
            gene_name = r['to']['genes'][0]['geneName']['value']
        except: 
            continue
        if gene_name != from_gene_name:
            continue
            
        # Protein ID
        to_protein_id = r['to']['primaryAccession']

        
        # Organism taxon number
        taxon_num = r['to']['organism']['taxonId']

        # Is the entry reviewed
        reviewed = ' reviewed ' in r['to']['entryType']

        # If you're not filtering by 'reviewed' or if the entry is reviewed
        if not filter_by_reviewed or reviewed == True:

            # If you're not filtering by taxa or if the taxa is your taxa
            if not filter_by_taxa or taxon_num in your_taxa:
                name2id_up.setdefault(from_gene_name, set()).add(to_protein_id)
                
    return name2id_up


def get_from2to_uniprot(results):
    '''
    FUNCTION:
    - Map the "from_id" to the "to_id" for all requests
    
    PARAMS:
    - results: nested dict/lists, results of UniProt API
    '''
    
    from2to, to2from = dict(), dict()

    for entry in results['results']:    
        from_id = entry['from']
        to_id = entry['to']['primaryAccession']
        from2to.setdefault(from_id, set()).add(to_id)
        to2from.setdefault(to_id, set()).add(from_id)
        
    return from2to, to2from


# Source for STRING API code: Alexander Pelletier
def get_filtered_string_interactors_STRINGAPI(string_ids, score_thresh):
    '''
    Get interactors for STRING IDs using the API, 
    only return those with combined_score >= score_thresh
    '''

    ### PPI API 
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv-no-header"
    method = "interaction_partners"
    request_url = "/".join([string_api_url, output_format, method])
    params = {
      "identifiers" : "%0d".join(string_ids), # your proteins
      "species" : 9606,                       # 9606 = humans 
      "required_score" : score_thresh*1000, # only report those with a score above 0.85
      "caller_identity" : "www.awesome_app.org" # your app name
    }
    response = requests.post(request_url, data=params)

    ### Save PPI results
    results = []
    for line in response.text.strip().split("\n"):
        l = line.strip().split("\t")

        query_ensp = l[0]
        query_name = l[2]
        partner_ensp = l[1]
        partner_name = l[3]
        combined_score = float(l[5])
        
        if combined_score >= score_thresh:
            results += [[query_ensp, query_name, partner_ensp, partner_name, combined_score]]
    return results



def batch_filtered_queries_STRINGAPI(string_ids, batch_size = 100, score_thresh = 0.9):
    '''
    API crashes if you put the whole list. batch the queries
    '''
    string_ids = list(string_ids)
    batched_results = []
    for i in range(0,len(string_ids),batch_size):
        i_end = min(i+batch_size, len(string_ids))
        ids = string_ids[i:i_end]
        # print(len(ids))
        results = get_filtered_string_interactors_STRINGAPI(ids,score_thresh=score_thresh)
        batched_results += results
        if i % 1000 == 0:
            print(i,'/',len(string_ids), end='\r')

    ### If proteins are found 
    if len(batched_results) > 0:
        result_df = pd.DataFrame(batched_results)
        result_df.columns = ['query_ensp', 'query_name', 
                            'partner_ensp', 'partner_name', 'combined_score']
        return result_df
    else:
        print("No results!!")
        return None
    
def k_hop_interactors_STRINGAPI(string_ids, k, score_thresh, debug=False, return_counts=False):
    '''
    Get interactors for STRING IDs using the API within k-hops, 
    only including interacting partners above score_threshold 
    '''
    all_proteins = set(string_ids)
    next_query_proteins = all_proteins
    return_df = None
    counts = []
    for i in range(k):
        # get interactors for the new proteins
        string_interactors_df = batch_filtered_queries_STRINGAPI(next_query_proteins, score_thresh=score_thresh)
        if string_interactors_df is not None:
            # append df
            if return_df is not None:
                return_df = return_df.append(string_interactors_df, ignore_index=True)
            else:
                return_df = string_interactors_df

            all_resulting_proteins = set(string_interactors_df['query_ensp']).union(set(string_interactors_df['partner_ensp']))

            # see how many we added
            new_all_proteins = all_proteins.union(all_resulting_proteins)
            next_query_proteins = new_all_proteins.difference(all_proteins)
            all_proteins = new_all_proteins
            if debug:
                print("%d total\t%d new proteins added for %d hop"%(len(all_proteins),len(next_query_proteins), i+1))
            counts += [(i+1,len(all_proteins),len(next_query_proteins))]
    if debug:
        print("%d total proteins"%(len(all_proteins)))
    if return_counts:
        return return_df, counts    
    return return_df

# k_hop_interactors_STRINGAPI(string_ids)  

