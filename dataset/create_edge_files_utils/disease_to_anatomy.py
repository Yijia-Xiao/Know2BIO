import json
import pandas as pd
import time
import os
try:
    from pubmed_api.pubmed import PubMedAPI
except:
    os.system('pip install pubmed-api')
    from pubmed_api.pubmed import PubMedAPI
from biomedkg_utils import switch_dictset_to_dictlist

def import_anatomy_terms():
    try:
        anatomy_name2id = json.load(open('output/anatomy2anatomy/meshterm-IS-meshid.json'))
    except:
        os.system('python anatomy_to_anatomy.py')
        anatomy_name2id = json.load(open('output/anatomy2anatomy/meshterm-IS-meshid.json'))
    anatomy_terms = list(anatomy_name2id.keys())
    
    return anatomy_terms
    
    
def import_disease_terms():
    try: 
        disease_name2id = json.load(open('output/disease2disease/meshterm-IS-meshid.json'))
    except:
        os.system('python disease_to_disease.py')
        disease_name2id = json.load(open('output/disease2disease/meshterm-IS-meshid.json'))
    disease_terms = list(disease_name2id.keys())

    return disease_terms   

def retrieve_topic_studying_pmids(terms, outpath):
    pa = PubMedAPI()
    mesh_to_pmid = {}
    NUM_TERMS = len(terms)
    for NTH_TERM, term in enumerate(terms):
        if NTH_TERM % 100 == 0:
            print(f'{NTH_TERM}/{NUM_TERMS}', end='\r')
        try:
            response = pa.extract(f'{term}[MeSH Terms]')
            time.sleep(1)
        except:
            print('Error with response for', term)#, 'Trying normal API')
            continue
            #try:
            #    response = call_eutils_api(term)
            #    print('Success')
            #except:
            #    print('Didnt work')
            #    continue
                
        pmids = response.pmids
        mesh_to_pmid[term] = pmids

    # Export
    with open(outpath,'w') as fout:
        json.dump(switch_dictset_to_dictlist(mesh_to_pmid), fout)

        

def retain_pmids_that_have_both_disease_and_anatomy_annotations():
    # Import
    anatomy_mesh_to_pmid = make_sets(json.load(open('output/anatomy2anatomy/anatomy_mesh_to_pmid.json')))
    disease_mesh_to_pmid = make_sets(json.load(open('output/disease2disease/disease_mesh_to_pmid.json')))
    num_anatomy_pmids = len(anatomy_pmids)
    num_disease_pmids = len(disease_pmids)
    print('Anatomy PMIDs', num_anatomy_pmids, 'Disease PMIDs', num_disease_pmids)

    disease_mesh_to_pmid = disease_mesh_to_pmid.copy()
    anatomy_mesh_to_pmid = anatomy_mesh_to_pmid.copy()

    # If any disease-studying PMID does not study an anatomy, remove the PMID (to make proportion calculation more faithful)
    print('#Disease PMIDs - #Removed PMIDs = #DiseaseAnatomy PMIDs')

    # For each disease
    for disease, d_pmids in disease_mesh_to_pmid.copy().items():

        # Remove disease-noanatomy PMIDs for each disease
        disease_mesh_to_pmid[disease] = set(d_pmids).intersection(set(anatomy_pmids))

        # If any PMIDs were removed, display the counts
        diff = len(d_pmids) - len(disease_mesh_to_pmid[disease])
        if diff != 0:
            remain_per = round(1-diff/len(d_pmids),2)*100
            print(len(d_pmids),'-',diff,'=',len(disease_mesh_to_pmid[disease]),'=', str(remain_per)+'%')

            # If the disease has no co-occurring anatomy MeSH terms, remove the disease
            if len(d_pmids) - diff == 0:
                disease_mesh_to_pmid.pop(disease)
                
    disease_mesh_to_pmid = switch_dictset_to_dictlist(disease_mesh_to_pmid)
    anatomy_mesh_to_pmid = switch_dictset_to_dictlist(anatomy_mesh_to_pmid)
    
    json.dump(disease_mesh_to_pmid, open('output/disease2disease/filtered_disease_mesh_to_pmid.json','w'))
    json.dump(anatomy_mesh_to_pmid, open('output/disease2disease/filtered_anatomy_mesh_to_pmid.json','w'))
                        
                
def calculate_relationships_between_disease_and_anatomy(disease_mesh_to_pmid, anatomy_mesh_to_pmid):

    disease_mesh_to_pmid = json.load(open('output/disease2disease/filtered_disease_mesh_to_pmid.json','w'))
    anatomy_mesh_to_pmid = json.load(open('output/disease2disease/filtered_anatomy_mesh_to_pmid.json','w'))
                     
    obs_disease_anatomy_prop = dict()
    tot = len(disease_mesh_to_pmid)

    # For each disease
    for i, (di, di_pmids) in enumerate(disease_mesh_to_pmid.items()):

        print(i,'/',tot, end='\r')

        # Initialize disease_i's dictionary
        obs_disease_anatomy_prop[di] = dict()

        # For each anatomy
        for aj, aj_pmids in anatomy_mesh_to_pmid.items():

            # Shared disease_i-anatomy_j PMIDs
            di_aj_inter_pmids = set(aj_pmids).intersection(set(di_pmids))

            # Observed ratio of anatomy_j-studying PMIDs that are disease_i-studying PMIDs out of disease_i-studying PMIDs
            # (I like thish better than Jaccard similairty)
            obs_disease_anatomy_prop[di][aj] = len(di_aj_inter_pmids)/len(di_pmids)

    with open('output/disease2anatomy/obs_disease_anatomy_prop.json','w') as fout:
        json.dump(obs_disease_anatomy_prop, fout)


    da_prop = json.load(open('output/disease2anatomy/obs_disease_anatomy_prop.json'))

    # Sort the Disease-Anatomy proportion relationships
    sort_da_prop = dict()
    for d in da_prop:
        sort_da_prop[d] = dict()
        tuplelist = sorted(da_prop[d].items(), key=lambda x:x[1], reverse=True)
        sortdict = dict(tuplelist)
        for k,v in sortdict.items():
            if v > 0.01: # this could change
                sort_da_prop[d][k] = v

    anatomy_name2id = json.load(open('output/anatomy2anatomy/meshterm-IS-meshid.json'))
    disease_name2id = json.load(open('output/disease2disease/meshterm-IS-meshid.json'))


    # Output to edges file
    with open('output/disease2anatomy/edges_disease2anatomy.csv','w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['Disease (MeSH)','Anatomy (MeSH)','Relationship','Weight'])

        for disease, anatomies in sort_da_prop.items():
            try: disease_id = disease_name2id[disease][0]
            except: print('bad disease', disease); continue

            for anatomy,weight in anatomies.items():
                try: anatomy_id = anatomy_name2id[anatomy][0]
                except: print('bad anatomy', anatomy); continue

                writer.writerow(['MeSH_Disease:'+disease_id, 'MeSH_Anatomy:'+anatomy_id, '-affects->', weight])

    df = pd.read_csv('output/disease2anatomy/edges_disease2anatomy.csv').drop_duplicates()
    df.to_csv('output/edges/edges_disease2anatomy.csv', index = False)
    df.to_csv('output/edges_to_use/Disease_(MeSH)_2_Anatomy_(MeSH).csv', index = False)


if __name__ == '__main__':
    anatomy_mesh_to_pmidoutpath = 'output/anatomy2anatomy/anatomy_mesh_to_pmid.json'
    disease_mesh_to_pmidoutpath = 'output/disease2disease/disease_mesh_to_pmid.json'

    anatomy_terms = import_anatomy_terms()
    disease_terms = import_disease_terms()

    retrieve_topic_studying_pmids(anatomy_terms, anatomy_mesh_to_pmidoutpath)
    retrieve_topic_studying_pmids(disease_terms, disease_mesh_to_pmidoutpath)
    
    retain_pmids_that_have_both_disease_and_anatomy_annotations()
    calculate_relationships_between_disease_and_anatomy()