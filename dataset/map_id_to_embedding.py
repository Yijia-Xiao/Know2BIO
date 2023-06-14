'''Get all names'''
import json

# Get all names
id2name_dicts = json.load(open('output/nodes/id2name_dicts.json'))
all_names = set()

for database in id2name_dicts:
    names = set(name for name_list in id2name_dicts[database].values() for name in name_list)
    all_names = all_names.union(names)
    
    
    
''' Add prefixes to the IDs '''
for database in id2name_dicts:
    prefix = database.replace(' ','_')
    if database != 'GO':
        id2name_dicts[database] = {str(prefix+':'+id_):names for id_,names in id2name_dicts[database].items() if ':' not in id_}
        
    
    
'''Create a master dictionary of IDs to names'''
id_to_names = dict()
for database in id2name_dicts:
    id_to_names.update(id2name_dicts[database])

# and for names to IDs
name_to_ids = dict()
for id_,names in id_to_names.items():
    for name in names:
        name_to_ids.setdefault(name, list()).append(id_)
        
        
        
'''Optional: Read in embeddings'''
name_to_embedding_dict = dict()
with open('output/nodes/node features/names_to_embeddings_dicts.txt','r') as fin:
    for i, line in enumerate(fin):
        if i % 1000 == 0:
            print(i, end='\r')
        
        
        try: kv = json.loads(line.strip())
        except:
            try: kv = json.loads('{'+line.strip())
            except: print(line); break
        
        name_to_embedding_dict.update(kv)
        
        if i > 1000:
            break
        #break # <-- remove this if you want, but it will crash the kernel likely
        
        
        
# NOTE: We could try to find which name is best for an ID or just pick the representative name at random
id_to_embedding = dict()
for name, embedding in name_to_embedding_dict.items():
    ids = name_to_ids[name]
    for id_ in ids:
        if id_ not in id_to_embedding:
            id_to_embedding[id_] = embedding
            
json.dump(id_to_embedding, open('output/node features/id_to_embedding.json','w'))