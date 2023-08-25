import urllib.request
import csv
import pandas as pd

def download_go_obo():
    # Download MeSH xml
    url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
    dest = 'input/go-basic.obo'
    urllib.request.urlretrieve(url, dest);

    
def go_obo_to_dict(go_obo_path):
    # Convert GO obo file to dict
    ID = ''
    go_dict = dict()
    with open(go_obo_path) as fin:
        for line in fin:
            if line.startswith('id: '):
                ID = line.split('id: ')[1].strip('\n')
                continue
            if ': ' in line and ID != '':
                k = line.split(': ')[0]
                v = line.split(': ')[1].strip('\n')
                go_dict.setdefault(ID,dict()).setdefault(k,[]).append(v)
    return go_dict
                
    
def get_go_ontology_names(go_dict):
    go_id_to_go_ont = dict()
    for go_id, values in go_dict.items():
        go_ont_name = values['namespace'][0]
        go_id_to_go_ont[go_id] = go_ont_name

    go_ontology_names = set(go_id_to_go_ont.values())
    go_ontology_names.remove('external')
    
    return go_ontology_names, go_id_to_go_ont


def export_go_edges(go_dict, go_id_to_go_ont):
    file_name = 'CC_MF_BP_(GO)_2_CC_MF_BP_(GO).csv'
    with open(f'output/go2go/{file_name}', 'w') as fout:  
        writer = csv.writer(fout)
        writer.writerow(['CC_MF_BP_(GO)', 'CC_MF_BP_(GO)', 'Relationship'])

        for go_id, values in go_dict.items():
            if not go_id.startswith('GO:'):
                continue
            go_id = go_id_to_go_ont[go_id]+':'+go_id.split('GO:')[1]

            # 'is_a' relationships
            try:
                is_a_gos = [go.split(' !')[0] for go in values['is_a']]
                for other_go in is_a_gos:
                    ont_type = go_id_to_go_ont[other_go]
                    other_go = ont_type+':'+other_go.split('GO:')[1]
                writer.writerow([go_id, other_go, '-is_a-'])
            except:
                pass

            # 'part_of', 'regulates', etc. relationships
            try:
                rel_gos = [r_g.split(' !')[0].split(' ') for r_g in values['relationship']]
                for rel, other_go in rel_gos:
                    other_go = ont_type+':'+other_go.split('GO:')[1]
                    writer.writerow([go_id, other_go, f'-{rel}->'])
            except:
                pass
    
    df = pd.read_csv(f'output/go2go/{file_name}')
    df.to_csv(f'output/edges/{file_name}', index=False)
    df.to_csv(f'output/edges_to_use/{file_name}', index=False)


def export_separate_ontology_go_edges(go_dict, go_ontology_names):
    for go2go_type in go_ontology_names:
        cap_go2go_type = '_'.join([w.capitalize() for w in go2go_type.split('_')])

        # Output GO-[rel]->GO
        with open(f'output/go2go/{cap_go2go_type}_(GO)_2_{cap_go2go_type}_(GO).csv','w') as fout:
            writer = csv.writer(fout)
            writer.writerow([f'{cap_go2go_type} (GO)', f'{cap_go2go_type} (GO)','Relationship'])
            go2go = dict()
            rel_counts = dict()

            for go, values in go_dict.items():
                if values['namespace'][0] == go2go_type:

                    for value in values:
                        ### Obsolete? ###
                        if value == 'is_obsolete' and go_dict[go][value] == ['true']:
                            continue
                        if value == 'replaced_by':
                            continue
                        if 'GO' not in go:
                            continue


                        ### Relationships ###
                        if value == 'is_a':
                            rel_type = '-'+value+'->'
                            start_node = go.split('GO:')[1]
                            start_node = go2go_type+':'+start_node
                            end_nodes = go_dict[go][value]
                            for end_node in end_nodes:
                                end_node = end_node.split('GO:')[1]
                                end_node = go2go_type+':'+end_node
                                go2go.setdefault(start_node, dict()).setdefault(rel_type,[]).append(end_node.split(' !')[0])
                                writer.writerow([start_node, end_node.split(' !')[0], rel_type])
                                rel_counts[rel_type] = rel_counts.setdefault(rel_type,0) + 1

                        elif value == 'relationship':
                            rels = go_dict[go][value]
                            for rel in rels:
                                rel_type = '-'+rel.split(' ')[0]+'->'
                                start_node = go.split('GO:')[1]
                                start_node = go2go_type+':'+start_node
                                end_node = rel.split(' ')[1]
                                end_node = end_node.split('GO:')[1]
                                end_node = go2go_type+':'+end_node
                                go2go.setdefault(start_node, dict()).setdefault(rel_type,[]).append(end_node)
                                writer.writerow([start_node, end_node, rel_type])
                                rel_counts[rel_type] = rel_counts.setdefault(rel_type,0) + 1    

        df = pd.read_csv(f'output/go2go/{cap_go2go_type}_(GO)_2_{cap_go2go_type}_(GO).csv')
        df.to_csv(f'output/go2go/{cap_go2go_type}_(GO)_2_{cap_go2go_type}_(GO).csv', index=False)

        df.to_csv(f'output/go2go/{cap_go2go_type}_(GO)_2_{cap_go2go_type}_(GO).csv', index=False)
        df.to_csv(f'output/edges/{cap_go2go_type}_(GO)_2_{cap_go2go_type}_(GO).csv', index=False)
        df.to_csv(f'output/edges_to_use/{cap_go2go_type}_(GO)_2_{cap_go2go_type}_(GO).csv', index=False)


        print(go2go_type)
        display(rel_counts) # Relationship counts
        display(df.head(3))
        
        
if __name__ == '__main__':
    download_go_obo()
    go_dict = go_obo_to_dict(go_obo_path='input/go-basic.obo')
    go_ontology_names, go_id_to_go_ont = get_go_ontology_names(go_dict)
    export_go_edges(go_dict, go_id_to_go_ont)
    export_separate_ontology_go_edges(go_dict, go_ontology_names)