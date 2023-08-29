import os
import csv
import pandas as pd
import json

def map_go_to_go_type():
    '''
    Function:
    - Map each go to its GO ontology: MF, BP, or CC
    Params:
    - go_obo_path (str): Path to the GO OBO file. 
    '''    
    os.system('wget -N -P input/ http://purl.obolibrary.org/obo/go/go-basic.obo')
    obo_path = 'input/go-basic.obo'
    
    go_to_go_type = dict()
    go_type_to_go = dict()
    ont_type_map = {'biological_process':'BP',
                   'molecular_function':'MF',
                   'cellular_component':'CC',
                   'external':'external'}

    with open(obo_path) as fin:
        for line in fin:

            # GO Term
            if line.startswith('id: '):
                go_id = line.split(' ')[1].strip()

            # Ontology (MF, BP, CC)
            elif line.startswith('namespace: '):
                ont_type = line.split(' ')[1].strip()
                #ont_type = ont_type_map[ont_type]
                go_to_go_type[go_id] = ont_type
                go_type_to_go.setdefault(ont_type, list()).append(go_id)

    go_to_go_type = go_to_go_type
    go_type_to_go = go_type_to_go
    
    return go_type_to_go, go_to_go_type


def download_go_annotation_file_gaf():
    go_gaf_url = 'http://geneontology.org/gene-associations/goa_human.gaf.gz'
    input_folder = 'output/'

    try:
        open(go_gaf_path)
    except:
        os.system(f'wget -N -P {input_folder} {go_gaf_url}')
        go_gaf_file = go_gaf_url.split('/')[-1]
        go_gaf_path = os.path.join(input_folder, go_gaf_file)
        os.system(f'gunzip {go_gaf_path}')


def map_protein_to_go():
    #os.system('wget -N -P input/ http://geneontology.org/gene-associations/goa_human.gaf.gz')
    #os.system('gunzip input/goa_human.gaf.gz')
    file = 'Protein_(UniProt)_2_GO_(GO).csv'
    bad_gos = set()

    with open(os.path.join('output/protein2go',file), 'w', newline='') as fin:
        writer = csv.writer(fin)
        writer.writerow(['Protein (UniProt)','GO (GO)','Relationship'])

        proteins_in_go, relations, go_terms, total = set(), dict(), set(), 0

        with open('input/goa_human.gaf') as f:
            for i, line in enumerate(f):
                if i > 40:
                    line = line.split('\t')
                    protein = line[1]
                    relation = line[3]
                    go_term = line[4]
                    try:
                        go_type = go_to_go_type[go_term]
                        go_term = go_type+':'+go_term.split(':')[1]

                        proteins_in_go.add(protein)
                        relations[relation] = relations.get(relation,0) + 1
                        go_terms.add(go_term)
                        total += 1

                        writer.writerow(['UniProt:'+protein, go_term, relation])
                    except:
                        bad_gos.add(go_term)

    print('Total Protein-GO Relationship:', total)
    print('Proteins w/GO Term:', len(proteins_in_go))
    print('GO Terms:', len(go_terms))
    print('Protein-GO Relationships:', len(relations))

    df = pd.read_csv(os.path.join('output/protein2go',file))
    df.to_csv(os.path.join('output/edges',file), index=False)
    df.to_csv(os.path.join('output/edges_to_use/',file), index=False)

    json.dump(list(proteins_in_go), open('output/protein2go/proteins_in_go.json','w'))    
    

if __name__ == '__main__':
    go_type_to_go, go_to_go_type = map_go_to_go_type()
    download_go_annotation_file_gaf()
    map_protein_to_go()