import json
import csv
import pandas as pd
import os
from compound_to_compound_alignment import extract_drugbank_xml, parse_drugbank_xml, map_drugbank_to_external_ids

def map_compound_to_drug_class():
    try:
        db_to_atc = json.load(open('output/compound2compound/db2atc.json'))
    except:
        extract_drugbank_xml()
        root = parse_drugbank_xml()
        dbs = map_drugbank_to_external_ids(root)  
        del root
        db_to_atc = json.load(open('output/compound2compound/db2atc.json'))
        
    atc_codes_in_db = [atc for atc_list in db_to_atc.values() for atc in atc_list]

    with open('output/compound2drugclass/atc_class_tree.csv','w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['ATC Child','ATC Parent','Relationship'])
        rel = 'is_subclass_of->'

        for atc in atc_codes_in_db:
            assert len(atc) == 7
            first_level = atc[:1]  # main anatomical group
            second_level = atc[:3] # therapeutic subgroup
            third_level = atc[:4]  # pharmacological subgroup
            fourth_level = atc[:5] # chemical subgroup
            fifth_level = atc      # chemical substance 
            atc_list = [first_level, second_level, third_level, 
                        fourth_level, fifth_level]

            writer.writerow(['ATC:'+second_level, 'ATC:'+first_level, rel])
            writer.writerow(['ATC:'+third_level, 'ATC:'+second_level, rel])
            writer.writerow(['ATC:'+fourth_level, 'ATC:'+third_level, rel])
            writer.writerow(['ATC:'+fifth_level, 'ATC:'+fourth_level, rel])

    atc_df = pd.read_csv('output/compound2drugclass/atc_class_tree.csv').drop_duplicates()
    atc_df.reset_index(inplace=True, drop=True)
    atc_df.to_csv('output/compound2drugclass/atc_class_tree.csv')

    print(len(set(atc_df['ATC Child']).union(set(atc_df['ATC Parent']))), 'ATC-ATC nodes')
    print(len(atc_df), 'ATC edges')
    atc_df.to_csv('output/compound2drugclass/atc_class_tree.csv',index=False)
    atc_df.to_csv('output/edges/edges_atc_to_atc.csv', index=False)
    atc_df.to_csv('output/edges_to_use/DrugClass_(ATC)_2_DrugClass_(ATC).csv',index=False)

if __name__ == '__main__':
    map_compound_to_drug_class()