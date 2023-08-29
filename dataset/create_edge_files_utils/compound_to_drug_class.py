import json
import csv
import pandas as pd
import os
from compound_to_compound_alignment import extract_drugbank_xml, parse_drugbank_xml, map_drugbank_to_external_ids
from biomedkg_utils import output_edgefile_onerel_noweight

def map_compound_to_drug_class():
    try:
        db_to_atc = json.load(open('output/compound2compound/db2atc.json'))
    except:
        extract_drugbank_xml()
        root = parse_drugbank_xml()
        dbs = map_drugbank_to_external_ids(root)  
        del root
        db_to_atc = json.load(open('output/compound2compound/db2atc.json'))
            
    output_edgefile_onerel_noweight(
        outpath='output/compound2drugclass/Compound_(DrugBank)_to_DrugClass_(ATC).csv',
        columns=['Compound (DrugBank)', 'Drug Class (ATC)', 'Relationship'],
        dictionary=db_to_atc,
        rel='-compound_classified_as_drug_class->',
        prefix_col1='DrugBank_Compound:',
        prefix_col2='ATC_Class:',
        )
    
    
def map_drug_class_to_drug_class():
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
        writer.writerow(['Drug Class (ATC)','Drug Class (ATC)','Relationship'])
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

            writer.writerow(['ATC_Class:'+second_level, 'ATC_Class:'+first_level, rel])
            writer.writerow(['ATC_Class:'+third_level, 'ATC_Class:'+second_level, rel])
            writer.writerow(['ATC_Class:'+fourth_level, 'ATC_Class:'+third_level, rel])
            writer.writerow(['ATC_Class:'+fifth_level, 'ATC_Class:'+fourth_level, rel])

    atc_df = pd.read_csv('output/compound2drugclass/atc_class_tree.csv').drop_duplicates()
    atc_df.reset_index(inplace=True, drop=True)
    atc_df.to_csv('output/compound2drugclass/atc_class_tree.csv')

    print(len(set(atc_df['Drug Class (ATC)']).union(set(atc_df['Drug Class (ATC).1']))), 'ATC-ATC nodes')
    print(len(atc_df), 'ATC edges')
    atc_df.to_csv('output/compound2drugclass/atc_class_tree.csv',index=False)
    atc_df.to_csv('output/edges/edges_atc_to_atc.csv', index=False)
    atc_df.to_csv('output/edges_to_use/DrugClass_(ATC)_2_DrugClass_(ATC).csv',index=False)

if __name__ == '__main__':
    map_compound_to_drug_class()
    map_drug_class_to_drug_class()