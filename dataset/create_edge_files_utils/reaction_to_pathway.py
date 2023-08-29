import requests
from biomedkg_utils import output_edgefile_onerel_noweight

def download_reactome_event_tree():
    url = 'https://reactome.org/ContentService/data/eventsHierarchy/9606?pathwaysOnly=false&resource=TOTAL&interactors=false'
    response = requests.get(url)
    response_dict = response.json()
    return response_dict


def map_reactome_event_tree(parent_entry, reaction_to_pathway):
    parent_id = parent_entry['stId']
    try:
        children_list = parent_entry['children']
    except:
        children_list = []
    for child_entry in children_list:
        child_id = child_entry['stId']
        
        # Map reactions to pathways
        child_is_reaction = ('reaction' in child_entry['type'].lower())
        parent_is_pathway = ('pathway' in parent_entry['type'].lower())
        if child_is_reaction and parent_is_pathway:
            reaction_to_pathway.setdefault(child_id, set()).add(parent_id)
        map_reactome_event_tree(child_entry, reaction_to_pathway)
        

def map_reactome_reaction_to_pathway():
    reaction_to_pathway = {}
    roots = download_reactome_event_tree()
    for root_entry in roots:
        map_reactome_event_tree(root_entry, reaction_to_pathway)
    output_edgefile_onerel_noweight(
        outpath='output/pathway2reaction/Pathway_(Reactome)_to_Reaction_(Reactome).csv',
        columns=['Pathway (Reactome)', 'Reaction (Reactome)', 'Relationship'],
        dictionary=reaction_to_pathway,
        rel='-reaction_involved_in_pathway->',
        prefix_col1='Reactome_Reaction:',
        prefix_col2='Reactome_Pathway:',
        )
    
if __name__ == '__main__':
    map_reactome_reaction_to_pathway()
    print('Finished mapping reaction to pathway')