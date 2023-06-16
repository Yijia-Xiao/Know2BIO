# Modified By Thai Tran and provided by Ping lab, modofied by me (Dylan Steinecke) https://github.com/CaseOLAP/covid-cvd-knowledgegraph/blob/main/01_DrugBankXML_CVDrugs_Extraction/parse_xml.py
   
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
    
    
    def getInteractingDrugs(ele):
        '''
        FUNCTION:
        - Find the DrugBank drugs that interact with this Drugbank drug

        PARAMS:
        - The xml node(?) of the main drug entry
        '''
        interacting_drugs = set()

        # Find the drug interactions section
        drug_interactions_section = ele.findall('{http://www.drugbank.ca}drug-interactions')
        for drug_int_section in drug_interactions_section:

            for drug_int in drug_int_section.iter():

                # Find the interacting drugs' sections
                for drug_int_el in drug_int.iter():

                    # Find the interacting drug IDs
                    drug_ids = drug_int_el.findall('{http://www.drugbank.ca}drugbank-id')
                    if len(drug_ids) > 0:
                        for drug_id in drug_ids:
                            interacting_drugs.add(drug_id.text)

        return interacting_drugs
