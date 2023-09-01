# Message to NeurIPS reviewers: We mistakenly submitted our latest revision to the "supplemental materials" instead of the main "pdf" location. 
Message to NeurIPS reviewers: We mistakenly submitted our latest revision to the "supplemental materials" instead of the main "pdf" location. 


**anatomy_to_anatomy**
The official xml file from MeSH was used to map anatomy MeSH IDs and MeSH tree numbers to each other, as well as MeSH tree numbers to each other to form the hierarchical relationships in the ontology. MeSH IDs were aligned to Uberon IDs via the official Uberon obo file (used in gene-to-anatomy) 

**compound_to_compound**
The compound_to_compound_alignment script aligned numerous compound identifiers in order to align DrugBank and MeSH IDs, two of the most prevalent IDs from data sources for different relationships in the scripts here. To produce this file, numerous resources were used to directly map DrugBank to MeSH IDs or to indirectly align the IDs (e.g., via DrugBank to UNII from DrugBank, then UNII to MeSH via MyChem.info). Resources include UMLS's MRCONSO.RRF file, DrugBank, MeSH, MyChem.info, the NIH's Inxight Drugs, KEGG, and TTD. In other scripts, DrugBank and MeSH compounds are mapped to one another via this mapping file.

Compound interactions were extracted from DrugBank. 

**compound_to_disease**
The majority of the compound-treats-disease and compound-biomarker_of-disease edges were from the Comparative Toxicogenomics Database. Additional edges were from PathFX (i.e., from repoDB) and Hetionet (reviewed by 3 physicians). 

**compound_to_drug_class**
Mappings from compounds to drug classes (ATC) were provided by DrugBank. 

**compound_to_gene**
Mapping compound to gene largely relies on CTD, though some relationships come from KEGG. Like many other compound mappings, this relies on the DrugBank-to-MeSH alignments from compound_to_compound_alignment.

**compound_to_pathway**
Mapping compounds to SMPDB pathways relies on DrugBank. Mapping compounds to Reactome pathways relies on Reactome, plus alignments to ChEBI compounds. Mapping compounds to KEGG pathways relies on KEGG. 

**compound_to_protein**
Most compound-to-protein relationships are from DrugBank. Some are taken from TTD, relying on mappings provided by TTD and aligning identifiers based on DrugBank- and TTD-provided identifiers.  

**disease_to_disease**
The official xml file from MeSH was used to map disease MeSH IDs and MeSH tree numbers to each other, as well as MeSH tree numbers to each other to form the hierarchical relationships in the ontology. 

To measure disease similarity, edges were obtained from DisGeNET's curated data. The UMLS-to-MeSH alignment was used (from compound_to_compound_alignment). 

Disease Ontology was used to align Disease Ontology to MeSH. Mondo and MyDisease.info were relied on to align Mondo to MeSH, DOID, OMIM, and UMLS. These alignments were used to align relationships from other scripts to the MeSH disease identifiers.

**compound_to_side_effect**
Mappings from compounds to the side effects they are associated with were provided by SIDER. This required alignments from PubChem to DrugBank (provided by DrugBank) and UMLS to MeSH (provided in compound_to_compound_alignment.py).

**disease_to_anatomy**
Disease and anatomy association mappings rely on MeSH for aligning the MeSH IDs and MeSH tree numbers and rely on the disease-anatomy coocurrences in PubMed articles' MeSH annotations.

**disease_to_pathway**
KEGG was used to map KEGG pathways to disease. Reactome was used to map Reactome pathways to diseases, relying on the DOID-to-MeSH alignments for disease.

**gene_to_anatomy**
Gene expression in anatomy was derived from Bgee. To align the Bgee-provided Ensembl gene IDs to Entrez, MyGene.info was used. To align the Bgee-provided Uberon anatomy IDs to MeSH, Uberon was used (see anatomy_to_anatomy)

**gene_to_disease**
Virtually all gene-disease associations were obtained from DisGeNET's entire dataset. Additional associations---many of which were already present in DisGeNET---were obtained from ClinVar, ClinGen, and PharmGKB. (Users may be interested in only using the curated evidence from DisGeNET or increasing the confidence score threshold for DisGeNET gene-disease association. We chose a threshold of 0.06 based on what a lead DisGeNET author mentioned to the Hetionet creator in a forum. 

**gene_to_protein**
We relied on UniProt and HGNC to map proteins to the genes that encode them. Notably, there is a very large overlap between these sources (~95%). HGNC currently broke, so only UniProt is being used. 

**go_to_go**
The source of the Gene Ontology ontologies is Gene Ontology itself. 

**go_to_protein**
The source of the mappings between proteins and their GO terms is Gene Ontology. 

**pathway_to_pathway**
The source of pathway hierarchy mappings for KEGG is KEGG and for Reactome is Reactome. (SMPDB does not have a hierarchy)

**protein_and_compound_to_reaction**
The source of mappings from proteins and compounds to reactions is Reactome. This file relies on alignments from ChEBI to DrugBank.

**protein_and_gene_to_pathway**
To map proteins and genes to pathways, KEGG was used for KEGG pathways (genes), Reactome for Reactome pathways (proteins and genes), and SMPDB for SMPDB pathways (proteins). 

**protein_to_gene_ie_transcription_factor_edges**
To map the proteins (i.e., transcription factors) to their targeted genes (i.e., the proteins that affect expression of particular genes), GRNdb's high confidence relationships virtually all derived from GTEx, were used. This also required aligning gene names to Entrez gene IDs through MyGene.info

**protein_to_protein**
Protein-protein interactions (i.e., functional associations) were derived from STRING. To map the STRING protein identifiers to UniProt, the UniProt API was used. A confidence threshold of 0.7 was used. (Users may adjust this in the script)

**reaction_to_pathway**
To map reactions to the pathways they participate in, Reactome was used.

**reaction_to_reaction**
To map reactions to reactions that precede them, Reactome was used. 
