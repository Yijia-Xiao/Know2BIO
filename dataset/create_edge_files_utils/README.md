**Anatomy-to-Anatomy**
The official xml file from MeSH was used to map anatomy MeSH IDs and MeSH tree numbers to each other, as well as MeSH tree numbers to each other to form the hierarchical relationships in the ontology. MeSH IDs were aligned to Uberon IDs via the official Uberon obo file (used in gene-to-anatomy) 

**Compound-to-Compound**
The compound_to_compound_alignment script aligned numerous compound identifiers in order to align DrugBank and MeSH IDs, two of the most prevalent IDs from data sources for different relationships in the scripts here. To produce this file, numerous resources were used to directly map DrugBank to MeSH IDs or to indirectly align the IDs (e.g., via DrugBank to UNII from DrugBank, then UNII to MeSH via MyChem.info). Resources include UMLS's MRCONSO.RRF file, DrugBank, MeSH, MyChem.info, the NIH's Inxight Drugs, KEGG, and TTD. In other scripts, DrugBank and MeSH compounds are mapped to one another via this mapping file.

Compound interactions were extracted from DrugBank. 

**Compound-to-Disease**
The majority of the compound-treats-disease and compound-biomarker_of-disease edges were from the Comparative Toxicogenomics Database. Additional edges were from PathFX (i.e., from repoDB) and Hetionet (reviewed by 3 physicians). 

**Disease-to-Disease**
The official xml file from MeSH was used to map disease MeSH IDs and MeSH tree numbers to each other, as well as MeSH tree numbers to each other to form the hierarchical relationships in the ontology. 
