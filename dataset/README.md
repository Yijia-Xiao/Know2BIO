# Know2BIO Dataset Construction

Scripts for constructing Know2BIO dataset. In order, run `compound_to_compound_alignment`, `compound_to_protein_compound_to_pathway`, `protein_to_pathway`, `protein_to_reaction`, `gene_to_protein`. The other scripts have no execution dependencies (files to download, intermediate files, etc.), they can be executed in any order.

Below is the datasheet for Know2BIO.

- [Motivation](#motivation)
- [Composition](#composition)
- [Collection process](#collection-process)
- [Preprocessing/cleaning/labeling](#preprocessingcleaninglabeling)
- [Uses](#uses)
- [Distribution](#distribution)
- [Maintenance](#maintenance)

## Motivation
Knowledge Graph Benchmark of Biomedical Instances and Ontologies (Know2Bio) is a comprehensive and evolving general-purpose heterogeneous knowledge graph encompassessing data from 29 diverse sources, representing 11 biomedical categories and capturing intricate relationships. The current version of Know2BIO consists of 216,000 nodes and 6,500,000 edges, and it can be automatically updated to ensure its relevance and currency. 

### For what purpose was the dataset created?
Know2Bio was created as a general-purpose biomedical knowledge graph. It is intended to be used as a real-world benchmark dataset for knowledge graph representtion learning models.

### Who created the dataset (e.g., which team, research group) and on behalf of which entity (e.g., company, institution, organization)?
Know2Bio was created by a team of researchers at the Scalable Analytics (ScAi) Laboratory at University of California, Los Angeles (UCLA) and Tsinghua University.

### Who funded the creation of the dataset?
This work was supported by National Science Foundation (NSF) 1829071, 2031187, 2106859, 2119643, 2200274, 2211557 to W.W., Research Awards from Amazon and NEC to W.W., National Institutes of Health (NIH) R35 HL135772 to P.P., NIH T32 HL13945 to A.R.P. and D.S., NIH T32 EB016640 to A.R.P., NSF Research Traineeship (NRT) 1829071 to A.R.P. and D.S., and the TC Laubisch Endowment to P.P. at UCLA.


## Composition
Know2Bio integrates data from 29 diverse sources, capturing intricate relationships across 11 biomedical categories. It currently consists of 216,000 nodes and 6,500,000 edges. Biomedical types are anatomy, biological process, cellular component, compounds/drugs, disease, drug class, genes, molecular function, pathways, proteins, and reactions.

### What do the instances that comprise the dataset represent (e.g., documents, photos, people, countries)?
The data set is available in three formats: 1) as raw input files (.csv) detailing individually extracted biomedical knowledge via API and downloads. These files also include intermediate files for mapping between ontologies as well as node features (e.g., text descriptions, sequence data, structure data), and edge weights which were not included in the combined dataset as they were not included in model evaluation. 2) a combined KG following the head-relation-tail (h,r,t) convention, as a comma-separated text file. These KG are released for the ontology view, instance view, and bridge view, as well as a combined whole KG. 3) To facilitate benchmark comparison between different KG embedding models, we also release the train, validation, and test split KGs. 

### How many instances are there in total (of each type, if appropriate)?
The current version of Know2BIO consists of 216,000 nodes and 6,500,000 edges across 11 biomedical categories. Tables describing detailed description of node and edge counts will be added soon. 

### Does the dataset contain all possible instances or is it a sample (not necessarily random) of instances from a larger set?
The dataset contains all possible instances. Individual input files, test/train/validation splits, and a subset of sampled dataset are also available.

### What data does each instance consist of?
Data are released as processed triples (head, relation, tail) as well as raw text description and sequence data.

### Is there a label or target associated with each instance?
There is no label or target associated with each instance.

### Is any information missing from individual instances?
Triples involving DrugBank are removed in this release due to licensing release issues. Mechanism to access full dataset will be released soon.

### Are relationships between individual instances made explicit (e.g., users’ movie ratings, social network links)?
Relationships between instances are explicit, as represented in the knowledge graph format.

### Are there recommended data splits (e.g., training, development/validation, testing)?
To facilitate benchmark comparison between different KG embedding models, we also release the train, validation, and test split KGs. 
These KG are released for the ontology view, instance view, and bridge view, as well as a combined whole KG. The resulting KG was split into a train and test KG. The largest strongly connected component from the main KG was used to generate the train and test split, designating as close to 30\% in the test set across all nodes as possible, for each view. The remaining portion is divided into train and validation sets with a ratio of 9:1. Table enumerating the dataset split to be added.

### Are there any errors, sources of noise, or redundancies in the dataset?
To construct our KG, we integrate data from 29 biomedical data sources spanning several biomedical disciplines. This data integration required careful selection of data sources from which we extracted data. It also required data identifiers (IDs) to be mapped to common IDs through various intermediary resources. This is needed to unify knowledge on each biomedical entity/concept because they (e.g., genes) are often represented by different IDs (e.g., gene name; IDs from Entrez, Ensembl, HGNC) in different data sources. However, this process can be circuitous. For example, to unify knowledge on compounds and the proteins they target (i.e., *Compound (DrugBank ID) -targets- Protein (UniProt ID)*) taken from the Therapeutic Target Database (TTD), the following relationships are aligned: *Compound (TTD ID) -targets- Protein (TTD ID)* from TTD, *Protein (TTD ID) -is- Protein (UniProt name)* from UniProt, and *Protein (UniProt name) -is- Protein (UniProt ID)*  from UniProt. This creates *Compound (TTD ID) -targets- Protein (UniProt ID)* edges. But to unify this with the same compounds represented by DrugBank IDs elsewhere in the KG, the following relationships are aligned: *Compound (DrugBank ID) -is- Compounds (old TTD, CAS, PubChem, and ChEBI IDs)* from DrugBank (4 relationships), and *Compounds (CAS, PubChem, and ChEBI) -is- Compound (new TTD)* from TTD (3 relationships).

Relationships are also backed by varying levels of evidence (e.g., confidence scores for protein-protein interactions from STRING, gene-disease associations from DisGeNET). To select appropriate thresholds for inclusion in our KG, we investigate how confidence scores are calculated, what past researchers have selected, KB author recommendations, and resulting data availability. 

### Is the dataset self-contained, or does it link to or otherwise rely on external resources (e.g., websites, tweets, other datasets)?
The dataset is self-contained. The source code used to generate the dataset link to external resources. Table describing websites and APIs used to be added.

### Does the dataset contain data that might be considered confidential (e.g., data that is protected by legal privilege or by doctor-patient confidentiality, data that includes the content of individuals’ non-public communications)?
No.

### Does the dataset contain data that, if viewed directly, might be offensive, insulting, threatening, or might otherwise cause anxiety?
No.

### Does the dataset relate to people? 
The dataset relates to biomedical research, human health, and disease.

### Does the dataset identify any subpopulations (e.g., by age, gender)?
No.

### Is it possible to identify individuals (i.e., one or more natural persons), either directly or indirectly (i.e., in combination with other data) from the dataset?
No.

### Does the dataset contain data that might be considered sensitive in any way (e.g., data that reveals racial or ethnic origins, sexual orientations, religious beliefs, political opinions or union memberships, or locations; financial or health data; biometric or genetic data; forms of government identification, such as social security numbers; criminal history)?
No.


## Collection process
To construct our KG, we integrate data from 29 biomedical data sources spanning several biomedical disciplines. This data integration required careful selection of data sources from which we extracted data. It also required data identifiers (IDs) to be mapped to common IDs through various intermediary resources. This is needed to unify knowledge on each biomedical entity/concept because they (e.g., genes) are often represented by different IDs (e.g., gene name; IDs from Entrez, Ensembl, HGNC) in different data sources. However, this process can be circuitous. For example, to unify knowledge on compounds and the proteins they target (i.e., *Compound (DrugBank ID) -targets- Protein (UniProt ID)*) taken from the Therapeutic Target Database (TTD), the following relationships are aligned: *Compound (TTD ID) -targets- Protein (TTD ID)* from TTD, *Protein (TTD ID) -is- Protein (UniProt name)* from UniProt, and *Protein (UniProt name) -is- Protein (UniProt ID)*  from UniProt. This creates *Compound (TTD ID) -targets- Protein (UniProt ID)* edges. But to unify this with the same compounds represented by DrugBank IDs elsewhere in the KG, the following relationships are aligned: *Compound (DrugBank ID) -is- Compounds (old TTD, CAS, PubChem, and ChEBI IDs)* from DrugBank (4 relationships), and *Compounds (CAS, PubChem, and ChEBI) -is- Compound (new TTD)* from TTD (3 relationships).

Relationships are also backed by varying levels of evidence (e.g., confidence scores for protein-protein interactions from STRING, gene-disease associations from DisGeNET). To select appropriate thresholds for inclusion in our KG, we investigate how confidence scores are calculated, what past researchers have selected, KB author recommendations, and resulting data availability.  

### How was the data associated with each instance acquired?
The data sources include various databases, knowledge bases, API services, and knowledge graphs: MyGene.info, MyChem.info, MyDisease.info, Bgee, KEGG, PubMed, MeSH, SIDER, UMLS, CTD, PathFX, DisGeNET, TTD, Hetionet, Uberon, Mondo, PharmGKB, DrugBank, Reactome, DO, ClinGen, ClinVar, UniProt, GO, STRING, InxightDrugs, SMPDB, HGNC, and GRNdb. 

### What mechanisms or procedures were used to collect the data (e.g., hardware apparatus or sensor, manual human curation, software program, software API)?
See above. Manual human inspection of API result was performed to detect errors and successful download.

### If the dataset is a sample from a larger set, what was the sampling strategy (e.g., deterministic, probabilistic with specific sampling probabilities)?

### Who was involved in the data collection process (e.g., students, crowdworkers, contractors) and how were they compensated (e.g., how much were crowdworkers paid)?
Data collection was performed by graduate students at UCLA.

### Over what timeframe was the data collected?
Creation of the data collection pipeline took place over a year and a half from publication and release; new dataset was generated a few months before release and will be periodically updated. 

### Were any ethical review processes conducted (e.g., by an institutional review board)?
No.

### Does the dataset relate to people?
No.

### Did you collect the data from the individuals in question directly, or obtain it via third parties or other sources (e.g., websites)?

### Were the individuals in question notified about the data collection?

### Did the individuals in question consent to the collection and use of their data?

### If consent was obtained, were the consenting individuals provided with a mechanism to revoke their consent in the future or for certain uses?

### Has an analysis of the potential impact of the dataset and its use on data subjects (e.g., a data protection impact analysis) been conducted?


## Preprocessing/cleaning/labeling
Raw data and intermediate files downloaded from API and other sources are included with the submission. Resulting files are represented as knowledge graph triples (i.e., head, relation, tail).

### Was any preprocessing/cleaning/labeling of the data done (e.g., discretization or bucketing, tokenization, part-of-speech tagging, SIFT feature extraction, removal of instances, processing of missing values)?
Namespace identifiers indicating which biomedical resource it originated from, were appended to node names (e.g., "Reactome:" was appended to all Reactome Pathways and Reactions). Relations with 0 weight were removed. Relation names were given to match knowledge graph triple convention.

### Was the “raw” data saved in addition to the preprocessed/cleaned/labeled data (e.g., to support unanticipated future uses)?
Yes. It is included.

### Is the software used to preprocess/clean/label the instances available?
Yes. It is included.


## Uses

### Has the dataset been used for any tasks already?
We have performed a benchmark for knowledge graph representation learning models.

### Is there a repository that links to any or all papers or systems that use the dataset?
If the accompanying paper is accepted, these can be found as papers which cite the publication.

### What (other) tasks could the dataset be used for?
This knowledge graph is a general-purpose knowledge graph for biomedical knowledge. It can be used for identifying drug targets and therapeutics, discovering biomarkers for disease, among other purposes.

### Is there anything about the composition of the dataset or the way it was collected and preprocessed/cleaned/labeled that might impact future uses?
This knowledge graph framework is easily extensible to other currently not included biomedical knowledge bases. In the preparation of this dataset, there is no unfair treatment of individuals or groups or other undesirable harms beyond what is already represented in the source biomedical knowledge bases. To our knowledge, we have not included knowledge bases of consern for these harms.

### Are there tasks for which the dataset should not be used?
This knowledge graph is intended for research purposes and should not be used for medical advise and clinical decision making without consulting a medical professional.


## Distribution
### Will the dataset be distributed to third parties outside of the entity (e.g., company, institution, organization) on behalf of which the dataset was created?
This dataset is available for download and reuse under the MIT license.

### How will the dataset will be distributed (e.g., tarball on website, API, GitHub)?
If this paper is accepted, we will release the DOI to access the dataset.

### When will the dataset be distributed?
This dataset will be distributed upon paper acceptabce.

### Will the dataset be distributed under a copyright or other intellectual property (IP) license, and/or under applicable terms of use (ToU)?
MIT license.

### Have any third parties imposed IP-based or other restrictions on the data associated with the instances?
A few data sources impose licensing restrictions and will not be included with our dataset release. Instructions on how to access license and retrieve the complete dataset will be released soon.

### Do any export controls or other regulatory restrictions apply to the dataset or to individual instances?
No.


## Maintenance

### Who is supporting/hosting/maintaining the dataset?
Our research group at the University of California, Los Angeles will support and maintain the dataset.

### How can the owner/curator/manager of the dataset be contacted (e.g., email address)?
Please contact Prof. Wei Wang (weiwang@cs.ucla.edu)

### Is there an erratum?
Not currently.

### Will the dataset be updated (e.g., to correct labeling errors, add new instances, delete instances)?
The dataset and accompanying code will be periodically updated and communicated through the project GitHub.

### If the dataset relates to people, are there applicable limits on the retention of the data associated with the instances (e.g., were individuals in question told that their data would be retained for a fixed period of time and then deleted)?
No.

### Will older versions of the dataset continue to be supported/hosted/maintained?
We will maintain versioning of new release of dataset.

### If others want to extend/augment/build on/contribute to the dataset, is there a mechanism for them to do so?
Users are welcome to extend, augment, build, and contribute to the dataset, per mechanisms on GitHub.
