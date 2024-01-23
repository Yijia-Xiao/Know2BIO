# Know2BIO Dataset

Know2BIO dataset is periodically updated and released. The latest release is from 2023-08-18 and consists of 3 versions:
1. Full dataset: The full Know2BIO Knowledge Graph consisting of 219,169 nodes, 6,181,160 edges, and their node features. This dataset requires licensing restrictions to access. Instructions on obtaining these license will be released soon. This dataset can be obtained by obtaining the necessary licenses and filling out the below webform.
2. Safe release: The Know2BIO Knowledge Graph without licensing-restricted information. The resulting KG is 152,845 nodes, 3,282,063 edges, and their node features. This dataset can be obtained by filling out the below webform.
3. Sampled safe release: A 1% sample of the Know2BIO safe release, representing roughly 1% proportional to all edges. This dataset is accesible under the 'sampled_know2bio_safe_release' folder. An example of their node features is included in this directory under 'sampled_node_features.json'. The full set of node features for this sampled dataset is accesible at this [link](https://drive.google.com/file/d/1t3EFiMhI5CNW6erm93f0eQt15F5uQCPC/view?usp=drive_link).

For accessing the full dataset of the safe release dataset, please fill out this webform: https://forms.gle/3HdKRtvW7ce9PKpw6

Summary:

Manually download this data into the dataset/create_edge_files_utils/input folder:
```
- UMLS thesaurus "MRCONSO.RRF" (make a UMLS account): https://www.nlm.nih.gov/research/umls/licensedcontent/umlsknowledgesources.html
- Drugbank full database (make a DrugBank account): https://go.drugbank.com/releases/5-1-9/downloads/all-full-database
- DrugBank compound structures (make a Drugbank account): https://go.drugbank.com/releases/5-1-9/downloads/all-structure-links
```

Execute these commands: 
```
python create_edge_files.py

bash
python ./prepare_kgs/prepare_kgs.py ./output/edges_to_use/ ./know2bio ./input_lists/all_kg_edges.txt
python ./prepare_kgs/split_dataset.py ./know2bio 0.8
python ./prepare_kgs/prepare_benchmark.py ./know2bio ../../benchmark/data
cd ../../benchmark/data/K2BIO
python n-n.py
```

# How to Construct the Know2BIO Knowledge Graph (KG)
## Create the edge files
Execute ```python create_edge_files_utils.py``` to create all edges. 

Note: This executes scripts in a suitable order. Note that some scripts must be run in this order (e.g., `compound_to_compound_alignment` and `gene_to_protein` should be run first) but others can be run in different orders or can be not run if you don't want to create the edges for the script's respective edge types. Alternatively, you can still create the edges but just choose to not use the edge file produced by the script. 

## Prepare the Knowledge Graph Files

Following construction of the Know2BIO dataset detailed above, the edge files must be assembled into a knowledge graph and prepared following a specific data format for use with our benchmark knowledge graph representation learning models.


## Assemble the Knowledge Graph from the Files
Scripts for constructing the knowledge graph for Know2BIO dataset are found in the `prepare_kgs` folder.

Set up the environment to run the script, using `pip install -r prepare_kgs_requirements.txt`

The script takes as input the directory where Know2BIO is downloaded, the output directory, and a list of all input files separated by which kg it is a part of. The kgs are separated into test/train/validation set following 80/10/10 split and are processed in the correct format for benchmarking.

Run `kg_sampler.py` to sample a subset of the knowledge graph for testing purposes.

### Step 0: (Optional) Only Assemble Part of the Knowledge Graph.
Including or excluding certain data files from Know2BIO allows the construction of a use-case specific knowledge graph. This is achieved by providing a tailored file list for the knowledge graph construction. Included within this repository is the `input_list` folder, specifying specific edge lists to be used for knowledge graph construction. The `ont_bridge_inst_list.txt` file specifies all edge_files generated from constructing Know2BIO. This text file specifies which files are considered as part of each view (i.e., instance view, ontology view, bridge view). To make a tailored file, we recommend copying this file and deleting the file names which should not be included.

An example of specific use cases is included within the `input_lists` folder and detailed below:
- A protein-protein interaction knowledge graph: This protein-centric knowledge graph includes only known protein-protein interactions, their relation to genes, and their relations to biological pathways. No disease or drug information is included, as their inclusion could hinder prediction of protein-protein interactions. An example input list for this knowledge graph is `protein_protein_interaction_kg_list.txt`
- A drug-target interaction knowledge graph: This knowledge graph focuses specifically on protein-drug interactions between DrugBank drugs, UniProt proteins, MeSH diseases, and Reactome biological pathways. This KG does not include data files from different knowledge bases which serve a similar purpose (e.g., using Reactome, instead of Reactome and KEGG; using DrugBank only, instead of DrugBank and MeSH compounds) to reduce complexity of the KG when integrating from multiple sources. An example input list for this knowledge graph is `drug_target_interaction_kg_list.txt`

Furthermore, these KGs can be constructed for testing predictive power of biologically relevant triples (e.g., train a model using entire Know2BIO and predict only on protein-drug edges) will be detailed soon.

### Step 1: Combine all individual edge files into a combined Knowledge Graph.
The `prepare_kgs.py` script prepares the input data from `output/edges_to_use` and splits the dataset into different views (i.e., instance view, ontology view, bridge view, and the aggregate view). These four combined knowledge graphs will be output to the `know2bio` folder. Run the below command to execute.

```bash
python ./prepare_kgs/prepare_kgs.py ./output/edges_to_use/ ./know2bio ./input_lists/ont_bridge_inst_list.txt
```

## Generate Train, Test, Validation datasets

### Step 2: Split the dataset for training and evaluation.
The `split_dataset.py` script splits the knowledge graphs constructed from the previous step, into separate knowledge graphs for training and evaluation. The input folder (with `kg1f_instances.txt`,`kg2f_ontologies.txt`, and `alignf_bridges.txt`) is provided as the first argument, with the proportion of the KG for training as the second argument. The train graph is required to span all nodes in the knowledge graph, for every connected component with greater than 10 nodes. 

```bash
python ./prepare_kgs/split_dataset.py ./know2bio 0.8
```

## Prepare Benchmark Dataset

### Step 3: Reformat the data for Benchmark modeling.
The `prepare_benchmark.py` script automatically reformats the resulting aggregate, instance, and ontology view data into the data folder of the benchmark section of the repository. This generates the `entity2id.txt`, `relation2id.txt`, `whole2id.txt`, `train2id.txt`, `test2id.txt`, and `valid2id.txt` files in the corresponding data folder. The `n-n.py` script must be executed to generate the remaining datafiles to sucessfully update the dataset for modeling.

```bash
python ./prepare_kgs/prepare_benchmark.py ./know2bio ../../benchmark/data
cd ../../benchmark/data/K2BIO
python n-n.py
```

## Prepare Multimodal Data (Optional)

### Step 4: Combine multimodal data.
The `prepare_multimodal_data.py` script automatically combines the resulting node features in the `multi-modal_data` folder as the `node_features.json` file. This file summarizes the node features for each node in the knowledge graph. To generate this file, make sure all parsed files are within the `multi-modal_data` folder and execute the below script.

```bash
python ./prepare_kgs/prepare_multimodal_data.py
```


# Datasheet for Know2BIO / Additional Details
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
Know2Bio was created as a general-purpose biomedical knowledge graph. It is intended to be used as a resource for biomedical discovery and a real-world benchmark dataset for knowledge graph representtion learning models.

### Who created the dataset (e.g., which team, research group) and on behalf of which entity (e.g., company, institution, organization)?
Due to review policy, author information will be released shortly.


### Who funded the creation of the dataset?
Due to review policy, funding details will be disclosed shortly.


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
Due to review policy, author information will be released shortly.

### Over what timeframe was the data collected?
Creation of the data collection pipeline took place over a year and a half from publication and release; new dataset was generated a few months before release and will be periodically updated. 

### Were any ethical review processes conducted (e.g., by an institutional review board)?
No.

### Does the dataset relate to people?
No.

### Did you collect the data from the individuals in question directly, or obtain it via third parties or other sources (e.g., websites)?
N.A.
### Were the individuals in question notified about the data collection?
N.A.
### Did the individuals in question consent to the collection and use of their data?
N.A.
### If consent was obtained, were the consenting individuals provided with a mechanism to revoke their consent in the future or for certain uses?
N.A.
### Has an analysis of the potential impact of the dataset and its use on data subjects (e.g., a data protection impact analysis) been conducted?
N.A.

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
Due to review policy, institute information will be released shortly.

### How can the owner/curator/manager of the dataset be contacted (e.g., email address)?
Due to review policy, author information will be released shortly.

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
