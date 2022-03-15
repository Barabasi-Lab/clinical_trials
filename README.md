# Clinical Trials 

## Python version

- **Version**: v1
- **Creator**: [Kishore Vasan] [viprankishore\@gmail.com](mailto:viprankishore@gmail.com) 
- **Last code update**: 03/14/2022 
- **Last data update**: 11/01/2020 (date the xml files were extracted)
- **Keywords**: Clinical Trials; Drug Innovation
- **Rights Statement**: Open Data


#### Generated files

- data/out/organized_ct_data -- clinical trials data 
-   ID \<chr>: Drug ID in DrugBank
-   Name \<chr>: Official drug name in DrugBank

- data/out/drug_mapped_ct_data -- clinical trials mapped drug data
- 


#### Functions 

#### `extract_ct_xml.py`

`organize_data(XML)` - parses the XML files to extract ct data and saves to organized_ct_data.csv
    - nct_id \<chr> : clinical trial id map of the trial
    - title \<chr> : title of the clinical trial
    - study_type \<chr> : type of clinical trials (interventional, behavioral etc.)
    - gender \<chr> : gender involved in the trail
    - min_age \<chr> : minimum age of the trial participants
    - max_age \<chr> : maximum age of the trial participants
    - status \<chr> : status of the trial (completed, recruiting etc.)
    - phase \<chr> : phase of the trial (phase 1, phase 2 etc.)
    - start_date \<chr> : start date of the trial 
    - location_countries \<chr> : location of the trials 
    - conditions \<chr> : disease conditions tested in the trial
    - keywords \<chr> : keywords involved in the trial
    - mesh_terms \<chr> : mesh terms of the trial 
    - results_pubs_pmid \<chr> : list of publications on the tirial
    - references_pmid \<chr> : references provided by the trial

` organize_funder_data()` - parses the funder data and saves it to `../data/out/funder_ct_data.csv`
   - nct_id \<chr>: clinical trial id map of the trial
   - funder_name \<chr>: name of the funder
   - funder_type \<chr>: type of funder (government, industry etc.)
   - funder_role \<chr> : role of the funder (lead or collaborator)

`organize_intervention_data()`- parses the intervention data and it to `../data/out/intervention_ct_data.csv`
   - nct_id \<chr> : clinical trial id map of the trial
    - intervention \<chr> : list of drugs/ products tested in the trial
    - intervention_type \<chr> : types of intervention (drug, product etc.)
   

#### `curate_drug_interventions_ct.py`

Input files : __ save the following files in the ../data/raw folder
1) intervention_ct_data.csv -- contains the xml parsed list of interventions in trials
2) all_drugbank_drugs.csv - parsed database containing drugbank drug id and names
3) drug_synonym.csv - containing drugbank id, name, and its corresponding synonym
4) products.csv -- containing drugbank id, name, and its corresponding products
5) drugs_external_identifiers.csv -- containings the drugbank id, name, and its corresponding external identifiers

Output: `drug_mapped_ct_data.csv`
   - nct_id \<chr> : clinical trial id map of the trial
   - Name \<chr> : official DrugBank name in lowercase 
   - intervention_type \<chr> : drug 


#### methodology: 

This file loads the intervention data and maps it through a five step process. 

1) Search for direct text matching with the drug name in DrugBank

2) Search for matching with synonyms of the drug names

3) Map the intervention names with the product names of the drugs

4) map the intervention names to external identifier (e.g. wikipedia)

5) fuzzy string match of the names with the drugbank names

   
#### Running the parser

- The latest XML data of all clinical trials can be downloaded from clinicaltrials.gov -- save it to /data/raw folder
- First, run the `extract_ct_xml.py` 
- Then run the `curate_drug_interventions_ct.py` files, 
- Finally, `read_data.py` curates drug data and clinical trials data 
- Data will be saved in the /data/out folder

### Data Stats:
```
Number of Trials: 356403
--
Number of Drug Trials: 127432
Proportion of drug trials mapped: 0.8709487813879738
Number of Interventions: 5694
--
drugbank...
Number of drugs: 6316
Number of targets: 3115
--
clinical trials...
Number of Targets: 2714
loading ppi network
Number of genes: 18508
Number of interactions: 332646
Name:
Type: Graph
Number of nodes: 18508
Number of edges: 326883
Average degree:  35.3234
--
Number of placebo trials: 1171
Number of placebo drugs: 590
--
Num druggable genes: 1327
--
N approved drugs: 1005
Number of drugs in CT mapped with approval dates: 956
Number of targets in CT: 1340
--
```
