# Clinical Trials 

## Python version

- **Version**: v1
- **Creator**: [Kishore Vasan] [viprankishore\@gmail.com](mailto:viprankishore@gmail.com) 
- **Last code update**: 03/14/2022 
- **Last data update**: 11/01/2020 
- **Keywords**: Clinical Trials; Drug Innovation
- **Rights Statement**: Open Data


#### Generated files

- data/out/organized_ct_data -- clinical trials data 
-   ID \<chr>: Drug ID in DrugBank
-   Name \<chr>: Official drug name in DrugBank

- data/out/drug_mapped_ct_data -- clinical trials mapped drug data
- 


#### Functions 




#### Running the parser
- The latest XML data of all clinical trials can be downloaded from clinicaltrials.gov
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
