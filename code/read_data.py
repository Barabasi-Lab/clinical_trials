"""
#!/usr/bin/env python
Author= Kishore Vasan

# Description:
## This script reads all the curated clinical trials data
## run the load_data() method to create global variables
## Usage: import using from read_data import *
"""

### List of data files needed:
# ct_data/clean_data/organized_ct_data.csv
# ct_data/clean_data/drug_mapped_ct_data_final_2.csv
# drug_data/all_drugbank_drugs.csv
# drug_data/PPI_covid_paper.csv
# ct_data/clean_data/placebo_trials.csv
# drug_data/druggable_genome.tsv
# drug_data/drug_approved_mapping.csv

# import packages
import pandas as pd
import numpy as np
from tqdm import tqdm_notebook
import networkx as nx
import time

# file to load all the ct data
def load_data():
    global df, drug_df, db_target_df, dt_trial_df, dt_start_year
    global ppi_g
    global placebo_trials, druggable_genome_df
    global drug_approval_dates

    df = pd.read_csv("ct_data/clean_data/organized_ct_data.csv")
    print("Number of Trials:", df.nct_id.nunique())
    print("--")

    drug_df = pd.read_csv("ct_data/clean_data/drug_mapped_ct_data.csv")
    print("Number of Drug Trials:", drug_df.nct_id.nunique())
    print("Proportion of drug trials mapped:", float(drug_df.nct_id.nunique())/float(df[~(df.intervention_types.isna()) & (df.intervention_types.str.contains('Drug'))].nct_id.nunique()))
    print("Number of Interventions:", drug_df.intervention.nunique())
    print("--")

    db_target_df = pd.read_csv("drug_data/all_drugbank_drugs.csv")
    db_target_df['Name'] = db_target_df.Name.str.lower()
    db_target_df = db_target_df[db_target_df.organism == "Humans"]
    print("drugbank...")
    print("Number of drugs:", db_target_df.Name.nunique())
    print("Number of targets:", db_target_df.Gene_Target.nunique())
    print("--")

    dt_trial_df = pd.merge(drug_df, db_target_df[['Name','Gene_Target','gene_type','known_action']], how='left',left_on='intervention',right_on='Name')
    print("clinical trials...")
    print("Number of Targets:", dt_trial_df.Gene_Target.nunique())
    dt_trial_df = pd.merge(dt_trial_df, df[['nct_id','start_date','phase','conditions']], how='inner')

    dt_start_year = dt_trial_df
    dt_start_year = dt_start_year[~dt_start_year.start_date.isna()]

    # get the year
    dt_start_year['start_date'] = [i.split(' ')[0]+' '+i.split(' ')[-1] if ',' in i else i for i in dt_start_year.start_date.tolist()]
    dt_start_year['month'] = [i.split(' ')[0] for i in dt_start_year.start_date.tolist()]
    dt_start_year['year'] = [i.split(' ')[-1] for i in dt_start_year.start_date.tolist()]
    dt_start_year['year'] = dt_start_year.year.astype(int)

    # add phase
    dt_start_year = pd.merge(dt_start_year, df[['nct_id','phase']], how ='left')
    dt_start_year['year'] = dt_start_year.year.astype(int)

    print("loading ppi network")
    ppi_data = pd.read_csv("drug_data/PPI_covid_paper.csv")
    ppi_data = ppi_data[['Symbol_A','Symbol_B']]
    ppi_data.columns = ['protein_1','protein_2']
    print("Number of genes:", len(set.union(set(ppi_data.protein_1), set(ppi_data.protein_2))))
    print("Number of interactions:", ppi_data.drop_duplicates().shape[0])

    ppi_g = nx.from_pandas_edgelist(ppi_data, "protein_1", "protein_2", create_using = nx.Graph())
    ppi_g.remove_edges_from(nx.selfloop_edges(ppi_g))
    print(nx.info(ppi_g))
    print("--")

    placebo_trials = pd.read_csv("ct_data/clean_data/placebo_trials.csv")
    placebo_trials = placebo_trials.drop_duplicates()
    placebo_trials['placebo'] = True
    print("Number of placebo trials:", placebo_trials.nct_id.nunique())
    print("Number of placebo drugs:", placebo_trials.drug_map.nunique())
    print("--")

    dt_start_year = pd.merge(dt_start_year, placebo_trials[['nct_id','drug_map','placebo']], how='left', left_on=['nct_id','intervention'], right_on=['nct_id','drug_map'])
    dt_start_year['placebo'] = dt_start_year.placebo.fillna(False)

    dt_trial_df = pd.merge(dt_trial_df, placebo_trials[['nct_id','drug_map','placebo']], how='left', left_on=['nct_id','intervention'], right_on=['nct_id','drug_map'])
    dt_trial_df['placebo'] = dt_trial_df.placebo.fillna(False)

    druggable_genome_df = pd.read_csv("drug_data/druggable_genome.tsv", sep='\t')
    #druggable_genome_df = druggable_genome_df[druggable_genome_df.category == 'DRUGGABLE GENOME']
    druggable_genome_set = set(druggable_genome_df.entrez_gene_symbol.tolist())
    druggable_genome_set = druggable_genome_set - set(list(ppi_g.nodes()))
    print("Num druggable genes:", len(druggable_genome_set))
    print("--")

    drug_approval_dates = pd.read_csv("drug_data/drug_approved_mapping.csv")
    print("N approved drugs:", drug_approval_dates.db_id.nunique())
    drug_approval_dates['Name'] = drug_approval_dates.Name.str.lower()
    drug_approval_dates.columns = ['db_id', 'product_name', 'intervention', 'labeler', 'start_marketting',
       'end_marketting', 'fda_app_num', 'approved', 'country', 'source',
       'Funder', 'approval_date', 'approval_year']
    drug_approval_dates = drug_approval_dates[drug_approval_dates.intervention.isin(dt_trial_df.intervention)]
    drug_approval_dates['approval_year'] = drug_approval_dates.approval_year.astype(int)

    #only consider FDA drugs
    drug_approval_dates = drug_approval_dates[drug_approval_dates.country == 'US']

    print("Number of drugs in CT mapped with approval dates:", drug_approval_dates.intervention.nunique())
    print("Number of targets in CT:", dt_trial_df[dt_trial_df.intervention.isin(drug_approval_dates.intervention)].Gene_Target.nunique())
    
    print("--")

    return df, drug_df, db_target_df, dt_trial_df, dt_start_year, ppi_g, placebo_trials, druggable_genome_df, drug_approval_dates
