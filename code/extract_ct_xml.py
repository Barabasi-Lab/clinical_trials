"""
#!/usr/bin/env python
Author= Kishore Vasan
Date = 03/15/2022
### Keywords: clinical trials;
# Description:
## This script loads the XML Files of the clinical trials
## and organizes them into readable csv files
### input: xml files of clinicaltrials.gov bulk download
### out: 1) organized_ct_data.csv 2) funder_ct_data.csv 3) intervention_ct_data.csv
"""

# import packages
import pandas as pd
import numpy as np
import glob
import xml.etree.ElementTree as ET
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
import time
import copy

import warnings
warnings.filterwarnings('ignore')

# Get ClinicalTrials.Gov data

files_list = glob.glob("../data/raw/NCT*/*.xml")
print("Number of Studies:", len(files_list))

# get the xml fields of a given clinical trial
def organize_data(ct_file_name):
    # create element tree object
    tree = ET.parse(ct_file_name)

    # get root element
    root = tree.getroot()

    # trial id
    nct_id = root.find('.//id_info//nct_id').text

    # trial title
    brief_title = root.find('.//brief_title').text

    # lead sponsor
    lead_sponsors = []
    lead_sponsors_types = []

    for sub_level in root.findall('.//sponsors//lead_sponsor'):
        lead_sponsors.append(sub_level.find('agency').text)
        if sub_level.find("agency_class") is not None:
            lead_sponsors_types.append(sub_level.find('agency_class').text)
        else:
            lead_sponsors_types.append("NA")
    lead_sponsors = ';'.join(lead_sponsors)
    lead_sponsors_types = ';'.join(lead_sponsors_types)

    # collaborators
    collaborators = []
    collaborators_types = []

    for sub_level in root.findall('.//sponsors//collaborator'):
        collaborators.append(sub_level.find('agency').text)
        if sub_level.find('agency_class') is not None:
            collaborators_types.append(sub_level.find("agency_class").text)
        else:
            collaborators_types.append("NA")
    collaborators = ';'.join(collaborators)
    collaborators_types = ';'.join(collaborators_types)

    ### eligibility
    gender = root.find(".//eligibility//gender")
    if gender is not None:
        gender = gender.text
    else:
        gender = ''

    min_age = root.find(".//eligibility//minimum_age")
    if min_age is not None:
        min_age = min_age.text.split(' ')[0]
    else:
        min_age = ''

    max_age = root.find(".//eligibility//maximum_age")
    if max_age is not None:
        max_age = max_age.text.split(' ')[0]
    else:
        max_age = ''

    ### study location
    location_countries = []

    for sub_level in root.findall("location_countries"):
        if sub_level.find("country") is not None:
            location_countries.append(sub_level.find("country").text)

    location_countries = ';'.join(location_countries)

    ### publications
    result_pubs = []

    for sub_level in root.findall('results_reference'):
        if sub_level.find("PMID") is not None:
            result_pubs.append(sub_level.find('PMID').text)

    result_pubs = ';'.join(result_pubs)

    references = []
    for sub_level in root.findall("reference"):
        if sub_level.find("PMID") is not None:
            references.append(sub_level.find('PMID').text)
    references = ';'.join(references)

    ### keywords
    keywords = []

    for sub_level in root.findall("keyword"):
        keywords.append(sub_level.text)
    keywords = ";".join(keywords)


    ### status, phase, and start date
    if root.find("overall_status") is not None:
        status = root.find("overall_status").text
    else:
        status = ''

    if root.find("phase") is not None:
        phase = root.find("phase").text
    else:
        phase = ""

    if root.find("start_date") is not None:
        start_date = root.find("start_date").text
    else:
        start_date = ''

    ### intervention
    interventions = []
    interventions_type = []

    for sub_level in root.findall('.//intervention'):
        interventions.append(sub_level.find("intervention_name").text)
        interventions_type.append(sub_level.find("intervention_type").text)

    interventions = ';'.join(interventions)
    intervention_types = ';'.join(interventions_type)

    ### condition a.k.a. disease
    conditions = []

    for sub_level in root.findall('.//condition'):
        conditions.append(sub_level.text)
    conditions = ';'.join(conditions)

    # mesh terms
    mesh_terms = root.findall('.//mesh_term')

    if len(mesh_terms) == 0:
        mesh_terms = ''
    else:
        mesh_terms = ';'.join([term.text for term in mesh_terms])

    # study type
    study_type = root.find("study_type").text

    return {'nct_id':nct_id, 'title':brief_title,'study_type':study_type,'gender':gender,'min_age':min_age,'max_age':max_age,
            'lead_sponsors':lead_sponsors,'lead_sponsor_type':lead_sponsors_types,
            'collaborators':collaborators, 'collaborator_types':collaborators_types,
            'interventions':interventions,'intervention_types':intervention_types,
            'status':status,'phase':phase, 'start_date':start_date,'location_countries':location_countries,
            'conditions':conditions,'keywords':keywords,'mesh_terms':mesh_terms,
           'result_pubs_pmid':result_pubs,'references_pmid':references}

    """# print stuff
    print(nct_id, study_type)
    print(lead_sponsor, lead_sponsor_type)
    print(keywords)
    print(location_countries)
    print(gender, min_age, max_age)
    print(interventions)
    print(intervention_types)
    print(status, phase, start_date)
    print(result_pubs)
    print(conditions)
    print(references)
    print("Lead sponsor:")
    print(lead_sponsors)
    print(lead_sponsors_type)
    print("collaborators")
    print(collaborators)
    print(collaborators_types)

    """

trials_list = []
df = pd.DataFrame()

for file_n in tqdm(files_list):
    trials_list.append(organize_data(file_n))

df = df.append(trials_list)
print(df.shape)
print(df.head())

df.to_csv("../data/out/organized_ct_data.csv", index = False)

## organize funder data

def organize_funder_data():
    trial_id = []
    funder_names = []
    funder_types = []
    funder_roles = []

    start_time = time.time()
    for _, r in df.iterrows():
        funder = [r['lead_sponsors']]
        funder_type = [r['lead_sponsor_type']]
        funder_role = ['lead']
        if type(r['collaborators']) != float:
            if len(r['collaborators'].split(';')) == len(r['collaborator_types'].split(';')):
                funder.extend(r['collaborators'].split(';'))
                funder_type.extend(r['collaborator_types'].split(';'))
            else:
                funder.extend(r['collaborators'].split(';'))
                funder_type.extend(['na']*len(r['collaborators'].split(';')))
            funder_role.extend(['collaborator']*len(r['collaborators'].split(';')))

        funder_names.extend(funder)
        funder_types.extend(funder_type)
        trial_id.extend([r['nct_id']]*len(funder))
        funder_roles.extend(funder_role)

    print("Time elapsed:", time.time() - start_time)
    return pd.DataFrame({'trial_id' : trial_id,'funder_name':funder_names,
                          'funder_type':funder_types,'funder_role':funder_roles})

funder_df = organize_funder_data()
funder_df = funder_df[funder_df.funder_name != '']
funder_df.head()

# save funder info
funder_df.to_csv("../data/out/funder_ct_data.csv", index = False)

### organize intervention data

def organize_intervention_data():
    intervention_names = []
    intervention_types = []
    nct_id = []

    counter = 0
    start_time = time.time()
    for _, r in df.iterrows():
        if type(r['interventions']) == float:
            counter+=1
            continue
        tmp_interventions = r['interventions'].split(';')
        tmp_interventions_types = r['intervention_types'].split(';')
        if len(tmp_interventions) != len(tmp_interventions_types):
            counter+=1
            continue

        intervention_names.extend(tmp_interventions)
        intervention_types.extend(tmp_interventions_types)
        nct_id.extend([r['nct_id']]*len(tmp_interventions))

    print("Time taken:", time.time() - start_time)
    print("Number of trials without interventions or intervention types:", counter)

    return pd.DataFrame({'nct_id':nct_id,'intervention':intervention_names,'intervention_type':intervention_types})

drug_df = organize_intervention_data()
drug_df['intervention']= drug_df.intervention.str.lower()
drug_df = drug_df.drop_duplicates()

remove_interventions = ['placebo','no intervention']
drug_df = drug_df[~drug_df.intervention.isin(remove_interventions)]
print(drug_df.shape)
print(drug_df.head())

# save file
drug_df.to_csv("../data/raw/intervention_ct_data.csv")
