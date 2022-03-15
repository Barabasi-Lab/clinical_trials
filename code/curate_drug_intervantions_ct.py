"""
#!/usr/bin/env python
Author= Kishore Vasan
Date = 03/15/2022

# Description:
## This script cleans the drug intervantion data in trials through a four step process
## and curates the drug intervention data

Input files : __ save the following files in the ../data/raw folder
1) intervention_ct_data.csv -- contains the xml parsed list of interventions in trials
1) all_drugbank_drugs.csv - parsed database containing drugbank drug id and names
2) drug_synonym.csv - containing drugbank id, name, and its corresponding synonym
3) products.csv -- containing drugbank id, name, and its corresponding products
4) drugs_external_identifiers.csv -- containings the drugbank id, name, and its corresponding external identifiers

## Specific steps

0.1 Ensure that the drug names from the intervention columns is already split
    i.e. ("A/B", "A or B", "A,B" becomes A,B) -- the `intervention_ct_data.csv` should already be split

1. map the drug intervention data to the DrugBank data -- this is done in the following steps :

    1.1) Search for direct text matching with the drug name in DrugBank

    1.2) Search for matching with synonyms of the drug names

    1.3) Map the intervention names with the product names of the drugs

    1.4) map the intervention names to external identifier (e.g. wikipedia)

    1.5) fuzzy string match of the names with the drugbank names

    1.6) save this data `../data/out/drug_mapped_ct_data.csv`

"""

# import files
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

###
## Step 1: Map drug interventions to Drug Bank (Direct Matching)
###

drug_df = pd.read_csv("../data/raw/intervention_ct_data.csv")
print("N trials:", drug_df.nct_id.nunique())
print("N interventions", drug_df.intervention.nunique())
print(drug_df.head())

db_target_df = pd.read_csv("../data/raw/all_drugbank_drugs.csv")
db_target_df['Name'] = db_target_df.Name.str.lower()
print("Number of drugs:", db_target_df.Name.nunique())
print("Number of targets:", db_target_df.Gene_Target.nunique())

print("N drugs with targets:", db_target_df[~db_target_df.Gene_Target.isna()].db_id.nunique())
print(db_target_df.head())

# data cleaning to get drug names from text
ids = []
new_int = []
int_type = []

all_drugs_list = set(db_target_df.Name.tolist())

start_time = time.time()
for _, r in drug_df.iterrows():

    if r['intervention_type'] != 'Drug':
        continue
        new_int.append(r['intervention'])
        ids.append(r['nct_id'])
        int_type.append(r['intervention_type'])

    int_names = r['intervention'].split('/')

    int_names_2 = []

    for i in int_names:
        if ',' in i:
            int_names_2.extend(i.split(','))
        else:
            int_names_2.append(i)

    tmp_l = set()

    for d_name in int_names_2:
        for j in all_drugs_list:
            if j in d_name:
                if d_name.index(j) == 0 or d_name[d_name.index(j)-1] == ' ':
                    tmp_l.add(j)

    new_int.extend(list(tmp_l))
    ids.extend([r['nct_id']]*len(tmp_l))
    int_type.extend([r['intervention_type']]*len(tmp_l))

print("time elapsed:", time.time() - start_time)
drug_df_2 = pd.DataFrame({'nct_id':ids,'intervention':new_int, 'intervention_type':int_type})
print(drug_df_2.head())

print("Number of trials:", drug_df_2[drug_df_2.intervention_type=='Drug'].nct_id.nunique())
print("Number of drugs:", drug_df_2[drug_df_2.intervention_type=='Drug'].intervention.nunique())
print("---")
print("original")
print("Number of trials:", drug_df[drug_df.intervention_type=='Drug'].nct_id.nunique())
print("Number of drugs:", drug_df[drug_df.intervention_type=='Drug'].intervention.nunique())

####
## Step 2: Map drug synonym
####

drug_synonym = pd.read_csv("../data/raw/drug_synonym.csv")
drug_synonym.columns = ['db_id','synonym']
drug_synonym = pd.merge(drug_synonym, db_target_df[['db_id','Name']], on='db_id')
drug_synonym['synonym'] = drug_synonym.synonym.str.lower()
drug_synonym = drug_synonym.drop_duplicates()

# remove synonyms with less than 5 letters
drug_synonym['synonym'] = drug_synonym.synonym.apply(lambda x: 1 if len(x) <5 else x)
drug_synonym = drug_synonym[drug_synonym.synonym != 1]
drug_synonym['synonym'] = drug_synonym['synonym'].astype(str)
drug_synonym['synonym'] = drug_synonym['synonym'].str.strip()

print("N drugs:", drug_synonym.Name.nunique())
print("N synonym:", drug_synonym.synonym.nunique())
print(drug_synonym.head())

#drug_synonym[drug_synonym.synonym.str.contains("interleukin-2")]

unmapped_interventions = drug_df[~drug_df.nct_id.isin(drug_df_2.nct_id)]
unmapped_interventions = unmapped_interventions[unmapped_interventions.intervention_type=='Drug']
print("N unmapped interventions:", unmapped_interventions.nct_id.nunique())
print("N unmapped_drugs:", unmapped_interventions.intervention.nunique())
print(unmapped_interventions.head())

clean_drugs = []
clean_nct_ids = []
mapped_synonym_list = []

n_synonym_error = 0
start_time = time.time()

# loop through all synonyms pairs
for _, r in drug_synonym.drop_duplicates().iterrows():

    # check for rows that have the synonym in them
    try:
        tmp_df = unmapped_interventions[unmapped_interventions.intervention.str.contains(r['synonym'])]
    except:
        n_synonym_error+=1
        #print(r['synonym'])
        # unable to match
        continue

    if len(tmp_df) >0 :
        clean_nct_ids.extend(list(tmp_df.nct_id.unique()))
        clean_drugs.extend([r['Name']]*tmp_df.nct_id.nunique())
        mapped_synonym_list.extend([r['synonym']]*tmp_df.nct_id.nunique())

print("Time elapsed:", time.time() - start_time)
print("N synonym unable to match", n_synonym_error)
synonym_mapped_drug_df = pd.DataFrame({"nct_id":clean_nct_ids,'intervention':clean_drugs,
                                       'synonym':mapped_synonym_list})
print("Number of drugs mapped:", synonym_mapped_drug_df.intervention.nunique())
print("Number of trials mapped:", synonym_mapped_drug_df.nct_id.nunique())
print(synonym_mapped_drug_df.head())

synonym_mapped_drug_df[synonym_mapped_drug_df.synonym=='interleukin-2']

t_df = unmapped_interventions[~unmapped_interventions.nct_id.isin(synonym_mapped_drug_df.nct_id)]

unmapped_names = list(t_df.intervention.unique())

###
## Step 3: Map drug products
###

drug_products = pd.read_csv("../data/raw/products.csv")
drug_products['product_name'] = drug_products.product_name.str.lower()
drug_products['Name'] = drug_products.Name.str.lower()
print("N drugs:", drug_products.Name.nunique())
print("N products:", drug_products.product_name.nunique())
print(drug_products.head())

print("Example:")
print(drug_products[drug_products.product_name=='statin'])

products_list = set(drug_products.product_name.tolist())
print("N products:", len(products_list))

clean_nct_ids = []
mapped_product_list = []

n_product_error = 0

start_time = time.time()

# loop through all products
for p in tqdm(products_list):
    # check for rows that have the product in them
    try:
        tmp_df = t_df[t_df.intervention.str.contains(p)]
    except:
        n_product_error+=1
        #print(r['product'])
        # unable to match
        continue

    if len(tmp_df) >0:
        clean_nct_ids.extend(list(tmp_df.nct_id.unique()))
        #clean_drugs.extend([r['Name']]*tmp_df.nct_id.nunique())
        mapped_product_list.extend([p]*tmp_df.nct_id.nunique())

print("Time elapsed:", time.time() - start_time)
print("N product unable to match", n_product_error)
product_mapped_drug_df = pd.DataFrame({"nct_id":clean_nct_ids,
                                       'product':mapped_product_list})
#print("Number of drugs mapped:", synonym_mapped_drug_df.intervention.nunique())
print("Number of trials mapped:", product_mapped_drug_df.nct_id.nunique())
print(product_mapped_drug_df.head())

unmapped_df = t_df[~t_df.nct_id.isin(product_mapped_drug_df.nct_id)]
print("N trials unmapped:", unmapped_df.nct_id.nunique())
print(unmapped_df.head())

####
## Step 4: Map to External Identifier
####

ext_ident_df = pd.read_csv("../data/raw/drugs_external_identifiers.csv")
ext_ident_df = ext_ident_df[ext_ident_df.identifier_resource=='Wikipedia']
ext_ident_df['identifier_name'] = ext_ident_df.identifier_name.str.lower()
ext_ident_df['Name'] = ext_ident_df.Name.str.lower()
print("N drugs:", ext_ident_df.nunique())
print(ext_ident_df.head())

product_drug_dict = {k: list(v) for k, v in drug_products.groupby('product_name')['Name']}

nct_id_list = []
db_id_list = []
drug_name_list = []
start_time = time.time()

for _,r in ext_ident_df.iterrows():
    tmp_df = unmapped_df[unmapped_df.intervention.str.contains(r['identifier_name'])]

    for _, r2 in tmp_df.iterrows():
        nct_id_list.append(r2['nct_id'])
        db_id_list.append(r['db_id'])
        drug_name_list.append(r['Name'])

print("Time elapsed:", time.time() - start_time)
wiki_mapped_df = pd.DataFrame({'nct_id':nct_id_list,
                              'intervention':drug_name_list})
wiki_mapped_df['intervention_type'] = 'Drug'
wiki_mapped_df = wiki_mapped_df.drop_duplicates()

print("N trials mapped:", wiki_mapped_df.nct_id.nunique())
print("N drugs mapped:", wiki_mapped_df.intervention.nunique())
print(wiki_mapped_df.head())

unmapped_df = unmapped_df[~unmapped_df.nct_id.isin(wiki_mapped_df.nct_id)]

print("N trials unmapped:", unmapped_df.nct_id.nunique())
print(unmapped_df.head())

####
### Step 5: Fuzzy matching drug names
####

drug_list = db_target_df.Name.unique()
drug_list = list(drug_list)
print("N drugs:", len(drug_list))

import Levenshtein
from scipy.spatial.distance import pdist

# given an intervention find the change needed to attain a drug
def get_closest_match(x):
    match_drug = []
    distance = []
    for d in drug_list:
        distance.append(Levenshtein.distance(x, d))
        match_drug.append(d)

    t_df = pd.DataFrame({'match_drug':match_drug,'distance':distance})
    t_df['intervention'] = x
    t_df = t_df[t_df.distance<5]

    return t_df

fuzzy_match = pd.DataFrame()

for u in tqdm(list(unmapped_df.intervention.unique())):
    fuzzy_match = fuzzy_match.append(get_closest_match(u), ignore_index = True)

print(fuzzy_match.head())
# remove unmatched interventions
fuzzy_match= fuzzy_match.dropna()

matched_dict = dict(zip(fuzzy_match[fuzzy_match.distance == 1].intervention,
                       fuzzy_match[fuzzy_match.distance == 1].match_drug))

map_list = []

for _, r in unmapped_df.iterrows():
    if r['intervention'] in matched_dict:
        map_list.append(matched_dict[r['intervention']])
    else:
        map_list.append('')

unmapped_df['map_intervention'] = map_list

print("N trials mapped:", unmapped_df[unmapped_df.map_intervention != ''].nct_id.nunique())
print("N interventions:", unmapped_df[unmapped_df.map_intervention != ''].map_intervention.nunique())

### combine all mappings
final_drug_map_df = copy.deepcopy(drug_df_2)
final_drug_map_df.index = range(len(final_drug_map_df))

# synonym mapping
for _, r in synonym_mapped_drug_df.iterrows():
    final_drug_map_df.loc[len(final_drug_map_df.index)] = [r['nct_id'],
                                                           r['intervention'],
                                                           'Drug']

# product mapping
for _, r in product_mapped_drug_df.iterrows():

    drugs_mapped_list_tmp = product_drug_dict[r['product']]

    for d in drugs_mapped_list_tmp:
        final_drug_map_df.loc[len(final_drug_map_df.index)] = [r['nct_id'],
                                                               d,
                                                               'Drug']
# wiki mapping
final_drug_map_df = final_drug_map_df.append(wiki_mapped_df, ignore_index = True)
final_drug_map_df = final_drug_map_df.drop_duplicates()

# fuzzy matching
final_drug_map_df = final_drug_map_df.append(unmapped_df[unmapped_df.map_intervention != ''][['nct_id',
                                                                                             'intervention',
                                                                                              'intervention_type']])

print("-----")
print("Final data counts:")
print("N Trials:", final_drug_map_df.nct_id.nunique())
print("N drugs:", final_drug_map_df.intervention.nunique())


print("Proportion of Trials mapped:", final_drug_map_df.nct_id.nunique()/drug_df[drug_df.intervention_type=='Drug'].nct_id.nunique())

# match columns with drugbank
final_drug_map_df.columns = ['nct_id','Name','intervention_type']

# save file
final_drug_map_df.to_csv("../data/out/drug_mapped_ct_data.csv", index=False)

### if you want to save the unmapped drug interventions
#unmapped_df = unmapped_df[unmapped_df.map_intervention == '']
#print("Example of unmapped intervention...")
#print(unmapped_df.head())
#unmapped_df.to_csv('ct_data/clean_data/unmapped_drug_trials.csv', index=False)

"""

####
## Extra: Get placebo trials and drugs
####

placebo_trials = drug_df[(drug_df.intervention.str.contains("placebo")) & (drug_df.intervention_type=='Drug')]
print("Number of trials with placebos:", placebo_trials.nct_id.nunique())
print("number of interventions:", placebo_trials.intervention.nunique())
print(placebo_trials.head())

# remove phrases
# placebo to, placebo for, placebo (for, placebo of
# begins with: placebo +
# ends with: and placebo, + placebo, or placebo

placebo_trials_2 = placebo_trials[~(placebo_trials.intervention.str.contains("placebo for")) & ~(placebo_trials.intervention.str.contains("placebo (for", regex=False))]
placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.str.contains("placebo to")]
placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.str.contains("placebo of")]

# starts with
placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.str.startswith("placebo+")]
placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.str.startswith("placebo +")]
placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.str.startswith("placebo plus")]
placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.str.startswith("plus placebo")]
placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.str.startswith("placebo and")]
placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.str.startswith("placebo or")]
placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.str.startswith("placebo matching")]
placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.str.startswith("placebo replacement")]
placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.str.startswith("placebo matched")]
placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.str.startswith("placebo /")]
placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.str.startswith("placebo/")]

# ends with
placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.str.endswith("+placebo")]
placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.str.endswith("+ placebo")]
placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.str.endswith("and placebo")]
placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.str.endswith("or placebo")]
placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.str.endswith("then placebo")]
placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.str.endswith("matching placebo")]
placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.str.endswith("/ placebo")]
placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.str.endswith("/placebo")]

# remove terms
remove_terms_1 = ['placebo oral tablet',
 'placebos',
 'placebo oral capsule',
 'matching placebo',
 'placebo comparator',
 'comparator: placebo',
 'placebo gel',
 'placebo tablet',
 'placebo capsule',
 'placebo patch',
 'placebo capsules',
 'placebo control',
 'placebo tablets',
 'placebo cream',
 'placebo pill',
 'oral placebo',
 'saline placebo',
 'placebo solution',
 'comparator: placebo (unspecified)',
 'placebo nasal spray',
 'matched placebo',
 'placebo injection',
 'inhaled placebo',
 'placebo group',
 'intranasal placebo',
 'placebo (saline)',
 'placebo administration',
 'placebo infusion',
 'intravenous placebo',
 'placebo iv',
 'placebo (normal saline)',
 'placebo ophthalmic solution',
 'placebo - cap',
 'placebo sc',
 'moxifloxacin placebo',
 'placebo treatment',
 'placebo - concentrate',
 'normal saline (placebo)',
 'placebo ointment',
 'placebo 2',
 'placebo 1',
 'iv placebo',
 'placebo dpi',
 'placebo drug',
 'placebo - sc',
 'placebo oral solution',
 'placebo - iv',
 'placebo arm']

placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.isin(remove_terms_1)]

remove_terms_2 = ['placebo (sugar pill)',
 'placebo single dose',
 'placebo eye drops',
 'administration of placebo',
 'placebo mouthwash',
 'placebo saline',
 'placebo pills',
 'vehicle (placebo)',
 'normal saline placebo',
 'placebo spray',
 'sugar pill (placebo)',
 'placebo sugar pill',
 'placebo inhaler',
 'placebo iv infusion',
 'placebo mdpi',
 'placebo oil',
 'placebo transdermal patch',
 'placebo film',
 'placebo (unspecified)',
 'placebo.',
 'placebo granules',
 'lactose placebo',
 'placebo (vehicle)',
 'placebo comparator: placebo',
 'placebo multiple doses',
 'placebo po',
 'sc placebo',
 'placebo medication',
 'placebo inhalation powder',
 'hec placebo gel',
 'placebo inhalation',
 'placebo bid',
 'placebo nasal aerosol',
 'matching placebo tablets',
 'b/f/taf placebo',
 'placebo (saline solution)',
 'placebo matched to atacicept',
 'placebo vaginal gel',
 'drug: placebo',
 'ftc/tdf placebo']

placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.isin(remove_terms_2)]

remove_terms_3 = ['placebo mdi',
 'placebo vaginal ring',
 'placebo sl tablet',
 'double-blind placebo',
 'placebo suspension',
 'pill placebo',
 'placebo intravenous',
 'placebo inhalation solution',
 'saline (placebo)',
 'placebo (high dose)',
 'sitagliptin placebo',
 'placebo intranasal spray',
 'inactive placebo',
 'placebo twice daily',
 'placebo peel',
 'placebo tablet to match 75 mg linzagolix tablet',
 'placebo drops',
 'placebo tablet to match 200 mg linzagolix tablet',
 'placebo vaginal tablet',
 'matching placebo patch',
 'placebo first',
 'placebo capsule to match add-back capsule',
 'placebo (artificial tears)',
 'placebo vehicle',
 'topical placebo',
 'srp plus placebo gel',
 'placebo - capsule',
 'e/c/f/tdf placebo',
 'placebo powder',
 'placebo 100 mg',
 'placebo 40 mg',
 'daily placebo',
 'matching placebo tablet',
 'group 1: placebo',
 'placebo (sc)',
 'matching placebo nasal spray',
 'transdermal placebo patch',
 'placebo multiple dose',
 'placebo: normal saline',
 'placebo 100mg',
 'placebo control group',
 'placebo collagen sponge',
 'bkm120 matching placebo']

placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.isin(remove_terms_3)]

remove_terms_4 = ['placebo matched to ivacaftor',
                  'placebo gas',
                  'topical placebo cream',
                  'placebo vaginal insert',
                  'asa placebo',
                  'placebo ointment (vehicle)',
                  'matched placebos',
                  'placebo (p)',
                 'placebo tablets bid',
                 'control: placebo',
 'placebo (id)','control placebo (c)']

placebo_trials_2 = placebo_trials_2[~placebo_trials_2.intervention.isin(remove_terms_4)]

placebo_trials_2 = placebo_trials_2.drop_duplicates()
print("Number of trials with placebos:", placebo_trials_2.nct_id.nunique())
print("number of interventions:", placebo_trials_2.intervention.nunique())

t = placebo_trials_2.intervention.value_counts().to_frame()

tagged_drug_names = []
synonym_mapped_names = []
mapping_drug_names = []

t = drug_synonym.drop_duplicates()

for j in tqdm(range(len(t))):

    try:
        tmp_df = placebo_trials_2[placebo_trials_2.intervention.str.contains(t.iloc[j]['synonym'], regex=False)]
    except:
        continue

    if len(tmp_df) > 0:
        mapping_drug_names.extend([t.iloc[j]['Name']]*len(tmp_df))
        synonym_mapped_names.extend([t.iloc[j]['synonym']]*len(tmp_df))
        tagged_drug_names.extend(tmp_df.intervention.tolist())

tmp_mapping_df = pd.DataFrame({'placebo_name':tagged_drug_names, 'drug_map':mapping_drug_names,'synonym_tagged':synonym_mapped_names})
print(tmp_mapping_df.head())

tmp_mapping_df = tmp_mapping_df.drop_duplicates(['placebo_name','drug_map'])

synonym_mapping_dict = dict(zip(tmp_mapping_df.placebo_name, tmp_mapping_df.drug_map))

final_placebo_trials = pd.merge(placebo_trials_2, tmp_mapping_df[['placebo_name', 'drug_map']], left_on='intervention',right_on='placebo_name')
final_placebo_trials = final_placebo_trials[['nct_id','drug_map','intervention_type']]
final_placebo_trials.head()

print("number of placebo trials:", final_placebo_trials.nct_id.nunique())
print("number of placebo drugs:", final_placebo_trials.drug_map.nunique())
final_placebo_trials.to_csv("ct_data/clean_data/placebo_trials.csv", index = False)
"""
