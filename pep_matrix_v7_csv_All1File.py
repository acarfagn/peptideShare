# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 19:32:15 2019

@author: amyk0
"""

#%%

import pandas as pd
import os
import re

#%%

# output file name
excel_file = 'peptide_matrix.xlsx'
# Set location where files are saved
directory = os.getcwd()
subfolder = directory + '\\input'

#%%
# Function to get unique list of peptides in Raw Data tab

def alpha_only(pep):
    p_ptm = re.compile('\([+-][0-9]{0,3}\.[0-9]{1,2}\)')
    p_mut = re.compile('\([a-z]{3,3}[\s]?[A-Z]?\)')
    pep2 = p_ptm.sub('',pep)
    pep_alpha = p_mut.sub('',pep2)
    return pep_alpha

def read_raw_peps(filename):
    '''
    Takes script output Excel file as input
    Reads Raw Data tab into dataframe
    Extracts list of unique peptides (including PTM indicators)
    and returns this list
    '''
    df = pd.read_excel(filename, sheet_name = 'Raw Data', skiprows = [0])
    peptides = df.Peptide.unique()
    peptides_long = [x for x in peptides if len(alpha_only(x)) > 7]
    
    accession_list = []
    for pep in peptides_long:
        s = df.loc[df.Peptide == pep]['Accession']
        acc = s.iloc[0]
        accession_list.append(acc)
    
    return peptides_long, accession_list

def read_pep_csv(filename):
    df0 = pd.read_csv(filename,usecols=['Peptide', 'Accession', '-10lgP'])
    df0['PepAlpha'] = df0['Peptide'].apply(alpha_only)
    criterion = df0['PepAlpha'].map(lambda x: len(x) > 7)
    df = df0[criterion]
    peptides_long = df.Peptide.to_list()
    accession_list = df.Accession.to_list()
    scores = df['-10lgP'].to_list()
        
    return peptides_long, accession_list, scores
    
    
def read_denovo_csv(filename):
    try:
        df = pd.read_csv(filename,usecols=['Peptide', 'Accession'])
        peptides = df.Peptide.unique()
        peptides_long = [x for x in peptides if len(alpha_only(x)) > 7]
        
        accession_list = []
        for pep in peptides_long:
            s = df.loc[df.Peptide == pep]['Accession']
            acc = s.iloc[0]
            accession_list.append(acc)
            
    except ValueError:
        df = pd.read_csv(filename,usecols=['Peptide'])
        peptides = df.Peptide.unique()
        peptides_long = [x for x in peptides if len(alpha_only(x)) > 7]
        
        accession_list = ['denovo'] * len(peptides_long)
        
    return peptides_long, accession_list 

#%%

def read_pep_csv_1file(filename):
    df0 = pd.read_csv(filename,usecols=['Peptide', 'Accession', '-10lgP'])
    df0['PepAlpha'] = df0['Peptide'].apply(alpha_only)
    criterion = df0['PepAlpha'].map(lambda x: len(x) > 7)
    df = df0[criterion]
    
    samples = df['xxx'].unique().tolist()
    
    peptides_long = df.Peptide.to_list()
    accession_list = df.Accession.to_list()
    scores = df['-10lgP'].to_list()
        
    return peptides_long, accession_list, scores

#%%

# Initialize dictionary where keys will be sample name and values
# will be unique peptides from Raw Data tab for each sample
peps_by_file = dict()
peps_by_acc = dict()
scores_by_file = dict()

#%%

# Iterate through the files
# Populate the dictionary
with os.scandir(subfolder) as it:
    for entry in it:
        if entry.name.endswith('.xls'):
            print('Processing file ', entry.name)
            x = entry.name.find('.')
            filename_abbrev = entry.name[:x]
            peptides, accessions = read_raw_peps(entry)
            peps_by_file[filename_abbrev] = peptides
            
            for i in range(len(peptides)):
                pep = peptides[i]
                if pep in peps_by_acc.keys():
                    if peps_by_acc[pep] == accessions[i]:
                        continue
                    else:
                        # Blank accession column gets 'nan'
                        # Need to backtrack to WHY accession column is blank for DB result!
                        if type(accessions[i]) == 'str':                            
                            peps_by_acc[pep] = peps_by_acc[pep] + ' | ' + accessions[i]
                else:
                    peps_by_acc[pep] = accessions[i]
            
        elif entry.name.endswith('.csv'):
            if 'DENOVO' in entry.name:
                print('Processing file ', entry.name)
                x = entry.name.find('.')
                filename_abbrev = entry.name[:x]
                peptides, accessions = read_denovo_csv(entry)
                peps_by_file[filename_abbrev] = peptides
                
                for i in range(len(peptides)):
                    pep = peptides[i]
                    if pep in peps_by_acc.keys():
                        if peps_by_acc[pep] == accessions[i]:
                            continue
                        else:
                            # Blank accession column gets 'nan'
                            # Need to backtrack to WHY accession column is blank for DB result!
                            if type(accessions[i]) == 'str':                            
                                peps_by_acc[pep] = peps_by_acc[pep] + ' | ' + accessions[i]
                    else:
                        peps_by_acc[pep] = accessions[i] 
            elif 'PEPTIDE' in entry.name:
                print('Processing file ', entry.name)
                x = entry.name.find('.')
                filename_abbrev = entry.name[:x]
                peptides, accessions, scores = read_pep_csv(entry)
                
                peps_by_file[filename_abbrev] = peptides
                scores_by_file[filename_abbrev] = scores
                
                for i in range(len(peptides)):
                    pep = peptides[i]
                    if pep in peps_by_acc.keys():
                        if peps_by_acc[pep] == accessions[i]:
                            continue
                        else:
                            # Blank accession column gets 'nan'
                            # Need to backtrack to WHY accession column is blank for DB result!
                            if type(accessions[i]) == 'str':                            
                                peps_by_acc[pep] = peps_by_acc[pep] + ' | ' + accessions[i]
                    else:
                        peps_by_acc[pep] = accessions[i]     

# Get parent list of peptides for all files being compared
parent_pep_list = []
parent_acc_list = []

for value in peps_by_file.values():
    for pep in value:
        if pep not in parent_pep_list:
            parent_pep_list.append(pep)

for j in range(len(parent_pep_list)):
    p = parent_pep_list[j]
    parent_acc_list.append(peps_by_acc[p])

#%%

# Initialize dictionary to store indicators of peptide presence
# Keys are sample names
# Each value is a list with binary indicator
# Iterate through the parent list of peptides;
# update value list: 1 means peptide is in sample, 0 means it is not
# ORDER is very important: order of value list should correspond
# to order of peptides in parent list; so, if element at index 3 of 
# value list is 1, that means peptide at index 3 in parent list is found
# in that sample

pep_count_dict = dict()

for sample in peps_by_file.keys():
    pep_count_dict[sample] = [0]*len(parent_pep_list)
    for i in range(len(parent_pep_list)):
        if parent_pep_list[i] in peps_by_file[sample]:
            ls = peps_by_file[sample]
            pidx = ls.index(parent_pep_list[i])
            sc = scores_by_file[sample][pidx]
            pep_count_dict[sample][i] = sc

# Prep for dataframe creation by adding parent peptide list as entry
# in dictionary
pep_count_dict['Peptide'] = parent_pep_list
pep_count_dict['Accession'] = parent_acc_list

# Make dataframe
pep_count_df = pd.DataFrame.from_dict(pep_count_dict, orient='columns')

#%%

# Reorder columns; Peptide list should be first column
# Subsequent columns should give sample names in alphabetical order
colnames = list(pep_count_df.columns)

colnames.sort()

colnames.remove('Peptide')
colnames.remove('Accession')
colnames.insert(0, 'Peptide')
colnames.insert(1, 'Accession')
#%%
pep_count_df = pep_count_df[colnames]
#%%
# Write the dataframe to Excel
writer = pd.ExcelWriter('output\\' + excel_file, engine='xlsxwriter')
workbook = writer.book
pep_count_df.to_excel(writer, index=False, header=True)
workbook.close()

