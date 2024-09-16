#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Title: Genomic Data Preprocessing and Merging Script
Filename: data_preprocessing.py

Description:
This script processes genomic datasets by merging coverage, heterozygosity, and GC
 content data from multiple sources. It calculates differential metrics between male 
 and female samples, cleans the data, and prepares it for further analysis. 
 The script also filters contigs belonging to multiple chromosomes.

Input:
- Coverage files (CSV format): Contains coverage data for multiple samples.
  Example:
    - cov_SRR11567976.csv
    - cov_SRR13996603.csv
    - cov_SRR19888834.csv
    - cov_SRR19888872.csv
    - cov_SRR19888877.csv
    - cov_SRR8893602.csv

- Heterozygosity files (CSV format): Contains heterozygosity data for the same set of samples.
  Example:
    - sample1.csv
    - sample2.csv
    - sample3.csv
    - sample4.csv
    - sample5.csv
    - sample6.csv

- GC Content files (Tab-separated format): Contains GC content data for each sample.
  Example:
    - srr_1156.txt
    - srr_1399.txt
    - srr_198834.txt
    - srr_198872.txt
    - srr_1988.txt
    - srr_8893.txt

- Classification file (CSV format): Contains classification information (chromosomes, contig IDs, etc.).
  Example:
    - clean_classified.csv

Output:
- Merged dataset including differential coverage, heterozygosity, and GC content metrics.
  Files generated:
    - xy_coverage.csv: Differential coverage between male and female samples.
    - xy_heterozygosity.csv: Differential heterozygosity.
    - xy_gcDif.csv: Differential GC content.
    - xy_data.csv: Final merged dataset combining all metrics.
    
- Statistics and summary outputs such as:
    - Number of contigs that belong to more than one chromosome.

Authors:
- Main Author: Erika Castaneda

Version:
- Python 3.11.5

Usage:
To run this script, ensure all input files are located in the correct directory. Then, execute it using Python:

"""

# Import necessary packages
import os
import pandas as pd
import numpy as np


# Load datasets
# Coverage data for multiple samples
SRR11567976 = pd.read_csv('cov_SRR11567976.csv', delimiter=',')
SRR13996603 = pd.read_csv('cov_SRR13996603.csv', delimiter=',')
SRR19888834 = pd.read_csv('cov_SRR19888834.csv', delimiter=',')
SRR19888872 = pd.read_csv('cov_SRR19888872.csv', delimiter=',')
SRR19888877 = pd.read_csv('cov_SRR19888877.csv', delimiter=',')
SRR8893602 = pd.read_csv('cov_SRR8893602.csv', delimiter=',')

# Heterozygosity data
h_SRR11567976 = pd.read_csv('sample1.csv', delimiter=',')
h_SRR13996603 = pd.read_csv('sample2.csv', delimiter=',')
h_SRR19888834 = pd.read_csv('sample3.csv', delimiter=',')
h_SRR19888872 = pd.read_csv('sample4.csv', delimiter=',')
h_SRR19888877 = pd.read_csv('sample5.csv', delimiter=',')
h_SRR8893602 = pd.read_csv('sample6.csv', delimiter=',')

# GC content data
gc_SRR11567976 = pd.read_csv('srr_1156.txt', delimiter='\t', header=None)
gc_SRR13996603 = pd.read_csv('srr_1399.txt', delimiter='\t', header=None)
gc_SRR19888834 = pd.read_csv('srr_198834.txt', delimiter='\t', header=None)
gc_SRR19888872 = pd.read_csv('srr_198872.txt', delimiter='\t', header=None)
gc_SRR19888877 = pd.read_csv('srr_1988.txt', delimiter='\t', header=None)
gc_SRR8893602 = pd.read_csv('srr_8893.txt', delimiter='\t', header=None)

# Classification data
classification = pd.read_csv('clean_classified.csv')

# Merging Coverage Data
# Create a dictionary to hold the coverage data for different samples
list_dataframes = {
    'SRR11567976': SRR11567976,
    'SRR13996603': SRR13996603,
    'SRR19888834': SRR19888834,
    'SRR19888872': SRR19888872,
    'SRR19888877': SRR19888877,
    'SRR8893602': SRR8893602
}

# Prepare dataset for coverage difference analysis
dataset_CC = classification[['Chromosomes', 'Contig_ID']]
merged_dataset = None

# Merge all sample dataframes on 'Contig_ID'
for name_sample, df in list_dataframes.items():
    if merged_dataset is None:
        merged_dataset = df.copy()
    else:
        merged_dataset = pd.merge(merged_dataset, df, on='Contig_ID', how='left')

# Combine with chromosome and contig data
xy_coverage = pd.merge(dataset_CC, merged_dataset, on='Contig_ID', how='left')

# Calculate differential coverage between male and female samples
male_columns = ['SRR13996603_no_duplicates TPM', 'SRR11567976_no_duplicates TPM', 'SRR8893602_no_duplicates TPM']
female_columns = ['SRR19888877_no_repeats TPM', 'SRR19888872_no_duplicates TPM', 'SRR19888834_no_duplicates TPM']

# Compute the coverage difference
xy_coverage['coverage_dif'] = (xy_coverage[female_columns].mean(axis=1) - xy_coverage[male_columns].mean(axis=1)) / 3
xy_coverage = xy_coverage[['Chromosomes', 'Contig_ID', 'coverage_dif']]

# Merging Heterozygosity Data
# Dictionary of heterozygosity dataframes for different samples
list_dataframes_h = {
    'SRR11567976': h_SRR11567976,
    'SRR13996603': h_SRR13996603,
    'SRR19888834': h_SRR19888834,
    'SRR19888872': h_SRR19888872,
    'SRR19888877': h_SRR19888877,
    'SRR8893602': h_SRR8893602
}

merged_dataset_h = None

# Merge all heterozygosity dataframes on 'Contig_ID'
for name_sample, df in list_dataframes_h.items():
    if merged_dataset_h is None:
        merged_dataset_h = df.copy()
    else:
        merged_dataset_h = pd.merge(merged_dataset_h, df, on='Contig_ID', how='left', suffixes=('', f'_{name_sample}'))

# Combine with chromosome and contig data
xy_heterozigosity = pd.merge(dataset_CC, merged_dataset_h, on='Contig_ID', how='left')

# Calculate differential heterozygosity between male and female samples
male_columns_h = ['Heterozigosity', 'Heterozigosity_SRR13996603', 'Heterozigosity_SRR8893602']
female_columns_h = ['Heterozigosity_SRR19888834', 'Heterozigosity_SRR19888872', 'Heterozigosity_SRR19888877']

# Compute heterozygosity difference
xy_heterozigosity['Heterozigosity_dif'] = (xy_heterozigosity[female_columns_h].sum(axis=1) - xy_heterozigosity[male_columns_h].sum(axis=1)) / 3
xy_heterozigosity = xy_heterozigosity[['Contig_ID', 'Heterozigosity_dif']].dropna()

# Merging GC Content Data
header = ["Contig_ID", "Start", "End", "GC"]

list_dataframes_gc = {
    'SRR11567976': gc_SRR11567976,
    'SRR13996603': gc_SRR13996603,
    'SRR19888834': gc_SRR19888834,
    'SRR19888872': gc_SRR19888872,
    'SRR19888877': gc_SRR19888877,
    'SRR8893602': gc_SRR8893602
}

merged_dataset_gc = None

# Merge all GC content dataframes on 'Contig_ID'
for name_sample, df in list_dataframes_gc.items():
    df.columns = header  # Assign headers
    if merged_dataset_gc is None:
        merged_dataset_gc = df.copy()
    else:
        merged_dataset_gc = pd.merge(merged_dataset_gc, df, on='Contig_ID', how='left', suffixes=('', f'_{name_sample}'))

# Remove '%' from GC content and convert to float
xy_gc_columns = [col for col in merged_dataset_gc.columns if col.startswith('GC')]
merged_dataset_gc[xy_gc_columns] = merged_dataset_gc[xy_gc_columns].apply(lambda x: x.str.replace('%', '').astype(float))

# Calculate differential GC content between male and female samples
male_columns_gc = ['GC', 'GC_SRR13996603', 'GC_SRR8893602']
female_columns_gc = ['GC_SRR19888834', 'GC_SRR19888872', 'GC_SRR19888877']

# Compute GC content difference
merged_dataset_gc['GC_dif'] = (merged_dataset_gc[female_columns_gc].sum(axis=1) - merged_dataset_gc[male_columns_gc].sum(axis=1)) / 3
xy_gcDif = merged_dataset_gc[['Contig_ID', 'GC_dif']]

# Merging Coverage, Heterozygosity, and GC Content into Final Dataset
xy_data = pd.merge(xy_coverage, xy_heterozigosity, on='Contig_ID', how='left')
xy_data = pd.merge(xy_data, xy_gcDif, on='Contig_ID', how='left')
xy_data = xy_data.dropna()

# Classification: Label contigs belonging to chromosomes X or Y
xy_data['class'] = xy_data['Chromosomes'].str.contains('X|Y').astype(int)

# Group by Contig_ID and compute GC mean
grouped_gc_xy = merged_dataset_gc.groupby('Contig_ID')['GC'].mean().reset_index()

# Filter contigs belonging to more than one chromosome
contig_chromosomes_counts = xy_data.groupby('Contig_ID')['Chromosomes'].nunique()
contigs_multiple_chromosomes = contig_chromosomes_counts[contig_chromosomes_counts > 1]
num_contigs_multiple_chromosomes = len(contigs_multiple_chromosomes)

print("Number of contigs belonging to more than one chromosome:", num_contigs_multiple_chromosomes)
