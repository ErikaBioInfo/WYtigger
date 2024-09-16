#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description
This script generates various plots to visualize differences in coverage, heterozygosity, 
and GC content between females and males based on processed genomic data. It includes 
scatter plots with normalized values and confidence intervals, as well as class proportions.

Inputs
    xy_data (DataFrame):
        Contains genomic data with columns for chromosomes, coverage difference, 
        heterozygosity difference, GC difference, and class labels.
        Required columns: 'Chromosomes', 'coverage_dif', 'Heterozygosity_dif', 'GC_dif', 'class'.
Outputs
    Graphs in .jpg and .pdf formats:
        scatter_plot_coverage.jpg: Scatter plot of normalized coverage differences between females and males.
        scatter_plot_coverage.pdf: PDF version of the scatter plot for coverage differences.
        scatter_plot_heterozygosity.pdf: Scatter plot of normalized heterozygosity differences between females and males.
        scatter_plot_GC.pdf: Scatter plot of normalized GC content differences between females and males.
Versions
    Version Python 3.11.5

Usage
    Prepare Data: Ensure xy_data is correctly processed and contains the required columns.
    Run Script: Execute the script in a Python environment with the necessary libraries installed.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import t
from sklearn.preprocessing import MinMaxScaler, Normalize
from matplotlib.cm import ScalarMappable

# Function to read and preprocess data
def read_and_preprocess_data(file_path):
    # Load data from CSV file
    data = pd.read_csv(file_path)
    
    # Define chromosome mapping
    mapping = {
        'chr1': '1', 'chr2': '2', 'chr3': '3', 'chr4': '4', 'chr5': '5', 'chr6': '6', 'chr7': '7', 'chr8': '8', 'chr9': '9',
        'chr10': '10', 'chr11': '11', 'chr12': '12', 'chr13': '13', 'chr14': '14', 'chr15': '15', 'chr16': '16', 'chr17': '17',
        'chr18': '18', 'chr19': '19', 'chr20': '20', 'chr21': '21', 'chr22': '22', 'chrY': 'Y', 'chrX': 'X', 'chrUn': 'Un'
    }
    
    # Replace values in the 'Chromosomes' column using the mapping
    data['Chromosomes'] = data['Chromosomes'].replace(mapping)
    
    return data

# Define custom order for chromosomes
custom_order = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15',
                '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'Un']

# Function to calculate confidence intervals
def calculate_confidence_intervals(data, column, confidence=0.95):
    ci_values = []
    for chromosome, group_data in data.groupby('Chromosomes'):
        n = len(group_data)
        mean_diff = group_data[column].mean()
        std_diff = group_data[column].std()
        margin_of_error = t.ppf((1 + confidence) / 2, n - 1) * (std_diff / np.sqrt(n))
        ci_values.append((chromosome, mean_diff, margin_of_error))
    
    ci_df = pd.DataFrame(ci_values, columns=['Chromosomes', 'Mean_diff', 'Margin_of_error'])
    ci_df['CI_lower'] = ci_df['Mean_diff'] - ci_df['Margin_of_error']
    ci_df['CI_upper'] = ci_df['Mean_diff'] + ci_df['Margin_of_error']
    return ci_df

# Function to plot normalized scatter plots
def plot_normalized_scatter(data, column, output_file, title, ylabel):
    # Filter out outliers
    Q1 = data[column].quantile(0.01)
    Q3 = data[column].quantile(0.99)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR
    filtered_data = data[(data[column] >= lower_bound) & (data[column] <= upper_bound)]
    
    # Normalize data
    scaler = MinMaxScaler(feature_range=(-1, 1))
    filtered_data[f'{column}_normalized'] = scaler.fit_transform(filtered_data[column].values.reshape(-1, 1))
    filtered_data['CI_lower_normalized'] = scaler.fit_transform(filtered_data['CI_lower'].values.reshape(-1, 1))
    filtered_data['CI_upper_normalized'] = scaler.fit_transform(filtered_data['CI_upper'].values.reshape(-1, 1))

    # Convert 'Chromosomes' to a categorical type with custom order
    filtered_data['Chromosomes'] = pd.Categorical(filtered_data['Chromosomes'], categories=custom_order, ordered=True)
    filtered_data = filtered_data.sort_values('Chromosomes')

    # Define color mapping
    cmap = sns.diverging_palette(240, 10, as_cmap=True)
    norm = Normalize(-1, 1)
    colors = cmap(norm(filtered_data[f'{column}_normalized']))

    # Plot scatter plot
    plt.figure(figsize=(30, 8))
    scatter = plt.scatter(filtered_data['Chromosomes'], filtered_data[f'{column}_normalized'], c=colors, s=40)
    plt.xticks(rotation=45, ha='center', fontsize=20)
    plt.xlabel('Chromosomes', fontsize=24)
    plt.ylabel(ylabel, fontsize=24)
    plt.title(title, fontsize=24)

    # Add colorbar
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=plt.gca())
    cbar.set_label('Normalized Value', fontsize=18)
    cbar.ax.tick_params(labelsize=18)

    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    plt.show()

# Main execution
if __name__ == "__main__":
    # Read and preprocess data
    file_path = 'data.csv'  # Replace with your file path
    xy_data = read_and_preprocess_data(file_path)

    # Calculate confidence intervals for each metric
    coverage_ci = calculate_confidence_intervals(xy_data, 'coverage_dif', confidence=0.95)
    heterozygosity_ci = calculate_confidence_intervals(xy_data, 'Heterozygosity_dif', confidence=0.95)
    gc_ci = calculate_confidence_intervals(xy_data, 'GC_dif', confidence=0.95)

    # Merge data with confidence intervals
    coverage_data = xy_data.merge(coverage_ci[['Chromosomes', 'CI_lower', 'CI_upper']], on='Chromosomes')
    heterozygosity_data = xy_data.merge(heterozygosity_ci[['Chromosomes', 'CI_lower', 'CI_upper']], on='Chromosomes')
    gc_data = xy_data.merge(gc_ci[['Chromosomes', 'CI_lower', 'CI_upper']], on='Chromosomes')

    # Plot graphs
    plot_normalized_scatter(coverage_data, 'coverage_dif', 'scatter_plot_coverage.pdf', 'Coverage Difference by Chromosome', 'Coverage Difference')
    plot_normalized_scatter(heterozygosity_data, 'Heterozygosity_dif', 'scatter_plot_heterozygosity.pdf', 'Heterozygosity Difference by Chromosome', 'Heterozygosity Difference')
    plot_normalized_scatter(gc_data, 'GC_dif', 'scatter_plot_GC.pdf', 'GC Content Difference by Chromosome', 'GC Difference')

    # Plot class proportion
    class_count = xy_data['class'].value_counts()
    print(f'Class 0: {class_count[0]}')
    print(f'Class 1: {class_count[1]}')
    print(f'Proportion: {round(class_count[0] / class_count[1], 2)}:1')

    plt.figure(figsize=(10, 6))
    class_count.plot(kind='bar', title='Count of Each Class')
    plt.xlabel('Class')
    plt.ylabel('Count')
    plt.show()
