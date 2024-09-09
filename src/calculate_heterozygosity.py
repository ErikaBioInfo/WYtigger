#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Calculate heterozygosity for each contig based on genotype data.

    Heterozygosity is defined as 1 minus the squared proportion of heterozygous genotypes.
    The input file should be in tab-separated format with columns for Contig_id, Position, and genotypes.

    Parameters:
    - input_file (str): Path to the input file in tab-separated format.
      The file should have columns: Contig_id, Position, and one or more genotype columns.
    - output_file (str): Path to the output CSV file where results will be saved.

    Outputs:
    - A CSV file containing two columns:
      - `Contig_ID`: Identifier for the contig.
      - `Heterozygosity`: Calculated heterozygosity for the contig.

    Author:
    - Erika Castaneda

    Version: Python 3.11.5
    
    Usage:
        python calculate_heterozygosity.py input_file.tsv -o output_file.csv

"""

import argparse
import pandas as pd

def calculate_heterozygosity(input_file, output_file):
    """
    Calculate heterozygosity for each contig based on genotype data.
    
    Parameters:
    - input_file (str): Path to the input file in tab-separated format.
    - output_file (str): Path to the output CSV file where results will be saved.
    """
    # Read the input file
    df = pd.read_csv(input_file, sep='\t', header=None)
    
    # Assign column names
    df.columns = ['Contig_id', 'Position'] + [f'Genotype_{i}' for i in range(1, len(df.columns) - 1)]

    # Initialize dictionaries to count homozygotes and heterozygotes per contig
    homozygote_count = {}
    heterozygote_count = {}

    # Iterate over the rows of the DataFrame
    for _, row in df.iterrows():
        contig = row['Contig_id']
        genotypes = row[2:]  # Adjusted to dynamically capture all genotype columns
        
        # Iterate over the genotypes in the row
        for genotype in genotypes:
            if isinstance(genotype, str) and genotype != './.':
                alleles = genotype.split('/')

                # Ensure there are at least two alleles before unpacking
                if len(alleles) >= 2:
                    allele_1, allele_2 = alleles

                    # Increment the count of homozygotes and heterozygotes per contig
                    homozygote_count[contig] = homozygote_count.get(contig, 0) + (allele_1 == allele_2)
                    heterozygote_count[contig] = heterozygote_count.get(contig, 0) + (allele_1 != allele_2)

    # Calculate heterozygosity and store results in a DataFrame
    result_data = []
    for contig in homozygote_count:
        total_genotypes = homozygote_count[contig] + heterozygote_count[contig]
        if total_genotypes > 0:
            heterozygosity = 1 - (heterozygote_count[contig] / total_genotypes) ** 2
            result_data.append([contig, heterozygosity])
        else:
            result_data.append([contig, None])  # Handle cases with no genotypes

    # Save the results to a CSV file
    result_df = pd.DataFrame(result_data, columns=['Contig_ID', 'Heterozygosity'])
    result_df.to_csv(output_file, index=False)

if __name__ == "__main__":
    # Configure the command-line argument parser
    parser = argparse.ArgumentParser(description='Calculate heterozygosity for contigs.')
    parser.add_argument('input_file', help='Input file (format: Contig_id Position Genotype).')
    parser.add_argument('-o', '--output_file', default='output.csv', help='Output CSV file.')

    # Parse the arguments
    args = parser.parse_args()

    # Call the function to calculate heterozygosity
    calculate_heterozygosity(args.input_file, args.output_file)
