# WYtigger
Identifying autosomal and sex-linked contigs using a machine learning-based classifier

# Materials & Methods

## Overview

This study utilized independent datasets from XY and WZ sex chromosome systems and four different machine learning (ML) models to classify contigs that belong to allosomes. Both systems were processed independently to account for their specific characteristics, leading to tailored data processing and modeling strategies.

## 1. Feature Selection

### Features

- **Coverage**: Coverage is used to identify sex-linked contigs by comparing sequencing read coverage between heterogametic and homogametic individuals.
- **Heterozygosity**: This measures the variation in heterozygous positions to differentiate between sex chromosomes and autosomes.
- **GC Content**: GC content differences between sex chromosomes and autosomes are used to identify allosome-linked contigs.

### Tools and Methods

- **SEX-DETector**, **DiscoverY**, **RADSex**, **FindZX**: Reviewed for feature effectiveness and biological relevance.
- **Coverage Calculation**: Based on WGS data and compared between males and females.
- **Heterozygosity Measurement**: Variants with quality ≥ 30 were used to calculate heterozygosity per contig.
- **GC Content Calculation**: Custom Python script using PySam library.

## 2. Data Collection and Processing

### Datasets

- **Species**: Homo sapiens (XY system) and Gallus domesticus (WZ system).
- **Data Components**: Contig ID, coverage, heterozygosity, GC content, and classification (1 or 0 for allosome).
- **Sample Quality**: Ensured by FastQC v0.12.1 and trimmed with Trimmomatic 0.39.

### Processing

- **Reference Genomes**: Used for aligning de novo assembly samples.
- **Assembly Quality**: Evaluated with QUAST 5.2.0 and aligned with NUCmer 3.1.
- **Dataset Merging**: Integrated using an in-house R script to associate contigs with chromosomes.

## 3. Alignment and Variant Calling

### Tools

- **RepeatMasker**: Masked annotated repeats.
- **Alignment**: Performed using bwa mem v0.7.17, sorted with samtools sort v1.7, and deduplicated with Picard tools V2.18.0.
- **Variant Calling**: Used bcftools 1.19 with mpileup and call options.

## 4. Feature Extraction

- **Coverage**: Calculated with coverM 0.6.1, normalized as TPM, and assessed differential coverage.
- **Heterozygosity**: Extracted from VCF files and calculated per contig.
- **GC Content**: Measured using a Python script with PySam library.

## 5. Model Training

### ML Models

- **Models**: Logistic Regression, Support Vector Machine (SVM), Random Forest, K-Nearest Neighbors.
- **Hyperparameters**: Tuned using grid search procedures.

### Training and Evaluation

- **Cross-Validation**: Used to train and evaluate models with balanced class representation.
- **Class Imbalance**: Addressed using class_weight='balanced' and SMOTE for resampling.

## 6. Performance Assessment

### Metrics

- **Evaluation**: AUC, precision, recall, F1 scores, accuracy, confusion matrices.
- **Technique**: Macro average for class balance in imbalanced datasets.
