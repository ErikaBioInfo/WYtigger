# Snakefile

import glob
import os

# Load configuration to get the input folder from config.yaml
configfile: "config.yaml"

# Define the directory with the input samples
INPUT_FOLDER = config["input_folder"]

# Automatically detect all sample files in the input folder
SAMPLES = [os.path.basename(f).replace("_1.fastq.gz", "") for f in glob.glob(os.path.join(INPUT_FOLDER, "*_1.fastq.gz"))]

# Rule to perform quality check with FastQC
rule quality_check:
    input:
        fasta=os.path.join(INPUT_FOLDER, "{sample}_1.fastq.gz")
    output:
        folder="fastqc/{sample}/"
    shell:
        "fastqc {input} -o {output}"

# Rule for quality trimming and adapter clipping with Trimmomatic
rule quality_trimming:
    input:
        sample1=os.path.join(INPUT_FOLDER, "{sample}_1.fastq.gz"),
        sample2=os.path.join(INPUT_FOLDER, "{sample}_2.fastq.gz")
    output:
        paired1="trimmed/{sample}_1_paired.fastq.gz",
        unpaired1="trimmed/{sample}_1_unpaired.fastq.gz",
        paired2="trimmed/{sample}_2_paired.fastq.gz",
        unpaired2="trimmed/{sample}_2_unpaired.fastq.gz"
    params:
        adapters="TruSeq3-PE.fa"
    threads: 10
    shell:
        """
        trimmomatic PE -phred33 -threads {threads} {input.sample1} {input.sample2} \
        {output.paired1} {output.unpaired1} {output.paired2} {output.unpaired2} \
        ILLUMINACLIP:{params.adapters}:2:30:10 LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:15 MINLEN:36
        """

# Rule for repeat removal using RepeatMasker
rule remove_repeats:
    input:
        fasta=os.path.join(INPUT_FOLDER, "de_novo_assambly_sample.fna")
    output:
        masked_fasta="masked/de_novo_assambly_sample_masked.fna"
    params:
        species="your_species"
    threads: 4
    shell:
        "RepeatMasker -species {params.species} -pa {threads} {input}"

# Rule to index the de novo assembly with BWA
rule bwa_index:
    input:
        fasta="masked/de_novo_assambly_sample_masked.fna"
    output:
        index_done="bwa_index.done"
    log:
        "bwa_index.log"
    shell:
        "bwa index {input} 2>&1 > {log} && touch {output.index_done}"

# Rule to map the samples against the de novo assembly using BWA
rule bwa_mapping:
    input:
        index_done="bwa_index.done",
        sample1="trimmed/{sample}_1_paired.fastq.gz",
        sample2="trimmed/{sample}_2_paired.fastq.gz",
        fasta="masked/de_novo_assambly_sample_masked.fna"
    output:
        sam="mapped/{sample}_mapped.sam"
    threads: 20
    shell:
        "bwa mem -t {threads} {input.fasta} {input.sample1} {input.sample2} > {output.sam}"

# Rule to sort SAM files with Samtools
rule samtools_sort:
    input:
        sam="mapped/{sample}_mapped.sam"
    output:
        bam="mapped/{sample}_mapped.bam"
    threads: 20
    shell:
        "samtools sort -@ {threads} -O bam -o {output.bam} {input.sam}"

# Rule to mark and remove duplicates using Picard
rule mark_duplicates:
    input:
        bam="mapped/{sample}_mapped.bam"
    output:
        bam="deduplicated/{sample}_no_duplicates.bam",
        metrics="deduplicated/{sample}_marked_duplicates_metrics.txt"
    shell:
        "java -jar picard.jar MarkDuplicates I={input.bam} O={output.bam} M={output.metrics} REMOVE_DUPLICATES=true"

# Rule to calculate coverage using CoverM
rule coverage_calculation:
    input:
        bam="deduplicated/{sample}_no_duplicates.bam"
    output:
        coverage="coverage/{sample}_coverage.csv"
    shell:
        "coverm contig --methods tpm --bam-files {input.bam} -o {output.coverage}"

# Rule for heterozygosity calculation using BCFtools and an in-house script
rule heterozygosity_calculation:
    input:
        fasta="masked/de_novo_assambly_sample_masked.fna",
        bams=expand("deduplicated/{sample}_no_duplicates.bam", sample=SAMPLES)
    output:
        variants="variants/variants.vcf",
        filtered_variants="variants/filtered_variants.vcf",
        genotypes="variants/genotypes.txt",
        heterozygosity="heterozygosity.txt"
    shell:
        """
        bcftools mpileup -Ou -f {input.fasta} {input.bams} | \
        bcftools call -mv -Ov -o {output.variants}
        bcftools filter -i 'QUAL >= 30' -Ov -o {output.filtered_variants} {output.variants}
        bcftools query -f '%CHROM\\t%POS\\t[%GT\\t]\\n' {output.filtered_variants} > {output.genotypes}
        python calculate_heterozygosity.py {output.genotypes} > {output.heterozygosity}
        """

# Rule for GC content calculation using an in-house Python script
rule gc_content_calculation:
    input:
        bam="deduplicated/{sample}_no_duplicates.bam",
        bed="masked/de_novo_assambly_sample_masked.fna.bed"
    output:
        gc_content="gc_content/{sample}_gc_content.txt"
    shell:
        "python gc_content.py {input.bam} {input.bed} > {output.gc_content}"


