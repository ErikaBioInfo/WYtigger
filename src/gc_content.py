#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This script calculates the GC content of sequences in a BAM file based on regions specified in a BED file. 
GC content is a measure of the proportion of guanine (G) and cytosine (C) bases in the sequences.

Inputs:
- BAM file: Contains alignment data with sequence reads.
- BED file: Specifies regions of interest with chromosome, start, and end positions.
- Output file: Where the GC content results are saved.

Outputs:
- A file with GC content percentages for each region specified in the BED file.
  - Columns:
    - `Region`: The original BED file line specifying the region (chromosome, start, end).
    - `GC_Content`: The percentage of GC bases in the region. If no reads are present, 'No Reads' is reported.

Author:
-Erika Castaneda

Version:
- Python 3.11.5

Usage:
- python gc_content.py <bam_file> <bed_file> <output_file>
"""

import pysam
import sys

def calculate_gc_content(bam_file, bed_file, output_file):
    """
    Calculate GC content for regions specified in a BED file based on reads in a BAM file.

    Parameters:
    - bam_file (str): Path to the BAM file containing alignment data.
    - bed_file (str): Path to the BED file specifying regions of interest.
    - output_file (str): Path to the output file where GC content results will be saved.

    Outputs:
    - A file containing GC content percentages for each region specified in the BED file.
      The output file will have the following columns:
      - `Region`: The original BED file line specifying the region.
      - `GC_Content`: The percentage of GC bases in the region. If no reads are present, 'No Reads' is reported.
    """
    try:
        bam = pysam.AlignmentFile(bam_file, 'rb')
        with open(bed_file) as bf, open(output_file, 'w') as out_file:
            for line in bf:
                line_parts = line.strip().split()
                chr = line_parts[0]
                start = int(line_parts[1])
                end = int(line_parts[2])

                read_data = bam.fetch(chr, start, end)
                total_bases = 0
                gc_bases = 0

                for read in read_data:
                    seq = read.query_sequence
                    total_bases += len(seq)
                    gc_bases += len([x for x in seq if x == 'C' or x == 'G'])

                if total_bases == 0:
                    gc_percent = 'No Reads'
                else:
                    gc_percent = '{0:.2f}%'.format(float(gc_bases) / total_bases * 100)

                out_file.write('{0}\t{1}\n'.format(line.strip(), gc_percent))

    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <bam_file> <bed_file> <output_file>")
        sys.exit(1)

    bam_file = sys.argv[1]
    bed_file = sys.argv[2]
    output_file = sys.argv[3]

    calculate_gc_content(bam_file, bed_file, output_file)
