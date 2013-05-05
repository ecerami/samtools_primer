#!/bin/sh

# First, make sure bowtie is installed
if ! type bowtie2 > /dev/null 2>&1; then
    echo "Could not find bowtie2.  Please adjust your path or download from:  http://bowtie-bio.sourceforge.net/bowtie2/index.shtml"
fi

# Align the simulated reads against the reference genome
bowtie2 -x indexes/e_coli -U simulated_reads/sim_reads.fq -S alignments/sim_reads_aligned.sam
