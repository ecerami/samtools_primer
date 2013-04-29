#!/bin/sh

# First, make sure bowtie is installed
if ! type bowtie > /dev/null 2>&1; then
    echo "Could not find bowtie.  Please adjust your path or download from:  http://bowtie-bio.sourceforge.net/index.shtml"
fi

# Align the simulated reads against the reference genome
bowtie -S e_coli simulated_reads/sim_reads.fq alignments/sim_reads_aligned.sam