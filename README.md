# SAMtools:  A Primer

by Ethan Cerami, Ph.D.

**keywords:**  samtools, next-gen, next-generation, sequencing, bowtie, sam, bam, primer, tutorial, how-to, introduction

## Revisons ##

* 0.1: April 23, 2013:  This document is currently a work in progress.  It is not ready for release.

## About ##

SAMtools is a popular open-source tool, commonly used in next-generation sequence analysis.  This short primer provides an introduction to using SAMtools, and is geared towards those new to next-generation sequence analysis.  The primer is also designed to be self-contained and hands-on, meaning that you only need to install SAMtools, and no other tools, and sample data sets are provided.  Terms in **bold** are also explained in the short glossary at the end of the document.

## Introduction to SAMtools and Next-Generation Sequence Analysis

Next-generation sequencing refers to new, cheaper, high-throughput technologies for sequencing the DNA or RNA of a single sample or a group of samples.  A typical work-flow for next-generation sequence analysis is usually composed of multiple steps:

1. DNA is extracted from a sample.
2. DNA is sequenced.
3. Raw sequencing reads are aligned to a reference genome.
4. The aligned reads are evaluated and visualized. 
5. Genomic variants, including single nucleotide polymorphisms (SNPs) and small insertions and deletions are identified.

SAMtools fits in at steps 4 and 5.  Most importantly, it can process aligned sequence reads, and manipulate them with ease.  For example, it can convert between the two most common file formats (SAM and BAM), sort and index files (for speedy retrieval later), and extract specific genomic regions of interest.  It also enables quality checking of reads, and automatic identification of genomic variants.

## Downloading and Installing SAMtools

To get started with the rest of this primer, download the SAMtools source code from [SourceForge.net](http://sourceforge.net/projects/samtools/files/).  Then, unzip and unpack the distribution to your preferred location &ndash;&ndash; for example, I have placed samtools at:  ~/libraries/samtools-0.1.19/.

Samtools is written in C, compiled with gcc and make, and has only two dependencies:  the [GNU curses library](http://www.gnu.org/software/ncurses/), and the [ZLib compression library](http://zlib.net/).  If you are using a modern variant of Linux or MacOS X, chances are good that you already have these libraries installed.

To build samtools, type:

	make
	
And, add the executables to your path.  For example, I have modified my `.bash_profile` like so:

	export SAMTOOLS_HOME=/Users/ecerami/libraries/samtools-0.1.19
	export PATH=$SAMTOOLS_HOME:$PATH
	export PATH=$SAMTOOLS_HOME/bcftools/:$PATH
	export PATH=$SAMTOOLS_HOME/misc/:$PATH

As an optional, but recommended step, copy the man page for `samtools.1` to one of your man page directories (1). 


## Tutorial:  A Complete Workflow in *E. coli*

To illustrate the use of SAMtools, the remainder of this document focuses on using SAMtools within a complete sample workflow for next-generation sequence analysis.  For simplicity, the tutorial uses a small set of simulated reads from *E. coli*.  I have chosen *E. coli* because its genome is quite small -- 4,649 kilobases, with 4,405 genes, and its entire genome fits into a single **FASTA** file of 4.8 megabytes.  The sample read file is also small -- just 790K -- enabling you to download both within a few minutes.  

If you are new to next-generation sequence analysis, you will soon find that one of the biggest obstacles is just finding and downloading sample data sets, and then downloading the very large reference genomes.  By focusing on *E. coli* and simulated data sets, you can start small, learn the tool sets, and then advance to other organisms and larger sample data sets.

### An Overview

The workflow below is organized into six steps:

1.  Generate a small set of simulated reads for *E. coli*.
2.  Align the reads to the reference *E. coli* genome.
3.  Convert the aligned reads from the **SAM** file format to **BAM**.
4.  Sort and index the **BAM** file.
5.  Identify genomic variants.
6.  Visualize the reads and genomic variants.

Steps 3-6 are focused on the use of SAMtools.  Steps 1-2 require the use of other tools.  Steps 1-2 do, however provide important context, may be helpful in the future, and are therefore described below.  That said, you are not required to actually perform any of steps yourself now, as intermediate files from these steps are available for download.  You can download these intermediate files directly and proceed with steps 3-6.

### Generate Simulated Reads for *E. Coli*

First, we need a small set of sample read data.  A number of tools, including ArtificialFastqGenerator, Mason, and SimSeq, will generate artificial or simulated sequence data for you.  For this tutorial, I chose to use the wgsim tool (created by Heng Li, also the creator of SAMtools).

  

## Footnotes

(1)  If you are uncertain where to copy the samtools man page, or are uncertain where your existing man pages are located, try typing:

	man -w

this will display the current set of man paths.  `samtools.1` is considered a "section 1" man page, reserved for user commands.  You therefore need to copy it into a directory named man1.  For example, on my system, I copied `samtools.1` to:  /usr/local/share/man/man1/.

## Glossary

**bcftools**:  a set of companion tools, currently bundled with samtools, for identifying and filtering genomics variants.

**bowtie**:  widely used, open source alignment software for aligning raw sequence reads to a reference genome.  

**BAM Format**:  binary, compressed format for storing SAM data.

**BCF Format**:  Binary call format.  Binary, compressed format for storing VCF data. 

**FASTA Format**:  text format for storing raw sequence data.  For example, the FASTA file at:  http://www.ncbi.nlm.nih.gov/nuccore/NC_008253 contains entire genome for Escherichia coli 536. 

**FASTQ Format**:  text format for storing raw sequence data along with quality scores for each base.

**Phred quality score**:  a score assigned to each base within a sequence, quantifying the accuracy of the base call.  Higher values correspond to higher accuracy. 

**SAM Flag**:  a single integer value (e.g. 16), which encodes multiple elements of meta-data regarding a read and its alignment.  Elements include: whether the read is one part of a paired-end read, whether the read aligns to the genome, and whether the read aligns to the forward or reverse strand of the genome.  A [useful online utility](http://picard.sourceforge.net/explain-flags.html) decodes a single SAM flag value into plain English.

**SAM Format**:  Text file format for storing sequence alignments against a reference genome.  See also BAM Format.

**samtools**:  widely used, open source command line tool for manipulating SAM/BAM files.  Includes options for converting, sorting, indexing and viewing SAM/BAM files.  The samtools distribution also includes bcftools, a set of command line tools for identifying and filtering genomics variants.  Created by Heng Li, currently of the Broad Institute.

**VCF Format**:  Variant call format.  Text file format for storing genomic variants, including single nucleotide polymorphisms, insertions, deletions and structural rearrangement.

## License
![Creative Commons](http://i.creativecommons.org/l/by/3.0/88x31.png)

SAMtools:  A Primer, by Ethan Cerami is licensed under a [Creative Commons Attribution 3.0 Unported License](http://creativecommons.org/licenses/by/3.0/deed.en_US).
