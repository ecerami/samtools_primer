# SAMtools:  A Primer

by Ethan Cerami, Ph.D.

**keywords:**  samtools, next-gen, next-generation, sequencing, bowtie, sam, bam, primer, tutorial, how-to, introduction

## Table of Contents

* [Revisons](#revisions)
* [About](#about)
* [Downloading the Sample Data Files](#downloading-the-sample-data-files)
* [Introduction to SAMtools and Next-Generation Sequence Analysis](#introduction-to-samtools-and-next-generation-sequence-analysis)
* [Downloading and Installing SAMtools](#downloading-and-installing-samtools)
* [Tutorial: A Complete Workflow for Identifying SNPs in E. coli](#tutorial--a-complete-workflow-for-identifying-snps-in-e-coli)
  * [An Overview](#an-overview)
  * [Step 1: Generate Simulated Reads](#step-1--generate-simulated-reads)
  * [Step 2: Align Reads to a Reference Genome](#step-2--align-reads-to-a-reference-genome)
  * [Understanding the SAM Format](#understanding-the-sam-format)
  * [Step 3: Convert the Aligned Reads from SAM to BAM](#step-3--convert-the-aligned-reads-from-sam-to-bam)
  * [Step 4: Sort and Index the BAM File](#step-4--sort-and-index-the-bam-file)
  * [Step 5: Identify Genomic Variants](#step-5--identify-genomic-variants)
  * [Understanding the VCF Format](#understanding-the-vcf-format)
  * [Step 6: Visualize Reads and Genomics Variants](#step-6--visualize-reads-and-genomics-variants)
* [References for Further Reading](#references-for-further-reading)
* [Footnotes](#footnotes)
* [Glossary](#glossary)
* [License](#license)

## Revisions ##

* 0.1: April 23, 2013:  This document is **currently a work in progress.  It is not ready for release.**

## About ##

SAMtools is a popular open-source tool, commonly used in next-generation sequence analysis.  This short primer provides an introduction to using SAMtools, and is geared towards those new to next-generation sequence analysis.  The primer is also designed to be self-contained and hands-on, meaning that you only need to install SAMtools, and no other tools, and sample data sets are provided.  Terms in **bold** are also explained in the short glossary at the end of the document.

## Downloading the Sample Data Files ##

You can access all the sample data files for this primer from the [companion github repository](https://github.com/ecerami/samtools_primer).

From the repository, you can:

* browse and download individual files.
* download a [complete zip file containing everything](https://github.com/ecerami/samtools_primer/archive/master.zip).
* clone the entire repository.  If you have never used Github before, there is a [comprehensive tutorial to get you started](https://help.github.com/articles/set-up-git).

## Introduction to SAMtools and Next-Generation Sequence Analysis

Next-generation sequencing refers to new, cheaper, high-throughput technologies for sequencing the DNA or RNA of a single sample or a group of samples.  A typical work-flow for next-generation sequence analysis is usually composed of multiple steps:

![Figure 1:  Typical Next-Generation Sequencing Workflow.](https://raw.github.com/ecerami/samtools_primer/master/figs/ngs_overview.png "Figure 1:  Typical Next-Generation Sequencing Workflow.")

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

## Tutorial:  A Complete Workflow for Identifying SNPs in E. coli

To illustrate the use of SAMtools, the remainder of this document focuses on using SAMtools within a complete sample workflow for next-generation sequence analysis.  For simplicity, the tutorial uses a small set of simulated reads from *E. coli*.  I have chosen *E. coli* because its genome is quite small -- 4,649 kilobases, with 4,405 genes, and its entire genome fits into a single **FASTA** file of 4.8 megabytes.  The sample read file is also small -- just 790K -- enabling you to download both within a few minutes.  

If you are new to next-generation sequence analysis, you will soon find that one of the biggest obstacles is just finding and downloading sample data sets, and then downloading the very large reference genomes.  By focusing on *E. coli* and simulated data sets, you can start small, learn the tool sets, and then advance to other organisms and larger sample data sets.

### An Overview

The workflow below is organized into six steps (see Figure 2 below).

![Figure 2:  Tutorial Workflow.](https://raw.github.com/ecerami/samtools_primer/master/figs/steps_overview.png "Figure 2:  Tutorial Workflow.")

Steps 3-6 are focused on the use of SAMtools.  Steps 1-2 require the use of other tools.  Steps 1-2 do, however provide important context, may be helpful in the future, and are therefore described below.  That said, you are not required to actually perform any of steps yourself now, as intermediate files from these steps are available for download.  You can download these intermediate files directly and proceed with steps 3-6.

### Step 1:  Generate Simulated Reads

First, we need a small set of sample read data.  A number of tools, including [ArtificialFastqGenerator](http://sourceforge.net/p/artfastqgen/wiki/Home/), [SimSeq](https://github.com/jstjohn/SimSeq), will generate artificial or simulated sequence data for you.  For this tutorial, I chose to use the [wgsim](https://github.com/lh3/wgsim) tool (created by Heng Li, also the creator of SAMtools).

The command line usage for wgsim is:

	wgsim [options] <in.ref.fa> <out.read1.fq> <out.read2.fq>
  
By default, wgsim therefore reads in a reference genome in the **FASTA** format, and generates simulated **paired-end** reads in the **FASTQ** format.

If you specify the full [*e. coli* genome](https://raw.github.com/ecerami/samtools_primer/master/tutorial/genomes/NC_008253.fna), wgsim will generate simulated reads across the entire genome.  However, for the tutorial, I chose to restrict the simulated reads to just the first 1,000 bases of the *e. coli* genome.  To do so, I  extracted the first 1K bases of the *e. coli* genome, placed these within a new file:  [NC_008253_truncated.fna](https://raw.github.com/ecerami/samtools_primer/master/tutorial/genomes/NC_008253_1K.fna), and ran:

	wgsim -N1000 -S1 genomes/NC_008253_truncated.fna output/sim_reads.fq /dev/null
	
Command line options are described below:

* `-N1000`:  directs wgsim to generate 1,000 simulated reads (default is set to:  1,000,000).
* `-S1`:  specifies a seed for the wgsim random number generator;  by specifying a seed value, one can reproducibily create simulated reads.
* `/dev/null`:  wgswim requires that you specify two output files (one for each set of pair-end reads).  However, if you set the second output to `/dev/null`, all data for the second set of reads will be discarded, and you can effectively generate single-end read data only.

wgsim will then output:

	[wgsim] seed = 1
	[wgsim_core] calculating the total length of the reference sequence...
	[wgsim_core] 1 sequences, total length: 1000
	gi|110640213|ref|NC_008253.1|	736	T	G	-

This indicates that wgsim has read in a reference sequence of 1K and has generated simulated reads to support exactly one artificial SNP, a T->G at position 736.  By default, all artificial reads have a read length of 70, and the reads are written in the FASTQ format. For illustrative purposes, the first read is shown below:

	@gi|110640213|ref|NC_008253.1|_418_952_1:0:0_1:0:0_0/1
	CCAGGCAGTGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCATCTGGTAGCGATGAT
	+
	2222222222222222222222222222222222222222222222222222222222222222222222

Note that in the FASTQ format, the first line specifies a unique sequence identifier (usually referred to as the QName), the second line specifies the sequence, and the fourth line specifies the **Phred quality score** for each base.  wgsim does not generate artificial quality scores, and all bases are simply set to 2, indicating that the bases have a 0.01995 probability of being called incorrectly (for additional details, refer to **Phred quality scores** in the glossary.)

You can download the [artificial reads from github](https://raw.github.com/ecerami/samtools_primer/master/tutorial/simulated_reads/sim_reads.fq) if you like, but this is not required for the rest of the tutorial.

### Step 2:  Align Reads to a Reference Genome

The next step is to align the artificial reads to the reference genome for *e. coli.*  For this, I have chosen to use the widely used [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) aligner, from Johns Hopkins University.  Again, you need not actually perform this step to use SAMtools, but it does provide important context, so I have included the details.

To align, use the command:

	bowtie2 -x indexes/e_coli -U simulated_reads/sim_reads.fq -S alignments/sim_reads_aligned.sam
	
This directs bowtie to align the simulated reads against the *e_coli* reference genome, and to output the alignments in the **SAM file format**.  Details regarding each command line option is provided below:

* `-x`:  the location of the reference genome index files.  Bowtie requires that reference genome indexes be downloaded or built, prior to alignment.  I have provided a [pre-built index for *e_coli* on the github repository](https://github.com/ecerami/samtools_primer/tree/master/tutorial/indexes) for this primer.

* `-U`:  list of unpaired reads to align.

* `-S`:  output alignments in the SAM format to the specified file.

For additional details on using bowtie2, refer to the [online manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml).

### Understanding the SAM Format

As SAMtools is primarily concerned with manipulating SAM files, it is useful to take a moment to examine the sample SAM file generated by bowtie, and to dive into the details of the SAM file format itself.  The first six lines from the bowtie SAM file are extracted below:

![Figure 2:  SAM File Format Example.](https://raw.github.com/ecerami/samtools_primer/master/figs/sam_format_example.png "Figure 3:  SAM File Format Example.")

SAM files consist of two types of lines:  headers and alignments.  Headers begin with @, and provide meta-data regarding the entire alignment file.  Alignments begin with any character except @, and describe a single alignment of a sequence read against the reference genome.

For now, we ignore the header-lines, auto-generated by bowtie, and focus instead on the alignments.  Two factors are most important to note.  First, each read in a FASTQ file may align to multiple regions within a reference genome, and an individual read can therefore result in multiple alignments.  In the SAM format, each of these alignments is reported on a separate line.  Second, each alignment has 11 mandatory fields, followed by a variable number of optional fields.  Each of the fields is described below in the table below:

<table>
	<tr>
		<th>Field Name</th>
		<th>Description</th>
		<th>Example from the e_coli SAM file</th>
	</tr>
	<tr>
		<td>QNAME</td>
		<td>Unique identifier of the read;  derived from the original FASTQ file.</td>
		<td>gi|110640213|ref|NC_008253.1|_418_952_1:0:0_1:0:0_0/1</td>
	</tr>
	<tr>
		<td>FLAG</td>
		<td>a single integer value (e.g. 16), which encodes multiple elements of meta-data regarding a read and its alignment.  Elements include: whether the read is one part of a paired-end read, whether the read aligns to the genome, and whether the read aligns to the forward or reverse strand of the genome.  A <a href="http://picard.sourceforge.net/explain-flags.html">useful online utility</a> decodes a single SAM flag value into plain English.
		</td>
		<td>16:  indicates that the read maps to the reverse strand of the genome.
		</td>
	</tr>
	<tr>
		<td>RNAME</td>
		<td>Reference genome identifier.  For organisms with multiple chromosomes, the RNAME is usually the chromosome number;  for example, in human, an RNAME of "chr3" indicates that the read aligns to chromosome 3.  For organisms with a single chromosome, the RNAME refers to the unique identifier associated with the full genome;  for example, in <i>e. coli</i>, the RNAME is set to:  gi|110640213|ref|NC_008253.1|.
		</td>
		<td>gi|110640213|ref|NC_008253.1|</td>
	</tr>
	<tr>
		<td>POS</td>
		<td>Left-most position within the reference genome where the alignment occurs.</td>
		<td>418</td>
	</tr>
	<tr>
		<td>MAPQ</td>
		<td>Quality of the genome mapping.  The MAPQ field uses a <B>phred-scaled probability</b> value to indicate the probability that the mapping to the genome is incorrect.  Higher values indicate increased confidence in the mapping.
		</td>
		<td>
		42:  indicates low probability (0.0000630957) that the mapping is incorrect.
		</td>
	<tr>
		<td>CIGAR</td>
		<td>A compact string that (partially) summarizes the alignment of the raw sequence read to the reference genome.  Three core abbreviations are used:  M for alignment match;  I for insertion;  and D for Deletion.  For example, a CIGAR string of 5M2I63M indicates that the first 5 base pairs of the read align to the reference, followed by 2 base pairs, which are unique to the read, and not in the reference genome, followed by an additional 63 base pairs of alignment.  Note that the CIGAR string is a partial representation of the alignment because M indicates an alignment match, but this could be indicative of an exact sequence match or a mismatch [2].  If you would like to determine match v. mismatch, you can consult the optional MD field, detailed below.
		</td>
		<td>
		70M:  70 base pairs within the read match the reference genome.
		</td>
	</tr>
	<tr>
		<td>RNEXT</td>
		<td>Reference genome identifier where the mate pair aligns.  Only applicable when processing paired-end sequencing data.  For example, an alignment with RNAME of "chr3" and RNEXT of "chr4" indicates that the mate and its pair span two chromosomes, indicating a possible structural rearrangement.  A value of * indicates that information is not available.
		</td>
		<td>
		the <i>e. coli</i> simulated data is single read sequencing data, and all RNEXT values are therefore set to * (no information available).
		</td>
	</tr>
	<tr>
		<td>PNEXT</td>
		<td>Position with the reference genome, where the second mate pair aligns.  As with RNEXT, this field is only applicable when processing paired-end sequencing data.  A value of 0 indicates that information is not available.
		</td>
		<td>
		the <i>e. coli</i> simulated data is single read sequencing data, and all RNEXT values are therefore set to 0 (no information available).
		</td>
	</tr>
	<tr>
		<td>TLEN</td>
		<td>Template Length.  Only applicable for paired-end sequencing data, TLEN is the size of the original DNA or RNA fragment, determined by examining both of the <b>paired-mates</b>, and counting bases from the left-most aligned base to the right-most aligned base.  A value of 0 indicates that TLEN information is not available.
		</td>
		<td>
		the <i>e. coli</i> simulated data is single read sequencing data, and all RNEXT values are therefore set to 0 (no information available).
		</td>
	</tr>
	<tr>
		<td>SEQ</td>
		<td>the raw sequence, as originally defined in the FASTQ file.</td>
		<td>CCAGGCAGTGGCAGGTGGCCACCG...</td>
	</tr>
	<tr>
		<td>QUAL</td>
		<td>The <b>Phred quality score</b> for each base, as originally defined in the FASTQ file.</td>
		<td>222222222222222222222222...:  as noted above, wgsim does not generate artificial quality scores, and all bases are simply set to 2, indicating that the bases have a 0.01995 probability of being called incorrectly.</td>
</table>

Following the 11 mandatory fields, each alignment can also have a variable number of optional fields.  For example, the bowtie generated SAM file for *e. coli* includes 8 additional fields for each alignment.  Most of these are specific to bowtie, and are beyond the scope of this document [3].  However, two fields are part of the official SAM format, and are worth noting.  These are detailed in the table below.  Note that optional SAM fields always use the format:  FIELD_NAME:DATA_TYPE:VALUE, where the DATA_TYPE of i indicates integer values and Z indicates string values [4].

<TABLE>
	<tr>
		<th>Field Name</th>
		<th>Description</th>
		<th>Examples from the e_coli SAM file</th>
	</tr>
	<tr>
		<td>NM</td>
		<td>The number of base pair changes required to change the aligned sequence to the reference sequence, known formally as the Edit Distance.
		</td>
		<td>
		NM:i:0 indicates a perfect match between the aligned sequence and the reference sequence.<BR>
		NM:i:1 indicates that the aligned sequence contain one mismatch.
		</td>
	</tr>
	<tr>
		<td>MD</td>
		<td>
		A string which summarizes the mismatch positions between the aligned read and the reference genome.
		</td>
		<td>
		MD:Z:8G61 indicates a single base pair mismatch.  Specifically, the aligned read matches the first 8 bases of the reference, after which it fails to match a G in the reference sequence, followed by 61 exact matches to the reference.
		</td>
</table>

### Step 3:  Convert the Aligned Reads from SAM to BAM

Now that you understand where alignments come from, know how to interpret SAM fields, and have a sample [SAM file to play with](tutorial/alignments/sim_reads_aligned.sam), you are finally read to tackle SAMtools.  As described above, our goal is to identify the set of genomic variants within the *e. coli* data set.  To do so, our first step is to convert the SAM file to BAM.  This is an important prerequisite, as all the downstream steps, including the identification of genomic variants and visualization of reads, require BAM as input.

To convert from SAM to BAM, use the samtools `view` command:

	samtools view -b -S -o alignments/sim_reads_aligned.bam alignments/sim_reads_aligned.sam

* `-b`: indicates that the output is BAM.
* `-S`:  indicates that the input is SAM.
* `-o`:  specified the name of the output file.

BAM files are stored in a compressed, binary format, and cannot be viewed directly.  However, you can use the same `view` command to display all aligments.  For example, running:

	samtools view alignments/sim_reads_aligned.bam | more

will display all your reads in the unix `more` paginated style.

You can also use `view` to only display reads which match your specific filtering criteria.  For example:

	samtools view -f 4 alignments/sim_reads_aligned.bam | more

* `-f INT`:  extracts only those reads which match the specified SAM flag.  In this case, we filter for only those reads with flag value of 4 = read fails to map to the reference genome.

or:

	samtools view -F 4 alignments/sim_reads_aligned.bam | more

* `-F INT`:  removes reads which match the specified SAM flag.  In this case, we remove reads with flag value of 4 = read fails to map to the reference genome, and display all other reads.

You can also try out the `-c` option, which does not output the reads, but rather outputs the number of reads which match your criteria.  For example:

	samtools view -c -f 4 alignments/sim_reads_aligned.bam

indicates that 34 of our artificial reads failed to align to the reference genome.

Finally, you can use the `-q` parameter to indicate a minimal quality mapping filter.  For example:

	samtools view -q 42 -c alignments/sim_reads_aligned.bam
	
outputs the total number of aligned reads (819) that have a mapping quality score of 42 or higher.

### Step 4:  Sort and Index the BAM File

The next step is to sort and index the BAM file.  There are two options for sorting BAM files:  by read name (`-n`), and by genomic location (default).  As our goal is to call genomic variants, and this requires that we "pile-up" all matching reads within a specific genomic location, we sort by location:

	samtools sort alignments/sim_reads_aligned.bam alignments/sim_reads_aligned.sorted

This reads the BAM file from `alignments/sim_reads_aligned.bam` and writes the sorted file to:  `alignments/sim_reads_aligned.sorted.bam`.

Once you have sorted your BAM file, you can then index it.  This enables tools, including SAMtools itself, and genomic viewers, such as the Broad Institute's Integrative Genomics Viewer (IGV), to perform efficient random access on the BAM file, resulting in greatly improved performance.  To do so, run:

	samtools index alignments/sim_reads_aligned.sorted.bam

This reads in the sorted BAM file, and creates a BAM index file with the `.bai` extension:  `sim_reads_aligned.sorted.bam.bai`.  Note that if you attempt to index a BAM file which has not been sorted, you will receive an error and indexing will fail.

<table width=100%>
	<tr bgcolor=#cccccc>
		<th align=left>
		A note about "Big Data"	
		</th>
	</tr>
	<tr>
		<td>
		If you are already conversant with command line tools and possess basic programming skills,
		you will quickly realize that there is nothing magical about next-generation sequence analysis.
		You just need to learn the terms, the file forms, and the tools, and tie them all together.
		However, there is one extra skill you will need to cultivate:  patience.  Once you move beyond
		toy examples, and enter the world of real sequence data, all the concepts remain the same, but
		all the files become much bigger, and everything takes longer -- usually a lot longer.  For
		example, sorting and indexing a single human exome BAM file can easily take 2-4 hours,
		depending on read depth and your computer set-up.  Worse, you may encounter an error after
		two hours of processing, then make an adjustment, and wait another two hours to see if the
		problem is fixed.  It's a bit like rolling giant boulders up a mountain, having the 
		occasional one slip through, and starting all over from the beginning.
		<br><br>
		Many bioinformaticians hedge their patience by running next-generation sequencing pipelines
		on large distributed compute clusters, so that they can at least move dozens of boulders at
		once -- however, this is a topic for another primer.  Until then, my best advice is to get
		the concepts right and start with small files.  Then, pack your patience and work your way
		up to larger data sets.
		</td>
	</tr>
</table>

### Step 5:  Identify Genomic Variants

Now, we are ready to identify the genomic variants from our reads.  Doing so requires two steps, and
while one can easily pipe these two steps together, I have broken them out into two distinct steps below for 
improved clarity.

The first step is to use the SAMtools `mpileup` command to calculate the **genotype likelihoods** supported by the aligned reads in our sample:

	samtools mpileup -g -f genomes/NC_008253.fna alignments/sim_reads_aligned.sorted.bam > variants/sim_variants.bcf

* `-g`:  directs SAMtools to output genotype likelihoods in the **binary call format (BCF)**.  This is a compressed binary format.
* `-f`:  directs SAMtools to use the specified reference genome.  A reference genome must be specified, and here we specify the reference genome for *e. coli*.

The `mpileup` command automatically scans every position supported by an aligned read, computes all the possible genotypes supported by these reads, and then computes the probability that each of these genotypes is truly present in our sample.

In our case, SAMtools scans the first 1000 bases in the *e. coli* genome, and the weighs the collective evidence of all reads to infer the true genotype.  For example, position 42 (reference = G) has 24 reads with a G and two reads with a T (total **read depth** = 26).  In this case, SAMtools concludes with high probability that the sample has a genotype of G, and that T reads are likely due to sequencing errors.  

By contrast, position 736 (reference = T) has 2 reads with a C and 66 reads with a G (total **read depth** = 68).  In this case, SAMtools concludes with high probability that the sample has a genotype of G.  For each position assayed, SAMtools computes all the possible genotypes, and then outputs all the results in the **binary call format (BCF)**.

The second step is to use `bcftools`, which is packaged with SAMtools, but located in the `bcftools` directory:

	bcftools view -c -v variants/sim_variants.bcf > variants/sim_variants.vcf

* `-c`:  directs bcftools to call SNPs.
* `-v`:  directs bcftools to only output potential variants

The `bcftools view` command uses the genotype likelihoods generated from the previous step to call SNPs and indels, and outputs the all identified variants in the **variant call format (VFC)**, the file format created for the [1000 Genomes Project](http://www.1000genomes.org/), and now widely used to represent genomic variants.

### Understanding the VCF Format

If you take a peak at our newly generated `sim_variants.vcf` file, and scroll down to the first line not starting with a # symbol, you will see our lone single nucleotide variant (SNV).  And, if all has gone well, it should match the SNV at position 736 output by `wgsim` way back in Step 1!  Congratulations!

Of course, to truly understand the variants identified in your sample, you must take a slight detour and understand the basics of the Variant Call Format (VCF).  We will not cover all the details here, but just enough to give you a basic understanding.

First, just like SAM files, VCF files contain header lines, which begin with # or ## symbols, and data lines, which do not begin with a # symbol.  The header contains important information, such as the name of the program which generated the file, the VCF format version number, the reference genome name, and information regarding individual columns within the file.

Each data line is required to have eight mandatory columns, including CHROM reference, POSition within the chromosome, REFerence sequence, and ALT sequence.  Directly after these eight is a FORMAT column used to describe the format for all subsequent sample columns.  VCF can contain information regarding multiple samples, and each sample gets its own column.

Within our VCF file, the FORMAT column is set to "PL".  If you search within the header, you can find that PL is defined as:

	##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
  
These are the same genotype likelihoods described above in step 5.

With this knowledge in hand, let's try to decipher the lone SNV we have identified:

![VCF Excerpt](https://raw.github.com/ecerami/samtools_primer/master/figs/vcf.png "VCF Excerpt")

According to the data columns, the SNV was identified in *e. coli* at position 736, where the reference is T, and our sequencing data supports two possible alternative sequences:  G and C.  INFO DP=68 indicates that the position has a **read depth** of 68, and the FORMAT column indicates that all samples (in this case, just our one) will include a list of genotype likelihoods.  In our case, this results in the cryptic string:  "64,178,0,66,147,60".

Based on the reference and the alternative sequences, SAMtools automatically calculates the following genotypes, in this exact order [5][6]:

<table>
<tr>
<th>Genotype</th>
<th>Phred Scaled Score</th>
<th>Likelihood p-value</th>
<tr>
<td>TT</td>
<td>64</td>
<td>3.981-07</td>
</tr>
<tr>
<td>TG</td>
<td>178</td>
<td>2.512e-18</td></tr>
<tr>
<td>GG</td>
<td>0</td>
<td>1.0</td>
</tr>
<tr>
<td>TC</td>
<td>66</td>
<td>2.512e-07</td></tr>
<tr>
<td>GC</td>
<td>147</td>
<td>1.995e-15</td></tr>
<tr>
<td>CC</td>
<td>60</td>
<td>1.000e-06</td></tr>
</table>
        
This indicates that the GG genotype is very likely the true genotype in your sample.  The reads in your sample therefore support a Single Nucleotide Variant (SNV) from T to G at position 736.

### Step 6:  Visualize Reads and Genomics Variants

## References for Further Reading

* Bowtie2 [Home](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), and [Online Manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml).

* Danecek P. et al.  **The variant call format and VCFtools**. Bioinformatics. 2011 Aug 1;27(15):2156-8.  [[PubMed](http://www.ncbi.nlm.nih.gov/pubmed/21653522)].

* Cock PJ et al.  **The Sanger FASTQ file format for sequences with quality scores, and the Solexa/Illumina FASTQ variants.**    Nucleic Acids Res. 2010 Apr;38(6)  [[PubMed](http://www.ncbi.nlm.nih.gov/pubmed/20015970)]

* Li H. et al.  **The Sequence Alignment/Map format and SAMtools.**  Bioinformatics. 2009 Aug 15;25(16):2078-9.  [[PubMed](http://www.ncbi.nlm.nih.gov/pubmed/19505943)]

* SAMTools [Home](http://samtools.sourceforge.net/), and [Manpage](http://samtools.sourceforge.net/samtools.shtml).

* The SAM Format Specification [[PDF](http://samtools.sourceforge.net/SAM1.pdf)].

* wgsim [repository on github](https://github.com/lh3/wgsim).

## Footnotes

[1]  If you are uncertain where to copy the samtools man page, or are uncertain where your existing man pages are located, try typing:

	man -w

this will display the current set of man paths.  `samtools.1` is considered a "section 1" man page, reserved for user commands.  You therefore need to copy it into a directory named man1.  For example, on my system, I copied `samtools.1` to:  /usr/local/share/man/man1/.

[2]  The ambiguous nature of M within the CIGAR string has caused some confusion within the community, and biostars.org has an [interesting thread on whether the CIGAR format itself should be redefined](http://www.biostars.org/p/17043/).

[3]  If you are curious, you can look up all the bowtie specific SAM fields in [bowtie2 manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) (refer to "optional fields" within the SAM output section.)

[4]  The complete list of data types and predefined optional fields are provided in the official [SAM Format Specification](http://samtools.sourceforge.net/SAM1.pdf).

[5]  How exactly is the order of genotypes determined?  This is formally specified in the [VCF version 4.1 specification](http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41).  According to the specification, each allele is assigned an index value, starting at 0.  For example, in a VCF row with reference T and alternative alleles G and C, we set T=0, G=1, and C=2.

With three alleles, and the assumption that we are working with a diploid organism, there are (3 choose 2 or P(3,2)) permutations for a total of 6 permutations.  Each permutation is described by their index values.  For example, TT (both alleles are T) is described as (0,0), whereas CC (both alleles are C) is described as (2,2).  The ordering of each permutations is then determined by the equation:  `F(j,k) = (k*(k+1)/2)+j`, where j,k refer to allele index values.  For example, TT (0,0) = 0, whereas CC (2,2) = (2*(2+1)/2) + 2 = 5.

[6] But, wait!  The *e. coli* genome is haploid, i.e. there is only a single set of genes.  And, there is therefore never the possibility of having two alleles at the same position.  SAMtools, however is designed to work with diploid genomes, and always calculates genotypes based on the assumption of a diploid genome.  To get around this assumption in our simulated case, we can assume that the *e. coli* reads come from a pool of samples, and the genotypes represent alternatives alleles present in the pool.

## Glossary

**alignment**:  the mapping of a raw sequence read to a location within a reference genome.  The mapping occurs because the sequences within the raw read match or align to sequences within the reference genome.  Alignment information is stored in the **SAM** or **BAM** file formats.

**bcftools**:  a set of companion tools, currently bundled with samtools, for identifying and filtering genomics variants.

**bowtie**:  widely used, open source alignment software for aligning raw sequence reads to a reference genome.  

**BAM Format**:  binary, compressed format for storing SAM data.

**BCF Format**:  Binary call format.  Binary, compressed format for storing VCF data. 

**CIGAR String**:  Compact Idiosyncratic Gapped Alignment Report.  A compact string that (partially) summarizes the alignment of the raw sequence read to the reference genome.  Three core abbreviations are used:  M for alignment match;  I for insertion;  and D for Deletion.  For example, a CIGAR string of 5M2I63M indicates that the first 5 base pairs of the read align to the reference, followed by 2 base pairs, which are unique to the read, and not in the reference genome, followed by an additional 63 base pairs of alignment.

**FASTA Format**:  text format for storing raw sequence data.  For example, the FASTA file at:  [http://www.ncbi.nlm.nih.gov/nuccore/NC_008253](http://www.ncbi.nlm.nih.gov/nuccore/NC_008253) contains entire genome for Escherichia coli 536. 

**FASTQ Format**:  text format for storing raw sequence data along with quality scores for each base.

**genotype likelihood**:  the probability that a specific genotype is present in the sample of interest.  Genotype likelihoods are usually expressed as a **Phred-scaled probability**, where P = 10 ^ (-Q/10).  For example, if the genotype TT (both alleles are T) at position 1,299,132 in human chromosome 12 (reference G) is 37, this translates to a probability of 10^(-37/10) = 0.0001995, meaning that there is very low probability that the reads in your sample support a TT genotype.  On the other hand, a genotype of GG at the same position with a score of 0 translates into a probability of 10^(-0) = 1, indicating extremely high probability that your sample contains a homozygous mutation of G to T.

**mate-pair**:  in paired-end sequencing, both ends of a single DNA or RNA fragment are sequenced, but the intermediate region is not sequenced.  The two ends which are sequenced form a pair, are are frequently referred to as mate-pairs.

**QNAME**:  unique identifier of a raw sequence read (also known as the Query Name).  Used in **FASTQ** and **SAM** files.

**paired-end sequencing**:  sequencing process where both ends of a single DNA or RNA fragment are sequenced, but the intermediate region is not sequenced.  Particularly useful for identifying structural rearrangements, including gene fusions.

**Phred-scaled probability**:  a scaled value (Q) used to compactly summarize a probability, where P = 10^(-Q/10).  For example, a Phred Q score of 10 translates to probability (P) = 10^(-10/10) = 0.1.  Phred-scaled probabilities are common in next-generation sequencing, and are used to represent multiple types of quality metrics, including quality of base calls and quality of mappings, and probabilities associated with specific genotypes.  The name Phred refers to the original Phred base-calling software, which first used and developed the scale.

**Phred quality score**:  a score assigned to each base within a sequence, quantifying the probability that the base was called incorrectly.  Scores use a **Phred-scaled probability** metric.  For example, a Phred Q score of 10 translates to P=10^(-10/10) = 0.1, indicating that the base has a 0.1 probability of being incorrect.  Higher Phred score  correspond to higher accuracy.  In the **FASTQ format**, Phred scores are represented as single ASCII letters.  For details on translating between Phred scores and ASCII values, refer to [Table 1 of this useful blog post from Damian Gregory Allis](http://www.somewhereville.com/?p=1508).  The name Phred refers to the original Phred base-calling software, originally developed at Washington University.

**read-length**:  the number of base pairs that are sequenced in an individual sequence read.

**read-depth**:  the number of sequence reads that pile up at the same genomic location.  For example, 30X read-depth coverage indicates that the genomic location is covered by 30 independent sequencing reads.  Increased read-depth translates into higher confidence for calling genomic variants.

**RNAME**:  reference genome identifier (also known as the Reference Name).  Within a SAM formatted file, the RNAME identifies the reference genome where the raw read aligns.

**SAM Flag**:  a single integer value (e.g. 16), which encodes multiple elements of meta-data regarding a read and its alignment.  Elements include: whether the read is one part of a paired-end read, whether the read aligns to the genome, and whether the read aligns to the forward or reverse strand of the genome.  A [useful online utility](http://picard.sourceforge.net/explain-flags.html) decodes a single SAM flag value into plain English.

**SAM Format**:  Text file format for storing sequence alignments against a reference genome.  See also BAM Format.

**samtools**:  widely used, open source command line tool for manipulating SAM/BAM files.  Includes options for converting, sorting, indexing and viewing SAM/BAM files.  The samtools distribution also includes bcftools, a set of command line tools for identifying and filtering genomics variants.  Created by Heng Li, currently of the Broad Institute.

**single-read sequencing**:  sequencing process where only one end of a single DNA or RNA fragment is sequenced.  Contrast with **paired-end** sequencing.

**VCF Format**:  Variant call format.  Text file format for storing genomic variants, including single nucleotide polymorphisms, insertions, deletions and structural rearrangement.

## License
![Creative Commons](http://i.creativecommons.org/l/by/3.0/88x31.png)

SAMtools:  A Primer, by Ethan Cerami is licensed under a [Creative Commons Attribution 3.0 Unported License](http://creativecommons.org/licenses/by/3.0/deed.en_US).
