#! /bin/bash

: <<DOC
Created by Sonya Erlandson, Peay Lab, Stanford University, January 2017.
Annotated by Patrick Ewing, Jordan Lab, University of Minnesota, February 16, 2017.

This is a USearch-based script to process illumina sequencing data for environmental samples. Unfortunately, the
maximum memory size is 4gb on the free version, which can use up to 24 threads. USEARCH documentation is at http://drive5.com/usearch/manual/uparse_pipeline.html

This is for USEARCH 9.2. MSI currently has USEARCH 8.1 installed, which won't work with this script. See "Setting Up USEARCH9.txt" for access.

The current iteration of this script is for the general fungal primers, ITS4 and 5.8SR primers, which target a 400bp fragment of ITS2. Sequences were generated
on a MiSeq using paired-end reads of 300bp. ~ 80 samples were multiplexed for this. DNA was from corn roots extracted using MoBio's PowerSoil kit.
DOC


# Start your interactive session with mesabi if you want to run this interactively.
#qsub -I -l nodes=1:ppn=8,mem=8gb,walltime=4:00:00

######################################################################################################
######################################################################################################
############                                                                      ####################
############                    INITIAL LOOK AT SEQUENCES                         ####################
############                                                                      ####################
######################################################################################################
######################################################################################################
: <<DOC
First, look at the quality of the sequences using fastqc. This is important for understanding how to set:
- Quality threshold for truncating the tails of reads
- Maximum number of differences allowed for merging read pairs.

Considerations are based on the amount of overlap based on the total fragment length. For example, I expect 400bp fragments
and have a 300bp read length. This means I can set quality thresholds quite high and still have lots of overlap, allowing for
more disagreements.
DOC

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

WORKINGDIR="~/bmgc_sequences/PZM_TEST"  # A copy of your raw reads should be in this directory, to make sure you don't overwrite them
FASTQC_OUT="FASTQC"

cd "${WORKINGDIR}"
mkdir "${FASTQC_OUT}"
module load fastqc  # UMN MSI command
fastqc "*.fastq" -o "${FASTQC_OUT}" -q

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

# Look at the html output using a web browser. Easiest is to ssh in your file manager. Open your file
# manager and in the address bar, type:
ssh://USER@login.msi.umn.edu/home/GROUP/USER/${WORKINGDIR}/${FASTQC_OUT}


######################################################################################################
######################################################################################################
############                                                                      ####################
############                            MERGE READS                               ####################
############                                                                      ####################
######################################################################################################
######################################################################################################
: <<DOC
The -fastq_mergepairs command does the following:
    - Merges reads into a single read
    - Trims primers and adapters

Important flags are below. See http://drive5.com/usearch/manual/merge_options.html
    -log <path>, a directory for pushing log files. Same as -report?
    -fastq_nostagger, since I have 250bp reads and a 400bp target length, sequences should not stagger
    -fastq_maxdiffs, maximum number of mismatches, defaults 5, higher is OK if you have a long overlap (but what does that mean?)
    -fastq_trunctail #, discards reads at the first base with a quality score <= #, default 2, higher often improves merging. Set it based on QC data from fastqc.
    -relabel ABC, for relabelling samples using something human-readable. The operator, @,  provides a nice default (-relabel @ ).

USEARCH will spit out some logs as it runs, as well. See the mediocre documentation at http://drive5.com/usearch/manual/merge_report.html

DOC


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

# Variables. Also be sure to check the id variable below. Default should work in most cases.
MAXDIFFS=7    # Maximum differences to accept. Higher increases number of merges. Error rates will be taken care of in the filtering steps, so feel free to raise this.
TRUNCTAIL=25  # Minimum quality score of a base at which to truncate. Higher increases

# Make directories for your new, merged reads
mkdir "./merged"
mkdir "./merged/stats"

#in raw reads file
#rename files to sample names and merge pairs
for f in "*_R1_*.fastq"; do
  # Pull identifying parts of each file name. See http://tldp.org/LDP/abs/html/string-manipulation.html for string manipulation
  id="${f%%_*}"    # unique part of each sample (excluding _*.fastq). This returns everything before the first "_". Current setup should work well for UMN BMGC illumina sequences.

  # See http://www.drive5.com/usearch/manual/cmd_fastq_mergepairs.html
  usearch9 -fastq_mergepairs "${id}*_R1*.fastq" \
    -fastqout "./merged/${id}_merged.fastq" \
    -relabel "${id}." \
    -log "./merged/stats/${id}_merge.log" \
    -fastq_nostagger \
    -fastq_maxdiffs "${MAXDIFFS}" \
    -fastq_trunctail "${TRUNCTAIL}" # the relabel flag changes fastq labels from jibberish to $id.readnumber

  done

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

: <<DOC
Example Output:

00:02 83Mb    100.0% 87.0% merged  # incidentally, this is ~ the percentage of sequences that have my primers.

Totals:
     43259  Pairs (43.3k)
     37620  Merged (37.6k, 86.96%)
      7692  Alignments with zero diffs (17.78%)
      5045  Too many diffs (> 7) (11.66%)
      7655  Fwd tails Q <= 25 trimmed (17.70%)
     41109  Rev tails Q <= 25 trimmed (95.03%)
         0  Fwd too short (< 64) after tail trimming (0.00%)
         3  Rev too short (< 64) after tail trimming (0.01%)
       569  No alignment found (1.32%)
         0  Alignment too short (< 16) (0.00%)
        22  Staggered pairs (0.05%) discarded
    162.79  Mean alignment length
    404.38  Mean merged length
      0.60  Mean fwd expected errors
      6.45  Mean rev expected errors
      0.12  Mean merged expected errors
DOC

######################################################################################################
######################################################################################################
############                                                                      ####################
############                            TRIM PRIMERS                              ####################
############                                                                      ####################
######################################################################################################
######################################################################################################
: <<DOC
Trim primers, which tend to remove variability within a sequence where they are. My primers are:
5.8SR (forward): TCGATGAAGAACGCAGCG (18 bp)
ITS4 (reverse):  TCCTCCGCTTATTGATATGC (20 bp)

Looking at the fastqc output, ~90% of sequences have primers as part of them. Truncate!

The -fastx_truncate command takes care of this. It truncates N positions from the 5' and 3' ends.
DOC

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

#Variables
LLENGTH=18  # 5' primer length
RLENGTH=20  # 3' primer length

mkdir trimmed

# trims 5' and 3' ends of each sequence and saves the resulting file as {ID}_trimmed.fastq in $WORKINGDIR/trimmed.
for f in merged/*.fastq; do
  id="${f#*/}"
  id="${id%_*}"
  usearch9 -fastx_truncate "$f" -stripleft ${LLENGTH} -stripright ${RLENGTH} \
    -fastqout "./trimmed/${id}_trimmed.fastq"
  done

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

######################################################################################################
######################################################################################################
############                                                                      ####################
############                          QUALITY CONTROL                             ####################
############                                                                      ####################
######################################################################################################
######################################################################################################
: <<DOC
Now filter sequences for quality. Uses the -fastq_filter command. The paper describing the algorithm is here:
<https://academic.oup.com/bioinformatics/article/31/21/3476/194979/Error-filtering-pair-assembly-and-error-correction>

Quality is based on overall sequence quality, recalculated after merging pairs. Basically, bases that agree get higher
quality scores than in the individual forward and reverse reads, and bases that disagree get lower quality scores, but
assigned as the base with the higher quality score.

Flags are:
    -fastq_maxee #, the maximum expected # of errors allowed in a sequence. 1.0 is recommended and should be sufficient.
    -fastq_truncqual #, a floor at which to truncate reads (from 5' end)
    and otheres that perform similar functions.

Documentation at <http://drive5.com/usearch/manual/cmd_fastq_filter.html>

DOC
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

# Variables
EMAX=1.0  # 1.0 is recommended.

# Make new directories for QC filtered
mkdir filtered
mkdir filtered/stats
mkdir filtered/discarded

# Filter for quality control using -fastq_filter command http://www.drive5.com/usearch/manual/cmd_fastq_filter.html
for f in trimmed/*trimmed.fastq; do
  id=${f#*/}
  id=${id%_*}
  usearch9 -fastq_filter $f -fastq_maxee ${EMAX} \
    -fastaout "./filtered/${id}_filtered.fasta" \
    -fastaout_discarded "./filtered/discarded/${id}_discarded.fasta" \
	-log "./filtered/stats/${id}_filter.log"
  done

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

: <<DOC
# Example Output
00:02 15Mb    100.0% Filtering, 98.8% passed
     37620  Reads (37.6k)
       467  Discarded reads with expected errs > 1.00
     37153  Filtered reads (37.2k, 98.8%)
DOC

######################################################################################################
######################################################################################################
############                                                                      ####################
############                   IDENTIFY AND COUNT UNIQUE SEQUENCES                ####################
############                                                                      ####################
######################################################################################################
######################################################################################################

: <<DOC
Now it's time to count unique sequences. USEARCH recommends doing this on the entire run of related samples. This facilitates downstream filtering.

Critical parameters:
  -sizeout, says that size annotations should be added
  -relable "ABC", relabel sequences as ABC#
  -minuniquesize, sets minimum abundance. Default 1; will be taken care of later.
  -tabbedout $PATH, produces a txt file with columns as output.
DOC

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

# Variables
$MIN_COUNT=1      # minimum number of occurrances allowed. 2 removes singletons.

##Concatenate all fasta files into a single file #####
mkdir reads
cat "./filtered/*.fasta" > "./reads/reads.fasta"  # places in $WORKINGDIR

## Run fastx_uniques
usearch9 -fastx_uniques "./reads/reads.fasta"  \
    -fastaout "./reads/uniques.fasta/"          \
    -sizeout        #important!                       \
    -relabel Uniq                          \

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

: <<DOC
Sample Output:

00:01 107Mb   100.0% DF                       
00:01 107Mb  37153 seqs, 3973 uniques, 2672 singletons (67.3%)
00:01 107Mb  Min size 1, median 1, max 18569, avg 9.35
00:01 102Mb   100.0% Writing reads/uniques.fasta

DOC

######################################################################################################
######################################################################################################
############                                                                      ####################
############        DENOISING APPROACH 1: REMOVE MISTAKES, SPURRIOUS OTUs         ####################
############                                                                      ####################
######################################################################################################
######################################################################################################
: <<DOC
THIS STEP IS ILLUMINA-SPECIFIC. 
See http://www.drive5.com/usearch/manual/cmd_unoise2.html
Now clean up the files and remove errors. The clustering algorithm (for forming "OTUs") works best after pooling, but this does decrease sensitivity to close variants. Recommends trying both before and after
pooling. The  uses the unoise2 algorithm and performs:
- sequence error ID and correction
- chimera removal
- PhiX removal
Paper: http://www.biorxiv.org/content/biorxiv/early/2016/10/15/081257.full.pdf
The basic approach is clustering sequences, then finding a centeroid based on abundance. Edgar calls
these clusters *ZOTUs* (zero-radius).

Possible flags:
  -ampout ABC, writes corrected amplicons (including chimeras). Not for downstream analysis
  -minampsize N, minimum abundnace sequences to feed the algorithm (for identifying clusters). Default 4. These won't be removed completely, however. Pooling means that most likely, only repeat singletons and PCR errors 
  -unoise_alpha N, critical parameter. Default 2.

FIRST
To make centeroids, the algorithm uses the skew = (seq count / cluster count) and the distance (sequence 
and centeroid, in differences including indels). These are related through the equation:
  beta(dist) = 1/2exp(alpha*dist + 1)
If skew < beta(dist), then more likely than not the sequence is an error of the centeroid, rather than
a separate cluster.
As ALPHA increases, beta(dist) decreases, meaning that more sequences are included in the cluster. This
increases the chance of getting a bad sequence, but also reduces sensitivity to overly small differences
among clusters. 

SECOND
The same matching criteria are used on the entire dataset (including reads with counts < -minampsize) to
reassign sequences to ZOTUs. 

THIRD
Chimera filtering. Takes parent sequences and query sequences. The parents must have abundances 
> 2x query, under assumption that they experience an extra PCR cycle (assumes perfect PCR efficiency).
Lowering this threshold increases 'chimeras' that are removed. This really only becomes an issue for
high-diversity datasets. Chimeras are detected by setting a max point error threshold introduced by
sequencing rather than PCR. If E(error) < 1 from the first step, then setting the threshold at 0 should be fine (and in fact, I don't think you can change it here).

DOC

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

# variables
ALPHA=2 # higher tends to increase number of clusters (and accepted sequences); lower decreases.
MINAMP=4 # lower increases number of clusters, but has an inconsistent effect on accepted sequences.

mkdir "./denoised"

usearch9 -unoise2 "./reads/uniques.fasta" \
  -unoise_alpha ${ALPHA} \
  -minampsize ${MINAMP} \
  -fastaout "./denoised/unoise_defaults.fasta"

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

: <<DOC
Example output: 
00:00 208Mb   100.0% Word stats
00:00 208Mb   100.0% Alloc rows
00:00 208Mb   100.0% Build index
00:01 179Mb   100.0% Reading ./reads/uniques.fasta
00:01 153Mb   100.0% 94 amplicons, 9984 bad (size >= 4)  # of 37000 seqs; after matching <4 abundances.
00:01 159Mb   100.0% 94 good, 0 chimeras               
00:01 159Mb   100.0% Writing amplicons  

DOC

######################################################################################################
######################################################################################################
############                                                                      ####################
############                  DENOISING APPROACH 2: OTU CLUSTERING                ####################
############                                                                      ####################
######################################################################################################
######################################################################################################
: <<DOC
Alternatively to the UNOISE approach is the 97% sequence similarity approach, which is more defensible. Uses
the `cluster_otus` command. The algorithm favors OTU centers around high-frequency reads and removes chimeras
in the process. 

Critical parameters:
  -minsize N, minimum read count to ID OTUs. Set 2 to discard singletons.
  -otu_radius_pct N, for setting sequence similarity = (1-N). Default=3 (for 97%)
  -otus ABC.fasta, for the output files (FASTA)
  -uparseout ABC.txt, for text output file
  -uparsealnout ABC.txt, for alignment of each query to the references.
Others are available. See documentation http://drive5.com/usearch/manual/cmd_cluster_otus.html
  - parsimony score options
  - alignment parameters
  - hueristics 
  
OTU clustering is followed by assigning reads to OTUs, using the -usearch_global command. This requires a database, such as the one from UNITE 
(2016-11-20 (ver. 7.1): https://unite.ut.ee/sh_files/utax_reference_dataset_20.11.2016.zip)


DOC

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

#variables
MINSIZE=2

mkdir "./denoised"

usearch9 -cluster_otus "./reads/uniques.fasta" \
  -minsize ${MINSIZE} \
  -otus "./denoised/uparse_defaults.fasta" \
  -uparseout "./denoised/uparse_defaults.txt"  
  
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

: <<DOC

00:01 48Mb    100.0% 56 OTUs, 34 chimeras

DOC

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

#Download the unite database from USEARCH (pre-sorted)
DB_URL="http://drive5.com/utax/data/utax_unite_v7.tar.gz"
TARGET="its2"
VERSION=7

mkdir "UNITE"
cd "UNITE"

wget "${DB_URL}"
tar -xvz *.tar.gz
rm *.tar.gz
cd ..
usearch9 -makeudb_utax ./UNITE/utaxref/unite_v${VERSION}/fasta/refdb.fa \
  -output ./UNITE/${TARGET}_ref.udb \
  -taxconfsin ./UNITE/utaxref/unite_v${VERSION}/taxconfs/${TARGET}.tc \
  -report ./UNITE/${TARGET}_report.txt

######################################################################################################
######################################################################################################
############                                                                      ####################
############                             ASSIGN TAXONOMY                          ####################
############                                                                      ####################
######################################################################################################
######################################################################################################

#assign taxonomy to denoised otus - I want to use my taxonomy-assigned otus as my database for mapping reads to otus.
#make taxonomy databases for sintax ITS and 16S
usearch -makeudb_sintax utax_reference_dataset_22.08.2016.fasta -output unite_ITS.udb
usearch -makeudb_sintax rdp_16s_v16_sp.fa -output rdp_16s.udb

#merge fastq reads from before filtering step (merged and renamed, but not filtered)
cat *fastq > reads.fastq


#run sintax algorithm
usearch -sintax reads.fastq -db unite_ITS.udb -tabbedout allsamples.sintax -strand both -sintax_cutoff 0.8


#make otu table
usearch -usearch_global reads.fastq -db denoised.fa -id 0.97 -strand both -otutabout otu_table.txt

