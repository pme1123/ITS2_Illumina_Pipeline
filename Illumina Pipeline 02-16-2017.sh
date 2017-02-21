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

WORKINGDIR=~/bmgc_sequences/PZM_TEST  # A copy of your raw reads should be in this directory, to make sure you don't overwrite them
FASTQC_OUT=FASTQC

cd ${WORKINGDIR}
mkdir ${FASTQC_OUT}
module load fastqc  # UMN MSI command
fastqc *.fastq -o ${FASTQC_OUT} -q

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
Documentation: http://www.drive5.com/usearch/manual/cmd_fastq_mergepairs.html

Important flags are below. See http://drive5.com/usearch/manual/merge_options.html
    -log <path>, a directory for pushing log files. Same as -report?
    -fastq_nostagger, since I have 250bp reads and a 400bp target length, sequences should not stagger
    -fastq_maxdiffs, maximum number of mismatches, defaults 5, higher is OK if you have a long overlap (but what does that mean?)
    -fastq_trunctail #, discards reads at the first base with a quality score <= #, default 2, higher often improves merging. Set it based on QC data from fastqc.
    -relabel ABC, for relabelling samples using something human-readable. The operator, @,  provides a nice default (-relabel @ ). It truncates at the first underscore adds a (.) to separate, and a number as read count.

USEARCH will spit out some logs as it runs, as well. See the mediocre documentation at http://drive5.com/usearch/manual/merge_report.html

In practice this step isn't so important for generating OTUs taxa as long as error rates are held low and >65% of reads can be merged. Maxing out read count increases discarding later and chimera flagging rates, anyway. 
DOC


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

# Variables. Also be sure to check the id variable below. Default should work in most cases.
MAXDIFFS=10    # Maximum differences to accept. Higher increases number of merges. Error rates will be taken care of in the filtering steps, so feel free to raise this.
TRUNCTAIL=20  # Minimum quality score of a base at which to truncate. Higher increases merges (be reducing mismatches). Not necessarily good... watch mean expected errors and alignment length
MERGED_OUT="PZM_TEST"

# Make directories for your new, merged reads
mkdir ./merged
mkdir ./merged/stats

usearch9 -fastq_mergepairs *_R1*.fastq \
  -fastqout ./merged/${MERGED_OUT}_merged.fastq \
  -relabel @ \
  -log ./merged/stats/${MERGED_OUT}_merge.log \
  -fastq_nostagger \
  -fastq_maxdiffs ${MAXDIFFS} \
  -fastq_trunctail ${TRUNCTAIL} # the relabel flag changes fastq labels from jibberish to $id.readnumber


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

: <<DOC
00:01 83Mb    100.0% 81.0% merged

Totals:
     43259  Pairs (43.3k)
     35035  Merged (35.0k, 80.99%)
      2464  Alignments with zero diffs (5.70%)
      7720  Too many diffs (> 10) (17.85%)
      3260  Fwd tails Q <= 20 trimmed (7.54%)
     35758  Rev tails Q <= 20 trimmed (82.66%)
       482  No alignment found (1.11%)
         0  Alignment too short (< 16) (0.00%)
        22  Staggered pairs (0.05%) discarded
    180.70  Mean alignment length
    404.31  Mean merged length
      0.52  Mean fwd expected errors
      6.16  Mean rev expected errors
      0.10  Mean merged expected errors  # awesome.

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
usearch9 -fastx_truncate ./merged/${MERGED_OUT}_merged.fastq \
  -stripleft ${LLENGTH} \
  -stripright ${RLENGTH} \
  -fastqout ./merged/${MERGED_OUT}_trimmed.fastq
  
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

#
# 00:01 37Mb    100.0% Processing, 0 (0.0%) too short
#

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
mkdir merged/stats
mkdir merged/discarded

# Filter for quality control using -fastq_filter command http://www.drive5.com/usearch/manual/cmd_fastq_filter.html

usearch9 -fastq_filter ./merged/${MERGED_OUT}_trimmed.fastq\
  -fastq_maxee ${EMAX} \
  -fastaout ./merged/${MERGED_OUT}_filtered.fasta \
  -fastaout_discarded ./merged/discarded/${MERGED_OUT}_discarded.fasta \
  -log ./merged/stats/${MERGED_OUT}_filter.log

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

: <<DOC
# Example Output
00:01 15Mb    100.0% Filtering, 99.3% passed
     35035  Reads (35.0k)                   
       232  Discarded reads with expected errs > 1.00
     34803  Filtered reads (34.8k, 99.3%)
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
$MIN_COUNT=1      # minimum number of occurrances allowed. 2 removes singletons. (not recommended)

## Run fastx_uniques
usearch9 -fastx_uniques ./merged/${MERGED_OUT}_filtered.fasta \
    -fastaout ./merged/${MERGED_OUT}_uniques.fasta \
    -sizeout \
    -relabel Uniq

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

: <<DOC
Sample Output:

00:00 33Mb    100.0% DF                                        
00:00 34Mb   34803 seqs, 3708 uniques, 2473 singletons (66.7%)
00:00 34Mb   Min size 1, median 1, max 17529, avg 9.39
00:00 33Mb    100.0% Writing ./filtered/PZM_TEST_uniques.fasta

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

usearch9 -unoise2 ./merged/${MERGED_OUT}_uniques.fasta \
  -unoise_alpha ${ALPHA} \
  -minampsize ${MINAMP} \
  -fastaout ./merged/${MERGED_OUT}_unoise_defaults.fasta

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

: <<DOC
Example output: 
00:00 208Mb   100.0% Word stats
00:00 208Mb   100.0% Alloc rows
00:00 208Mb   100.0% Build index
00:00 179Mb   100.0% Reading ./filtered/PZM_TEST_uniques.fasta
00:00 153Mb   100.0% 93 amplicons, 9201 bad (size >= 4)       
00:00 159Mb   100.0% 93 good, 0 chimeras               
00:00 159Mb   100.0% Writing amplicons    

DOC

######################################################################################################
######################################################################################################
############                                                                      ####################
############                  DENOISING APPROACH 2: OTU CLUSTERING                ####################
############                         AND READ ASSIGNMENT                          ####################
############                                                                      ####################
######################################################################################################
######################################################################################################
: <<DOC
Alternatively to the UNOISE approach is the 97% sequence similarity approach, which is more defensible. Uses
the `-cluster_otus` command. The algorithm favors OTU centers around high-frequency reads and removes 
chimeras in the process. 

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


DOC
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

#variables for cluster_otus
MINSIZE=2

mkdir ./denoised

usearch9 -cluster_otus ./merged/${MERGED_OUT}_uniques.fasta \
  -minsize ${MINSIZE} \
  -otus ./merged/${MERGED_OUT}_uparse_defaults.fasta \
  -uparseout ./merged/${MERGED_OUT}_uparse_defaults.txt  
  
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#
: <<DOC

00:01 48Mb    100.0% 54 OTUs, 26 chimeras

DOC
######################################################################################################
######################################################################################################
############                                                                      ####################
############                          ASSIGN TAXONOMY                             ####################
############                                                                      ####################
######################################################################################################
######################################################################################################
: <<DOC
Then predict taxonomy using `-utax`. This is a similar algorithm to RDP (naive baysian). It tends to have 
lower errors, but also lower sensitivity, than other algorithms. According to 
(Richardson et al, 2016; DOI:10.1111/1755-0998.12628), optimizing the cutoff can increase sensitivity
without increasing error rates; to do this effectively would require mock datasets. Parameters:
  -utax_cutoff [0, 1], minimum confidence for keeping a taxa. Default 0.9 (generally recommended; 
    however, decreasing to 0.65 might give better results).
  -db ABC.udb, database performed by the makeudb_utax command, see Make ITS2 Unite Database.sh
  -strand {}, options "plus" or "both". Both searches both forward and reverse-compliments. Plus is 
    only forward.
  -id [0, 1], matching OTUs. 97% (0.97) recommended.
  -rdpout PATH, for an output file in RDP format
DOC
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

#variables for -utax
CUTOFF=0.90 #0.9 default; recommended. 
STRAND="both"
MATCH_ID=0.97

PATH_TO_DB=./UNITE/UNITEv7_its2_ref.udb
#PATH_TO_CLUSTERS=./denoised/unoise_defaults.fasta
NAME_OUT=${MERGED_OUT}"_uparse_"${CUTOFF}  # string identifying the output .utax and .txt files
METHOD="uparse"

usearch9 -utax ./merged/${MERGED_OUT}_${METHOD}_defaults.fasta \
  -db ${PATH_TO_DB} \
  -strand ${STRAND} \
  -id ${MATCH_ID} \
  -utax_cutoff ${CUTOFF} \
  -utaxout ./merged/${NAME_OUT}_OTUs.utax \
  -alnout ./merged/${NAME_OUT}_alignment.txt
  
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#  
: <<DOC

00:00 212Mb   100.0% 54 seqs, 64.8% at phylum, 9.3% genus (P > 0.90)  # uparse method
00:00 212Mb   100.0% 93 seqs, 30.1% at phylum, 4.3% genus (P > 0.90)  # unoise method 

00:00 212Mb   100.0% 54 seqs, 64.8% at phylum, 33.3% genus (P > 0.65)  # uparse method
00:00 212Mb   100.0% 93 seqs, 30.1% at phylum, 17.2% genus (P > 0.65)  # unoise method

Overall, uparse matches a higher percentage but also has fewer OTUs:
  UNOISE has 93 OTUs
  UPARSE has 54 OTUs
Therefore, in terms of absolute matches, the methods are equivalent. Plus, regardless of cutoff or method,
I find 9 matches to Glomeromycota

DOC

######################################################################################################
######################################################################################################
############                                                                      ####################
############                               MAKE OTU TABLE                         ####################
############                                                                      ####################
######################################################################################################
######################################################################################################
: <<DOC
Reassign all reads to taxa. This includes singletons, which Edgar argues are often point mutations of legitimate sequences. Reference: http://www.drive5.com/usearch/manual/cmd_usearch_global.html

This uses the `-usearch_global` command. Critical parameters are:
  -strand {}, options plus or both, for forward or both, respectively.
  -id NN, for matching sequences. Recommend 0.97.
  -otutabout PATH.json, for txt file (QIIME format)
  -biomout PATH.json, for BIOM format which is becoming dominant
  -mothur_shared_out PATH, for mothur
  -maxaccepts NN, to increase maximum accepted OTU assignments per sequence beyond 1.

The input file should be the trimmed sequences, no quality filtering. If they don't match, they're no good
anyway. If they're close, they're probably a legitimate error. 
DOC

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

#variables
STRAND="both"
MATCH_ID=0.97

PATH_TO_OTUS=./merged/${MERGED_OUT}_${METHOD}_defaults.fasta  # from clustering OTUs (UNOISE or UPARSE)
NAME_OUT="unoise${CUTOFF}"  # string identifying the output .utax and .txt files

mkdir ./merged/OTU_TABLE

usearch9 -usearch_global ./merged/${MERGED_OUT}_trimmed.fast* \
  -db ${PATH_TO_OTUS} \
  -strand ${STRAND} \
  -id ${MATCH_ID} \
  -log ./merged/OTU_TABLE/${MERGED_OUT}_table_log.txt \
  -otutabout ./merged/OTU_TABLE/${MERGED_OUT}_table.txt \
  -biomout ./merged/OTU_TABLE/${MERGED_OUT}_table.json

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#
  
: <<DOC
00:00 131Mb   100.0% Searching PZM_TEST_trimmed.fastq, 99.5% matched

DOC
