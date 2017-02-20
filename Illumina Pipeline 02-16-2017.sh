#! /bin/bash

: <<DOC
Created by Sonya Erlandson, Peay Lab, Stanford University, January 2017.
Annotated by Patrick Ewing, Jordan Lab, University of Minnesota, February 16, 2017.

This is a USearch-based script to process illumina sequencing data for environmental samples. Unfortunately, the
maximum memory size is 4gb on the free version, which can use up to 24 threads. USEARCH documentation is at http://drive5.com/usearch/manual/uparse_pipeline.html

This is for USEARCH 9.1. MSI currently has USEARCH 8.1 installed, which won't work with this script. See "Setting Up USEARCH9.txt" for access.

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
mkdir "merged"
mkdir "merged/stats"

#in raw reads file
#rename files to sample names and merge pairs
for f in "*_R1_*.fastq"; do
  # Pull identifying parts of each file name. See http://tldp.org/LDP/abs/html/string-manipulation.html for string manipulation
  id="${f%%_*}"    # unique part of each sample (excluding _*.fastq). This returns everything before the first "_". Current setup should work well for UMN BMGC illumina sequences.

  # See http://www.drive5.com/usearch/manual/cmd_fastq_mergepairs.html
  usearch9 -fastq_mergepairs "${id}*_R1*.fastq" \
    -fastqout "merged/${id}_merged.fastq" \
    -relabel "${id}." \
    -log "merged/stats/${id}_merge.log" \
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
    -fastqout "trimmed/${id}_trimmed.fastq"
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
    -fastaout filtered/${id}_filtered.fasta \
    -fastaout_discarded filtered/discarded/${id}_discarded.fasta \
	-log filtered/stats/${id}_filter.log
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
Now it's time to count unique sequences. USEARCH recommends doing this on the entire run of related samples. This
facilitates downstream filtering:
- singletons
- chimeras
-
DOC

# Variables
$MIN_COUNT=1      # minimum number of occurrances allowed. 2 removes singletons.

##### Concatenate all fasta files into a single file #####
mkdir reads

cat filtered/*.fasta > reads/reads.fasta  # places in $WORKINGDIR

usearch9 -fastx_uniques reads/reads.fasta  \
    -fastaout reads/uniques.fasta          \
    -sizeout                               \
    -relabel Uniq                          \


##### Find unique sequences and tabulate abundance for each unique sequence #####
usearch -fastx_uniques reads.fasta -fastaout uniques.fasta -sizeout
# Documentation: http://www.drive5.com/usearch/manual/cmd_fastx_uniques.html

#stats for cedar creek fungi
#00:14 1.2Gb   100.0% Reading reads.fasta
#00:23 1.3Gb   100.0% DF
#00:23 1.4Gb  4171086 seqs, 890625 uniques, 715745 singletons (80.4%)
#00:23 1.4Gb  Min size 1, median 1, max 108294, avg 4.68
#00:43 1.3Gb   100.0% Writing uniques.fasta

######################################################################################################
######################################################################################################
############                                                                      ####################
############                     REMOVE MISTAKES, SPURRIOUS OTUs                  ####################
############                                                                      ####################
######################################################################################################
######################################################################################################

# denoise samples - remove sequencing error, chimeras, PhiX. Correct abundances.
# UNOISE is Illumina-specific. See http://www.drive5.com/usearch/manual/unoise_algo.html
usearch -unoise2 uniques.fasta -tabbedout out.txt -fastaout denoised.fasta

#stats

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

