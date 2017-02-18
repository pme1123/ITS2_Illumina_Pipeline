#! /bin/bash

# Created by Sonya Erlandson, Peay Lab, Stanford University, January 2017.
# Annotated by Patrick Ewing, Jordan Lab, University of Minnesota, February 16, 2017.

# This is a USearch-based script to process illumina sequencing data for environmental samples. Unfortunately, the
# maximum memory size is 4gb on the free version, which can use up to 24 threads. USEARCH documentation is at http://drive5.com/usearch/manual/uparse_pipeline.html

# Start your interactive session with mesabi for most of this.
qsub -I -l nodes=1:ppn=4,mem=8gb,walltime=4:00:00  #4

# Currently, MSI has USEARCH version 8.1.1756. This script is for version 9.2, which you can add to your home directory.
# Alternatively, add /home/jordannr/ewing069/Programs to your $PATH. usearch is called usearch9 (to distinguish from MSI's version).

# From your terminal (unix, requires `rsync`)
rsync -e ssh -avz local_path_to/usearch username@login.msi.umn.edu:/home/group/username

# Then, switch to MSI
# Add the path to usearch to your $PATH so you can execute it
nano ~/.profile  # edit your .profile configuration file
# Add the following to the end of the file
if [ -d "$HOME/path_to_usearch" ] ; then
        PATH="$PATH:$HOME/path_to_usearch"
fi

# make it executable
cd ~/path_to_usearch
chown 770 usearch  # gives entire group permissions to read/write/access.

######################################################################################################
######################################################################################################
############                                                                      ####################
############                            MERGE READS                               ####################
############                                                                      ####################
######################################################################################################
######################################################################################################

# the -fastq_mergepairs command does the following:
#   - Merges reads into a single read
#   - Trims primers and adapters
#   

# Make directories for your new, merged reads
WORKINGDIR="V1_PZM"  # A copy of your raw reads should be in this directory

mkdir $WORKINGDIR/merged
mkdir $WORKINGDIR/merged/stats

#in raw reads file
#rename files to sample names and merge pairs
for f in *_R1_*.fastq; do
    # Pull identifying parts of each file name. See http://tldp.org/LDP/abs/html/string-manipulation.html for string manipulation
	id=${f%%_*}    # unique part of each sample (excluding _*.fastq). This returns everything before the first "_". Current setup should work well for UMN BMGC illumina sequences. 
	
	# See http://www.drive5.com/usearch/manual/cmd_fastq_mergepairs.html
	usearch9 -fastq_mergepairs ${id}*_R1*.fastq -fastqout merged/${id}_merged.fastq \
	    -relabel "$id." -log merged/stats/${id}_merge.log -fastq_nostagger
    # the relabel flag changes fastq labels from jibberish to $id.readnumber
	
	done

#00:00 175Mb   100.0% 33.5% merged

#Totals:
#     43259  Pairs (43.3k)
#     14473  Merged (14.5k, 33.46%)
#       411  Alignments with zero diffs (0.95%)
#     28349  Too many diffs (> 5) (65.53%)
#         0  Fwd tails Q <= 2 trimmed (0.00%)
#        37  Rev tails Q <= 2 trimmed (0.09%)
#       437  No alignment found (1.01%)
#         0  Alignment too short (< 16) (0.00%)
#        22  Staggered pairs (0.05%) merged & trimmed
#    197.33  Mean alignment length
#    404.24  Mean merged length
#      0.31  Mean fwd expected errors
#      4.70  Mean rev expected errors
#      0.07  Mean merged expected errors


######################################################################################################
######################################################################################################
############                                                                      ####################
############                          QUALITY CONTROL                             ####################
############                                                                      ####################
######################################################################################################
######################################################################################################

# Make new directories for QC filtered
mkdir filtered
mkdir filtered/stats
cd merged

# Filter for quality control using -fastq_filter command http://www.drive5.com/usearch/manual/cmd_fastq_filter.html
for f in *_merged.fastq; do
	id=${f%%_*}
	
	usearch9 -fastq_filter ${id}_merged.fastq -fastq_maxee 1.0 -fastaout ../filtered/${id}_filtered.fasta -log ../filtered/stats/${id}_filter.log
	
	done

######################################################################################################
######################################################################################################
############                                                                      ####################
############                   IDENTIFY AND COUNT UNIQUE SEQUENCES                ####################
############                                                                      ####################
######################################################################################################
######################################################################################################

##### Concatenate all fasta files into a single file #####

# with working directory = $WORKINGDIR/filtered
cat *.fa > ../reads.fasta  # places in $WORKINGDIR

# change current directory back to $WORKINGDIR
cd ..

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

