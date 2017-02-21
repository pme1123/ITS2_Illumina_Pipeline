: <<DOC

To run -utax to assign taxa to OTUs, you need an input database of OTUs to match against. This script uses the UNITE database and usearch methods to make 
the database. Note that different methods are required to run `-sintax`, which is an alternative (less tested) way of doing the same thing. See 
http://drive5.com/usearch/manual/utax_or_sintax.html As of writing, neither algorithm is published. Independent analysis 
(Richardson et al, 2016; DOI:10.1111/1755-0998.12628) shows UTAX is less sensitive than RDP, but also has a lower error rate, for plant ID at the genus 
level. Larger databases increase match rate for both. This observation agrees withvUSEARCH's documentation: http://drive5.com/usearch/manual/tax_se_its2.html

OTU clustering is followed by assigning reads to OTUs, using the -usearch_global command. This requires a database, such as the one from UNITE.
(2016-11-20 (ver. 7.1): https://unite.ut.ee/sh_files/utax_reference_dataset_20.11.2016.zip). I recommend downloading the pre-trained dataset from 
USEARCH's website: http://drive5.com/usearch/manual/utax_downloads.html

The main -makeudb_utax documentation is at http://drive5.com/usearch/manual/cmd_makeudb_utax.html.

The -taxconfsin file contains parameters for taxonomy assignments so you don't need to retrain the algorithm.
You can train the alogrithm if you want. See http://drive5.com/usearch/manual/taxconfs_file.html

DOC

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#

# variables
DB_URL="http://drive5.com/utax/data/utax_unite_v7.tar.gz"
TARGET="its2"
VERSION=7
WORKINGDIR="~/bmgc_sequences/PZM_TEST"  # one level above where you want to place the database. 

# make working directory
mkdir "UNITE"
cd "UNITE"

# download, extract, clean up
wget "${DB_URL}"
tar -xvz utax*.tar.gz
rm utax*.tar.gz

# make the database
usearch9 -makeudb_utax "./utaxref/unite_v${VERSION}/fasta/refdb.fa" \
  -output "./UNITEv${VERSION}_${TARGET}_ref.udb" \
  -taxconfsin "./utaxref/unite_v${VERSION}/taxconfs/${TARGET}.tc" \
  -report "./UNITEv${VERSION}_${TARGET}_report.txt"  # probably doesn't show anything useful for UNITE as input

# return to working directory
cd ..

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#
  
: <<DOC
Output:

00:01 76Mb    100.0% Reading ./UNITE/utaxref/unite_v7/fasta/refdb.fa
00:01 42Mb    100.0% Converting to upper case                       
00:02 43Mb    100.0% Word stats              
00:02 43Mb    100.0% Alloc rows
00:04 163Mb   100.0% Build index
00:04 166Mb   100.0% Initialize taxonomy data
00:04 167Mb   100.0% Building name table     
00:04 167Mb  12840 names, tax levels min 1, avg 5.4, max 7
00:04 167Mb  Buffers (42004 seqs)
00:04 184Mb   100.0% Seqs

DOC
