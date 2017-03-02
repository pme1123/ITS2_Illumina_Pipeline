# ITS2_Illumina_Pipeline
Pipeline for analyzing ITS2 soil fungal community data. Optimized for output from the University of Minnesota Genomics Center using the MSI resources. 

All steps are tested for copy-paste into the terminal using an interactive session (see code at top of Illumina Pipeline to start an interactive session). If you set your variables it might run through a batch command, but copy-paste only takes ~ 10 minutes for a complete run. 

# Steps
1. Set up access to the latest USEARCH version (currently v9.1). See <a href="https://github.com/pme1123/ITS2_Illumina_Pipeline/blob/master/Setting%20Up%20USEARCH9.txt">Setting up USEARCH9.txt</a>. 
2. Download and prepare your ITS2 database. See <a href="https://github.com/pme1123/ITS2_Illumina_Pipeline/blob/master/Make%20ITS2%20Unite%20Database.sh">Make ITS2 Unite Database.sh</a>.
3. Run through the pipeline. See <a href="https://github.com/pme1123/ITS2_Illumina_Pipeline/blob/master/Illumina%20Pipeline%2002-16-2017.sh">Illumina Pipeline.sh</a>. 