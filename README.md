# Temperature-dependent small RNA expression depends on wild genetic backgrounds of *Caenorhabditis briggsae*

Scripts for analyzing small RNA sequencing data in C. briggsae from Fusca et al. (2022)

# Description of scripts
### demultiplex.pl
This Perl script was used to demultiplex our raw sequencing files and remove the 4 nucleotide 5' barcodes. This script was written by Wei Wang.


### Get_Genes_That_Lifted_To_HK104.R
After using flo to transfer gene annotation coordinates from the AF16 reference genome to our HK104 pseudo-reference genome, this script was used to filter the results to keep only the genes for which all constituent parts (e.g. exons, UTRs) were successfully lifted from AF16 to HK104. 
