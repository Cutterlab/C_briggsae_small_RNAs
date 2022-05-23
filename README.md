# Temperature-dependent small RNA expression depends on wild genetic backgrounds of *Caenorhabditis briggsae*

Scripts for analyzing small RNA sequencing data in C. briggsae from Fusca et al. (2022), as well as the pseudo-reference genome sequence for strain HK104 used in these analyses. 

## HK104 pseudo-reference genome sequence
HK104 Pseudo-reference Genome.zip contains the files for the pseudo-reference genome sequence for HK104, based off of the reference sequence for AF16. Contains the genome FASTA file (c_briggsae_HK104_NGM_reference.fa), GFF3 files of annotated repeat regions (c_briggsae_HK104_NGM_Alignments_repeat_annotations.gff3) and transposable elements (c_briggsae_HK104_NGM_Alignments_transposable_element_annotations.gff3), and the GTF file of gene annotations (lifted_c_briggsae_HK104_NGM_Alignments_canonical_geneset.gtf).

## Description of scripts
### demultiplex.pl
This Perl script was used to demultiplex our raw sequencing files and remove the 4 nucleotide 5' barcodes. This script was written by Wei Wang.

### genSpecRef.py
This Python script was used to apply SNPs and indels from the HK104 genome to the AF16 reference genome sequence, in order to make a genotype-specific reference genome sequence for HK104 in FASTA format. This script was written by Santiago Sánchez-Ramírez.

### Get_Genes_That_Lifted_To_HK104.R
After using flo to transfer gene annotation coordinates from the AF16 reference genome to our HK104 pseudo-reference genome, this R script was used to filter the results to keep only the genes for which all constituent parts (e.g. exons, UTRs) were successfully lifted from AF16 to HK104. 

### Quantification_Annotation_Script_AF16.R and Quantification_Annotation_Script_HK104.R
These R scripts create genome annotation files to be used with Small_RNA_Counting.R, using the relevant reference genome. Quantification_Annotation_Script_AF16.R is used for the AF16 reference genome, and Quantification_Annotation_Script_HK104.R is used for the HK104 pseudo-reference genome. Both of these scripts are modified from the script [c_elegans_reference_prep.R](https://github.com/ClaycombLab/Charlesworth_2020/blob/master/c_elegans_reference_prep.R), used in [Charlesworth et al. (2021)](https://academic.oup.com/nar/article/49/15/8836/6331683). 

### Script to count all small RNAs of every type

### 22G-RNA_Differential_Expression.R
This R script uses the results from the previous script to count the total number of 22G-RNAs aligning antisense to all transposons, repeats, pseudogenes, and protein-coding genes that are present in both genome annotations, and normalizes this expression data. It then identifies features that are differentially expressed between genotypes and temperatures, as well as those showing a genotype-temperature interaction on 22G-RNA expression. This script also fits linear models to 22G-RNA expression values, for individual features as well as entire chromosomes.

### 22G-RNA_mRNA_Expression_Comparisons.R
This R script performs comparisons between 22G-RNA expression (this study, as calculated in 22G-RNA_Differential_Expression.R) and mRNA expression ([Mark et al. 2019](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.15185)) for 12,623 genes with both types of data available, including separate comparisons for genes targeted by either CSR-1 or WAGO-1. 
