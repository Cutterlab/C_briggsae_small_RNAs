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
After using flo to transfer gene annotation coordinates from the AF16 reference genome to our HK104 pseudo-reference genome, this script was used to filter the results to keep only the genes for which all constituent parts (e.g. exons, UTRs) were successfully lifted from AF16 to HK104. 

### Script to prepare objects for the counting script (for both AF16 and HK104)

### Script to count all small RNAs of every type

### 22G-RNA_Differential_Expression.R
This script uses the results from the previous script to count the total number of 22G-RNAs aligning antisense to all transposons, repeats, pseudogenes, and protein-coding genes that are present in both genome annotations, and normalizes this expression data. It then identifies genes that are differentially expressed between genotypes and temperatures, as well as those showing a genotype-temperature interaction on 22G-RNA expression.
