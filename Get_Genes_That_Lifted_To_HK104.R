# Load in the lifted-over HK104 gene GTF (c_briggsae_HK104_NGM_Alignments_canonical_geneset.gtf)
c_briggsae_HK104_NGM_Alignments_canonical_geneset <- 
  read.delim("~/c_briggsae_HK104_NGM_Alignments_canonical_geneset.gtf", header=FALSE, comment.char="#", quote = "")

# Load in all of the features that didn't lift over from the AF16 GTF to this HK104 GTF (unlifted.gff3)
unlifted <- read.delim("~/unlifted.gff3", header=FALSE, comment.char="#", quote = "")

# Add a column to the HK104 GTF dataframe that gives the WormBase ID of the gene that each row belongs to
c_briggsae_HK104_NGM_Alignments_canonical_geneset$V10 <- 
  unlist(lapply(c_briggsae_HK104_NGM_Alignments_canonical_geneset$V9, function(x) {return(unlist(strsplit(as.character(x), ";"))[1])}))
c_briggsae_HK104_NGM_Alignments_canonical_geneset$V10 <- 
  unlist(lapply(c_briggsae_HK104_NGM_Alignments_canonical_geneset$V10, function(x) {return(unlist(strsplit(as.character(x), " "))[2])}))

# Add a column to the set of unlifted features that gives the WormBase ID of the gene that each row belongs to
unlifted$V10 <- unlist(lapply(unlifted$V9, function(x) {return(unlist(strsplit(as.character(x), ";"))[1])}))
unlifted$V10 <- unlist(lapply(unlifted$V10, function(x) {return(unlist(strsplit(as.character(x), " "))[2])}))

# Get all of the genes that didn't lift over completely (i.e. where either the gene itself, or any of its child features like exons or UTRs, 
# didn't lift over
unlifted_HK104_genes <- unique(unlifted$V10)


# Get all of the genes that lifted over completely to HK104 coordinates 
lifted_HK104_genes <- 
  unique(c_briggsae_HK104_NGM_Alignments_canonical_geneset$V10[!(c_briggsae_HK104_NGM_Alignments_canonical_geneset$V10 %in% unlifted_HK104_genes)])

# Filter the lifted-over HK104 GTF to get only the rows belonging to a gene that lifted over completely
lifted_c_briggsae_HK104_NGM_Alignments_canonical_geneset <- 
  c_briggsae_HK104_NGM_Alignments_canonical_geneset[c_briggsae_HK104_NGM_Alignments_canonical_geneset$V10 %in% lifted_HK104_genes, ]

# Save this new dataframe of completely lifted genes (excluding the new column of gene names I added) as a GTF file
write.table(lifted_c_briggsae_HK104_NGM_Alignments_canonical_geneset[ , 1:9], 
            file = "lifted_c_briggsae_HK104_NGM_Alignments_canonical_geneset.gtf", 
            col.names = F, row.names = F, quote = F, sep = "\t")