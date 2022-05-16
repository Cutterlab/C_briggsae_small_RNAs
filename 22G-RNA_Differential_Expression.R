# Attach packages needed for this pipeline
library(ggplot2)
library(gridExtra)
library(bitops)
library(GenomicAlignments)
library(GenomicRanges)
library(rtracklayer)
library(stringr)
library(seqinr)
library(dplyr)
library(scales)
library(ggpubr)
library(ComplexHeatmap)
library(limma)
library(edgeR)
library(VennDiagram)
library(ggcorrplot)
require(reshape2)
require(data.table)


# Load in previously-calculated results
load("~/results.RData")
load("~/resultsSizes.RData")
load("~/WS272.HK104_gene.whole.transposons.repeats.RData")
load("~/WS272.HK104_gene.parts.transposons.repeats.RData")


# The first section of this script calculates small RNA expression of each gene
# using every type of small RNA (of ANY length or first base), instead of just 
# 22G-RNAs. While we did not use these overall small RNA expression counts for 
# further analysis, this code is presented here since some objects created here 
# are reused in our analysis of 22G-RNA expression

# For each sample, add up all of the counts of all small RNAs of any type
# for each gene and then store the results in allSampleGeneCounts
allSampleGeneCounts <- list()
for (name in names(results)) {
  print(name)
  pooledResults <- data.frame()
  for (rnaCategory in names(results[[name]])) {
    if (rnaCategory != "no_feature")
      pooledResults <- rbind(pooledResults, results[[name]][[rnaCategory]])
  }
  
  geneCounts <- aggregate(pooledResults$count, 
                          by = list(pooledResults$gene), FUN = sum)
  colnames(geneCounts) <- c("Gene", "Count")

  
  allSampleGeneCounts[[name]] <- geneCounts
}

# Using allSampleGeneCounts, create a dataframe where every row is a gene, 
# every column is a different sample, and every cell contains the raw
# count of total small RNA expression of any type for that gene and sample
counted_genes <- unique(bind_rows(allSampleGeneCounts)$Gene)
allRawGeneCounts <- data.frame(Gene = counted_genes)

getRawCountFromSample <- function(gene, sample) {
  return(allSampleGeneCounts[[sample]]$Count[allSampleGeneCounts[[sample]]$Gene == gene])
}

listOfRawCounts <- list()

# Go through every sample and for each gene, get the raw expression value of this
# gene for every type of small RNA, or 0 if this gene has no detected expression in this sample
for (sample in names(results)) {
  listOfRawCounts[[sample]] <- as.numeric(lapply(allRawGeneCounts$Gene, getRawCountFromSample, sample))
  listOfRawCounts[[sample]][is.na(listOfRawCounts[[sample]])] <- 0
}
allRawGeneCounts <- cbind(allRawGeneCounts, as.data.frame(listOfRawCounts))

# Remove the quotation marks around each WormBase ID in the set of unlifted genes
# ("unlifted_HK104_genes" was created in the R script Get_Genes_That_Lifted_To_HK104.R)
unlifted_HK104_genes <- gsub("\"", "", unlifted_HK104_genes)

# Subset this dataframe of genes to remove those that did not lift over completely
# from AF16 to HK104 coordinates
lifted_allRawGeneCounts <- allRawGeneCounts[!(allRawGeneCounts$Gene %in% unlifted_HK104_genes), ]
lifted_allRawGeneCounts <- lifted_allRawGeneCounts[lifted_allRawGeneCounts$Gene %in% gene.whole$gene_id, ]


# Add a column to this dataframe classifying the type of gene in each row
get_gene_biotype_from_annotations <- function(gene){
  return(unique(gene.whole$gene_biotype[gene.whole$gene_id == gene])[1])
}
gene_biotypes <- lapply(lifted_allRawGeneCounts$Gene, 
                        get_gene_biotype_from_annotations)
lifted_allRawGeneCounts$Gene_Biotype <- unlist(gene_biotypes)


# There are 4 genes with a biotype of "NA" that are not in the set of lifted-over
# genes, but also not in the set of unlifted genes. To be safe, remove these
# 4 genes
lifted_allRawGeneCounts <- lifted_allRawGeneCounts[!(is.na(lifted_allRawGeneCounts$Gene_Biotype)), ]


# There are also 13 genes that aren't in the annotation object for AF16, but are
# in the annotation object for HK104, so these genes all have 0 counts in all 9
# AF16 replicates. Again, to be safe, remove these 13 genes from the dataset of
# gene expression values

# Load in the annotations for AF16
load("~/WS272.gene.whole.transposons.repeats.RData")
load("~/WS272.gene.parts.transposons.repeats.RData")

lifted_allRawGeneCounts <- lifted_allRawGeneCounts[lifted_allRawGeneCounts$Gene %in% gene.whole$gene_id, ]

# Set rownames of the dataframe to be the name of each gene
row.names(lifted_allRawGeneCounts) <- lifted_allRawGeneCounts$Gene




# The remainder of this script is used to calculate the expression of each gene
# based ONLY on 22G-RNAs, rather than every type of small RNA as was done above.


# For each sample, add up all of the 22G-RNAs that aligned antisense transposons, 
# repeats, pseudogenes, and protein-coding genes, and then store the
# results, classifying small RNAs between 21-23 bases long and beginning with a
# G as "22G" RNAs
load("~/WS272.HK104_gene.whole.transposons.repeats.RData")
load("~/WS272.HK104_gene.parts.transposons.repeats.RData")

allSample22GCounts <- list()
for (name in names(results)) {
  pooledResults <- data.frame()
  for (rnaCategory in c("transposable_element_AS", "repeat_region_AS", 
                        "pseudogene_AS", "protein_coding_AS")) {
      pooledResults <- rbind(pooledResults, results[[name]][[rnaCategory]])
  }
  
  pooledResults <- pooledResults[pooledResults$first == "G" & pooledResults$length %in% c(21, 22, 23), ]
  geneCounts <- aggregate(pooledResults$count, 
                          by = list(pooledResults$gene), FUN = sum)
  colnames(geneCounts) <- c("Gene", "Count")
  
  allSample22GCounts[[name]] <- geneCounts
}


# Using the dataset of 22G-RNA expression, create a dataframe where every row is a gene, 
# every column is a different sample, and every cell contains the raw
# count of 22G RNAs for that gene and sample
counted_22G_genes <- unique(bind_rows(allSample22GCounts)$Gene)
allRaw22GCounts <- data.frame(Gene = counted_22G_genes)


getRaw22GCountFromSample <- function(gene, sample) {
  return(allSample22GCounts[[sample]]$Count[allSample22GCounts[[sample]]$Gene == gene])
}


listOfRaw22GCounts <- list()
# Go through all 18 samples and for each gene, get the raw expression value of this
# gene, or 0 if this gene has no detected expression in this sample, for 22G RNAs
for (sample in names(results)) {
  listOfRaw22GCounts[[sample]] <- as.numeric(lapply(allRaw22GCounts$Gene, getRaw22GCountFromSample, sample))
  listOfRaw22GCounts[[sample]][is.na(listOfRaw22GCounts[[sample]])] <- 0
}
allRaw22GCounts <- cbind(allRaw22GCounts, as.data.frame(listOfRaw22GCounts))


# Subset this dataframe of raw 22G counts to keep only those genes that
# are in both the HK104 and the AF16 annotation file (and are thus possible targets
# for alignment in this pipeline) - these genes that are in both annotations are
# stored in lifted_allRawGeneCounts, which was created in the first section of
# this script
lifted_allRaw22GCounts <- allRaw22GCounts[allRaw22GCounts$Gene %in% rownames(lifted_allRawGeneCounts), ]
rownames(lifted_allRaw22GCounts) <- lifted_allRaw22GCounts$Gene


# Add a column to this dataframe classifying the type of gene in each row
get_gene_biotype_from_annotations <- function(gene){
  return(unique(gene.whole$gene_biotype[gene.whole$gene_id == gene])[1])
}
lifted_allRaw22GCounts$Gene_Biotype <- unlist(lapply(lifted_allRaw22GCounts$Gene, 
                                                      get_gene_biotype_from_annotations))


# Filter this dataframe to remove lowly-expressed genes (i.e. genes without at
# least 1 Read Per Million of expression in at least 3 of the 18 replicates, using 
# the # of perfectly-aligned 22G or 26G reads in each replicate as the total for RPM normalization)
RPM_22G_counts <- cpm(lifted_allRaw22GCounts[, 2:19])


# Get just the genes that have at least 1 RPM of expression in at least 3 replicates
expressed_raw_22G_counts <- lifted_allRaw22GCounts[rowSums(RPM_22G_counts >= 1) >= 3, ]


# Create an edgeR object to store this gene expression data, and annotate each
# replicate based on genotype and temperature. Also reorder the columns so that
# all genotype-temperature combinations are grouped together, and exclude the
# HK30-1 replicate that had issues being sequenced
edgeR_22G_expression_data <- DGEList(expressed_raw_22G_counts[ , c(2, 6, 7, 3, 8, 9, 10, 4, 5, 11:16, 18, 19)])
edgeR_22G_expression_data$samples$group <- c(rep(c("AF"), 9), rep(c("HK"), 8))
edgeR_22G_expression_data$samples$temp <- c(rep(c("14", "20", "30", "14", "20"), each = 3), "30", "30")


# Set the levels of the temperature annotations so that 20 is the baseline, by
# making 20 the first level
edgeR_22G_expression_data$samples$temp <- factor(edgeR_22G_expression_data$samples$temp, 
                                                  levels = c("20", "14", "30"))

# Set library scaling factors to 1
edgeR_22G_expression_data_no_TMM <- calcNormFactors(edgeR_22G_expression_data, method = "none")

# Make a linear model for differential-expression analysis for 22G-RNA expression. 
# Model is: expression = genotype + temperature + genotype X temperature interaction
differential_expression_22G_model <- model.matrix(~ edgeR_22G_expression_data$samples$group + 
                                                    edgeR_22G_expression_data$samples$temp + 
                                                    edgeR_22G_expression_data$samples$group*edgeR_22G_expression_data$samples$temp)

# Rename column names in the model so that they're easier to read
colnames(differential_expression_22G_model) <- gsub("edgeR_22G_expression_data\\$samples\\$", 
                                                    "", colnames(differential_expression_22G_model))

# Normalize 22G-RNA expression data to log2-CPM values using voom
normalized_22G_expression_no_TMM <- voom(edgeR_22G_expression_data_no_TMM, 
                                  design = differential_expression_22G_model)

# Fit normalized 22G-RNA expression data to the linear models using limma, using 
# "trend = T" to specify that this is RNA-seq data
fit_22G_no_TMM <- lmFit(normalized_22G_expression_no_TMM, differential_expression_22G_model)
bayes_fit_22G_no_TMM <- eBayes(fit_22G_no_TMM, trend = T) 


# Get the genes that have a significant genotype effect on 22G expression
genotype_differences_22G_no_TMM <- topTable(bayes_fit_22G_no_TMM, coef = 2, 
                                            adjust.method = "BH", number = nrow(normalized_22G_expression_no_TMM$E))
genotype_diff_genes_22G_no_TMM <- 
  rownames(genotype_differences_22G_no_TMM)[genotype_differences_22G_no_TMM$adj.P.Val <= 0.05]

# Get the genes that have a significant temperature effect on 22G expression
temperature_differences_22G_no_TMM <- topTable(bayes_fit_22G_no_TMM, coef = 3:4, 
                                               adjust.method = "BH", number = nrow(normalized_22G_expression_no_TMM$E))
temperature_diff_genes_22G_no_TMM <- 
  rownames(temperature_differences_22G_no_TMM)[temperature_differences_22G_no_TMM$adj.P.Val <= 0.05]

# Get the genes that have a significant genotype X temperature interaction on 22G expression
interaction_differences_22G_no_TMM <- topTable(bayes_fit_22G_no_TMM, coef = 5:6, 
                                               adjust.method = "BH", number = nrow(normalized_22G_expression_no_TMM$E))
interaction_diff_genes_22G_no_TMM <- 
  rownames(interaction_differences_22G_no_TMM)[interaction_differences_22G_no_TMM$adj.P.Val <= 0.05]


# For 22G expression, split genes into 5 groups: those
# that aren't differentially-expressed at all, those with a genotype x
# temperature interaction, those with a genotype effect only, those with a
# temperature effect only, and those with an additive genotype & temperature
# effect but not an interaction
additive_diff_genes_22G_no_TMM <- 
  intersect(genotype_diff_genes_22G_no_TMM, temperature_diff_genes_22G_no_TMM)
additive_diff_genes_22G_no_TMM <- 
  additive_diff_genes_22G_no_TMM[!(additive_diff_genes_22G_no_TMM %in% interaction_diff_genes_22G_no_TMM)]

genotype_only_diff_genes_22G_no_TMM <- 
  genotype_diff_genes_22G_no_TMM[!(genotype_diff_genes_22G_no_TMM %in% temperature_diff_genes_22G_no_TMM) 
                                 & !(genotype_diff_genes_22G_no_TMM %in% interaction_diff_genes_22G_no_TMM)]

temperature_only_diff_genes_22G_no_TMM <- 
  temperature_diff_genes_22G_no_TMM[!(temperature_diff_genes_22G_no_TMM %in% genotype_diff_genes_22G_no_TMM) 
                                    & !(temperature_diff_genes_22G_no_TMM %in% interaction_diff_genes_22G_no_TMM)]

non_diff_genes_22G_no_TMM <- rownames(interaction_differences_22G_no_TMM)[!(
  rownames(interaction_differences_22G_no_TMM) %in% interaction_diff_genes_22G_no_TMM| 
    rownames(interaction_differences_22G_no_TMM) %in% genotype_diff_genes_22G_no_TMM | 
    rownames(interaction_differences_22G_no_TMM) %in% temperature_diff_genes_22G_no_TMM)]
