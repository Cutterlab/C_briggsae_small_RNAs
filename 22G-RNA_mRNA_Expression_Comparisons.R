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


# log2-CPM normalize the dataframe of all 16,916 genes with 22G-RNA expression (rather than just
# the 14,283 genes with enough 22G expression to use for differential-expression analysis)
normalized_all_22G_expression <- cpm(lifted_allRaw22GCounts[,  c(2, 6, 7, 3, 8, 9, 10, 4, 5, 11:16, 18, 19)], log = T)


# Load in mRNA expression data from  Mark et al. 2019 - the file being loaded in 
# here is Supplementary File S1 from that paper
Mark_mRNA_expression_data <- read.csv("C:/Users/Daniel/Downloads/mec15185-sup-0002-files1 (1).csv", 
                                      stringsAsFactors=FALSE)

# Load in mappings between C. briggsae WormBase IDs and gene names from WormBase 
# release WS272 (https://downloads.wormbase.org/species/c_briggsae/PRJNA10731/annotation/geneIDs/c_briggsae.PRJNA10731.WS272.geneIDs.txt.gz)
Cbriggsae_WS272_geneIDs <- read.csv("C:/Users/Daniel/Desktop/c_briggsae.PRJNA10731.WS272.geneIDs.txt", 
                                    header=FALSE, stringsAsFactors=FALSE)

# Add WormBase IDs to the dataframe of mRNA expression data from Mark et al., based
# on the list of gene names from WormBase
Mark_mRNA_expression_data <- merge(Mark_mRNA_expression_data, Cbriggsae_WS272_geneIDs[ , c(2, 4)], 
                                   by.x = "GeneID", by.y = "V4", all.x = T)


# Get just the genes in the mRNA expression data that are also in the dataframe
# of 16,916 genes with 22G expression data, and just the genes in the 22G-RNA
# dataframe that are also in that new mRNA expression dataframe
Mark_mRNA_expression_data_22G_genes <-
  Mark_mRNA_expression_data[Mark_mRNA_expression_data$V2 %in%
                              rownames(normalized_all_22G_expression), ]

twentyTwoG_expression_data_mRNA_genes <-
  normalized_all_22G_expression[rownames(normalized_all_22G_expression) %in%
                                  Mark_mRNA_expression_data_22G_genes$V2, ]

# Reorder the rows of the mRNA dataframe so that they're in the same order as the 22G dataframe
rownames(Mark_mRNA_expression_data_22G_genes) <- Mark_mRNA_expression_data_22G_genes$V2
Mark_mRNA_expression_data_22G_genes <- 
  Mark_mRNA_expression_data_22G_genes[rownames(twentyTwoG_expression_data_mRNA_genes), ]


# Get the average log2-CPM 22G expression for the 12,623 genes with 22G-RNA expression, 
# for all 6 genotype-temperature combinations
AF14_22G_sRNA <- apply(twentyTwoG_expression_data_mRNA_genes[, 1:3], 1, mean)
AF20_22G_sRNA <- apply(twentyTwoG_expression_data_mRNA_genes[, 4:6], 1, mean)
AF30_22G_sRNA <- apply(twentyTwoG_expression_data_mRNA_genes[, 7:9], 1, mean)
HK14_22G_sRNA <- apply(twentyTwoG_expression_data_mRNA_genes[, 10:12], 1, mean)
HK20_22G_sRNA <- apply(twentyTwoG_expression_data_mRNA_genes[, 13:15], 1, mean)
HK30_22G_sRNA <- apply(twentyTwoG_expression_data_mRNA_genes[, 16:17], 1, mean)

# Get the average log2-CPM mRNA expression for the 12,623 genes with 22G-RNA expression, 
# for all 6 genotype-temperature combinations
AF14_22G_mRNA <- apply(Mark_mRNA_expression_data_22G_genes[, 2:4], 1, mean)
AF20_22G_mRNA <- apply(Mark_mRNA_expression_data_22G_genes[, 5:7], 1, mean)
AF30_22G_mRNA <- apply(Mark_mRNA_expression_data_22G_genes[, 8:10], 1, mean)
HK14_22G_mRNA <- apply(Mark_mRNA_expression_data_22G_genes[, 11:13], 1, mean)
HK20_22G_mRNA <- apply(Mark_mRNA_expression_data_22G_genes[, 14:16], 1, mean)
HK30_22G_mRNA <- apply(Mark_mRNA_expression_data_22G_genes[, 17:19], 1, mean)

# Get the genes in the top 10% of 22G-RNA expression for each of the 6 genotype-
# temperature combinations
top_AF14_22G_sRNA <- AF14_22G_sRNA[AF14_22G_sRNA >= quantile(AF14_22G_sRNA, probs = seq(0, 1, 0.1))["90%"]]
top_AF20_22G_sRNA <- AF20_22G_sRNA[AF20_22G_sRNA >= quantile(AF20_22G_sRNA, probs = seq(0, 1, 0.1))["90%"]]
top_AF30_22G_sRNA <- AF30_22G_sRNA[AF30_22G_sRNA >= quantile(AF30_22G_sRNA, probs = seq(0, 1, 0.1))["90%"]]
top_HK14_22G_sRNA <- HK14_22G_sRNA[HK14_22G_sRNA >= quantile(HK14_22G_sRNA, probs = seq(0, 1, 0.1))["90%"]]
top_HK20_22G_sRNA <- HK20_22G_sRNA[HK20_22G_sRNA >= quantile(HK20_22G_sRNA, probs = seq(0, 1, 0.1))["90%"]]
top_HK30_22G_sRNA <- HK30_22G_sRNA[HK30_22G_sRNA >= quantile(HK30_22G_sRNA, probs = seq(0, 1, 0.1))["90%"]]

top_AF14_22G_mRNA <- AF14_22G_mRNA[names(top_AF14_22G_sRNA)]
top_AF20_22G_mRNA <- AF20_22G_mRNA[names(top_AF20_22G_sRNA)]
top_AF30_22G_mRNA <- AF30_22G_mRNA[names(top_AF30_22G_sRNA)]
top_HK14_22G_mRNA <- HK14_22G_mRNA[names(top_HK14_22G_sRNA)]
top_HK20_22G_mRNA <- HK20_22G_mRNA[names(top_HK20_22G_sRNA)]
top_HK30_22G_mRNA <- HK30_22G_mRNA[names(top_HK30_22G_sRNA)]


# Combine the average 22G-RNA and mRNA expression values for all 12,623 genes into 1 dataframe
combined_22G_mRNA_expression <- cbind(AF14_22G_sRNA, AF20_22G_sRNA, AF30_22G_sRNA, HK14_22G_sRNA, 
                                      HK20_22G_sRNA, HK30_22G_sRNA, AF14_22G_mRNA, AF20_22G_mRNA, 
                                      AF30_22G_mRNA, HK14_22G_mRNA, HK20_22G_mRNA, HK30_22G_mRNA)

# Combine the average 22G-RNA and mRNA expression values for just the genes in the top 
# 10% of 22G-RNA expression into 1 dataframe
combined_top_22G_mRNA_expression <- cbind(top_AF14_22G_sRNA, top_AF20_22G_sRNA, 
                                          top_AF30_22G_sRNA, top_HK14_22G_sRNA, 
                                          top_HK20_22G_sRNA, top_HK30_22G_sRNA, 
                                          top_AF14_22G_mRNA, top_AF20_22G_mRNA, 
                                          top_AF30_22G_mRNA, top_HK14_22G_mRNA, 
                                          top_HK20_22G_mRNA, top_HK30_22G_mRNA)

# Get just the correlations between 22G and mRNA expression within the same 
# genotype-temperature replicate, for just the genes in the top 10% 
# of 22G expression
top_22G_mRNA_correlations <- diag(cor(combined_top_22G_mRNA_expression, 
                                      method = "spearman")[1:6, 7:12])
names(top_22G_mRNA_correlations) <- c("AF14", "AF20", "AF30", "HK14", "HK20", "HK30")


# Get just the correlations between 22G and mRNA expression within the same 
# genotype-temperature replicate, for all 12,623 genes with 22G-RNA
# expression data available
all_22G_mRNA_correlations <- diag(cor(combined_22G_mRNA_expression, 
                                      method = "spearman")[1:6, 7:12])
names(all_22G_mRNA_correlations) <- c("AF14", "AF20", "AF30", "HK14", "HK20", "HK30")

# Calculate the significance of the correlation between 22G-RNA and mRNA expression
# in each genotype-temperature combination, for both the set of all 12,623 genes
# and only the genes in the top 10% of 22G-RNA expression
cor.test(AF14_22G_mRNA, AF14_22G_sRNA, method = "spearman")
cor.test(AF20_22G_mRNA, AF20_22G_sRNA, method = "spearman")
cor.test(AF30_22G_mRNA, AF30_22G_sRNA, method = "spearman")
cor.test(HK14_22G_mRNA, HK14_22G_sRNA, method = "spearman")
cor.test(HK20_22G_mRNA, HK20_22G_sRNA, method = "spearman")
cor.test(HK30_22G_mRNA, HK30_22G_sRNA, method = "spearman")
cor.test(top_AF14_22G_mRNA, top_AF14_22G_sRNA, method = "spearman")
cor.test(top_AF20_22G_mRNA, top_AF20_22G_sRNA, method = "spearman")
cor.test(top_AF30_22G_mRNA, top_AF30_22G_sRNA, method = "spearman")
cor.test(top_HK14_22G_mRNA, top_HK14_22G_sRNA, method = "spearman")
cor.test(top_HK20_22G_mRNA, top_HK20_22G_sRNA, method = "spearman")
cor.test(top_HK30_22G_mRNA, top_HK30_22G_sRNA, method = "spearman")


# Get the genes showing different patterns of mRNA expression (genotype-only, temperature-
# only, additive genotype & temperature, genotype-temperature inteaction, and no effect),
# for comparison with the genes showing those same effects on 22G-RNA expression
mRNA_interaction_genes <- 
  Mark_mRNA_expression_data$V2[!(is.na(Mark_mRNA_expression_data$V2)) & Mark_mRNA_expression_data$DEcategory == "gxe_genes"]

mRNA_genotype_only_genes <- 
  Mark_mRNA_expression_data$V2[!(is.na(Mark_mRNA_expression_data$V2)) & Mark_mRNA_expression_data$DEcategory == "g_genes"]

mRNA_temperature_only_genes <- 
  Mark_mRNA_expression_data$V2[!(is.na(Mark_mRNA_expression_data$V2)) & Mark_mRNA_expression_data$DEcategory == "e_genes"]

mRNA_additive_genes <- 
  Mark_mRNA_expression_data$V2[!(is.na(Mark_mRNA_expression_data$V2)) & Mark_mRNA_expression_data$DEcategory == "gne_genes"]

mRNA_no_effect_genes <- 
  Mark_mRNA_expression_data$V2[!(is.na(Mark_mRNA_expression_data$V2)) & Mark_mRNA_expression_data$DEcategory == "no_de_genes"]


# Function to compare two Spearman's correlation coefficients from 22G-RNA : mRNA
# comparisons and see if the correlation coefficients are significantly different
# from each other, using a z-test after Fisher's r-to-z transformation:
zTestCorrelationcoefficients <- function(mRNA_1, sRNA_1, mRNA_2, sRNA_2) {
  
  corr_1 <- cor(mRNA_1, sRNA_1, method = "spearman")
  corr_2 <- cor(mRNA_2, sRNA_2, method = "spearman")
  
  z_corr_1 <- 0.5 * log((1 + corr_1) / (1 - corr_1))
  z_corr_2 <- 0.5 * log((1 + corr_2) / (1 - corr_2))
  
  z_score <- (z_corr_1 - z_corr_2) / (((length(mRNA_1) - 3) ^ -1 + (length(mRNA_2) - 3) ^ -1)^0.5)
  
  p_value <- 2 * pnorm(abs(z_score), mean = 0, sd = 1, lower.tail = F)
  
  return(p_value)
}

# Example of running the above function
zTestCorrelationcoefficients(HK20_22G_mRNA, HK20_22G_sRNA, HK14_22G_mRNA, HK14_22G_sRNA)


# Get just the 22G-RNA expression and mRNA expression for the CSR-1 targets with
# data available
AF14_22G_sRNA_CSR1 <- AF14_22G_sRNA[names(AF14_22G_sRNA) %in% briggsae_CSR1_targets]
AF20_22G_sRNA_CSR1 <- AF20_22G_sRNA[names(AF20_22G_sRNA) %in% briggsae_CSR1_targets]
AF30_22G_sRNA_CSR1 <- AF30_22G_sRNA[names(AF30_22G_sRNA) %in% briggsae_CSR1_targets]
HK14_22G_sRNA_CSR1 <- HK14_22G_sRNA[names(HK14_22G_sRNA) %in% briggsae_CSR1_targets]
HK20_22G_sRNA_CSR1 <- HK20_22G_sRNA[names(HK20_22G_sRNA) %in% briggsae_CSR1_targets]
HK30_22G_sRNA_CSR1 <- HK30_22G_sRNA[names(HK30_22G_sRNA) %in% briggsae_CSR1_targets]

AF14_22G_mRNA_CSR1 <- AF14_22G_mRNA[names(AF14_22G_mRNA) %in% briggsae_CSR1_targets]
AF20_22G_mRNA_CSR1 <- AF20_22G_mRNA[names(AF20_22G_mRNA) %in% briggsae_CSR1_targets]
AF30_22G_mRNA_CSR1 <- AF30_22G_mRNA[names(AF30_22G_mRNA) %in% briggsae_CSR1_targets]
HK14_22G_mRNA_CSR1 <- HK14_22G_mRNA[names(HK14_22G_mRNA) %in% briggsae_CSR1_targets]
HK20_22G_mRNA_CSR1 <- HK20_22G_mRNA[names(HK20_22G_mRNA) %in% briggsae_CSR1_targets]
HK30_22G_mRNA_CSR1 <- HK30_22G_mRNA[names(HK30_22G_mRNA) %in% briggsae_CSR1_targets]


# For each genotype-temperature combination, get the 22G-RNA expression and mRNA expression 
# for just the CSR-1 targets in the top 10% of 22G-RNA expression among CSR-1 targets 
# in that genotype-temperature combination
top_AF14_22G_sRNA_CSR1 <- AF14_22G_sRNA_CSR1[AF14_22G_sRNA_CSR1 >= quantile(AF14_22G_sRNA_CSR1, probs = seq(0, 1, 0.1))["90%"]]
top_AF20_22G_sRNA_CSR1 <- AF20_22G_sRNA_CSR1[AF20_22G_sRNA_CSR1 >= quantile(AF20_22G_sRNA_CSR1, probs = seq(0, 1, 0.1))["90%"]]
top_AF30_22G_sRNA_CSR1 <- AF30_22G_sRNA_CSR1[AF30_22G_sRNA_CSR1 >= quantile(AF30_22G_sRNA_CSR1, probs = seq(0, 1, 0.1))["90%"]]
top_HK14_22G_sRNA_CSR1 <- HK14_22G_sRNA_CSR1[HK14_22G_sRNA_CSR1 >= quantile(HK14_22G_sRNA_CSR1, probs = seq(0, 1, 0.1))["90%"]]
top_HK20_22G_sRNA_CSR1 <- HK20_22G_sRNA_CSR1[HK20_22G_sRNA_CSR1 >= quantile(HK20_22G_sRNA_CSR1, probs = seq(0, 1, 0.1))["90%"]]
top_HK30_22G_sRNA_CSR1 <- HK30_22G_sRNA_CSR1[HK30_22G_sRNA_CSR1 >= quantile(HK30_22G_sRNA_CSR1, probs = seq(0, 1, 0.1))["90%"]]

top_AF14_22G_mRNA_CSR1 <- AF14_22G_mRNA_CSR1[names(top_AF14_22G_sRNA_CSR1)]
top_AF20_22G_mRNA_CSR1 <- AF20_22G_mRNA_CSR1[names(top_AF20_22G_sRNA_CSR1)]
top_AF30_22G_mRNA_CSR1 <- AF30_22G_mRNA_CSR1[names(top_AF30_22G_sRNA_CSR1)]
top_HK14_22G_mRNA_CSR1 <- HK14_22G_mRNA_CSR1[names(top_HK14_22G_sRNA_CSR1)]
top_HK20_22G_mRNA_CSR1 <- HK20_22G_mRNA_CSR1[names(top_HK20_22G_sRNA_CSR1)]
top_HK30_22G_mRNA_CSR1 <- HK30_22G_mRNA_CSR1[names(top_HK30_22G_sRNA_CSR1)]


# Combine the average 22G-RNA and mRNA expression values for all CSR-1 targets into 1 dataframe
combined_22G_mRNA_expression_CSR1 <- cbind(AF14_22G_sRNA_CSR1, AF20_22G_sRNA_CSR1, AF30_22G_sRNA_CSR1, HK14_22G_sRNA_CSR1, 
                                           HK20_22G_sRNA_CSR1, HK30_22G_sRNA_CSR1, AF14_22G_mRNA_CSR1, AF20_22G_mRNA_CSR1, 
                                           AF30_22G_mRNA_CSR1, HK14_22G_mRNA_CSR1, HK20_22G_mRNA_CSR1, HK30_22G_mRNA_CSR1)

# Combine the average 22G-RNA and mRNA expression values for just the CSR-1 targets in the top 10% of 22G expression into 1 dataframe
combined_top_22G_mRNA_expression_CSR1 <- cbind(top_AF14_22G_sRNA_CSR1, top_AF20_22G_sRNA_CSR1, 
                                               top_AF30_22G_sRNA_CSR1, top_HK14_22G_sRNA_CSR1, 
                                               top_HK20_22G_sRNA_CSR1, top_HK30_22G_sRNA_CSR1, 
                                               top_AF14_22G_mRNA_CSR1, top_AF20_22G_mRNA_CSR1, 
                                               top_AF30_22G_mRNA_CSR1, top_HK14_22G_mRNA_CSR1, 
                                               top_HK20_22G_mRNA_CSR1, top_HK30_22G_mRNA_CSR1)

# Get just the correlations between 22G and mRNA expression within the same genotype-temperature replicate, for just the genes in the top 10% 
# of 22G expression
top_22G_mRNA_correlations_CSR1 <- diag(cor(combined_top_22G_mRNA_expression_CSR1, 
                                           method = "spearman")[1:6, 7:12])
names(top_22G_mRNA_correlations_CSR1) <- c("AF14", "AF20", "AF30", "HK14", "HK20", "HK30")


# Get just the correlations between 22G and mRNA expression within the same genotype-temperature replicate, for all CSR-1 targets
all_22G_mRNA_correlations_CSR1 <- diag(cor(combined_22G_mRNA_expression_CSR1, 
                                           method = "spearman")[1:6, 7:12])
names(all_22G_mRNA_correlations_CSR1) <- c("AF14", "AF20", "AF30", "HK14", "HK20", "HK30")


# Get just the 22G-RNA expression and mRNA expression for the WAGO-1 targets with
# data available
AF14_22G_sRNA_WAGO1 <- AF14_22G_sRNA[names(AF14_22G_sRNA) %in% briggsae_orthologs_elegans_WAGO1_targets]
AF20_22G_sRNA_WAGO1 <- AF20_22G_sRNA[names(AF20_22G_sRNA) %in% briggsae_orthologs_elegans_WAGO1_targets]
AF30_22G_sRNA_WAGO1 <- AF30_22G_sRNA[names(AF30_22G_sRNA) %in% briggsae_orthologs_elegans_WAGO1_targets]
HK14_22G_sRNA_WAGO1 <- HK14_22G_sRNA[names(HK14_22G_sRNA) %in% briggsae_orthologs_elegans_WAGO1_targets]
HK20_22G_sRNA_WAGO1 <- HK20_22G_sRNA[names(HK20_22G_sRNA) %in% briggsae_orthologs_elegans_WAGO1_targets]
HK30_22G_sRNA_WAGO1 <- HK30_22G_sRNA[names(HK30_22G_sRNA) %in% briggsae_orthologs_elegans_WAGO1_targets]

AF14_22G_mRNA_WAGO1 <- AF14_22G_mRNA[names(AF14_22G_mRNA) %in% briggsae_orthologs_elegans_WAGO1_targets]
AF20_22G_mRNA_WAGO1 <- AF20_22G_mRNA[names(AF20_22G_mRNA) %in% briggsae_orthologs_elegans_WAGO1_targets]
AF30_22G_mRNA_WAGO1 <- AF30_22G_mRNA[names(AF30_22G_mRNA) %in% briggsae_orthologs_elegans_WAGO1_targets]
HK14_22G_mRNA_WAGO1 <- HK14_22G_mRNA[names(HK14_22G_mRNA) %in% briggsae_orthologs_elegans_WAGO1_targets]
HK20_22G_mRNA_WAGO1 <- HK20_22G_mRNA[names(HK20_22G_mRNA) %in% briggsae_orthologs_elegans_WAGO1_targets]
HK30_22G_mRNA_WAGO1 <- HK30_22G_mRNA[names(HK30_22G_mRNA) %in% briggsae_orthologs_elegans_WAGO1_targets]


# For each genotype-temperature combination, get the 22G-RNA expression and mRNA expression for just the WAGO-1 targets 
# in the top 10% of 22G-RNA expression among WAGO-1 targets in that genotype-temperature combination
top_AF14_22G_sRNA_WAGO1 <- AF14_22G_sRNA_WAGO1[AF14_22G_sRNA_WAGO1 >= quantile(AF14_22G_sRNA_WAGO1, probs = seq(0, 1, 0.1))["90%"]]
top_AF20_22G_sRNA_WAGO1 <- AF20_22G_sRNA_WAGO1[AF20_22G_sRNA_WAGO1 >= quantile(AF20_22G_sRNA_WAGO1, probs = seq(0, 1, 0.1))["90%"]]
top_AF30_22G_sRNA_WAGO1 <- AF30_22G_sRNA_WAGO1[AF30_22G_sRNA_WAGO1 >= quantile(AF30_22G_sRNA_WAGO1, probs = seq(0, 1, 0.1))["90%"]]
top_HK14_22G_sRNA_WAGO1 <- HK14_22G_sRNA_WAGO1[HK14_22G_sRNA_WAGO1 >= quantile(HK14_22G_sRNA_WAGO1, probs = seq(0, 1, 0.1))["90%"]]
top_HK20_22G_sRNA_WAGO1 <- HK20_22G_sRNA_WAGO1[HK20_22G_sRNA_WAGO1 >= quantile(HK20_22G_sRNA_WAGO1, probs = seq(0, 1, 0.1))["90%"]]
top_HK30_22G_sRNA_WAGO1 <- HK30_22G_sRNA_WAGO1[HK30_22G_sRNA_WAGO1 >= quantile(HK30_22G_sRNA_WAGO1, probs = seq(0, 1, 0.1))["90%"]]

top_AF14_22G_mRNA_WAGO1 <- AF14_22G_mRNA_WAGO1[names(top_AF14_22G_sRNA_WAGO1)]
top_AF20_22G_mRNA_WAGO1 <- AF20_22G_mRNA_WAGO1[names(top_AF20_22G_sRNA_WAGO1)]
top_AF30_22G_mRNA_WAGO1 <- AF30_22G_mRNA_WAGO1[names(top_AF30_22G_sRNA_WAGO1)]
top_HK14_22G_mRNA_WAGO1 <- HK14_22G_mRNA_WAGO1[names(top_HK14_22G_sRNA_WAGO1)]
top_HK20_22G_mRNA_WAGO1 <- HK20_22G_mRNA_WAGO1[names(top_HK20_22G_sRNA_WAGO1)]
top_HK30_22G_mRNA_WAGO1 <- HK30_22G_mRNA_WAGO1[names(top_HK30_22G_sRNA_WAGO1)]


# Combine the average 22G-RNA and mRNA expression values for all WAGO-1 targets into 1 dataframe
combined_22G_mRNA_expression_WAGO1 <- cbind(AF14_22G_sRNA_WAGO1, AF20_22G_sRNA_WAGO1, AF30_22G_sRNA_WAGO1, HK14_22G_sRNA_WAGO1, 
                                            HK20_22G_sRNA_WAGO1, HK30_22G_sRNA_WAGO1, AF14_22G_mRNA_WAGO1, AF20_22G_mRNA_WAGO1, 
                                            AF30_22G_mRNA_WAGO1, HK14_22G_mRNA_WAGO1, HK20_22G_mRNA_WAGO1, HK30_22G_mRNA_WAGO1)

# Combine the average 22G and mRNA expression values for just the WAGO-1 targets in the top 10% of 22G-RNA expression into 1 dataframe
combined_top_22G_mRNA_expression_WAGO1 <- cbind(top_AF14_22G_sRNA_WAGO1, top_AF20_22G_sRNA_WAGO1, 
                                                top_AF30_22G_sRNA_WAGO1, top_HK14_22G_sRNA_WAGO1, 
                                                top_HK20_22G_sRNA_WAGO1, top_HK30_22G_sRNA_WAGO1, 
                                                top_AF14_22G_mRNA_WAGO1, top_AF20_22G_mRNA_WAGO1, 
                                                top_AF30_22G_mRNA_WAGO1, top_HK14_22G_mRNA_WAGO1, 
                                                top_HK20_22G_mRNA_WAGO1, top_HK30_22G_mRNA_WAGO1)

# Get just the correlations between 22G-RNA and mRNA expression within the same genotype-temperature replicate, for just the genes in the top 10% 
# of 22G expression
top_22G_mRNA_correlations_WAGO1 <- diag(cor(combined_top_22G_mRNA_expression_WAGO1, 
                                            method = "spearman")[1:6, 7:12])
names(top_22G_mRNA_correlations_WAGO1) <- c("AF14", "AF20", "AF30", "HK14", "HK20", "HK30")


# Get just the correlations between 22G and mRNA expression within the same genotype-temperature replicate, for all WAGO-1 targets
all_22G_mRNA_correlations_WAGO1 <- diag(cor(combined_22G_mRNA_expression_WAGO1, 
                                            method = "spearman")[1:6, 7:12])
names(all_22G_mRNA_correlations_WAGO1) <- c("AF14", "AF20", "AF30", "HK14", "HK20", "HK30")


# Calculate the significance of the correlation between 22G-RNA and mRNA expression
# in each genotype-temperature combination for CSR-1 targets and WAGO-1 targets, 
# for both the set of all genes and only the genes in the top 10% of 22G-RNA expression
cor.test(AF14_22G_mRNA_CSR1, AF14_22G_sRNA_CSR1, method = "spearman")
cor.test(AF20_22G_mRNA_CSR1, AF20_22G_sRNA_CSR1, method = "spearman")
cor.test(AF30_22G_mRNA_CSR1, AF30_22G_sRNA_CSR1, method = "spearman")
cor.test(HK14_22G_mRNA_CSR1, HK14_22G_sRNA_CSR1, method = "spearman")
cor.test(HK20_22G_mRNA_CSR1, HK20_22G_sRNA_CSR1, method = "spearman")
cor.test(HK30_22G_mRNA_CSR1, HK30_22G_sRNA_CSR1, method = "spearman")

cor.test(top_AF14_22G_mRNA_CSR1, top_AF14_22G_sRNA_CSR1, method = "spearman")
cor.test(top_AF20_22G_mRNA_CSR1, top_AF20_22G_sRNA_CSR1, method = "spearman")
cor.test(top_AF30_22G_mRNA_CSR1, top_AF30_22G_sRNA_CSR1, method = "spearman")
cor.test(top_HK14_22G_mRNA_CSR1, top_HK14_22G_sRNA_CSR1, method = "spearman")
cor.test(top_HK20_22G_mRNA_CSR1, top_HK20_22G_sRNA_CSR1, method = "spearman")
cor.test(top_HK30_22G_mRNA_CSR1, top_HK30_22G_sRNA_CSR1, method = "spearman")

cor.test(AF20_22G_mRNA_WAGO1, AF20_22G_sRNA_WAGO1, method = "spearman")
cor.test(AF30_22G_mRNA_WAGO1, AF30_22G_sRNA_WAGO1, method = "spearman")
cor.test(HK14_22G_mRNA_WAGO1, HK14_22G_sRNA_WAGO1, method = "spearman")
cor.test(HK20_22G_mRNA_WAGO1, HK20_22G_sRNA_WAGO1, method = "spearman")
cor.test(HK30_22G_mRNA_WAGO1, HK30_22G_sRNA_WAGO1, method = "spearman")

cor.test(top_AF14_22G_mRNA_WAGO1, top_AF14_22G_sRNA_WAGO1, method = "spearman")
cor.test(top_AF20_22G_mRNA_WAGO1, top_AF20_22G_sRNA_WAGO1, method = "spearman")
cor.test(top_AF30_22G_mRNA_WAGO1, top_AF30_22G_sRNA_WAGO1, method = "spearman")
cor.test(top_HK14_22G_mRNA_WAGO1, top_HK14_22G_sRNA_WAGO1, method = "spearman")
cor.test(top_HK20_22G_mRNA_WAGO1, top_HK20_22G_sRNA_WAGO1, method = "spearman")
cor.test(top_HK30_22G_mRNA_WAGO1, top_HK30_22G_sRNA_WAGO1, method = "spearman")
