# Load C. briggsae cSR-1, and C. elegans CSR-1 and WAGO-1 targets from nar-02785-a-2014-File012.xlsx, which
# is Supplementary Table S2 from Tu et al. 2015 (https://academic.oup.com/nar/article/43/1/208/1008699)

library(readxl)
nar_02785_a_2014_File012 <- read_excel("nar-02785-a-2014-File012.xlsx", 
                                       skip = 6) 
briggsae_CSR1_targets <- 
  unique(nar_02785_a_2014_File012$`C. briggsae genes`[nar_02785_a_2014_File012$`C. briggsae CSR-1 targets` == "TRUE"])

# Get the C. briggsae orthologs of C. elegans CSR-1 targets (to determine the overlap
# with C. briggsae CSR-1 targets)
briggsae_orthologs_elegans_CSR1_targets <- 
  unique(nar_02785_a_2014_File012$`C. briggsae genes`[nar_02785_a_2014_File012$`C. elegans CSR-1 targets` == "TRUE"])

# Get the C. briggsae orthologs of C. elegans WAGO-1 targets
briggsae_orthologs_elegans_WAGO1_targets <- 
  unique(nar_02785_a_2014_File012$`C. briggsae genes`[nar_02785_a_2014_File012$`C. elegans WAGO-1 targets` == "TRUE" & 
                                                        !(nar_02785_a_2014_File012$`C. briggsae genes` == "NA")])

# Load C. elegans EGO-1 targets from 1-s2.0-S009286740901174X-mmc5 (1).xls, which is
# Supplementary Table S5 from Claycomb et al. 2009 (https://www.cell.com/cell/fulltext/S0092-8674(09)01174-X)
X1_s2_0_S009286740901174X_mmc5_1_ <- read_excel("1-s2.0-S009286740901174X-mmc5 (1).xls", 
                                                skip = 6)

# To convert cosmid IDs to WormBase IDs, load in C. elegans gene IDs from
# WormBase release WS277, which can be found at 
# https://downloads.wormbase.org/releases/WS277/species/c_elegans/PRJNA13758/annotation/c_elegans.PRJNA13758.WS277.geneIDs.txt.gz
c_elegans.PRJNA13758.WS277.geneIDs <- read.csv("c_elegans.PRJNA13758.WS277.geneIDs.txt", header=FALSE)

# Get WormBase IDs of the EGO-1 targets in C. elegans, then get the C. briggsae orthologs of those
elegans_ego1_22G_targets <- 
  as.character(c_elegans.PRJNA13758.WS277.geneIDs$V2[c_elegans.PRJNA13758.WS277.geneIDs$V4 %in% 
                                                       X1_s2_0_S009286740901174X_mmc5_1_$`cosmid ID`])
briggsae_orthologs_elegans_ego1_22G_targets <- 
  unique(nar_02785_a_2014_File012$`C. briggsae genes`[nar_02785_a_2014_File012$`C. elegans genes` 
                                                      %in% elegans_ego1_22G_targets & 
                                                        nar_02785_a_2014_File012$`C. briggsae genes` != "NA"])


# Load C. elegans RRF-1 targets from elegans_rrf1_22G_targets.txt, which contains
# the RRF-1 targets from Supplementary Table S4 from Vasale et al. 2010
# (https://www.pnas.org/doi/10.1073/pnas.0911908107)
elegans_rrf1_22G_targets <- read.table("elegans_rrf1_22G_targets.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)

elegans_rrf1_22G_targets <- as.character(c_elegans.PRJNA13758.WS277.geneIDs$V2[c_elegans.PRJNA13758.WS277.geneIDs$V4 %in% elegans_rrf1_22G_targets$V1])

briggsae_orthologs_elegans_rrf1_22G_targets <- 
  unique(nar_02785_a_2014_File012$`C. briggsae genes`[nar_02785_a_2014_File012$`C. elegans genes` %in% 
                                                        elegans_rrf1_22G_targets & nar_02785_a_2014_File012$`C. briggsae genes` != "NA"])


# Load C. elegans HRDE-1 targets from nature11352-s2.xlsx, which is Supplementary 
# Table 2 from Buckley et al. 2012 (https://www.nature.com/articles/nature11352#Sec7). 
# The HRDE-1 targets are the members of "GeneList_FINN1top500.txt" and "Sheet2"
# with at least 50 RPKM of expression. Then, get the C. briggsae orthologs of these
# C. elegans targets. 
hrde1_protein_coding_frame <- read_excel("nature11352-s2.xlsx", sheet = "GeneList_FINN1top500.txt")
hrde1_pseudogene_frame <- read_excel("nature11352-s2.xlsx", sheet = "Sheet2")

elegans_HRDE1_targets <- union(hrde1_protein_coding_frame$geneName[hrde1_protein_coding_frame$`FINN1coIPsiRNA(rpkm)` >= 50], 
                               hrde1_pseudogene_frame$pseudogene[hrde1_pseudogene_frame$rpkm >= 50])

elegans_HRDE1_targets <- as.character(union(c_elegans.PRJNA13758.WS277.geneIDs$V2[c_elegans.PRJNA13758.WS277.geneIDs$V4 %in% 
                                                                                    elegans_HRDE1_targets], 
                                            c_elegans.PRJNA13758.WS277.geneIDs$V2[c_elegans.PRJNA13758.WS277.geneIDs$V3 %in% elegans_HRDE1_targets]))

briggsae_orthologs_elegans_HRDE1_targets <- 
  unique(nar_02785_a_2014_File012$`C. briggsae genes`[nar_02785_a_2014_File012$`C. elegans genes` %in% 
                                                        elegans_HRDE1_targets & nar_02785_a_2014_File012$`C. briggsae genes` != "NA"])


# Load in C. elegans targets of CSR-1a and CSR-1b from media-3 (1).xlsx, which are in 
# "IP Enrichment (22G reads)" from Supplementary Table S2 from Charlesworth et al. 2021
# (https://academic.oup.com/nar/article/49/15/8836/6331683). Then, get the C. briggsae
# orthologs of these C. elegans targets
new_elegans_csr1_isoform_targets <- read_excel("media-3 (1).xlsx", 
                                               sheet = "IP Enrichment (22G reads)", skip = 1)

new_elegans_csr1a_targets <- unique(new_elegans_csr1_isoform_targets$Gene[new_elegans_csr1_isoform_targets$Argonaute == "CSR-1a" & new_elegans_csr1_isoform_targets$`Enriched?` == TRUE & new_elegans_csr1_isoform_targets$Biotype == "protein_coding_AS"])
briggsae_orthologs_elegans_new_csr1a_targets <- unique(nar_02785_a_2014_File012$`C. briggsae genes`[nar_02785_a_2014_File012$`C. elegans genes` %in% new_elegans_csr1a_targets])
briggsae_orthologs_elegans_new_csr1a_targets <- briggsae_orthologs_elegans_new_csr1a_targets[briggsae_orthologs_elegans_new_csr1a_targets != "NA"]

new_elegans_csr1b_targets <- unique(new_elegans_csr1_isoform_targets$Gene[new_elegans_csr1_isoform_targets$Argonaute == "CSR-1b" & new_elegans_csr1_isoform_targets$`Enriched?` == TRUE & new_elegans_csr1_isoform_targets$Biotype == "protein_coding_AS"])
new_elegans_csr1b_targets <- new_elegans_csr1b_targets[new_elegans_csr1b_targets != "GFP"]
briggsae_orthologs_elegans_new_csr1b_targets <- unique(nar_02785_a_2014_File012$`C. briggsae genes`[nar_02785_a_2014_File012$`C. elegans genes` %in% new_elegans_csr1b_targets])
briggsae_orthologs_elegans_new_csr1b_targets <- briggsae_orthologs_elegans_new_csr1b_targets[briggsae_orthologs_elegans_new_csr1b_targets != "NA"]
