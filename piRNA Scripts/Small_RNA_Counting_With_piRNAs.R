# This script is modified from the script "sRNA_quant.R", found at 
# https://github.com/ClaycombLab/Charlesworth_2020

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
require(reshape2)
require(data.table)

# NOTE: This script was run using AF16 reads aligned to the AF16 reference genome, 
# and using HK104 reads aligned to the HK104 pseudo-reference genome

# The code between here and the save() statements at the end of the script needs to be run 
# once for each replicate (so 18 times in total). The variable "name" defines which
# replicate is used for this particular run of this part of the script (so
# running the code below will only perform the counting for the AF14-1 replicate). 
# The BAM files used here (1 for each replicate) were the output from running STAR.

name <- "HK30-3"
bam_file <- paste("~/../Desktop/briggsae_small_RNA_alignments/", name, "Aligned.sortedByCoord.out.bam",
                  sep = "")

# Load in the BAM file and filter out reads that did not align perfectly to
# to genome
param <- ScanBamParam(what = c("seq", "flag", "qname"), reverseComplement=TRUE)
bam <- readGAlignments(bam_file, param=param)
total.count <- length(unique(mcols(bam)$qname))
perfect <- grepl("(^\\d+M$)|(^\\d+M\\d+N\\d+M$)", cigar(bam))
bam.unperfect <- bam[!perfect]
bam <- bam[perfect]
bam.unperfect <- bam.unperfect[!(mcols(bam.unperfect)$qname %in% mcols(bam)$qname)]
perfect.count <- length(unique(mcols(bam)$qname))
unperfect.count <- length(unique(mcols(bam.unperfect)$qname))


# Load in annotation data for read counting - these are either AF16 annotations
# or HK104 annotations, depending on which sample we're trying to count.
# These files were created in Quantification_Annotation_Script_AF16_With_piRNAs.R and
# Quantification_Annotation_Script_HK104_With_piRNAs.R
if (grepl("AF", name)) {
  load("~/WS272.gene.whole.transposons.repeats_with_piRNAs.rdata")
  load("~/WS272.gene.parts.transposons.repeats_with_piRNAs.rdata")
} else if(grepl("HK", name)) {
  load("~/WS272.HK104_gene.whole.transposons.repeats_with_piRNAs.RData")
  load("~/WS272.HK104_gene.parts.transposons.repeats_with_piRNAs.RData")
}


# Define functions and vectors to be used for read counting and annotation
count_reads <- function(ov, bam) {
  # fix names and merge to remove ov.parts edge cases
  colnames(ov)[2] <- "gene"
  colnames(ov.parts)[2] <- "feature"
  ov.merged <- merge(ov, ov.parts, by=c("gene", "queryHits", "reads"))
  # determine which reads map to more than one feature
  ov.merged[, shared := .N, by=reads]
  # determine distribution of unique reads
  ov.merged[, prior := sum(shared == 1), by=list(gene, feature)]
  ov.merged[, sum.prior := sum(prior), by=reads]
  # partition reads between genes
  ov.merged[, fraction := prior / sum.prior ]
  ov.merged[is.nan(fraction), fraction := 1 / shared]  # can this be optimized
  # add sequence properties of reads
  ov.merged[, seq := as.character(mcols(bam[queryHits])$seq)]
  ov.merged[, c("length", "first") := list(nchar(seq), substr(seq, 1, 1)),]
  # count reads
  count <- ov.merged[,
                     .(count=sum(fraction)),
                     by=.(gene, feature, first, length)
                     ]
  return(count)
}

prepare_biotype_reference <- function(reference, biotype, sense){
  reference.biotype <- reference[reference$gene_biotype == biotype]

  if (biotype %in% biotypes.antisense) {
    if (sense == "expected") {
      reference.biotype <- invertStrand(reference.biotype)
      reference.biotype$gene_biotype <- paste(biotype, "_AS", sep="")
    }
  } else {
    if (sense == "unexpected") {
      reference.biotype <- invertStrand(reference.biotype)
      reference.biotype$gene_biotype <- paste(biotype, "_AS", sep="")
    }
  }
  # update whole.ref for parts references
  return(reference.biotype)
}


biotypes <- c(
  "transposable_element", "repeat_region", "rRNA", "snoRNA", "snRNA", "tRNA",  
  "miRNA", "ncRNA", "protein_coding", "pseudogene", "lincRNA", "piRNA"
)

biotypes.antisense <-
  c("transposable_element", "repeat_region", "ncRNA", "protein_coding", "pseudogene")



# Count how many reads map onto each type of feature, for both sense and antisense
results_with_piRNAs[[name]] <- list()
seen <- c()
counted.count <- 0
print(paste("biotype", "#_features", "counts"))
for(sense in c("expected", "unexpected")) {
  for (biotype in biotypes) {
    # prepare biotype specifc gene.whole and find overlaps


    gene.whole.biotype <- prepare_biotype_reference(gene.whole, biotype, sense)
    gene.parts.biotype <- prepare_biotype_reference(gene.parts, biotype, sense)


    biotype <- as.character(gene.whole.biotype[1]$gene_biotype)
    ov <- findOverlaps(bam, gene.whole.biotype, type = "within")
    ov.parts <- findOverlaps(bam, gene.parts.biotype)
    ov.parts <- as.data.table(ov.parts)
    ov.parts[, reads := mcols(bam[ov.parts[, queryHits]])$qname]


    ov.parts$gene <- gene.whole[gene.parts.biotype[ov.parts$subjectHits]$whole.ref]$parts.ref
    ov.parts$subjectHits <- gene.parts.biotype[ov.parts$subjectHits]$type

    # little sanity check
    stopifnot(sum(gene.parts.biotype$gene_id != gene.whole.biotype[gene.whole[gene.parts.biotype$whole.ref]$parts.ref]$gene_id) == 0)


    # map alignments (queryHits) to reads
    ov <- as.data.table(ov)
    ov[, reads := mcols(bam[ov[, queryHits]])$qname]

    # remove previously seen alignments and update seen set
    ov <- ov[!(reads %in% seen), ]
    toDistribute <- length(unique(ov[, reads]))
    potentials <- unique(ov[, subjectHits])
    seen <- union(seen, unique(ov[, reads]))  # update seen

    # count reads
    counted <- 0
    if (nrow(ov)) {
      counts <- count_reads(ov, bam)
      counts[, gene := gene.whole.biotype[gene]$gene_id]
      results_with_piRNAs[[name]][[biotype]] <- counts
      counted <- counts[, sum(count)]
    }
    counted.count <- counted.count + counted

    print(paste(biotype, length(gene.whole.biotype), counted))
  }
}

# Get count information for reads that didn't map onto any feature
results_with_piRNAs[[name]][["no_feature"]] <- data.table(
  gene = "no_feature",
  length = nchar(mcols(bam[!(mcols(bam)$qname %in% seen)])$seq),
  first = substr(mcols(bam[!(mcols(bam)$qname %in% seen)])$seq, 1, 1),
  count = 1 / as.numeric(table(mcols(bam[!(mcols(bam)$qname %in% seen)])$qname)[mcols(bam[!(mcols(bam)$qname %in% seen)])$qname])
)[,
  .(count=sum(count)),
  by=.(gene, first, length)
  ]


# Count the number of reads that were loaded in, that aligned perfectly to the
# reference genome, and that mapped onto an annotated feature
results.sizes_with_piRNAs[[name]] <- c(
  mapped=total.count,
  perfect=perfect.count,
  counted=counted.count
)

# Filter the BAM file to keep only the primary alignments, and update the
# resulting dataframe to also contain the first base of each RNA
primary_bam_frame <- as.data.frame(bam[!bitAnd(0x100, mcols(bam)$flag)])
get_first_base <- function(RNA){return(unlist(strsplit(RNA, ""))[1])}

primary_bam_frame$first <- unlist(lapply(primary_bam_frame$seq, get_first_base))
primary_bam_frame <- primary_bam_frame[, c(1, 5, 7, 12)]

bam_frames[[name]] <- primary_bam_frame



# After the above code is run for all 18 replicates, save the results
save(results_with_piRNAs, file = "results_with_piRNAs.RData")
save(results.sizes_with_piRNAs, file = "resultsSizes_with_piRNAs.RData")
save(bam_frames, file = "bam_frames.RData")


