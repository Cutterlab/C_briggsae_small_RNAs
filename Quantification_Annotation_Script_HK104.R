# This script is modified from the script "c_elegans_reference_prep.R", found at 
# https://github.com/ClaycombLab/Charlesworth_2020

# Attach packages needed for this pipeline
library(rtracklayer)
library(GenomicAlignments)
library(data.table)
library(BiocParallel)
library(bitops)

# This script uses gene annotations for genes 
# (lifted_c_briggsae_HK104_NGM_Alignments_canonical_geneset.gtf), repeats 
# (c_briggsae_HK104_NGM_Alignments_repeat_annotations.gff3), and transposons
# (c_briggsae_HK104_NGM_Alignments_transposable_element_annotations.gff3), all
# of which are available at 
# https://github.com/Cutterlab/C_briggsae_small_RNAs/blob/main/HK104%20Pseudo-reference%20Genome.zip


base.chromosomes <- c("I", "II", "III", "IV", "V", "X")

# Prepare combined reference at gene level
# read in genes file
genes <- "lifted_c_briggsae_HK104_NGM_Alignments_canonical_geneset.gtf"
genes <- import(genes)
genes <- genes[seqnames(genes) %in% base.chromosomes]
genes <- genes[genes$type == "gene" & grepl("WBGene", genes$gene_id)]
mcols(genes) <- mcols(genes)[, c("gene_id", "gene_biotype")]

# readin transposable elements from transposons
transposons <- "c_briggsae_HK104_NGM_Alignments_transposable_element_annotations.gff3"
transposons <- import(transposons)
transposons <- transposons[seqnames(transposons) %in% base.chromosomes]
transposons <- transposons[transposons$type == "transposable_element"]
transposons$gene_id <- transposons$Family
transposons$gene_id[is.na(transposons$gene_id) | transposons$gene_id == "UNKNOWN"] <- 
  transposons$Name[is.na(transposons$gene_id) | transposons$gene_id == "UNKNOWN"]
transposons$gene_biotype <- transposons$type
mcols(transposons) <- mcols(transposons)[, c("gene_id", "gene_biotype")]

# readin repeat elements from transposons
repeats <- "c_briggsae_HK104_NGM_Alignments_repeat_annotations.gff3"
repeats <- import(repeats)
repeats <- repeats[seqnames(repeats) %in% base.chromosomes]
repeats <- repeats[grepl("_CB", repeats$Target) & !grepl("tRNA", repeats$Target) & repeats$type == "repeat_region"]
repeats$gene_id <- gsub("^([^ ]+).*", "\\1", repeats$Target)
repeats$gene_biotype <- repeats$type
mcols(repeats) <- mcols(repeats)[, c("gene_id", "gene_biotype")]

gene.whole <- c(genes, transposons, repeats)
gene.whole$parts.ref <- data.table(gene_biotype=gene.whole$gene_biotype)[, count := seq(.N), by=gene_biotype]$count

save(gene.whole, file="WS272.HK104_gene.whole.transposons.repeats.RData")

# read in canonical geneset and collapse transcripts at gene level into
# disjoint set with type of each element set to the union of overlapping
# features
gtf_file <- "lifted_c_briggsae_HK104_NGM_Alignments_canonical_geneset.gtf"
gtf <- import(gtf_file)
gtf <- gtf[seqnames(gtf) %in% base.chromosomes]
gtf <- c(gtf)


gtf.genes <- bplapply(unique(grep("(WBGene)", gtf$gene_id, value=TRUE)), function(gene, gtf){

  gene.parts <- gtf[
    gtf$type %in% c("exon", "three_prime_utr", "five_prime_utr") & 
      (grepl("WBGene", gtf$gene_id)) &
      gtf$gene_id == gene
    ]
  # convert string to binary encoding (introns=8, full_feature=16)
  gene.parts$type <- as.integer(c("five_prime_utr"=1, "exon"=2, "three_prime_utr"=4)[as.character(gene.parts$type)])
  gene.biotype <- gtf[gtf$type == "gene" & gtf$gene_id == gene]$gene_biotype
  parts.by.transcripts <- split(gene.parts, gene.parts$transcript_id)
  parts.by.transcripts.disjoint <- lapply(parts.by.transcripts, disjoin, with.revmap=TRUE)
  
  for (name in names(parts.by.transcripts.disjoint)) {
    parts.by.transcripts.disjoint[[name]]$type <- sapply(parts.by.transcripts.disjoint[[name]]$revmap, function(x){
      # NaN because UTRs should never overlap
      c(1, 2, 1, 4, NaN, 4, NaN)[Reduce(bitops::bitOr, parts.by.transcripts[[name]][x]$type)]
    })
    parts.by.transcripts.disjoint[[name]]$revmap <- NULL
    introns <- setdiff(range(parts.by.transcripts.disjoint[[name]]), parts.by.transcripts.disjoint[[name]])
    introns$type <- rep(as.integer(8), length(introns))
    parts.by.transcripts.disjoint[[name]] <- c(
      parts.by.transcripts.disjoint[[name]],
      introns
    )
  }
  
  parts.by.transcripts.disjoint <- unlist(GRangesList(parts.by.transcripts.disjoint))
  parts.disjoint <- disjoin(parts.by.transcripts.disjoint, with.revmap=TRUE)
  parts.disjoint$type <- sapply(parts.disjoint$revmap, function(x){
    Reduce(bitops::bitOr, parts.by.transcripts.disjoint[x]$type)
  })
  parts.disjoint$revmap <- NULL
  parts.disjoint$gene_id <- gene
  parts.disjoint$gene_biotype <- gene.biotype
  parts.disjoint
}, gtf)
gtf.genes <- unlist(GRangesList(gtf.genes))
gtf.genes$whole.ref <- structure(1:(length(genes)), names=c(genes$gene_id))[gtf.genes$gene_id]
save(gtf.genes, file="WS272.HK104_gene.parts.RData")

# read in transposable elements, use entire transposable element as desired feature
transposons <-
  "c_briggsae_HK104_NGM_Alignments_transposable_element_annotations.gff3"
transposons <- import(transposons)
transposons <- transposons[seqnames(transposons) %in% base.chromosomes]
transposons <- transposons[transposons$type == "transposable_element"]
transposons$gene_id <- transposons$Family
transposons$gene_id[is.na(transposons$gene_id) | transposons$gene_id == "UNKNOWN"] <- 
  transposons$Name[is.na(transposons$gene_id) | transposons$gene_id == "UNKNOWN"]
transposons$gene_biotype <- transposons$type
mcols(transposons) <- mcols(transposons)[, c("type", "gene_id", "gene_biotype")]
transposons$type <- as.integer(16)  # represents full feature in current encoding scheme
transposons$whole.ref <- length(genes) + 1:length(transposons)

# readin repeat elements from gff3
repeats <- "c_briggsae_HK104_NGM_Alignments_repeat_annotations.gff3"
repeats <- import(repeats)
repeats <- repeats[seqnames(repeats) %in% base.chromosomes]
repeats <- repeats[grepl("_CB", repeats$Target) & !grepl("tRNA", repeats$Target) & repeats$type == "repeat_region"]
repeats$gene_id <- gsub("^([^ ]+).*", "\\1", repeats$Target)
repeats$gene_biotype <- repeats$type
mcols(repeats) <- mcols(repeats)[, c("type", "gene_id", "gene_biotype")]
repeats$type <- as.integer(16)  # represents full gene in current encoding scheme
repeats$whole.ref <- length(genes)  + length(transposons) + 1:length(repeats)

# combine gene parts and transposon annotations and save
gene.parts <- c(gtf.genes, transposons, repeats)
stopifnot(sum(gene.parts$gene_id != gene.whole[gene.parts$whole.ref]$gene_id) == 0)
save(gene.parts, file="WS272.HK104_gene.parts.transposons.repeats.RData")
