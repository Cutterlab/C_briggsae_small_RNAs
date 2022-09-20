# This script is modified from the script "c_elegans_reference_prep.R", found at 
# https://github.com/ClaycombLab/Charlesworth_2020

# Attach packages needed for this pipeline
library(rtracklayer)
library(GenomicAlignments)
library(data.table)
library(BiocParallel)
library(bitops)

# This script uses gene annotations for genes 
# (c_briggsae.PRJNA10731.WS272.canonical_geneset_with_piRNAs.gtf), repeats 
# (c_briggsae.PRJNA10731.WS272.annotations.repeat.gff3), and transposons
# (c_briggsae.PRJNA10731.WS272.annotations.transposable_element.gff3).  
# c_briggsae.PRJNA10731.WS272.canonical_geneset_with_piRNAs.gtf was made by 
# combining the gene annotations from c_briggsae.PRJNA10731.WS272.canonical_geneset.gtf
# (available at https://downloads.wormbase.org/releases/WS272/species/c_briggsae/PRJNA10731/)
# and the annotated piRNA loci in cbg_piRNAs_motifs_and_21Us.bed from
# Supplemental Data S1 of Beltran et al. 2019 (https://doi.org/10.1016/j.devcel.2018.12.026)
# c_briggsae.PRJNA10731.WS272.annotations.repeat.gff3 and 
# c_briggsae.PRJNA10731.WS272.annotations.transposable_element.gff3
# were created by extracting lines containing "repeat" and "transposable_element", 
# respectively, from c_briggsae.PRJNA10731.WS272.annotations.gff3, which is also 
# available at https://downloads.wormbase.org/releases/WS272/species/c_briggsae/PRJNA10731/ 



base.chromosomes <- c("I", "II", "III", "IV", "V", "X")

# Prepare combined reference at gene level
# read in genes file
genes <- "c_briggsae.PRJNA10731.WS272.canonical_geneset_with_piRNAs.gtf"
genes <- import(genes)
genes <- genes[seqnames(genes) %in% base.chromosomes]
genes <- genes[genes$type == "gene" & grepl("WBGene", genes$gene_id)]
mcols(genes) <- mcols(genes)[, c("gene_id", "gene_biotype")]

# readin transposable elements from transposons
transposons <- "c_briggsae.PRJNA10731.WS272.annotations.transposable_element.gff3"
transposons <- import(transposons)
transposons <- transposons[seqnames(transposons) %in% base.chromosomes]
transposons <- transposons[transposons$type == "transposable_element"]
transposons$gene_id <- transposons$Family
transposons$gene_id[is.na(transposons$gene_id) | transposons$gene_id == "UNKNOWN"] <- 
  transposons$Name[is.na(transposons$gene_id) | transposons$gene_id == "UNKNOWN"]
transposons$gene_biotype <- transposons$type
mcols(transposons) <- mcols(transposons)[, c("gene_id", "gene_biotype")]

# readin repeat elements from transposons
repeats <- "c_briggsae.PRJNA10731.WS272.annotations.repeat.gff3"
repeats <- import(repeats)
repeats <- repeats[seqnames(repeats) %in% base.chromosomes]
repeats <- repeats[grepl("_CB", repeats$Target) & !grepl("tRNA", repeats$Target) & repeats$type == "repeat_region"]
repeats$gene_id <- gsub("^([^ ]+).*", "\\1", repeats$Target)
repeats$gene_biotype <- repeats$type
mcols(repeats) <- mcols(repeats)[, c("gene_id", "gene_biotype")]

gene.whole <- c(genes, transposons, repeats)
gene.whole$parts.ref <- data.table(gene_biotype=gene.whole$gene_biotype)[, count := seq(.N), by=gene_biotype]$count

save(gene.whole, file="WS272.gene.whole.transposons.repeats_with_piRNAs.RData")

# read in canonical geneset and collapse transcripts at gene level into
# disjoint set with type of each element set to the union of overlapping
# features
gtf_file <- "c_briggsae.PRJNA10731.WS272.canonical_geneset_with_piRNAs.gtf"
gtf <- import(gtf_file)
gtf <- gtf[seqnames(gtf) %in% base.chromosomes]
gtf <- c(gtf)


gtf.genes <- lapply(unique(grep("(WBGene)", gtf$gene_id, value=TRUE)), function(gene, gtf, counter){

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
}, gtf, counter)
gtf.genes <- unlist(GRangesList(gtf.genes))
gtf.genes$whole.ref <- structure(1:(length(genes)), names=c(genes$gene_id))[gtf.genes$gene_id]
save(gtf.genes, file="WS272.gene.parts_with_piRNAs.RData")

# read in transposable elements, use entire transposable element as desired feature
transposons <-
  "c_briggsae.PRJNA10731.WS272.annotations.transposable_element.gff3"
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
repeats <- "c_briggsae.PRJNA10731.WS272.annotations.repeat.gff3"
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
save(gene.parts, file="WS272.gene.parts.transposons.repeats_with_piRNAs.RData")
