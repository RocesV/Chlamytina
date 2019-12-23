### Annotation files ChromHMM

library(Biostrings)
library(seqRFLP)
library(GenomicFeatures)

## 1 Chromsizes
Nucleus <- readDNAStringSet(filepath = "../Chromatin States/all_together/inputs/Assembly only chr (just in case)/Creinhardtii_281_v5.0.fa")
Organelles <- readDNAStringSet(filepath = "../Chromatin States/all_together/inputs/Assembly only chr (just in case)/C.reinhardtii_v5.3_genomic_scaffold_plastids.fasta")
class(Nucleus)
filter <- grep("scaffold", names(Nucleus))
Nucleus <- Nucleus[-filter] 
names(Nucleus) <- gsub("chromosome_", "chr" ,names(Nucleus))
names(Organelles)
Organelles <- Organelles[c(55,56)]
names(Organelles) <- c("chrC", "chrM")
Chlamytina_genome <- c(Nucleus, Organelles) 

chromsizes <- data.frame(names(Chlamytina_genome), width(Chlamytina_genome))
write.table(chromsizes, "../Chromatin States/all_together/inputs/cellfilemarkstable_chromsizes/Cre5_5.txt", quote = FALSE,
            col.names = F, row.names = F, sep = "\t")

## 2 COORDS
Nuclear <- makeTxDbFromGFF("../Chromatin States/all_together/inputs/Annotation_Enrichments/Creinhardtii_281_v5.5.gene_exons.gff3")
Chloroplast <- makeTxDbFromGFF("../Chromatin States/all_together/inputs/Annotation_Enrichments/C.reinhardtii_chloroplast.gff3")
Mitochondrion <- makeTxDbFromGFF("../Chromatin States/all_together/inputs/Annotation_Enrichments/C.reinhardtii_mitochondrion.gff3")
All.Anno <- list(Nuclear = Nuclear, Chl = Chloroplast, Mit = Mitochondrion)

Exons.bed <- list()
Introns.bed <- list()
genes.bed <- list()
Intergenic.bed <- list()
III.UTR.bed <- list()
V.UTR.bed <- list()
TSS.bed <- list()
TSS.2kb.bed <- list()
TES.bed <- list()


# Exons
Exons <- exons(All.Anno$Nuclear)
Exons.bed[["Nuclear"]] <- data.frame(chr = as.character(Exons@seqnames), start = Exons@ranges@start, end = (Exons@ranges@start + Exons@ranges@width - 1))
Exons.bed$Nuclear$chr <- gsub("chromosome_", "chr", Exons.bed$Nuclear$chr)
filter <- grep("scaffold", Exons.bed$Nuclear$chr)
Exons.bed$Nuclear <- Exons.bed$Nuclear[-filter,]
Exons.bed$Nuclear <- Exons.bed$Nuclear[!duplicated(Exons.bed$Nuclear),]  
Exons <- exons(All.Anno$Chl)  
Exons.bed[["Chl"]] <- data.frame(chr = as.character(Exons@seqnames), start = Exons@ranges@start, end = (Exons@ranges@start + Exons@ranges@width - 1))
Exons.bed$Chl$chr <- rep("chrC", nrow(Exons.bed$Chl))
Exons <- exons(All.Anno$Mit)
Exons.bed[["Mit"]] <- data.frame(chr = as.character(Exons@seqnames), start = Exons@ranges@start, end = (Exons@ranges@start + Exons@ranges@width - 1))
Exons.bed$Mit$chr <- rep("chrM", nrow(Exons.bed$Mit))

Exons.bed.f <- rbind(Exons.bed$Nuclear, Exons.bed$Chl, Exons.bed$Mit)
levels <- c(paste0("chr", c(1:17)), "chrC", "chrM")
Exons.bed.f$chr <- factor(Exons.bed.f$chr, levels, ordered = TRUE)
Exons.bed.f <- Exons.bed.f[order(Exons.bed.f$chr, Exons.bed.f$start),]
write.table(Exons.bed.f, file = "../Chromatin States/all_together/inputs/Annotation_Enrichments/COORDS/ChlamytinaExon.Cre55.bed", col.names = F, row.names = F, sep = "\t", 
            quote = F)

# Introns
Introns <- intronsByTranscript(All.Anno$Nuclear)
Introns <- unlist(Introns, use.names = F)
Introns.bed[['Nuclear']] <- data.frame(chr = as.character(Introns@seqnames), start = Introns@ranges@start, end = (Introns@ranges@start + Introns@ranges@width - 1))
Introns.bed$Nuclear$chr <- gsub("chromosome_", "chr", Introns.bed$Nuclear$chr)
filter <- grep("scaffold", Introns.bed$Nuclear$chr)
Introns.bed$Nuclear <- Introns.bed$Nuclear[-filter,]
Introns.bed$Nuclear <- Introns.bed$Nuclear[!duplicated(Introns.bed$Nuclear),]
Introns <- intronsByTranscript(All.Anno$Chl)
Introns <- unlist(Introns, use.names = F)
Introns.bed[['Chl']] <- data.frame(chr = as.character(Introns@seqnames), start = Introns@ranges@start, end = (Introns@ranges@start + Introns@ranges@width - 1))
Introns.bed$Chl$chr <- rep("chrC", nrow(Introns.bed$Chl))
 # Mit no introns

Introns.bed.f <- rbind(Introns.bed$Nuclear, Introns.bed$Chl)
Introns.bed.f$chr <- factor(Introns.bed.f$chr, levels, ordered = TRUE)
Introns.bed.f <- Introns.bed.f[order(Introns.bed.f$chr, Introns.bed.f$start),]
write.table(Introns.bed.f, file = "../Chromatin States/all_together/inputs/Annotation_Enrichments/COORDS/ChlamytinaIntron.Cre55.bed", col.names = F, row.names = F, sep = "\t", 
            quote = F)

# Genes and Intergenic
Genes <- genes(All.Anno$Nuclear)
Intergenic <- gaps(Genes)
genes.bed[['Nuclear']] <- data.frame(chr = as.character(Genes@seqnames), start = Genes@ranges@start, end = (Genes@ranges@start + Genes@ranges@width - 1))
Intergenic.bed[['Nuclear']] <- data.frame(chr = as.character(Intergenic@seqnames), start = Intergenic@ranges@start, end = (Intergenic@ranges@start + Intergenic@ranges@width - 1))
genes.bed$Nuclear$chr <- gsub("chromosome_", "chr", genes.bed$Nuclear$chr)
Intergenic.bed$Nuclear$chr <- gsub("chromosome_", "chr", Intergenic.bed$Nuclear$chr)
filter1 <- grep("scaffold", genes.bed$Nuclear$chr)
filter2 <- grep("scaffold", Intergenic.bed$Nuclear$chr)
genes.bed$Nuclear <- genes.bed$Nuclear[-filter1,]
Intergenic.bed$Nuclear <- Intergenic.bed$Nuclear[-filter2,]

Genes <- genes(All.Anno$Chl)
Intergenic <- gaps(Genes)
genes.bed[['Chl']] <- data.frame(chr = as.character(Genes@seqnames), start = Genes@ranges@start, end = (Genes@ranges@start + Genes@ranges@width - 1))
Intergenic.bed[['Chl']] <- data.frame(chr = as.character(Intergenic@seqnames), start = Intergenic@ranges@start, end = (Intergenic@ranges@start + Intergenic@ranges@width - 1))
genes.bed$Chl$chr <- rep("chrC", nrow(genes.bed$Chl))
Intergenic.bed$Chl$chr <- rep("chrC", nrow(Intergenic.bed$Chl))

Genes <- genes(All.Anno$Mit)
Intergenic <- gaps(Genes)
genes.bed[['Mit']] <- data.frame(chr = as.character(Genes@seqnames), start = Genes@ranges@start, end = (Genes@ranges@start + Genes@ranges@width - 1))
Intergenic.bed[['Mit']] <- data.frame(chr = as.character(Intergenic@seqnames), start = Intergenic@ranges@start, end = (Intergenic@ranges@start + Intergenic@ranges@width - 1))
genes.bed$Mit$chr <- rep("chrM", nrow(genes.bed$Mit))
Intergenic.bed$Mit$chr <- rep("chrM", nrow(Intergenic.bed$Mit))

genes.bed.f <- rbind(genes.bed$Nuclear, genes.bed$Chl, genes.bed$Mit)
genes.bed.f$chr <- factor(genes.bed.f$chr, levels, ordered = TRUE)
genes.bed.f <- genes.bed.f[order(genes.bed.f$chr, genes.bed.f$start),]
write.table(genes.bed.f, file = "../Chromatin States/all_together/inputs/Annotation_Enrichments/COORDS/ChlamytinaGenes.Cre55.bed", col.names = F, row.names = F, sep = "\t", 
            quote = F)
Intergenic.bed.f <- rbind(Intergenic.bed$Nuclear, Intergenic.bed$Chl, Intergenic.bed$Mit)
Intergenic.bed.f$chr <- factor(Intergenic.bed.f$chr, levels, ordered = TRUE)
Intergenic.bed.f <- Intergenic.bed.f[order(Intergenic.bed.f$chr, Intergenic.bed.f$start),]
write.table(Intergenic.bed.f, file = "../Chromatin States/all_together/inputs/Annotation_Enrichments/COORDS/ChlamytinaIntergenic.Cre55.bed", col.names = F, row.names = F, sep = "\t", 
            quote = F)

# UTRs/TSS/TSS2kb/TES/TSS Anchor/TES Anchoer

III.UTR <- unlist(threeUTRsByTranscript(All.Anno$Nuclear), use.names = F)
III.UTR.bed[['Nuclear']] <- data.frame(chr = as.character(III.UTR@seqnames), start = III.UTR@ranges@start, end = (III.UTR@ranges@start + III.UTR@ranges@width - 1))
III.UTR.bed$Nuclear$chr <- gsub("chromosome_", "chr", III.UTR.bed$Nuclear$chr)
filter <- grep("scaffold", III.UTR.bed$Nuclear$chr)
III.UTR.bed$Nuclear <- III.UTR.bed$Nuclear[-filter,]
III.UTR.bed$Nuclear <-  III.UTR.bed$Nuclear[!duplicated(III.UTR.bed$Nuclear),]
# no III.UTRs Chl and Mit

III.UTR.bed$Nuclear$chr <- factor(III.UTR.bed$Nuclear$chr, levels, ordered = TRUE)
III.UTR.bed$Nuclear <- III.UTR.bed$Nuclear[order(III.UTR.bed$Nuclear$chr, III.UTR.bed$Nuclear$start),]
write.table(III.UTR.bed$Nuclear, file = "../Chromatin States/all_together/inputs/Annotation_Enrichments/COORDS/Chlamytina3UTR.Cre55.bed", col.names = F, row.names = F, sep = "\t", 
            quote = F)

V.UTR <- unlist(fiveUTRsByTranscript(All.Anno$Nuclear), use.names = F)
V.UTR.bed[['Nuclear']] <- data.frame(chr = as.character(V.UTR@seqnames), start = V.UTR@ranges@start, end = (V.UTR@ranges@start + V.UTR@ranges@width - 1))
V.UTR.bed$Nuclear$chr <- gsub("chromosome_", "chr", V.UTR.bed$Nuclear$chr)
filter <- grep("scaffold", V.UTR.bed$Nuclear$chr)
V.UTR.bed$Nuclear <- V.UTR.bed$Nuclear[-filter,]
V.UTR.bed$Nuclear <-  V.UTR.bed$Nuclear[!duplicated(V.UTR.bed$Nuclear),]
# no III.UTRs Chl and Mit

V.UTR.bed$Nuclear$chr <- factor(V.UTR.bed$Nuclear$chr, levels, ordered = TRUE)
V.UTR.bed$Nuclear <- V.UTR.bed$Nuclear[order(V.UTR.bed$Nuclear$chr, V.UTR.bed$Nuclear$start),]
write.table(V.UTR.bed$Nuclear, file = "../Chromatin States/all_together/inputs/Annotation_Enrichments/COORDS/ChlamytinaVUTR.Cre55.bed", col.names = F, row.names = F, sep = "\t", 
            quote = F)

# TES

V.UTR <- fiveUTRsByTranscript(All.Anno$Nuclear)
III.UTR <- threeUTRsByTranscript(All.Anno$Nuclear)

TES.chr <- list()
TES.start <- list()
TES.end <- list()
TES.strand <- list()
for(i in 1:length(III.UTR)){
  TES <- III.UTR[i]
  TES.chr[[i]] <- as.character(TES@unlistData@seqnames@values)
  TES.strand[[i]] <- TES@unlistData@strand@values
  if(TES@unlistData@strand@values != "-"){
    TES.end[[i]] <- max(TES@unlistData@ranges@start + TES@unlistData@ranges@width - 1)
    TES.start[[i]] <- TES.end[[i]] - 1
  }else if(TES@unlistData@strand@values == "-"){
    TES.start[[i]] <- min(TES@unlistData@ranges@start + TES@unlistData@ranges@width - 1)
    TES.end[[i]] <- TES.start[[i]] + 1
  }
}

TES.COORD <- data.frame(chr = unlist(TES.chr), start = unlist(TES.start), end = unlist(TES.end), strad = unlist(TES.strand))
TES.COORD$chr <- gsub("chromosome_", "chr", TES.COORD$chr)
filter <- grep("scaffold", TES.COORD$chr)
TES.COORD <- TES.COORD[-filter,]
TES.COORD <-  TES.COORD[!duplicated(TES.COORD),]
TES.COORD$chr <- factor(TES.COORD$chr, levels, ordered = TRUE)
TES.COORD <- TES.COORD[order(TES.COORD$chr, TES.COORD$start),]
write.table(TES.COORD[,1:3], file = "../Chromatin States/all_together/inputs/Annotation_Enrichments/COORDS/ChlamytinaTES.Cre55.bed", col.names = F, row.names = F, sep = "\t", 
            quote = F)

TES.ANCHOR.start <- rep(NA, nrow(TES.COORD)) 
TES.ANCHOR.start[which(TES.COORD$strad != "-")] <- TES.COORD$end[which(TES.COORD$strad != "-")]
TES.ANCHOR.start[which(TES.COORD$strad == "-")] <- TES.COORD$start[which(TES.COORD$strad == "-")]
TES.ANCHOR <- data.frame(chr = TES.COORD$chr, start = TES.ANCHOR.start, strand = TES.COORD$strad)
write.table(TES.ANCHOR, file = "../Chromatin States/all_together/inputs/Annotation_Enrichments/ANCHOR/ChlamytinaTES.Cre55.bed", col.names = F, row.names = F, sep = "\t", 
            quote = F)

# TSS

TSS.chr <- list()
TSS.start <- list()
TSS.end <- list()
TSS.strand <- list()
for(i in 1:length(V.UTR)){
  TSS <- V.UTR[i]
  TSS.chr[[i]] <- as.character(TSS@unlistData@seqnames@values)
  TSS.strand[[i]] <- TSS@unlistData@strand@values
  if(TSS@unlistData@strand@values != "-"){
    TSS.start[[i]] <- min(TSS@unlistData@ranges@start + TSS@unlistData@ranges@width - 1)
    TSS.end[[i]] <- TSS.start[[i]] + 1
  }else if(TSS@unlistData@strand@values == "-"){
    TSS.end[[i]] <- max(TSS@unlistData@ranges@start + TSS@unlistData@ranges@width - 1)
    TSS.start[[i]] <- TSS.end[[i]] - 1
  }
}

TSS.COORD <- data.frame(chr = unlist(TSS.chr), start = unlist(TSS.start), end = unlist(TSS.end), strand = unlist(TSS.strand))
TSS.COORD$chr <- gsub("chromosome_", "chr", TSS.COORD$chr)
filter <- grep("scaffold", TSS.COORD$chr)
TSS.COORD <- TSS.COORD[-filter,]
TSS.COORD <-  TSS.COORD[!duplicated(TSS.COORD),]
TSS.COORD$chr <- factor(TSS.COORD$chr, levels, ordered = TRUE)
TSS.COORD <- TSS.COORD[order(TSS.COORD$chr, TSS.COORD$start),]
write.table(TSS.COORD[,1:3], file = "../Chromatin States/all_together/inputs/Annotation_Enrichments/COORDS/ChlamytinaTSS.Cre55.bed", col.names = F, row.names = F, sep = "\t", 
            quote = F)

TSS.ANCHOR.start <- rep(NA, nrow(TSS.COORD)) 
TSS.ANCHOR.start[which(TSS.COORD$strand != "-")] <- TSS.COORD$start[which(TSS.COORD$strand != "-")]
TSS.ANCHOR.start[which(TSS.COORD$strand == "-")] <- TSS.COORD$end[which(TSS.COORD$strand == "-")]
TSS.ANCHOR <- data.frame(chr = TSS.COORD$chr, start = TSS.ANCHOR.start, strand = TSS.COORD$strand)
write.table(TSS.ANCHOR, file = "../Chromatin States/all_together/inputs/Annotation_Enrichments/ANCHOR/ChlamytinaTSS.Cre55.bed", col.names = F, row.names = F, sep = "\t", 
            quote = F)

TSS.2kb.COORD <- TSS.COORD
TSS.2kb.COORD$start <- TSS.2kb.COORD$start - 2000
TSS.2kb.COORD$end <- TSS.2kb.COORD$end + 2000
write.table(TSS.2kb.COORD[,1:3], file = "../Chromatin States/all_together/inputs/Annotation_Enrichments/COORDS/ChlamytinaTSS2kb.Cre55.bed", col.names = F, row.names = F, sep = "\t", 
            quote = F)

# CpG Islands: EMBOSS CpGplot
for(i in 1:length(Chlamytina_genome)){
  Chr.fasta <- data.frame(names = names(Chlamytina_genome[i]), sequences = as.character(Chlamytina_genome[i]))
  tmp <- dataframe2fas(Chr.fasta, file = paste0("./", names(Chlamytina_genome[i]), ".fasta"))
}
# DONE

# Repeated and low comple regions
Cre.norm <- read.fasta("../Chromatin States/all_together/inputs/Assembly only chr (just in case)/Chlamydomonas_reinhardtii.Chlamydomonas_reinhardtii_v5.5.dna.toplevel.fa")
Cre.sm <- read.fasta("../Chromatin States/all_together/inputs/Assembly only chr (just in case)/Chlamydomonas_reinhardtii.Chlamydomonas_reinhardtii_v5.5.dna_sm.toplevel.fa")
Norm.splitt <- strsplit(x = Cre.norm, split = "")
Sm.splitt <- strsplit(x = Cre.sm, split = "")

Norm.splitt.seqs <- Norm.splitt[which((1:length(Norm.splitt) %% 2) != 1)]
Sm.split.seqs <- Sm.splitt[which((1:length(Sm.splitt) %% 2) != 1)]

Repeated_lowcomplex.bed <- list()
for(i in 1:17){
  logicalis <- Norm.splitt.seqs[[i]] == Sm.split.seqs[[i]]
  differences <- which(logicalis == FALSE)
  intervals <- list()
  for(j in 1:length(differences)){
    intervals[[j]] <- differences[j] == (differences[j - 1] + 1)
    if(length(differences[j] == (differences[j - 1] + 1)) == 0){ intervals[[j]] <- FALSE}
  }
  start <- differences[!unlist(intervals, use.names = F)]
  end <- differences[(which(unlist(intervals, use.names = F) == FALSE) - 1)[2:length(which(unlist(intervals, use.names = F) == FALSE))]]
  if(length(start) - length(end) == 1){
    end <- c(end, differences[length(differences)])
  }else if(length(start - length(end) > 1)){stop("\n WTF is this \n")}
  chr <- rep(paste0("chr",i), length(start))
  cat("\n Quality Check \n")
  if(length(which((start < end) == FALSE)) == 0){ cat("Passed")
    }else {stop("Check Start-End Coordinates")}
  bed <- data.frame(chr = chr, start = start, end = end) 
  Repeated_lowcomplex.bed[[i]] <- bed
}

RepLow <- do.call(rbind, Repeated_lowcomplex.bed)
write.table(RepLow, file = "../Chromatin States/all_together/inputs/Annotation_Enrichments/COORDS/ChlamytinaRepsLC.Cre55.bed", col.names = F, row.names = F, sep = "\t", 
            quote = F)

# Segmentation Wei2015 CS: One bed per State

CS_C <- read.delim("clipboard", header = F)
CS_N <- read.delim("clipboard", header = F)
CS_S <- read.delim("clipboard", header = F)

CS_Wao <- list(C = CS_C, N = CS_N, S = CS_S)

CS.beds <- lapply(CS_Wao, function(x){
  x[,1] <- gsub("chr_", "chr", x[,1])
  x[,1] <- gsub("chromosome_", "chr", x[,1])
  filter <- grep("scaffold", x[,1])
  x <- x[-filter,1:4]
  SState <- list()
  for(i in levels(as.factor(x[,4]))){
    SState[[i]] <- x[which(x[,4] == i), 1:3]
  }
  SState
})

for(j in 1:length(CS.beds)){
  for(i in 1:length(CS.beds[[j]])){
    write.table(CS.beds[[j]][[i]], file= paste0(names(CS.beds[j]), sep = "_CS", i, ".bed"), row.names = F, col.names = F, quote = F, sep = "\t")
  }
}
  








