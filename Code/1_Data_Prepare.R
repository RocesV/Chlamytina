#!/usr/bin/env Rscript

### Future: Add diff expression and whole proteome and whole differential backgrounds
### 0. Load/install packages

mypkgs <- c("utils", "nVennR", "readxl", "optparse", "GenomicFeatures", "R.utils")
logicals <- is.element(mypkgs, installed.packages()[,1])
base::sapply(mypkgs[logicals], FUN = function(x){  suppressPackageStartupMessages(library(x, character.only = TRUE))})
base::sapply(mypkgs[!logicals], FUN = function(x){
  if(x != "GenomicFeatures")
  install.packages(x)
  if(x == "GenomicFeatures"){
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    
    BiocManager::install("GenomicFeatures")
  }
})
base::sapply(mypkgs[!logicals], FUN = function(x){ suppressPackageStartupMessages(library(x, character.only = TRUE))})

### 1. Define args

option_list = list(
  make_option(c("-f", "--file1"), type = "character", default = NULL, help = "Dataset1 file name", metavar = "character"),
  make_option(c("-f", "--file2"), type = "character", default = NULL, help = "Dataset2 file name", metavar = "character"),
  make_option(c("-f", "--file3"), type = "character", default = NULL, help = "Dataset3 file name", metavar = "character"),
  make_option(c("-i", "--intersect"), type = "logical", default = TRUE, help = "CreIDs intra-inter group specific discrimination [default = %default]", metavar = "character"),
  make_option(c("-o", "--out"), type = "character", default = "./Data/BED_Input/", help = "Output directory [default = %default]", metavar = "character"),
  make_option(c("-c", "--chromosome"), type = "logical", default = TRUE, help = "If true, non-chromosome mapped (scaffolds ...) proteins are not taked into account [default = %default]", metavar = "character")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file1)){
  print_help(opt_parser)
  stop("At least one Dataset is needed", call.=FALSE)
}

if(is.null(opt$file2)){ opt$intersect = FALSE}

### 2. Program
unzip(zipfile = "./Data/DB/Phytozome_download.zip", exdir = "./Data/DB/unzip")

inputs <- opt[grep("f",names(opt))]
inputs <- base::lapply(inputs, FUN = function(x){
  if(is.null(x)){pass}
  if(!is.null(x)){x}
})
inputs <- base::lapply(inputs, FUN = function(x){
  if(length(grep(pattern = ".txt", x = x, fixed = TRUE)) == 1){
    input <- read.table(x, header = TRUE, sep = " ")
  }else if(length(grep(pattern = ".xls", x = x, fixed = TRUE)) == 1){
    input <- read_excel(x, col_names = TRUE)
  }else if(length(grep(pattern = ".xlsx", x = x, fixed = TRUE)) == 1){
    input <- read_excel(x, col_names = TRUE)
  }else(stop("Non supported format. Please try txt or excel tab separated with headers"))
  input
})

## pass to 5.5 version
DB <- read.table("./Data/DB/ChlamydomonasTranscriptNameConversionBetweenReleases.Mch12b.txt", header = F)
colnames(DB) <- c("5.5", "3.1", "Genbank","4", "4.3", "u5", "u9", "5.3.1")

inputs <- lapply(inputs, FUN = function(x){
  x <- x[,grep("Cre*.", x)[1]]
  versions <- list()
  for(i in 1:ncol(DB)){
    versions[i] <- length(which((as.data.frame(x)[,1] %in% DB[,i]) == TRUE))
  }
  version <- colnames(DB)[which(unlist(versions,use.names = F) == max(unlist(versions,use.names = F)))]
  if(as.numeric(version) == 5.5){
    x <- as.data.frame(x)[,1]
    x <- x[!is.na(x)] 
    x
  }else if(as.numeric(version) != 5.5){
    newIDs <- list()
    for(i in 1:nrow(as.data.frame(x))){
      if(length(as.character(DB[which(DB[,version] == as.data.frame(x)[i,1]),1])[1]) != 0){
        newIDs[[i]] <- as.character(DB[which(DB[,version] == as.data.frame(x)[i,1]),1])[1] 
      }else if(length(as.character(DB[which(DB[,version] == as.data.frame(x)[i,1]),1])[1]) == 0){
        pass
      }}
    newIDs <- unlist(newIDs, use.names = F)
    newIDs <- newIDs[!is.na(newIDs)]
    newIDs
  }
})


## intersection
if(opt$intersect){
  cat("\n Intersect is defined as TRUE \n")
  myV <- nVennR::plotVenn(inputs, nCycles = 7000, opacity = 0.2, borderWidth = 3, systemShow = T, fontScale = 2)
  names <- lapply(myV$reg, FUN = function(x){
    file1 <- length(which((x %in% myV$orig$file1) == TRUE))
    file2 <- length(which((x %in% myV$orig$file2) == TRUE))
    file3 <- length(which((x %in% myV$orig$file3) == TRUE))
    if(file1 > 0 & file2 > 0 & file3 > 0){ "file1_file2_file3"
      }else if(file1 > 0 & file2 > 0 & file3 == 0){ "file1_file2"
      }else if(file1 > 0 & file2 == 0 & file3 > 0){ "file1_file3"
      }else if(file1 == 0 & file2 > 0 & file3 > 0){ "file2_file3"
      }else if(file1 == 0 & file2 == 0 & file3 > 0){ "file3"
      }else if(file1 == 0 & file2 > 0 & file3 == 0){ "file2"
      }else if(file1 > 0 & file2 == 0 & file3 == 0){ "file1"
      }
  })
  names(myV$reg) <- unlist(names, use.names = F)
  inputs <- myV$reg
} else{cat("\n Intersect is defined as FALSE \n")}

## convert to BED
R.utils::gunzip(filename = "./Data/DB/unzip/Phytozome/PhytozomeV12_unrestricted/Creinhardtii/annotation/Creinhardtii_281_v5.5.gene_exons.gff3.gz", destname = "./Data/DB/unzip/Phytozome/PhytozomeV12_unrestricted/Creinhardtii/annotation/Creinhardtii_281_v5.5.gene_exons.gff3")
TxDb.Cre <- makeTxDbFromGFF(file = "./Data/DB/unzip/Phytozome/PhytozomeV12_unrestricted/Creinhardtii/annotation/Creinhardtii_281_v5.5.gene_exons.gff3")
Transcripts <- transcripts(TxDb.Cre)

#QC
cat("\n Computing CreIDs coherence ... \n")
QC <- lapply(inputs, function(x){
  x %in% Transcripts$tx_name
})
if(length(which(unlist(QC,use.names = F) == FALSE)) > 0){stop("\n Something go wrong with CreIDs. Check \n")}
if(length(which(unlist(QC,use.names = F) == FALSE)) == 0){cat("\n Passed \n")}

AllIDs <- unlist(inputs, use.names = F)
chr <- list()
start <- list()
end <- list()
for(j in 1:length(AllIDs)){
    chr[j] <- as.character(Transcripts@seqnames[which(Transcripts$tx_name == AllIDs[j])])
    start[j] <- Transcripts@ranges@start[which(Transcripts$tx_name == AllIDs[j])]
    end[j] <- Transcripts@ranges@start[which(Transcripts$tx_name == AllIDs[j])] + Transcripts@ranges@width[which(Transcripts$tx_name == AllIDs[j])] - 1
  }

All.bed <- data.frame(row.names = AllIDs, unlist(chr, use.names = F), unlist(start, use.names = F), unlist(end, use.names = F))

inputs.2 <- lapply(inputs, FUN = function(x){
  input <- All.bed[x,]
  input[,1] <- gsub("chromosome_", "chr", input[,1])
  if(opt$chromosome){
    cat("\n Chromosome filtering is defines as TRUE \n")
    filter <- grep("scaffold", input[,1])
    if(length(filter) == 0){
      input
    }else if(length(filter) > 0){
      input <- input[-filter,]
      input
      }
  }else if(opt$chromosome){ input}
})

levels <- paste0("chr", c(1:17))

for(i in 1:length(inputs.2)){
  colnames(inputs.2[[i]]) <- c("chr", "start", "end")
  inputs.2[[i]]$chr <- factor(inputs.2[[i]]$chr, levels, ordered = TRUE)
  inputs.2[[i]] <- inputs.2[[i]][order(inputs.2[[i]]$chr, inputs.2[[i]]$start),] 
  write.table(inputs.2[i], file = paste0(opt$out,names(inputs.2[i]), ".bed"), row.names = FALSE, col.names = F, sep = "\t", quote = FALSE)
}

