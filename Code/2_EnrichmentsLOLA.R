#!/usr/bin/env Rscript

##### 0. Load/primary pkgs and args #####

# Primary pkgs
mypkgs1 <- c("cowsay", "optparse")
logicals1 <- is.element(mypkgs1, installed.packages()[,1])
tmp <- base::sapply(mypkgs1[logicals1], FUN = function(x){  suppressPackageStartupMessages(library(x, character.only = TRUE))})
tmp <- base::sapply(mypkgs1[!logicals1], FUN = function(x){
  install.packages(x, repos = "https://cloud.r-project.org/")
})
tmp <- base::sapply(mypkgs1[!logicals1], FUN = function(x){ suppressPackageStartupMessages(library(x, character.only = TRUE))})

# Arguments
option_list = list(
  make_option(c("-A", "--file1"), type = "character", default = NULL, help = "First BED file path. Any BED file with chr, start and end or DataPrepare output", metavar = "character"),
  make_option(c("-B", "--file2"), type = "character", default = NULL, help = "Second BED file path", metavar = "character"),
  make_option(c("-C", "--file3"), type = "character", default = NULL, help = "Third BED file path", metavar = "character"),
  make_option(c("-D", "--file4"), type = "character", default = NULL, help = "Fourth BED file path", metavar = "character"),
  make_option(c("-E", "--file5"), type = "character", default = NULL, help = "Fifth file path", metavar = "character"),
  make_option(c("-F", "--file6"), type = "character", default = NULL, help = "Sixth file path", metavar = "character"),
  make_option(c("-G", "--file7"), type = "character", default = NULL, help = "Seventh file path", metavar = "character"),
  make_option(c("-H", "--file8"), type = "character", default = NULL, help = "Eighth file path", metavar = "character"),
  make_option(c("-I", "--file9"), type = "character", default = NULL, help = "Nineth file path", metavar = "character"),
  make_option(c("-J", "--file10"), type = "character", default = NULL, help = "Tenth file path", metavar = "character"),
  make_option(c("-b", "--background"), type = "character", default = NULL, help = "Background BED file path. The set of regions tested for enrichments", metavar = "character"),
  make_option(c("-l", "--list"), type = "logical", default = FALSE, help = "If true, the rest of args are ignored and list all the possible files for one regionDB", metavar = "character"),
  make_option(c("-o", "--out"), type = "character", default = "./Outputs/", help = "Output directory [default = %default]", metavar = "character"),
  make_option(c("-r", "--database"), type = "character", default = "Marks", help = "regionDB used. Options: Marks (epigenetic marks) or States (chromatin states) [default = %default]", metavar = "character"),
  make_option(c("-c", "--cores"), type = "character", default = "1", help = "Number of cores [default = %default]", metavar = "character")
  );

opt_parser = OptionParser(option_list = option_list, usage = "2_EnrichmentsLOLA.R [file] [file] [file] [background] [database] ... [options]");
opt = parse_args(opt_parser);

if(as.logical(opt$list)){
  if(opt$database == "Marks"){
    cat("\n All:", list.files("./Data/regionDB/Chlamytina/Marks/regions/"))
  }else if(opt$database == "States"){
    cat("\n Waiting for chromHMM... \n")
  }
  stop(say(what = "STOP: Displaying file lists ", by = "poop"), call.=FALSE)
}

if (is.null(opt$file1)){
  print_help(opt_parser)
  stop(say(what = "STOP: At least one file is needed", by = "poop"), call.=FALSE)
}

if (is.null(opt$background)){
  print_help(opt_parser)
  stop(say(what = "STOP: Please select a background ", by = "poop"), call.=FALSE)
}

##### 1. Secondary pkgs: checks and install #####

say(paste0("Welcome to EnrichmentsLOLA ! ", Sys.time()), by = "rabbit", what_color = "white", by_color = "yellow")
cat("\n Checking-installing-loading needed libs and packages ... \n")

mypkgs <- c("simpleCache", "LOLA", "GenomicRanges", "dplyr", "data.table", "ggplot2", "reshape2", "pheatmap", "RColorBrewer", "scales")
logicals <- is.element(mypkgs, installed.packages()[,1])
tmp <- base::sapply(mypkgs[logicals], FUN = function(x){  suppressPackageStartupMessages(library(x, character.only = TRUE))})
tmp <- base::sapply(mypkgs[!logicals], FUN = function(x){
  if(x != "LOLA" & x != "GenomicRanges")
    install.packages(x, repos = "https://cloud.r-project.org/")
  if(x == "LOLA"){
    BiocManager::install("LOLA")
  }
  if(x == "GenomicRanges"){
    BiocManager::install("GenomicRanges")
  }
})
tmp <- base::sapply(mypkgs[!logicals], FUN = function(x){ suppressPackageStartupMessages(library(x, character.only = TRUE))})

##### 2. Import files and args #####

cat("\n Importing files and arguments ... \n")

if(length(grep(pattern = ".txt", x = opt$background, fixed = TRUE)) == 0 & length(grep(pattern = ".bed", x = opt$background, fixed = TRUE)) == 0){
  stop(say(what = "STOP: Non supported format. Please try .txt or .bed tab separated", by = "poop"), call.=FALSE)
}
background <- read.table(opt$background)
colnames(background) <- c("chr", "start", "end")
background <- background[,c(1,2,3)]
Universe <- makeGRangesFromDataFrame(background, ignore.strand = TRUE, keep.extra.columns = TRUE, seqnames.field = "chr", start.field = "start", end.field = "end")

inputs <- opt[grep("file", names(opt))]
inputs <- base::lapply(inputs, FUN = function(x){
  if(length(grep(pattern = ".txt", x = x, fixed = TRUE)) == 0 & length(grep(pattern = ".bed", x = x, fixed = TRUE)) == 0){
    stop(say(what = "STOP: Non supported format. Please try .txt or .bed tab separated", by = "poop"), call.=FALSE)
  }
  input <- read.table(x)
  colnames(input) <- c("chr", "start", "end")
  input <- input[,c(1,2,3)]
  input <- makeGRangesFromDataFrame(input, ignore.strand = TRUE, keep.extra.columns = TRUE, seqnames.field = "chr", start.field = "start", end.field = "end")
})
UserSets <- GRangesList(inputs)
opt$cores <- as.numeric(opt$cores)

##### 3. Load regionDB and enrichment #####

cat("\n Loading regionDB and running enrichments ... \n")
regionDB <- loadRegionDB(dbLocation = "./Data/regionDB/Chlamytina/", collections = opt$database)

regionResults <- lapply(UserSets, FUN = function(x){
  Results <- runLOLA(x, Universe, regionDB, cores = opt$cores)
  Results
})

##### 4. Plot results and export #####

cat("\n Plotting significative results (0.05) ... \n")
## HEATMAP

if(opt$database == "Marks"){
  for(i in 1:length(regionResults)){
    name <- strsplit(opt[[names(regionResults[i])]], split = ".", fixed = T)[[1]][1]
    if(name == ""){ name <- strsplit(opt[[names(regionResults[i])]], split = ".", fixed = T)[[1]][2]}
    name <- strsplit(name, split = "/")[[1]][length(strsplit(name, split = "/")[[1]])]
    regionResults[[i]]$file <- rep(name, nrow(regionResults[[i]]))
    regionResults[[i]]$oddsRatio[which(regionResults[[i]]$pValueLog < 1.30103)] <- 1
  }
  
  regionResults <- do.call(rbind, regionResults)
  regionResults <- regionResults[which(regionResults$oddsRatio > 1.0),]
  df <- data.frame(id = regionResults$file, condition = regionResults$description, variable = paste(regionResults$cellType, regionResults$tissue, regionResults$description, sep = " "), value = regionResults$oddsRatio)
  df <- reshape(df[,-2], idvar = "id", v.names = "value",  timevar = "variable", direction = "wide")
  df[is.na(df)] <- 1
  Conditions <- data.frame(Conditions = colnames(df)[-1])
  Conditions$Conditions <- as.character(Conditions$Conditions)
  for(i in 1:nrow(Conditions)){
    Conditions$Conditions[i] <- strsplit(as.character(Conditions$Conditions[i]), split = " ")[[1]][3]
  }
  colnames(df) <- gsub("value.", "", colnames(df), fixed = T)
  colnames(df) <- gsub("control", "", colnames(df), fixed = T)
  colnames(df) <- gsub("Nitrogen", "", colnames(df), fixed = T)
  colnames(df) <- gsub("Sulphur", "", colnames(df), fixed = T)
  rownames(df) <- df$id
  df <- df[,-1]
  rownames(Conditions) <- colnames(df)
  col <- colorRampPalette(brewer.pal(7, "Greens"))(100)
  df <- t(as.matrix(df))
  title <- strsplit(opt$background, split = ".", fixed = T)[[1]][1]
  if(title == ""){title <- strsplit(opt$background, split = ".", fixed = T)[[1]][2]}
  title <- strsplit(title, split = "/")[[1]][length(strsplit(title, split = "/")[[1]])]
  Annotation_colors <- list(Conditions = c(control = "darkgreen",Nitrogen = "cyan3", Sulphur = "gold3"))
  pdf(file = paste0(opt$out,title,".pdf"), paper = "a4r", height = 21, width = 28, onefile = T)
  pheatmap(df, scale = "none", color = col,  cluster_rows = F, cluster_cols = T,  clustering_method = "ward.D2", cellwidth = 75, annotation_row = Conditions, annotation_names_row = FALSE, main = title, annotation_colors = Annotation_colors)
  dev.off()
}

say(paste0("EnrichmentsLOLA has finished ! ", Sys.time()), by = "rabbit", what_color = "white", by_color = "yellow")




