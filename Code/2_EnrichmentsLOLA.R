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
  make_option(c("-l", "--list"), type = "logical", default = FALSE, help = "If true, the rest of args are ignored and list all the files for one regionDB", metavar = "character"),
  make_option(c("-o", "--out"), type = "character", default = "./Outputs/", help = "Output directory [default = %default]", metavar = "character"),
  make_option(c("-r", "--database"), type = "character", default = "MMarks", help = "regionDB used. Options: Marks (epigenetic marks by original conditions), MMarks (merged Marks wo conditions) or CS_Control, CS_N, CS_S (Ngan et al., Nat.Plants 2015, Chromatin States !Nitrogen !Sulfur) or CS_Chlamytina (Updated Chromatin states with 5mC, 6mA and MNase) [default = %default]", metavar = "character"),
  make_option(c("-c", "--cores"), type = "character", default = "1", help = "Number of cores [default = %default]", metavar = "character")
  );

opt_parser = OptionParser(option_list = option_list, usage = "2_EnrichmentsLOLA.R [file] [file] [file] [background] [database] ... [options]");
opt = parse_args(opt_parser);

if(as.logical(opt$list)){
  if(opt$database == "Marks"){
    cat("\n All files:", list.files("./Data/regionDB/Chlamytina/Marks/regions/"))
  }else if(opt$database == "MMarks"){
    cat("\n All files:", list.files("./Data/regionDB/Chlamytina/MMarks/regions/"))
  }else if(opt$database == "CS_Control"){
    cat("\n All files:", list.files("./Data/regionDB/Chlamytina/CS_Control/regions/"))
  }else if(opt$database == "CS_N"){
    cat("\n All files:", list.files("./Data/regionDB/Chlamytina/CS_N/regions/"))
  }else if(opt$database == "CS_S"){
    cat("\n All files:", list.files("./Data/regionDB/Chlamytina/CS_S/regions/"))
  }else if(opt$database == "CS_Chlamytina"){
    cat("\n All files:", list.files("./Data/regionDB/Chlamytina/CS_Chlamytina/regions/"))
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
    regionResults[[i]]$oddsRatio[which(regionResults[[i]]$pValueLog < 1.30103)] <- 0.5
  }
  
  regionResults <- do.call(rbind, regionResults)
  regionResults <- regionResults[which(regionResults$oddsRatio > 0.5),]
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
  col <- colorRampPalette(brewer.pal(9, "Purples"))(250)
  df <- t(as.matrix(df))
  df[is.infinite(df)] <- 25
  title <- strsplit(opt$background, split = ".", fixed = T)[[1]][1]
  if(title == ""){title <- strsplit(opt$background, split = ".", fixed = T)[[1]][2]}
  title <- strsplit(title, split = "/")[[1]][length(strsplit(title, split = "/")[[1]])]
  Annotation_colors <- list(Conditions = c(control = "darkgreen",Nitrogen = "cyan3", Sulphur = "gold3", light = "white", dark = "black"))
  pdf(file = paste0(opt$out,title,opt$database,".pdf"), paper = "a4r", height = 21, width = 28, onefile = T)
  pheatmap(df, scale = "none", color = col,  cluster_rows = F, cluster_cols = T,  clustering_method = "ward.D2", cellwidth = 75, annotation_row = Conditions, annotation_names_row = FALSE, main = title, annotation_colors = Annotation_colors)
  dev.off()
}

if(opt$database == "MMarks"){
  for(i in 1:length(regionResults)){
    name <- strsplit(opt[[names(regionResults[i])]], split = ".", fixed = T)[[1]][1]
    if(name == ""){ name <- strsplit(opt[[names(regionResults[i])]], split = ".", fixed = T)[[1]][2]}
    name <- strsplit(name, split = "/")[[1]][length(strsplit(name, split = "/")[[1]])]
    regionResults[[i]]$file <- rep(name, nrow(regionResults[[i]]))
    regionResults[[i]]$oddsRatio[which(regionResults[[i]]$pValueLog < 1.30103)] <- 0.5
  }
  
  regionResults <- do.call(rbind, regionResults)
  df <- data.frame(id = regionResults$file, condition = regionResults$description, variable = paste(regionResults$cellType, regionResults$tissue, regionResults$description, sep = " "), value = regionResults$oddsRatio)
  df <- reshape(df[,-2], idvar = "id", v.names = "value",  timevar = "variable", direction = "wide")
  Conditions <- data.frame(Conditions = colnames(df)[-1])
  Conditions$Conditions <- as.character(Conditions$Conditions)
  for(i in 1:nrow(Conditions)){
    Conditions$Conditions[i] <- strsplit(as.character(Conditions$Conditions[i]), split = " ")[[1]][3]
  }
  colnames(df) <- gsub("value.", "", colnames(df), fixed = T)
  colnames(df) <- gsub("NA", "", colnames(df), fixed = T)
  colnames(df) <- gsub("active", "", colnames(df), fixed = T)
  colnames(df) <- gsub("repressive", "", colnames(df), fixed = T)
  rownames(df) <- df$id
  df <- df[,-1]
  colnames(df)[grep("nucleosome", colnames(df))] <- "nucleosome"
  rownames(Conditions) <- colnames(df)
  col <- colorRampPalette(brewer.pal(9, "Greens"))(250)
  df <- t(as.matrix(df))
  df[is.infinite(df)] <- 25
  title <- strsplit(opt$background, split = ".", fixed = T)[[1]][1]
  if(title == ""){title <- strsplit(opt$background, split = ".", fixed = T)[[1]][2]}
  title <- strsplit(title, split = "/")[[1]][length(strsplit(title, split = "/")[[1]])]
  Annotation_colors <- list(Conditions = c(active = "lightsteelblue1", repressive = "darksalmon", nucleosome = "gold3"))
  pdf(file = paste0(opt$out,title,opt$database,".pdf"), paper = "a4r", height = 21, width = 28, onefile = T)
  pheatmap(df, scale = "none", color = col,  cluster_rows = F, cluster_cols = T,  clustering_method = "ward.D2", cellwidth = 75, annotation_row = Conditions, annotation_names_row = FALSE, main = title, annotation_colors = Annotation_colors)
  dev.off()
}

if(opt$database == "CS_Control" | opt$database == "CS_N" | opt$database == "CS_S"){
  for(i in 1:length(regionResults)){
    name <- strsplit(opt[[names(regionResults[i])]], split = ".", fixed = T)[[1]][1]
    if(name == ""){ name <- strsplit(opt[[names(regionResults[i])]], split = ".", fixed = T)[[1]][2]}
    name <- strsplit(name, split = "/")[[1]][length(strsplit(name, split = "/")[[1]])]
    regionResults[[i]]$file <- rep(name, nrow(regionResults[[i]]))
    regionResults[[i]]$oddsRatio[which(regionResults[[i]]$pValueLog < 1.30103)] <- 0.5
  }
  
  regionResults <- do.call(rbind, regionResults)
  df <- data.frame(id = regionResults$file, functions = regionResults$description, location = regionResults$cellType , name = rep(NA, 16), value = regionResults$oddsRatio)
  for(i in 1:nrow(regionResults)){
    name <- strsplit(regionResults$filename[i], split = ".", fixed = T)[[1]][1]
    df$name[i] <- name 
    regionResults$filename[i] <- name
  }
  df2 <- reshape(df[,-c(2,3)], idvar = "id", v.names = "value",  timevar = "name", direction = "wide")
  colnames(df2) <- gsub("value.", "", colnames(df2), fixed = T)
  colnames(df2) <- gsub("NA", "", colnames(df2), fixed = T)
  rownames(df2) <- df2$id
  df2 <- df2[,-1]
  if(opt$database == "CS_Control"){
    Conditions <- data.frame(row.names = paste0("C_CS",c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)), functions = rep(NA,16), locations = rep(NA, 16), evolution =rep(NA,16)) 
  }else if(opt$database == "CS_N"){
    Conditions <- data.frame(row.names = paste0("N_CS",c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)), functions = rep(NA,16), locations = rep(NA, 16), evolution =rep(NA,16))
  }else if(opt$database == "CS_S"){
    Conditions <- data.frame(row.names = paste0("S_CS",c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)), functions = rep(NA,16), locations = rep(NA, 16), evolution =rep(NA,16))
  }
  x <- 0
  for(i in rownames(Conditions)){
    x <- x+1
    Conditions$functions[x] <- regionResults$description[which(regionResults$filename == i)[1]]
    Conditions$locations[x] <- regionResults$cellType[which(regionResults$filename == i)[1]]
    Conditions$evolution[x] <- regionResults$tissue[which(regionResults$filename == i)[1]]
  }
  col <- colorRampPalette(brewer.pal(9, "OrRd"))(250)
  df2 <- t(as.matrix(df2))
  df2[is.infinite(df2)] <- 25
  title <- strsplit(opt$background, split = ".", fixed = T)[[1]][1]
  if(title == ""){title <- strsplit(opt$background, split = ".", fixed = T)[[1]][2]}
  title <- strsplit(title, split = "/")[[1]][length(strsplit(title, split = "/")[[1]])]
  Annotation_colors <- list(functions = c(Transcribed = "lightsteelblue1", Repressed = "darksalmon", Bivalent = "gold3", Promoter = "lightsteelblue3", Heterochromatin = "black", 'no info' = "grey"),
                            locations = c('no info' = "grey", Intragenic = "darkseagreen", "3'gene" = "darkseagreen1", "5'gene" = "darkseagreen2"),
                            evolution = c(Conserved = "burlywood", Algal = "aquamarine"))
  pdf(file = paste0(opt$out,title,opt$database,".pdf"), paper = "a4r", height = 21, width = 28, onefile = T)
  pheatmap(df2, scale = "none", color = col,  cluster_rows = F, cluster_cols = T,  clustering_method = "ward.D2", cellwidth = 75, annotation_row = Conditions, annotation_names_row = FALSE, main = title, annotation_colors = Annotation_colors)
  
  dev.off()
}

if(opt$database == "CS_Chlamytina" ){
  for(i in 1:length(regionResults)){
    name <- strsplit(opt[[names(regionResults[i])]], split = ".", fixed = T)[[1]][1]
    if(name == ""){ name <- strsplit(opt[[names(regionResults[i])]], split = ".", fixed = T)[[1]][2]}
    name <- strsplit(name, split = "/")[[1]][length(strsplit(name, split = "/")[[1]])]
    regionResults[[i]]$file <- rep(name, nrow(regionResults[[i]]))
    regionResults[[i]]$oddsRatio[which(regionResults[[i]]$pValueLog < 1.30103)] <- 0.5
  }
  
  regionResults <- do.call(rbind, regionResults)
  df <- data.frame(id = regionResults$file, pref.marks = regionResults$description, location = regionResults$cellType , name = rep(NA, nrow(regionResults)), value = regionResults$oddsRatio)
  for(i in 1:nrow(regionResults)){
    name <- strsplit(regionResults$filename[i], split = ".", fixed = T)[[1]][1]
    name <- gsub("E", "CS", x = name)
    df$name[i] <- name 
    regionResults$filename[i] <- name
  }
  df2 <- reshape(df[,-c(2,3)], idvar = "id", v.names = "value",  timevar = "name", direction = "wide")
  colnames(df2) <- gsub("value.", "", colnames(df2), fixed = T)
  colnames(df2) <- gsub("NA", "", colnames(df2), fixed = T)
  rownames(df2) <- df2$id
  df2 <- df2[,-1]
  Conditions <- data.frame(row.names = paste0("CS",c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)), functions = rep(NA,23), locations = rep(NA, 23)) 
  x <- 0
  for(i in rownames(Conditions)){
    x <- x+1
    Conditions$functions[x] <- regionResults$description[which(regionResults$filename == i)[1]]
    Conditions$locations[x] <- regionResults$cellType[which(regionResults$filename == i)[1]]
  }
  Conditions$functions <- gsub("Nucleosomes", "Nucl", Conditions$functions)
  Conditions$functions <- gsub("RNA", "R", Conditions$functions)
  Conditions$functions <- gsub("X", "", Conditions$functions)
  Conditions$locations <- gsub("CpGIslands", "CpGI", Conditions$locations)
  Conditions$locations <- gsub("Intron", "Intr", Conditions$locations)
  col <- colorRampPalette(brewer.pal(9, "OrRd"))(250)
  df2 <- t(as.matrix(df2))
  df2[is.infinite(df2)] <- 25
  title <- strsplit(opt$background, split = ".", fixed = T)[[1]][1]
  if(title == ""){title <- strsplit(opt$background, split = ".", fixed = T)[[1]][2]}
  title <- strsplit(title, split = "/")[[1]][length(strsplit(title, split = "/")[[1]])]
  pref.marks <- Conditions$functions
  locations <- levels(as.factor(Conditions$locations))
  col.marks <- c(brewer.pal(12, "Set3"), brewer.pal(11, "Paired"))
  names(col.marks) <- pref.marks
  col.locations <- c(brewer.pal(11, "BrBG"), brewer.pal(10,"PuOr"))
  names(col.locations) <- locations
  Annotation_colors <- list(functions = col.marks, locations = col.locations)
  pdf(file = paste0(opt$out,title,opt$database,".pdf"), paper = "a4r", height = 21, width = 28, onefile = T)
  pheatmap(df2, scale = "none", color = col,  cluster_rows = F, cluster_cols = T,  clustering_method = "ward.D2", cellwidth = 50, annotation_row = Conditions, annotation_names_row = FALSE, main = title, annotation_colors = Annotation_colors, fontsize = 7, fontsize_row = 10, fontsize_col = 10)
  
  dev.off()
}

say(paste0("EnrichmentsLOLA has finished ! ", Sys.time()), by = "rabbit", what_color = "white", by_color = "yellow")




