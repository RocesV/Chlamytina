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
options(warn = -1)
chlamytina_say <- function(text) {
  msg <- cowsay::say(
    text,
    by = "rabbit",
    type = "string",
    what_color = "white",
    by_color = "yellow",
    width = 60
  )
  msg <- gsub("\\[nosig\\]", "[Roces, V]", msg)
  cat(msg, "\n", sep = "")
}

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
  make_option(c("-m", "--significance_metric"), type = "character", default = "pValueLog", help = "Metric for target LOLA significance. Options: qValue, pValueLog [default = %default]", metavar = "character"),
  make_option(c("-t", "--significance_threshold"), type = "double", default = NA, help = "Threshold for --significance_metric. Defaults: 0.05 for qValue; 1.30103 for pValueLog", metavar = "numeric"),
  make_option(c("-p", "--permutation_path"), type = "character", default = NULL, help = "Folder with expression-matched permutation BEDs from 1_DataPrepare.R. Useful when only one contrast is available", metavar = "character"),
  make_option(c("-e", "--empirical_threshold"), type = "double", default = 0.05, help = "Empirical FDR threshold for permutation mode [default = %default]", metavar = "numeric"),
  make_option(c("-r", "--database"), type = "character", default = "MMarks", help = "regionDB used. Options: Marks (epigenetic marks by original conditions), MMarks (merged Marks wo conditions) or CS_Control, CS_N, CS_S (Ngan et al., Nat.Plants 2015, Chromatin States !Nitrogen !Sulfur) or CS_Chlamytina (Updated Chromatin states with 5mC, 6mA and MNase) [default = %default]", metavar = "character"),
  make_option(c("-V", "--version"), type = "character", default = "v5", help = "Genome/regionDB version. Options: v5, v6 [default = %default]", metavar = "character"),
  make_option(c("-c", "--cores"), type = "integer", default = 1, help = "Number of cores [default = %default]", metavar = "integer")
  );

opt_parser = OptionParser(option_list = option_list, usage = "2_EnrichmentsLOLA.R [file] [file] [file] [background] [database] ... [options]");
opt = parse_args(opt_parser);
normalize_genome_version <- function(version) {
  version <- tolower(version)
  if(version %in% c("v5", "5", "v5.5", "5.5", "v5.6", "5.6")) return("v5")
  if(version %in% c("v6", "6", "v6.1", "6.1")) return("v6")
  stop(say("STOP: Unsupported genome version. Use v5 or v6", by = "poop"), call. = FALSE)
}

opt$version <- normalize_genome_version(opt$version)
region_db_root <- if(opt$version == "v5") "./Data/regionDB/Chlamytina_v5" else "./Data/regionDB/Chlamytina_v6"
normalize_significance_metric <- function(metric) {
  metric <- tolower(metric)
  if(metric == "qvalue") return("qValue")
  if(metric == "pvaluelog") return("pValueLog")
  stop(say("STOP: Unsupported significance metric. Use qValue or pValueLog", by = "poop"), call. = FALSE)
}
opt$significance_metric <- normalize_significance_metric(opt$significance_metric)
if(is.na(opt$significance_threshold)){
  opt$significance_threshold <- if(opt$significance_metric == "qValue") 0.05 else -log10(0.05)
}
if(is.na(opt$significance_threshold) || opt$significance_threshold <= 0){
  stop(say("STOP: --significance_threshold must be numeric and > 0", by = "poop"), call. = FALSE)
}
if(is.na(opt$empirical_threshold) || opt$empirical_threshold <= 0 || opt$empirical_threshold > 1){
  stop(say("STOP: --empirical_threshold must be numeric, > 0 and <= 1", by = "poop"), call. = FALSE)
}
if(!is.null(opt$permutation_path) && opt$permutation_path == ""){
  opt$permutation_path <- NULL
}

if(as.logical(opt$list)){
  collection_dir <- file.path(region_db_root, opt$database, "regions")
  if(!dir.exists(collection_dir)){
    stop(say(paste0("STOP: regionDB collection not found: ", collection_dir), by = "poop"), call. = FALSE)
  }
  cat("\n All files:", list.files(collection_dir))
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

if(!grepl("[/\\\\]$", opt$out)){
  opt$out <- paste0(opt$out, .Platform$file.sep)
}
if(!dir.exists(opt$out)){
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
}
if(!dir.exists(opt$out)){
  stop(say(paste0("STOP: Output directory could not be created: ", opt$out), by = "poop"), call. = FALSE)
}

##### 1. Secondary pkgs: checks and install #####

chlamytina_say(paste0("Welcome to EnrichmentsLOLA ! ", Sys.time()))
cat("\n Checking-installing-loading needed libs and packages ... \n")

mypkgs <- c("BiocManager", "simpleCache", "LOLA", "GenomicRanges", "qvalue", "dplyr", "data.table", "ggplot2", "reshape2", "pheatmap", "RColorBrewer", "scales")
logicals <- is.element(mypkgs, installed.packages()[,1])
tmp <- base::sapply(mypkgs[logicals], FUN = function(x){  suppressPackageStartupMessages(library(x, character.only = TRUE))})
bioc_pkgs <- c("LOLA", "GenomicRanges", "qvalue")
tmp <- base::sapply(mypkgs[!logicals], FUN = function(x){
  if(x %in% bioc_pkgs){
    BiocManager::install(x, update = FALSE, ask = FALSE)
  }else{
    install.packages(x, repos = "https://cloud.r-project.org/")
  }
})
tmp <- base::sapply(mypkgs[!logicals], FUN = function(x){ suppressPackageStartupMessages(library(x, character.only = TRUE))})

##### 1.5 Some helper functions #####

quiet_lola <- function(expr) {
  invisible(utils::capture.output(value <- suppressMessages(force(expr))))
  value
}
get_file_label <- function(path) {
  name <- basename(path)
  sub("\\.[^.]*$", "", name)
}
read_bed_granges <- function(path) {
  if(length(grep(pattern = ".txt", x = path, fixed = TRUE)) == 0 & length(grep(pattern = ".bed", x = path, fixed = TRUE)) == 0){
    stop(say(what = "STOP: Non supported format. Please try .txt or .bed tab separated", by = "poop"), call.=FALSE)
  }
  bed <- read.table(path)
  if(ncol(bed) < 3){
    stop(say(paste0("STOP: BED file has fewer than 3 columns: ", path), by = "poop"), call. = FALSE)
  }
  colnames(bed)[1:3] <- c("chr", "start", "end")
  bed <- bed[,c(1,2,3)]
  GenomicRanges::makeGRangesFromDataFrame(bed, ignore.strand = TRUE, keep.extra.columns = TRUE, seqnames.field = "chr", start.field = "start", end.field = "end")
}
lola_is_significant <- function(results, metric, threshold) {
  if(metric == "qValue"){
    if(!("qValue" %in% colnames(results))){
      stop(say("STOP: qValue column is not available in LOLA results. Install qvalue or use --significance_metric pValueLog", by = "poop"), call. = FALSE)
    }
    return(!is.na(results$qValue) & results$qValue <= threshold)
  }
  !is.na(results$pValueLog) & results$pValueLog >= threshold
}
apply_lola_significance_floor <- function(results, metric, threshold) {
  significant <- lola_is_significant(results, metric, threshold)
  results$oddsRatio[!significant] <- 0.5
  results
}
make_heatmap_breaks <- function(mat, n_colors) {
  vals <- as.numeric(mat)
  vals <- vals[is.finite(vals)]
  if(length(vals) == 0) return(seq(0, 1, length.out = n_colors + 1))
  r <- range(vals)
  if(r[1] == r[2]){
    pad <- max(abs(r[1]) * 0.01, 0.01)
    r <- r + c(-pad, pad)
  }
  seq(r[1], r[2], length.out = n_colors + 1)
}
add_lola_file_names <- function(regionResults, opt) {
  for(i in seq_along(regionResults)){
    name <- get_file_label(opt[[names(regionResults[i])]])
    regionResults[[i]]$file <- rep(name, nrow(regionResults[[i]]))
  }
  regionResults
}
write_lola_results_table <- function(results, path) {
  write.table(as.data.frame(results), file = path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}
read_permutation_sets <- function(permutation_path) {
  if(!dir.exists(permutation_path)){
    stop(say(paste0("STOP: permutation folder not found: ", permutation_path), by = "poop"), call. = FALSE)
  }
  permutation_files <- sort(list.files(permutation_path, pattern = "\\.bed$", full.names = TRUE))
  if(length(permutation_files) == 0){
    stop(say(paste0("STOP: no permutation BED files found in: ", permutation_path), by = "poop"), call. = FALSE)
  }
  permutation_sets <- vector("list", length(permutation_files))
  for(i in seq_along(permutation_files)){
    if(i == 1 || i %% 25 == 0 || i == length(permutation_files)){
      cat(paste0("\n \t Loading permutation ", i, "/", length(permutation_files), ": ", basename(permutation_files[i]), " \n"))
    }
    permutation_sets[[i]] <- read_bed_granges(permutation_files[i])
  }
  names(permutation_sets) <- sub("\\.bed$", "", basename(permutation_files))
  GenomicRanges::GRangesList(permutation_sets)
}
calculate_permutation_results <- function(combinedResults, target_user_set, metric, threshold, empirical_threshold) {
  combinedResults$userSet <- as.character(combinedResults$userSet)
  targetResults <- combinedResults[combinedResults$userSet == target_user_set, , drop = FALSE]
  permutationResults <- combinedResults[combinedResults$userSet != target_user_set, , drop = FALSE]
  targetResults$targetSignificant <- lola_is_significant(targetResults, metric, threshold)
  n_perm <- integer(nrow(targetResults))
  n_perm_ge <- integer(nrow(targetResults))
  empirical_p <- numeric(nrow(targetResults))
  for(i in seq_len(nrow(targetResults))){
    p <- permutationResults$pValueLog[permutationResults$dbSet == targetResults$dbSet[i]]
    n_perm[i] <- length(p)
    n_perm_ge[i] <- sum(p >= targetResults$pValueLog[i], na.rm = TRUE)
    empirical_p[i] <- (n_perm_ge[i] + 1) / (n_perm[i] + 1)
  }
  targetResults$nPermutations <- n_perm
  targetResults$nPermutationGreaterEqual <- n_perm_ge
  targetResults$empiricalPValue <- empirical_p
  targetResults$empiricalQValue <- p.adjust(empirical_p, method = "BH")
  targetResults$permutationSignificant <- !is.na(targetResults$empiricalQValue) & targetResults$empiricalQValue <= empirical_threshold
  targetResults$passedPermutationFilter <- targetResults$targetSignificant & targetResults$permutationSignificant
  targetResults
}

##### 2. Import files and args #####

cat("\n Importing files and arguments ... \n")

Universe <- read_bed_granges(opt$background)
inputs <- opt[grep("file", names(opt))]
inputs <- base::lapply(inputs, FUN = function(x){
  read_bed_granges(x)
})
UserSets <- GenomicRanges::GRangesList(inputs)
opt$cores <- as.numeric(opt$cores)
if(is.na(opt$cores) || opt$cores < 1){
  stop(say("STOP: --cores must be numeric and >= 1", by = "poop"), call. = FALSE)
}

##### 3. Load regionDB and enrichment #####

cat("\n Loading regionDB and running enrichments ... \n")
regionDB <- quiet_lola(LOLA::loadRegionDB(dbLocation = region_db_root, collections = opt$database))
background_title <- get_file_label(opt$background)

if(is.null(opt$permutation_path)){
  regionResults <- lapply(UserSets, FUN = function(x){
    quiet_lola(LOLA::runLOLA(x, Universe, regionDB, cores = opt$cores))
  })
  regionResults <- add_lola_file_names(regionResults, opt)
  allRegionResults <- do.call(rbind, regionResults)
  if(opt$significance_metric == "qValue" && !("qValue" %in% colnames(allRegionResults))){
    stop(say("STOP: qValue column is not available in LOLA results. Install qvalue or use --significance_metric pValueLog", by = "poop"), call. = FALSE)
  }
  write_lola_results_table(allRegionResults, paste0(opt$out, background_title, opt$database, "_LOLA_results.txt"))
}else{
  if(length(UserSets) != 1){
    stop(say("STOP: --permutation_path can only be used with one target BED file", by = "poop"), call. = FALSE)
  }
  permutationSets <- read_permutation_sets(opt$permutation_path)
  target_option_name <- names(UserSets)[1]
  target_file_label <- get_file_label(opt[[target_option_name]])
  combinedSets <- GenomicRanges::GRangesList(c(as.list(UserSets), as.list(permutationSets)))
  cat(paste0("\n \t Running LOLA for target and ", length(permutationSets), " expression-matched permutation sets using ", opt$cores, " core(s) ... \n"))
  combinedResults <- as.data.frame(quiet_lola(LOLA::runLOLA(combinedSets, Universe, regionDB, cores = opt$cores)))
  target_user_set <- target_option_name
  if(!(target_user_set %in% as.character(combinedResults$userSet))){
    if("1" %in% as.character(combinedResults$userSet)){
      target_user_set <- "1"
    }else{
      stop(say("STOP: Could not identify target userSet in LOLA permutation results", by = "poop"), call. = FALSE)
    }
  }
  targetResults <- combinedResults[as.character(combinedResults$userSet) == target_user_set, , drop = FALSE]
  targetResults$file <- rep(target_file_label, nrow(targetResults))
  write_lola_results_table(targetResults, paste0(opt$out, background_title, opt$database, "_LOLA_results.txt"))
  permutationAwareResults <- calculate_permutation_results(
    combinedResults = combinedResults,
    target_user_set = target_user_set,
    metric = opt$significance_metric,
    threshold = opt$significance_threshold,
    empirical_threshold = opt$empirical_threshold
  )
  permutationAwareResults$file <- rep(target_file_label, nrow(permutationAwareResults))
  write_lola_results_table(permutationAwareResults, paste0(opt$out, background_title, opt$database, "_LOLA_permutation_results.txt"))
  permutationPlotResults <- permutationAwareResults
  permutationPlotResults$oddsRatio[!permutationPlotResults$passedPermutationFilter] <- 0.5
  cat(paste0("\n \t Permutation LOLA finished: ", length(permutationSets), " permutation sets processed \n"))
  regionResults <- list()
  regionResults[[target_option_name]] <- permutationPlotResults
  if(sum(permutationAwareResults$passedPermutationFilter, na.rm = TRUE) == 0){
    cat("\n \t WARNING: No enrichments passed target and permutation filters. Heatmap will show all tested marks/states at baseline. \n")
  }
}

##### 4. Plot results and export #####

cat(paste0("\n Plotting significant results with ", opt$significance_metric, " threshold ", signif(opt$significance_threshold, 4), " ... \n"))
## HEATMAP

if(opt$database == "Marks"){
  for(i in 1:length(regionResults)){
    name <- strsplit(opt[[names(regionResults[i])]], split = ".", fixed = T)[[1]][1]
    if(name == ""){ name <- strsplit(opt[[names(regionResults[i])]], split = ".", fixed = T)[[1]][2]}
    name <- strsplit(name, split = "/")[[1]][length(strsplit(name, split = "/")[[1]])]
    regionResults[[i]]$file <- rep(name, nrow(regionResults[[i]]))
    regionResults[[i]] <- apply_lola_significance_floor(regionResults[[i]], opt$significance_metric, opt$significance_threshold)
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
  cluster_cols <- ncol(df) > 1
  title <- strsplit(opt$background, split = ".", fixed = T)[[1]][1]
  if(title == ""){title <- strsplit(opt$background, split = ".", fixed = T)[[1]][2]}
  title <- strsplit(title, split = "/")[[1]][length(strsplit(title, split = "/")[[1]])]
  Annotation_colors <- list(Conditions = c(control = "darkgreen",Nitrogen = "cyan3", Sulphur = "gold3", light = "white", dark = "black"))
  pdf(file = paste0(opt$out,title,opt$database,".pdf"), paper = "a4r", height = 21, width = 28, onefile = T)
  pheatmap(df, scale = "none", color = col, breaks = make_heatmap_breaks(df, length(col)), cluster_rows = F, cluster_cols = cluster_cols, clustering_method = "ward.D2", cellwidth = 75, annotation_row = Conditions, annotation_names_row = FALSE, main = title, annotation_colors = Annotation_colors)
  invisible(dev.off())
}

if(opt$database == "MMarks"){
  for(i in 1:length(regionResults)){
    name <- strsplit(opt[[names(regionResults[i])]], split = ".", fixed = T)[[1]][1]
    if(name == ""){ name <- strsplit(opt[[names(regionResults[i])]], split = ".", fixed = T)[[1]][2]}
    name <- strsplit(name, split = "/")[[1]][length(strsplit(name, split = "/")[[1]])]
    regionResults[[i]]$file <- rep(name, nrow(regionResults[[i]]))
    regionResults[[i]] <- apply_lola_significance_floor(regionResults[[i]], opt$significance_metric, opt$significance_threshold)
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
  cluster_cols <- ncol(df) > 1
  title <- strsplit(opt$background, split = ".", fixed = T)[[1]][1]
  if(title == ""){title <- strsplit(opt$background, split = ".", fixed = T)[[1]][2]}
  title <- strsplit(title, split = "/")[[1]][length(strsplit(title, split = "/")[[1]])]
  Annotation_colors <- list(Conditions = c(active = "lightsteelblue1", repressive = "darksalmon", nucleosome = "gold3"))
  pdf(file = paste0(opt$out,title,opt$database,".pdf"), paper = "a4r", height = 21, width = 28, onefile = T)
  pheatmap(df, scale = "none", color = col, breaks = make_heatmap_breaks(df, length(col)), cluster_rows = F, cluster_cols = cluster_cols, clustering_method = "ward.D2", cellwidth = 75, annotation_row = Conditions, annotation_names_row = FALSE, main = title, annotation_colors = Annotation_colors)
  invisible(dev.off())
}

if(opt$database == "CS_Control" | opt$database == "CS_N" | opt$database == "CS_S"){
  for(i in 1:length(regionResults)){
    name <- strsplit(opt[[names(regionResults[i])]], split = ".", fixed = T)[[1]][1]
    if(name == ""){ name <- strsplit(opt[[names(regionResults[i])]], split = ".", fixed = T)[[1]][2]}
    name <- strsplit(name, split = "/")[[1]][length(strsplit(name, split = "/")[[1]])]
    regionResults[[i]]$file <- rep(name, nrow(regionResults[[i]]))
    regionResults[[i]] <- apply_lola_significance_floor(regionResults[[i]], opt$significance_metric, opt$significance_threshold)
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
  cluster_cols <- ncol(df2) > 1
  title <- strsplit(opt$background, split = ".", fixed = T)[[1]][1]
  if(title == ""){title <- strsplit(opt$background, split = ".", fixed = T)[[1]][2]}
  title <- strsplit(title, split = "/")[[1]][length(strsplit(title, split = "/")[[1]])]
  Annotation_colors <- list(functions = c(Transcribed = "lightsteelblue1", Repressed = "darksalmon", Bivalent = "gold3", Promoter = "lightsteelblue3", Heterochromatin = "black", 'no info' = "grey"),
                            locations = c('no info' = "grey", Intragenic = "darkseagreen", "3'gene" = "darkseagreen1", "5'gene" = "darkseagreen2"),
                            evolution = c(Conserved = "burlywood", Algal = "aquamarine"))
  pdf(file = paste0(opt$out,title,opt$database,".pdf"), paper = "a4r", height = 21, width = 28, onefile = T)
  pheatmap(df2, scale = "none", color = col, breaks = make_heatmap_breaks(df2, length(col)), cluster_rows = F, cluster_cols = cluster_cols, clustering_method = "ward.D2", cellwidth = 75, annotation_row = Conditions, annotation_names_row = FALSE, main = title, annotation_colors = Annotation_colors)
  invisible(dev.off())
}

if(opt$database == "CS_Chlamytina" ){
  for(i in 1:length(regionResults)){
    name <- strsplit(opt[[names(regionResults[i])]], split = ".", fixed = T)[[1]][1]
    if(name == ""){ name <- strsplit(opt[[names(regionResults[i])]], split = ".", fixed = T)[[1]][2]}
    name <- strsplit(name, split = "/")[[1]][length(strsplit(name, split = "/")[[1]])]
    regionResults[[i]]$file <- rep(name, nrow(regionResults[[i]]))
    regionResults[[i]] <- apply_lola_significance_floor(regionResults[[i]], opt$significance_metric, opt$significance_threshold)
  }
  
  regionResults <- do.call(rbind, regionResults)
  df <- data.frame(id = regionResults$file, pref.marks = regionResults$description, location = regionResults$cellType , name = rep(NA, nrow(regionResults)), value = regionResults$oddsRatio)
  for(i in 1:nrow(regionResults)){
    name <- strsplit(regionResults$filename[i], split = ".", fixed = T)[[1]][1]
    name <- gsub("^E", "CS", x = name)
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
  cluster_cols <- ncol(df2) > 1
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
  pheatmap(df2, scale = "none", color = col, breaks = make_heatmap_breaks(df2, length(col)), cluster_rows = F, cluster_cols = cluster_cols, clustering_method = "ward.D2", cellwidth = 50, annotation_row = Conditions, annotation_names_row = FALSE, main = title, annotation_colors = Annotation_colors, fontsize = 7, fontsize_row = 10, fontsize_col = 10)
  invisible(dev.off())
}

chlamytina_say(paste0("EnrichmentsLOLA has finished ! ", Sys.time()))




