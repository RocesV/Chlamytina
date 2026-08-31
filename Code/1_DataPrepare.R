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
  make_option(c("-A", "--file1"), type = "character", default = NULL, help = "Dataset1 file path. First column CreIDs. Other columns quantification data.", metavar = "character"),
  make_option(c("-a", "--condition1"), type = "character", default = NULL, help = "Dataset1 Condition vector. It representes replicates for each treatment, separated by - \n \t Condition vector must contain all your replicates. For example (9 samples): 1) 3-6 will set a contrast between the first three replicates and the last six \n \t 2) 3-3-3 will set all possible two by two contrasts between the three treatments \n", metavar = "character"),
  make_option(c("-B", "--file2"), type = "character", default = NULL, help = "Dataset2 file path", metavar = "character"),
  make_option(c("-b", "--condition2"), type = "character", default = NULL, help = "Dataset2 Condition vector", metavar = "character"),
  make_option(c("-C", "--file3"), type = "character", default = NULL, help = "Dataset3 file path", metavar = "character"),
  make_option(c("-q", "--condition3"), type = "character", default = NULL, help = "Dataset3 Condition vector", metavar = "character"),
  make_option(c("-D", "--file4"), type = "character", default = NULL, help = "Dataset4 file path", metavar = "character"),
  make_option(c("-x", "--condition4"), type = "character", default = NULL, help = "Dataset4 Condition vector", metavar = "character"),
  make_option(c("-E", "--file5"), type = "character", default = NULL, help = "Dataset5 file path", metavar = "character"),
  make_option(c("-y", "--condition5"), type = "character", default = NULL, help = "Dataset5 Condition vector", metavar = "character"),
  make_option(c("-d", "--differential"), type = "logical", default = TRUE, help = "If true, differential expression is performed [default = %default]"),
  make_option(c("-s", "--sva"), type = "logical", default = FALSE, help = "If true, sva removing unwanted variation is performed. Only for n>10-15 samples datasets. [default = %default]"),
  make_option(c("-i", "--intersect"), type = "logical", default = TRUE, help = "CreIDs intra-inter group specific discrimination [default = %default]", metavar = "character"),
  make_option(c("-o", "--out"), type = "character", default = "./Outputs/", help = "Output directory [default = %default]", metavar = "character"),
  make_option(c("-g", "--chromosome"), type = "logical", default = TRUE, help = "If true, non-chromosome mapped (scaffolds ...) proteins are not taked into account [default = %default]", metavar = "character"),
  make_option(c("-n", "--normalization"), type = "character", default = "none", help = "Normalization metric used (limma-only). Options: normalizeQuantiles, none \n \t It is advisable to set this argument as none and preprocess the data with other pkgs like Processomics [default = %default]", metavar = "character"),
  make_option(c("-m", "--de_method"), type = "character", default = "limma", help = "Differential expression method. Options: limma, DESeq2 [default = %default]", metavar = "character"),
  make_option(c("-f", "--diff_background"), type = "character", default = "differential", help = "Differential background mode. Options: differential, permutations. permutations is useful when only one contrast is available [default = %default]", metavar = "character"),
  make_option(c("-p", "--n_permutations"), type = "integer", default = 500, help = "Number of expression-matched permutation sets when --diff_background permutations [default = %default]", metavar = "integer"),
  make_option(c("-S", "--permutation_seed"), type = "integer", default = 123, help = "Random seed for permutation background generation [default = %default]", metavar = "integer"),
  make_option(c("-c", "--cores"), type = "integer", default = 1, help = "Number of cores for permutation background generation [default = %default]", metavar = "integer"),
  make_option(c("-V", "--version"), type = "character", default = "v5", help = "Genome annotation version. Options: v5, v6 [default = %default]", metavar = "character")
);

opt_parser = OptionParser(option_list = option_list, usage = "1_DataPrepare.R [file] [condition] [file] [condition] ... [options]");
opt = parse_args(opt_parser);
normalize_genome_version <- function(version) {
  version <- tolower(version)
  if(version %in% c("v5", "5", "v5.5", "5.5", "v5.6", "5.6")) return("v5")
  if(version %in% c("v6", "6", "v6.1", "6.1")) return("v6")
  stop(say("STOP: Unsupported genome version. Use v5 or v6", by = "poop"), call. = FALSE)
}
opt$version <- normalize_genome_version(opt$version)
normalize_de_method <- function(method) {
  method <- tolower(method)
  if(method == "limma") return("limma")
  if(method %in% c("deseq2", "deseq")) return("deseq2")
  stop(say("STOP: Unsupported DE method. Use limma or DESeq2", by = "poop"), call. = FALSE)
}
normalize_normalization <- function(normalization) {
  normalization <- tolower(normalization)
  if(normalization == "none") return("none")
  if(normalization == "normalizequantiles") return("normalizeQuantiles")
  stop(say("STOP: Unsupported normalization. Use none or normalizeQuantiles", by = "poop"), call. = FALSE)
}
normalize_diff_background <- function(mode) {
  mode <- tolower(mode)
  if(mode %in% c("differential", "default")) return("differential")
  if(mode %in% c("permutation", "permutations")) return("permutations")
  stop(say("STOP: Unsupported differential background mode. Use differential or permutations", by = "poop"), call. = FALSE)
}
opt$de_method <- normalize_de_method(opt$de_method)
opt$normalization <- normalize_normalization(opt$normalization)
opt$diff_background <- normalize_diff_background(opt$diff_background)
opt$n_permutations <- as.integer(opt$n_permutations)
opt$permutation_seed <- as.integer(opt$permutation_seed)
opt$cores <- as.integer(opt$cores)
de_alpha <- 0.05
if(opt$de_method == "deseq2" && opt$normalization != "none"){
  stop(say("STOP: normalizeQuantiles cannot be used with DESeq2. Use --normalization none", by = "poop"), call. = FALSE)
}
if(is.na(opt$n_permutations) || opt$n_permutations < 1){
  stop(say("STOP: --n_permutations must be an integer >= 1", by = "poop"), call. = FALSE)
}
if(is.na(opt$permutation_seed)){
  stop(say("STOP: --permutation_seed must be an integer", by = "poop"), call. = FALSE)
}
if(is.na(opt$cores) || opt$cores < 1){
  stop(say("STOP: --cores must be an integer >= 1", by = "poop"), call. = FALSE)
}
if(opt$diff_background == "permutations" && !opt$differential){
  stop(say("STOP: --diff_background permutations requires --differential TRUE", by = "poop"), call. = FALSE)
}
genome_config <- list(
  v5 = list(
    label = "v5.5/v5.6",
    gff = "./Data/DB/v5/Creinhardtii_281_v5.5.gene_exons.gff3.gz",
    universe = "./Data/DB/v5/Universe.bed",
    conversion_table = "./Data/DB/v5/ChlamydomonasTranscriptNameConversionBetweenReleases.Mch12b.txt",
    conversion_cols = c("5.5", "3.1", "Genbank", "4", "4.3", "u5", "u9", "5.3.1"),
    target_col = "5.5",
    convert_ids = TRUE,
    chromosomes = paste0("chr", 1:17)
  ),
  v6 = list(
    label = "v6.1",
    gff = "./Data/DB/v6/CreinhardtiiCC_4532_707_v6.1.gene_exons.gff3.gz",
    universe = "./Data/DB/v6/Universe.bed",
    conversion_table = NULL,
    conversion_cols = NULL,
    target_col = NULL,
    convert_ids = FALSE,
    chromosomes = sprintf("chromosome_%02d", 1:17)
  )
)
cfg <- genome_config[[opt$version]]
normalize_seqnames <- function(x, genome_version) {
  if(genome_version == "v5") {
    x <- sub("^chromosome_0?([0-9]+)$", "chr\\1", x)
  } else if(genome_version == "v6") {
    x <- sub("^plastome$", "cpDNA", x)
    x <- sub("^mitogenome$", "mtDNA", x)
    x <- sub("^mating_type_plus$", "MT_plus_R", x)
  }
  x
}

if (is.null(opt$file1)){
  print_help(opt_parser)
  stop(say(what = "STOP: At least one dataset is needed", by = "poop"), call.=FALSE)
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

chlamytina_say(paste0("Welcome to DataPrepare ! ", Sys.time()))
cat("\n Checking-installing-loading needed libs and packages ... \n")

mypkgs <- c("utils", "nVennR", "readxl", "R.utils", "BiocManager", "GenomicFeatures", "GenomicRanges", "GenomeInfoDb", "IRanges", "txdbmaker", "limma", "sva", "DESeq2", "SummarizedExperiment")
logicals <- is.element(mypkgs, installed.packages()[,1])
tmp <- base::sapply(mypkgs[logicals], FUN = function(x){  suppressPackageStartupMessages(library(x, character.only = TRUE))})
bioc_pkgs <- c("GenomicFeatures", "GenomicRanges", "GenomeInfoDb", "IRanges", "txdbmaker", "limma", "sva", "DESeq2", "SummarizedExperiment")
tmp <- base::sapply(mypkgs[!logicals], FUN = function(x){
  if(x %in% bioc_pkgs){
    BiocManager::install(x, update = FALSE, ask = FALSE)
  }else{
    install.packages(x, repos = "https://cloud.r-project.org/")
  }
})
tmp <- base::sapply(mypkgs[!logicals], FUN = function(x){ suppressPackageStartupMessages(library(x, character.only = TRUE))})
quiet_make_txdb_from_gff <- function(gff) {
  suppressWarnings(suppressMessages({
    invisible(utils::capture.output(txdb <- txdbmaker::makeTxDbFromGFF(file = gff)))
    txdb
  }))
}
strip_version_suffix <- function(x) {
  sub("\\.v[0-9]+(\\.[0-9]+)?$", "", as.character(x))
}
transcript_to_gene_id <- function(x) {
  x <- strip_version_suffix(x)
  x <- sub("\\.t[0-9]+(\\.[0-9]+)*$", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x
}
make_feature_bed <- function(gr, ids, feature_level) {
  ids <- strip_version_suffix(ids)
  keep <- !is.na(ids) & ids != ""
  data.frame(
    id = ids[keep],
    chr = as.character(GenomeInfoDb::seqnames(gr))[keep],
    start = IRanges::start(gr)[keep],
    end = IRanges::end(gr)[keep],
    feature_level = feature_level,
    stringsAsFactors = FALSE
  )
}
get_feature_bed <- local({
  feature_bed <- NULL
  function() {
    if(is.null(feature_bed)){
      txdb <- quiet_make_txdb_from_gff(cfg$gff)
      transcripts_gr <- GenomicFeatures::transcripts(txdb)
      transcript_ids <- strip_version_suffix(transcripts_gr$tx_name)
      transcript_bed <- make_feature_bed(transcripts_gr, transcript_ids, "transcript")
      genes_gr <- GenomicFeatures::genes(txdb)
      gene_ids <- names(genes_gr)
      if(is.null(gene_ids) || length(gene_ids) == 0){
        gene_ids <- genes_gr$gene_id
      }
      gene_ids <- transcript_to_gene_id(gene_ids)
      gene_bed <- make_feature_bed(genes_gr, gene_ids, "gene")
      feature_bed_tmp <- rbind(transcript_bed, gene_bed)
      feature_bed_tmp <- feature_bed_tmp[!duplicated(feature_bed_tmp$id), , drop = FALSE]
      feature_bed_tmp$chr <- normalize_seqnames(feature_bed_tmp$chr, opt$version)
      rownames(feature_bed_tmp) <- feature_bed_tmp$id
      feature_bed <<- feature_bed_tmp
    }
    feature_bed
  }
})
get_genome_ids <- local({
  ids <- NULL
  function() {
    if(is.null(ids)){
      ids <<- rownames(get_feature_bed())
    }
    ids
  }
})
get_bed_eligible_ids <- local({
  bed_ids <- NULL
  function() {
    if(is.null(bed_ids)){
      feature_bed <- get_feature_bed()
      if(opt$chromosome){
        feature_bed <- feature_bed[feature_bed$chr %in% cfg$chromosomes, , drop = FALSE]
      }
      bed_ids <<- rownames(feature_bed)
    }
    bed_ids
  }
})
get_v5_id_map <- local({
  id_map <- NULL
  function() {
    if(is.null(id_map)){
      feature_ids <- rownames(get_feature_bed())
      id_map <<- setNames(feature_ids, feature_ids)
      DB <- read.table(cfg$conversion_table, header = FALSE, sep = "\t", quote = "", comment.char = "#", fill = TRUE, stringsAsFactors = FALSE)
      colnames(DB) <- cfg$conversion_cols
      DB[DB == "--"] <- NA
      add_aliases <- function(source, target) {
        source <- strip_version_suffix(source)
        target <- strip_version_suffix(target)
        keep <- !is.na(source) & source != "" & !is.na(target) & target %in% feature_ids
        id_map <<- c(id_map, setNames(target[keep], source[keep]))
      }
      target_tx <- as.character(DB[[cfg$target_col]])
      for(col in cfg$conversion_cols){
        source_tx <- as.character(DB[[col]])
        add_aliases(source_tx, target_tx)
        add_aliases(transcript_to_gene_id(source_tx), transcript_to_gene_id(target_tx))
      }
      id_map <<- id_map[!duplicated(names(id_map))]
    }
    id_map
  }
})
clean_input_ids <- function(input, file_name) {
  ids <- as.character(input[,1])
  missing <- is.na(ids) | trimws(ids) == ""
  if(any(missing)){
    cat(paste0("\n \t WARNING: ", sum(missing), " rows with empty IDs were dropped from ", file_name, ". \n"))
    input <- input[!missing, , drop = FALSE]
    ids <- as.character(input[,1])
  }
  if(cfg$convert_ids){
    id_map <- get_v5_id_map()
    mapped <- unname(id_map[ids])
    unresolved <- is.na(mapped)
    mapped[unresolved] <- unname(id_map[strip_version_suffix(ids[unresolved])])
    keep <- !is.na(mapped)
    if(any(!keep)){
      cat(paste0("\n \t WARNING: ", sum(!keep), " IDs not found in the v5 conversion table or genome annotation were dropped from ", file_name, ". \n"))
    }
    input <- input[keep, , drop = FALSE]
    input[,1] <- mapped[keep]
  }else{
    genome_ids <- get_genome_ids()
    mapped <- ids
    unresolved <- !(mapped %in% genome_ids)
    mapped[unresolved] <- strip_version_suffix(mapped[unresolved])
    keep <- mapped %in% genome_ids
    if(any(!keep)){
      cat(paste0("\n \t WARNING: ", sum(!keep), " IDs not found in the ", cfg$label, " genome annotation were dropped from ", file_name, ". \n"))
    }
    input <- input[keep, , drop = FALSE]
    input[,1] <- mapped[keep]
  }
  feature_bed <- get_feature_bed()
  levels <- feature_bed$feature_level[match(as.character(input[,1]), rownames(feature_bed))]
  cat(paste0("\n \t Detected ", sum(levels == "gene", na.rm = TRUE), " gene-level IDs and ", sum(levels == "transcript", na.rm = TRUE), " transcript-level IDs in ", file_name, ". \n"))
  ids <- as.character(input[,1])
  duplicates <- duplicated(ids)
  if(any(duplicates)){
    cat(paste0("\n \t WARNING: ", sum(duplicates), " duplicated IDs were dropped from ", file_name, ". First occurrence was kept. \n"))
    input <- input[!duplicates, , drop = FALSE]
  }
  if(nrow(input) == 0){
    stop(say(paste0("STOP: No valid genome-mapped IDs remain in ", file_name), by = "poop"), call. = FALSE)
  }
  input
}

##### 2. Import tables and args:OK #####

cat("\n Importing tables and arguments ... \n")
conditions <- opt[grep("condition",names(opt))]
inputs <- opt[grep("file",names(opt))]
inputs <- base::lapply(inputs, FUN = function(x){
  if(length(grep(pattern = ".txt", x = x, fixed = TRUE)) == 1){
    input <- read.table(x, header = TRUE, sep = "", check.names = FALSE)
  }else if(length(grep(pattern = ".xls", x = x, fixed = TRUE)) == 1){
    input <- suppressMessages(read_excel(x, col_names = TRUE))
  }else if(length(grep(pattern = ".xlsx", x = x, fixed = TRUE)) == 1){
    input <- suppressMessages(read_excel(x, col_names = TRUE))
  } else{stop(say("STOP: Non supported format. Please try txt or excel tab separated with headers", by = "poop"))}
  input <- as.data.frame(input)
  input <- clean_input_ids(input, x)
  if(opt$differential && opt$de_method == "limma" && opt$normalization == "normalizeQuantiles"){
    input[,-1] <- limma::normalizeQuantiles(as.matrix(input[,-1]), ties = TRUE)
  }
  input
})

if(length(inputs) == 1 & opt$intersect & !opt$differential){
  cat("\n Because only one file is detected and intersect is defined as TRUE, differential is forced as TRUE \n")
  opt$differential <- TRUE
}

##### 2.5 Some helpers for DESeq2/permutations:OK #####

prepare_deseq_counts <- function(input, samples) {
  ids <- as.character(input[,1])
  if(any(is.na(ids) | ids == "")){
    stop(say("STOP: DESeq2 input contains empty IDs", by = "poop"), call. = FALSE)
  }
  if(anyDuplicated(ids)){
    stop(say("STOP: DESeq2 input contains duplicated IDs", by = "poop"), call. = FALSE)
  }
  counts <- as.matrix(input[, samples, drop = FALSE])
  suppressWarnings(storage.mode(counts) <- "numeric")
  if(any(is.na(counts))){
    stop(say("STOP: DESeq2 input contains NA/non-numeric values", by = "poop"), call. = FALSE)
  }
  if(any(counts < 0)){
    stop(say("STOP: DESeq2 input contains negative values", by = "poop"), call. = FALSE)
  }
  non_integer <- abs(counts - round(counts)) > .Machine$double.eps^0.5
  if(any(non_integer)){
    cat(paste0("\n \t WARNING: DESeq2 input contains ", sum(non_integer), " non-integer count-like values. Values are rounded before DESeq2 analysis. \n"))
  }
  counts <- round(counts)
  storage.mode(counts) <- "integer"
  rownames(counts) <- ids
  colnames(counts) <- colnames(input)[samples]
  counts
}

parallel_lapply <- function(X, FUN, cores) {
  if(cores > 1 && .Platform$OS.type != "windows"){
    return(parallel::mclapply(X, FUN, mc.cores = cores))
  }
  if(cores > 1 && .Platform$OS.type == "windows"){
    cat("\n \t WARNING: permutation generation uses serial mode on Windows \n")
  }
  lapply(X, FUN)
}

count_differential_contrasts <- function(inputs, conditions) {
  if(length(inputs) != length(conditions)){
    stop(say("STOP: Same number of condition vectors and file inputs is required", by = "poop"), call. = FALSE)
  }
  sum(vapply(seq_along(inputs), function(l){
    groups <- as.numeric(strsplit(conditions[[l]], split = "-")[[1]])
    if(any(is.na(groups)) || length(groups) < 2){
      stop(say("STOP: Each differential condition vector must define at least two treatments", by = "poop"), call. = FALSE)
    }
    choose(length(groups), 2)
  }, numeric(1)))
}

make_expression_matched_permutations <- function(input, samples, target_ids, n_permutations, cores, seed, bins = 20, eligible_ids = NULL) {
  ids <- as.character(input[,1])
  abundance_matrix <- as.matrix(input[, samples, drop = FALSE])
  suppressWarnings(storage.mode(abundance_matrix) <- "numeric")
  abundance <- rowMeans(abundance_matrix, na.rm = TRUE)
  names(abundance) <- ids
  usable <- !is.na(ids) & ids != "" & is.finite(abundance)
  if(sum(!usable) > 0){
    cat(paste0("\n \t WARNING: ", sum(!usable), " genes/proteins with non-finite abundance ignored for permutation matching \n"))
  }
  ids <- ids[usable]
  abundance <- abundance[ids]
  if(!is.null(eligible_ids)){
    eligible_ids <- unique(as.character(eligible_ids))
    usable_target_ids <- intersect(unique(as.character(target_ids)), ids)
    lost_target_ids <- setdiff(usable_target_ids, eligible_ids)
    if(length(lost_target_ids) > 0){
      cat(paste0("\n \t WARNING: ", length(lost_target_ids), " target DE IDs are outside the BED-exportable chromosome set and will be excluded from permutation matching. \n"))
    }
    ids <- intersect(ids, eligible_ids)
    abundance <- abundance[ids]
  }
  target_ids <- intersect(unique(as.character(target_ids)), ids)
  candidate_ids <- setdiff(ids, target_ids)
  if(length(target_ids) == 0){
    stop(say("STOP: permutation background cannot be generated because the target DE set is empty", by = "poop"), call. = FALSE)
  }
  if(length(candidate_ids) < length(target_ids)){
    stop(say("STOP: Not enough non-target genes/proteins for permutation background. Use Universe or --diff_background differential", by = "poop"), call. = FALSE)
  }
  n_bins <- min(bins, max(1, length(unique(abundance))))
  abundance_rank <- rank(abundance, ties.method = "first")
  abundance_bin <- ceiling(abundance_rank / max(abundance_rank) * n_bins)
  abundance_bin[abundance_bin < 1] <- 1
  abundance_bin[abundance_bin > n_bins] <- n_bins
  names(abundance_bin) <- ids
  target_bin_counts <- table(abundance_bin[target_ids])
  insufficient_bins <- names(target_bin_counts)[vapply(names(target_bin_counts), function(bin_id){
    sum(abundance_bin[candidate_ids] == as.integer(bin_id)) < as.integer(target_bin_counts[[bin_id]])
  }, logical(1))]
  if(length(insufficient_bins) > 0){
    cat("\n \t WARNING: Some expression bins do not have enough non-target genes/proteins. Neighboring expression bins will be used; consider Universe or --diff_background differential if results look unstable. \n")
  }
  sample_one <- function(i) {
    set.seed(seed + i)
    selected <- character(0)
    for(bin_id in names(target_bin_counts)){
      needed <- as.integer(target_bin_counts[[bin_id]])
      bin_num <- as.integer(bin_id)
      candidate_pool <- character(0)
      for(radius in 0:n_bins){
        candidate_pool <- candidate_ids[abs(abundance_bin[candidate_ids] - bin_num) <= radius]
        candidate_pool <- setdiff(candidate_pool, selected)
        if(length(candidate_pool) >= needed) break
      }
      if(length(candidate_pool) < needed){
        stop("Not enough genes/proteins to sample a matched permutation set")
      }
      selected <- c(selected, sample(candidate_pool, needed, replace = FALSE))
    }
    selected
  }
  permutations <- parallel_lapply(seq_len(n_permutations), sample_one, cores)
  names(permutations) <- sprintf("Permutation_%03d", seq_along(permutations))
  permutations
}
n_diff_contrasts <- if(opt$differential) count_differential_contrasts(inputs, conditions) else 0
if(opt$diff_background == "permutations" && n_diff_contrasts != 1){
  stop(say("STOP: --diff_background permutations is only supported when exactly one differential contrast is available", by = "poop"), call. = FALSE)
}
get_nvenn_cycles <- function(n_sets) {
  if(n_sets > 5) return(1000)
  7000
}

##### 3. Diff expression / sva:OK #####

if(opt$differential){
  cat(paste0("\n Performing differential expression with ", opt$de_method, " ... \n"))
  Differential <- list()
  Permutation_background <- list()
  for(l in 1:length(inputs)){
    # format, args and qc
    if(length(inputs) != length(conditions)){ stop(say("STOP: Same number of condition vectors and file inputs is required", by = "poop"))}
    groups <- as.numeric(strsplit(conditions[[l]], split = "-")[[1]])
    if(sum(groups) != ncol(inputs[[l]][,-1])){stop(say(paste0("STOP: Your condition vector dont match the number of replicates in", names(inputs[l])),by = "poop"))}
    comb.contrast <- combn(1:length(groups), 2)
    cat(paste0("\n" ,names(inputs[l]), sep = ": ", length(groups), " treatments and ", ncol(comb.contrast), " two-by-two contrasts \n"))
    groups.list <- list()
    for(z in 1:length(groups)){
      if(z == 1){ groups.list[[z]] <- 1:groups[z]
      } else if(z != 1){groups.list[[z]] <- (max(groups.list[[z - 1]]) + 1):(max(groups.list[[z - 1]] + groups[z]))}
    }
    for(i in 1:ncol(comb.contrast)){
      treatment1 <- paste0("Treatment", comb.contrast[1,i])
      treatment2 <- paste0("Treatment", comb.contrast[2,i])
      contrastss <- paste0(treatment2, "-", treatment1)
      samples <- c(groups.list[[comb.contrast[1,i]]], groups.list[[comb.contrast[2,i]]]) + 1
      Treatment <- factor(
        c(
          rep(treatment1, groups[comb.contrast[1,i]]),
          rep(treatment2, groups[comb.contrast[2,i]])
        ),
        levels = c(treatment1, treatment2)
      )
      PhenoData <- data.frame(samples = colnames(inputs[[l]][,samples]), Treatment = Treatment)
      use_sva <- opt$sva
      if(length(samples) < 12){
        use_sva <- FALSE
        cat("\n \t sva turned OFF because contrast nsamples < 12 \n")
      }
      if(opt$de_method == "limma"){
        mod <- model.matrix(~0+Treatment, data=PhenoData)
        colnames(mod) <- levels(PhenoData$Treatment)
        mod0 <- model.matrix(~1,data=PhenoData)
        colnames(mod0) <- "samples"
        if(use_sva){
          cat("\n \t sva: ON  \n")
          n.sv <- sva::num.sv(as.matrix(inputs[[l]][,samples]), mod, method="leek")
          max_sv <- max(0, nrow(mod) - qr(mod)$rank - 1)
          n.sv <- min(n.sv, max_sv)
          cat(paste0("\n \t sva: ", n.sv, " unknown batch effects founded \n"))
          if(n.sv > 0){
            svobj <- sva::sva(as.matrix(inputs[[l]][,samples]), mod, mod0, n.sv=n.sv)
            colnames(svobj$sv) <- paste(rep("col",ncol(svobj$sv)),c(1:ncol(svobj$sv)),sep="")
            mod <- cbind(mod, svobj$sv)
          }
        }else{
          cat("\n \t sva: OFF \n")
        }
        fit <- limma::lmFit(inputs[[l]][,samples], mod)
        contrast.matrix <- limma::makeContrasts(contrastss, levels=mod)
        fit2 <- limma::contrasts.fit(fit,contrast.matrix)
        fit2 <- limma::eBayes(fit2)
        DE_filtered <- limma::topTable(fit2, number = Inf, coef = 1, p.value = de_alpha, sort.by = "logFC")
        DE_filtered <- cbind(Accession = inputs[[l]][rownames(DE_filtered),1], DE_filtered)
      }else if(opt$de_method == "deseq2"){
        count_data <- prepare_deseq_counts(inputs[[l]], samples)
        coldata <- data.frame(row.names = colnames(count_data), Treatment = Treatment)
        dds <- DESeq2::DESeqDataSetFromMatrix(
          countData = count_data,
          colData = coldata,
          design = ~ Treatment
        )
        dds <- dds[rowSums(DESeq2::counts(dds)) > 0,]
        if(nrow(dds) == 0){
          stop(say("STOP: DESeq2 has no genes/proteins with counts > 0", by = "poop"), call. = FALSE)
        }
        if(use_sva){
          cat("\n \t svaseq: ON  \n")
          dds <- DESeq2::estimateSizeFactors(dds)
          dat <- DESeq2::counts(dds, normalized = TRUE)
          dat <- dat[rowMeans(dat) > 1, , drop = FALSE]
          coldata_df <- as.data.frame(SummarizedExperiment::colData(dds))
          mod <- model.matrix(~ Treatment, data = coldata_df)
          mod0 <- model.matrix(~ 1, data = coldata_df)
          max_sv <- max(0, nrow(mod) - qr(mod)$rank - 1)
          n.sv <- if(nrow(dat) > 1 && max_sv > 0) sva::num.sv(dat, mod, method="leek") else 0
          n.sv <- min(n.sv, max_sv)
          cat(paste0("\n \t svaseq: ", n.sv, " unknown batch effects founded \n"))
          if(n.sv > 0){
            svseq <- sva::svaseq(dat, mod, mod0, n.sv = n.sv)
            for(k in 1:n.sv){
              dds[[paste0("SV", k)]] <- svseq$sv[,k]
            }
            design(dds) <- as.formula(paste("~", paste(c(paste0("SV", 1:n.sv), "Treatment"), collapse = " + ")))
          }
        }else{
          cat("\n \t svaseq: OFF \n")
        }
        dds <- DESeq2::DESeq(dds, quiet = TRUE)
        res <- DESeq2::results(dds, contrast = c("Treatment", treatment2, treatment1), alpha = de_alpha)
        res <- as.data.frame(res)
        res <- res[order(res$log2FoldChange, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
        DE_filtered <- res[!is.na(res$padj) & res$padj < de_alpha, , drop = FALSE]
        DE_filtered <- cbind(Accession = rownames(DE_filtered), DE_filtered)
      }
      file.name <- paste0(names(inputs[l]), contrastss)
      if(opt$diff_background == "permutations"){
        cat(paste0("\n \t Generating ", opt$n_permutations, " expression-matched permutation sets for ", contrastss, ". This can take a while ... \n"))
        Permutation_background[[file.name]] <- make_expression_matched_permutations(
          input = inputs[[l]],
          samples = samples,
          target_ids = DE_filtered[,1],
          n_permutations = opt$n_permutations,
          cores = opt$cores,
          seed = opt$permutation_seed,
          eligible_ids = get_bed_eligible_ids()
        )
        cat(paste0("\n \t Permutation background for ", contrastss, " generated \n"))
      }
      cat(paste0("\n \t ", contrastss, " | ", nrow(inputs[[l]]), " total proteins/genes , ", nrow(DE_filtered), " differential proteins/genes (FDR ", de_alpha, ") \n"))
      write.table(DE_filtered, file = paste0(opt$out, file.name, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = F)
      Differential[[file.name]] <- DE_filtered[,1]
    }
  }}

##### 4. Version ID liftover:OK #####

inputs <- lapply(inputs, FUN = function(x){ x[,1] })
if(opt$differential){
  VictorgoestoBED <- list(inputs = inputs, differential = Differential)
  if(opt$diff_background == "permutations"){
    VictorgoestoBED$permutations <- Permutation_background
  }
}else if(!opt$differential){
  VictorgoestoBED <- list(inputs = inputs)
}
map_nested_id_list <- function(x, FUN) {
  if(is.list(x)){
    return(lapply(x, function(y){ map_nested_id_list(y, FUN) }))
  }
  FUN(x)
}
valid_feature_ids <- rownames(get_feature_bed())
if(cfg$convert_ids){
  cat(paste0("\n Input IDs were normalized to ", cfg$label, " during import. No additional CreID conversion is performed. \n"))
}else{
  cat(paste0("\n Using ", cfg$label, " IDs directly. No CreID conversion is performed. \n"))
}
VictorgoestoBED <- map_nested_id_list(VictorgoestoBED, function(x){
  x <- strip_version_suffix(as.character(x[!is.na(x)]))
  x[x %in% valid_feature_ids]
})

##### 5. Background selection:OK #####

cat(paste0("\n Universe background for ", cfg$label, " can be found at ", cfg$universe, " \n"))
cat("\n File background may be a good reference for your enrichments! \n")
if(length(inputs) > 1 && opt$intersect){
  Global_background <- unique(unlist(VictorgoestoBED$inputs, use.names = FALSE))
  VictorgoestoBED$Global_background <- Global_background
  cat("\n For the args defined, Global background may be a good reference for your enrichments! \n")
}
if(opt$differential && opt$intersect){
  Diff_background <- unique(unlist(VictorgoestoBED$differential, use.names = FALSE))
  VictorgoestoBED$Diff_background <- Diff_background
  if(n_diff_contrasts == 1){
    cat("\n \t WARNING: Only one differential contrast was detected. Diff_background is identical to the query set and is usually not recommended as LOLA background. Use Universe, file background, or --diff_background permutations. \n")
  }else{
    cat("\n For the args defined, Differential background may be interesting for your enrichments! \n")
  }
}
if(opt$differential && n_diff_contrasts == 1 && !opt$intersect){
  cat("\n \t Only one differential contrast was detected. Intersections cannot be computed with one contrast; consider Universe, file background, or --diff_background permutations for LOLA. \n")
}
if(opt$diff_background == "permutations"){
  cat("\n \t Expression-matched permutation BEDs will be exported under the permutations folder and should be passed to 2_EnrichmentsLOLA.R with --permutation_path. \n")
}

##### 6. Intersection:OK #####

if(length(inputs) > 1){ 
  if(opt$intersect){
    cat("\n Intersect is defined as TRUE \n")
    myV <- nVennR::plotVenn(VictorgoestoBED$inputs, nCycles = get_nvenn_cycles(length(VictorgoestoBED$inputs)), opacity = 0.2, borderWidth = 3, systemShow = F, fontScale = 2, outFile = paste0(opt$out, "intersection.svg"))
    myV2 <- nVennR::listVennRegions(myV)
    nonames <- names(myV2)
    nonames <- sapply(nonames, USE.NAMES = F,FUN = function(x){
      strsplit(x, split = "(", fixed = T)[[1]][2]
    })
    names(myV2) <- paste0("Uniq_", nonames)
    names(myV2) <- gsub(")", "", names(myV2), fixed = T)
    names(myV2) <- gsub(", ", "", names(myV2), fixed = T)
    VictorgoestoBED$intersect_files <- myV2
    if(opt$differential){
      cat("\n Differential is defined as TRUE \n")
      if(length(VictorgoestoBED$differential) > 10){
        cat(paste0("\n \t WARNING: Differential nVennR skipped because ", length(VictorgoestoBED$differential), " DE contrast sets were detected (>10). diff_intersection.svg and Diff_uniq_*.bed files will not be generated. \n"))
        cat("\n \t Useful LOLA inputs are still computed: individual contrast BEDs, file*.bed inputs, Global_background.bed, and Diff_background.bed. To obtain differential intersection BEDs, reduce the number of contrasts or redesign the input groups. \n")
      }else{
        myV <- nVennR::plotVenn(VictorgoestoBED$differential, nCycles = get_nvenn_cycles(length(VictorgoestoBED$differential)), opacity = 0.2, borderWidth = 3, systemShow = F, fontScale = 2, outFile = paste0(opt$out, "diff_intersection.svg"))
        myV2 <- nVennR::listVennRegions(myV)
        nonames <- names(myV2)
        nonames <- sapply(nonames, USE.NAMES = F,FUN = function(x){
          strsplit(x, split = "(", fixed = T)[[1]][2]
        })
        names(myV2) <- paste0("Diff_uniq_", nonames)
        names(myV2) <- gsub(")", "", names(myV2), fixed = T)
        names(myV2) <- gsub(", ", "", names(myV2), fixed = T)
        VictorgoestoBED$intersect_diff <- myV2
      }
    }
  } else{
    cat("\n Intersect is defined as FALSE \n")
  }
}else if(length(inputs) == 1){
  if(!opt$intersect){
    cat("\n Intersect is defined as FALSE \n")
  }else if(!opt$differential){
    cat("\n Intersect is defined as TRUE, but only one input set is available and differential is FALSE. Intersections cannot be computed. \n")
  }else if(n_diff_contrasts == 1){
    cat("\n Intersect is defined as TRUE, but only one differential contrast is available. Differential intersections cannot be computed. \n")
  }else if(n_diff_contrasts > 1){
    cat("\n Intersect is defined as TRUE \n")
    cat("\n Differential is defined as TRUE \n")
    if(length(VictorgoestoBED$differential) > 10){
      cat(paste0("\n \t WARNING: Differential nVennR skipped because ", length(VictorgoestoBED$differential), " DE contrast sets were detected (>10). diff_intersection.svg and Diff_uniq_*.bed files will not be generated. \n"))
      cat("\n \t Useful LOLA inputs are still computed: individual contrast BEDs, file*.bed inputs, and Diff_background.bed. To obtain differential intersection BEDs, reduce the number of contrasts or redesign the input groups. \n")
    }else{
      myV <- nVennR::plotVenn(VictorgoestoBED$differential, nCycles = get_nvenn_cycles(length(VictorgoestoBED$differential)), opacity = 0.2, borderWidth = 3, systemShow = F, fontScale = 2, outFile = paste0(opt$out, "diff_intersection.svg"))
      myV2 <- nVennR::listVennRegions(myV)
      nonames <- names(myV2)
      nonames <- sapply(nonames, USE.NAMES = F,FUN = function(x){
        strsplit(x, split = "(", fixed = T)[[1]][2]
      })
      names(myV2) <- paste0("Diff_uniq_", nonames)
      names(myV2) <- gsub(")", "", names(myV2), fixed = T)
      names(myV2) <- gsub(", ", "", names(myV2), fixed = T)
      VictorgoestoBED$intersect_diff <- myV2
    }
  }
}

##### 7. From ID table to bed #####

if(!file.exists(cfg$gff)){
  stop(say(paste0("STOP: GFF file not found: ", cfg$gff), by = "poop"), call. = FALSE)
}
Feature.bed <- get_feature_bed()
cat("\n Computing CreIDs coherence ... \n")
QC <- lapply(VictorgoestoBED$inputs, function(x){
  x %in% rownames(Feature.bed)
})
if(length(which(unlist(QC, use.names = FALSE) == FALSE)) > 0){
  stop(say("STOP: Something go wrong with CreIDs. Check", by = "poop"))
}
if(length(which(unlist(QC, use.names = FALSE) == FALSE)) == 0){
  cat("\n Passed \n")
}
cat("\n Converting to .bed ... \n")
AllIDs <- unique(unlist(VictorgoestoBED$inputs, use.names = FALSE))
if(length(AllIDs) == 0){
  stop(say("STOP: No genome-mapped IDs remain for BED conversion after ID normalization", by = "poop"), call. = FALSE)
}
All.bed <- Feature.bed[AllIDs, c("chr", "start", "end"), drop = FALSE]
if(opt$chromosome){
  cat("\n Chromosome filtering is defined as TRUE \n")
  All.bed <- All.bed[All.bed[,1] %in% cfg$chromosomes, , drop = FALSE]
  if(nrow(All.bed) == 0){
    stop(say("STOP: No BED intervals remain after chromosome filtering. Check genome version or use --chromosome FALSE", by = "poop"), call. = FALSE)
  }
}
colnames(All.bed) <- c("chr", "start", "end")
if(opt$chromosome){
  All.bed$chr <- factor(All.bed$chr, cfg$chromosomes, ordered = TRUE)
  All.bed <- All.bed[order(All.bed$chr, All.bed$start),]
  All.bed$chr <- as.character(All.bed$chr)
}else{
  All.bed <- All.bed[order(All.bed$chr, All.bed$start),]
}

##### 8. Export beds #####

Permutation_beds <- NULL
if("permutations" %in% names(VictorgoestoBED)){
  Permutation_beds <- VictorgoestoBED$permutations
}
for(x in 1:length(VictorgoestoBED)){
  if(names(VictorgoestoBED)[x] == "permutations") next
  
  if(class(VictorgoestoBED[[x]]) == "list"){
    for(y in 1:length(VictorgoestoBED[[x]])){
      VictorgoestoBED[[x]][[y]] <- unique(VictorgoestoBED[[x]][[y]])
      VictorgoestoBED[[x]][[y]] <- VictorgoestoBED[[x]][[y]][VictorgoestoBED[[x]][[y]] %in% row.names(All.bed)]
      VictorgoestoBED[[x]][[y]] <- VictorgoestoBED[[x]][[y]][order(VictorgoestoBED[[x]][[y]])]
      VictorgoestoBED[[x]][[y]] <- All.bed[VictorgoestoBED[[x]][[y]],]
      write.table(VictorgoestoBED[[x]][[y]], file = paste0(opt$out,names(VictorgoestoBED[[x]][y]), ".bed"), row.names = FALSE, col.names = F, sep = "\t", quote = FALSE)
    }
  }else if(class(VictorgoestoBED[[x]]) != "list"){
    VictorgoestoBED[[x]] <- unique(VictorgoestoBED[[x]])
    VictorgoestoBED[[x]] <- VictorgoestoBED[[x]][VictorgoestoBED[[x]] %in% row.names(All.bed)]
    VictorgoestoBED[[x]] <- VictorgoestoBED[[x]][order(VictorgoestoBED[[x]])]
    VictorgoestoBED[[x]] <- All.bed[VictorgoestoBED[[x]],]
    write.table(VictorgoestoBED[[x]], file = paste0(opt$out,names(VictorgoestoBED[x]), ".bed"), row.names = FALSE, col.names = F, sep = "\t", quote = FALSE)
  }
}
if(!is.null(Permutation_beds)){
  for(contrast in names(Permutation_beds)){
    permutation_dir <- file.path(opt$out, "permutations", contrast)
    dir.create(permutation_dir, recursive = TRUE, showWarnings = FALSE)
    permutation_index <- data.frame(permutation = names(Permutation_beds[[contrast]]), n_regions = NA_integer_)
    for(y in seq_along(Permutation_beds[[contrast]])){
      original_n <- length(unique(Permutation_beds[[contrast]][[y]]))
      ids <- unique(Permutation_beds[[contrast]][[y]])
      ids <- ids[ids %in% row.names(All.bed)]
      ids <- ids[order(ids)]
      perm_bed <- All.bed[ids,]
      permutation_index$n_regions[y] <- nrow(perm_bed)
      if(nrow(perm_bed) < original_n){
        cat(paste0("\n \t WARNING: ", names(Permutation_beds[[contrast]])[y], " lost IDs during BED conversion/filtering \n"))
      }
      write.table(perm_bed, file = file.path(permutation_dir, paste0(names(Permutation_beds[[contrast]])[y], ".bed")), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
    }
    write.table(permutation_index, file = file.path(permutation_dir, "permutation_index.txt"), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    cat(paste0("\n \t Permutation BEDs written to ", permutation_dir, " \n"))
  }
}

chlamytina_say(paste0("DataPrepare has finished ! ", Sys.time()))

