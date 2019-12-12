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
  make_option(c("-A", "--file1"), type = "character", default = NULL, help = "Dataset1 file path. First column CreIDs. Other columns quantification data.", metavar = "character"),
  make_option(c("--condition1"), type = "character", default = NULL, help = "Dataset1 Condition vector. It representes replicates for each treatment, separated by - \n \t Condition vector must contain all your replicates. For example (9 samples): 1) 3-6 will set a contrast between the first three replicates and the last six \n \t 2) 3-3-3 will set all possible two by two contrasts between the three treatments \n", metavar = "character"),
  make_option(c("-B", "--file2"), type = "character", default = NULL, help = "Dataset2 file path", metavar = "character"),
  make_option(c("--condition2"), type = "character", default = NULL, help = "Dataset2 Condition vector", metavar = "character"),
  make_option(c("-C", "--file3"), type = "character", default = NULL, help = "Dataset3 file path", metavar = "character"),
  make_option(c("--condition3"), type = "character", default = NULL, help = "Dataset3 Condition vector", metavar = "character"),
  make_option(c("-D", "--file4"), type = "character", default = NULL, help = "Dataset4 file path", metavar = "character"),
  make_option(c("--condition4"), type = "character", default = NULL, help = "Dataset4 Condition vector", metavar = "character"),
  make_option(c("-E", "--file5"), type = "character", default = NULL, help = "Dataset5 file path", metavar = "character"),
  make_option(c("--condition5"), type = "character", default = NULL, help = "Dataset5 Condition vector", metavar = "character"),
  make_option(c("-d", "--differential"), type = "logical", default = TRUE, help = "If true, differential expression limma based test is performed [default = %default]"),
  make_option(c("-s", "--sva"), type = "logical", default = FALSE, help = "If true, sva removing unwanted variation is performed. Only for n>10-15 samples datasets. [default = %default]"),
  make_option(c("-i", "--intersect"), type = "logical", default = TRUE, help = "CreIDs intra-inter group specific discrimination [default = %default]", metavar = "character"),
  make_option(c("-o", "--out"), type = "character", default = "./Outputs/", help = "Output directory [default = %default]", metavar = "character"),
  make_option(c("-c", "--chromosome"), type = "logical", default = TRUE, help = "If true, non-chromosome mapped (scaffolds ...) proteins are not taked into account [default = %default]", metavar = "character"),
  make_option(c("-n", "--normalization"), type = "character", default = "none", help = "Normalization metric used. Options: normalizeQuantiles (limma), none \n \t It is advisable to set this argument as none and preprocess the data with other pkgs like Processomics [default = %default]", metavar = "character")
);

opt_parser = OptionParser(option_list = option_list, usage = "1_DataPrepare.R [file] [condition] [file] [condition] ... [options]");
opt = parse_args(opt_parser);

if (is.null(opt$file1)){
  print_help(opt_parser)
  stop(say(what = "STOP: At least one dataset is needed", by = "poop"), call.=FALSE)
}

##### 1. Secondary pkgs: checks and install #####

say(paste0("Welcome to DataPrepare ! ", Sys.time()), by = "rabbit", what_color = "white", by_color = "yellow")
cat("\n Checking-installing-loading needed libs and packages ... \n")

mypkgs <- c("utils", "nVennR", "readxl", "R.utils","BiocManager","GenomicFeatures", "limma", "sva")
logicals <- is.element(mypkgs, installed.packages()[,1])
tmp <- base::sapply(mypkgs[logicals], FUN = function(x){  suppressPackageStartupMessages(library(x, character.only = TRUE))})
tmp <- base::sapply(mypkgs[!logicals], FUN = function(x){
  if(x != "GenomicFeatures" & x != "limma" & x != "sva")
  install.packages(x, repos = "https://cloud.r-project.org/")
  if(x == "GenomicFeatures"){
    BiocManager::install("GenomicFeatures")
  }
  if(x == "limma"){
    BiocManager::install("limma")
  }
  if(x == "sva"){
    BiocManager::install("sva")
  }
})
tmp <- base::sapply(mypkgs[!logicals], FUN = function(x){ suppressPackageStartupMessages(library(x, character.only = TRUE))})

##### 2. Import tables and args:OK #####

cat("\n Importing tables and arguments ... \n")
conditions <- opt[grep("condition",names(opt))]
inputs <- opt[grep("file",names(opt))]
inputs <- base::lapply(inputs, FUN = function(x){
  if(length(grep(pattern = ".txt", x = x, fixed = TRUE)) == 1){
    input <- read.table(x, header = TRUE, sep = " ")
  }else if(length(grep(pattern = ".xls", x = x, fixed = TRUE)) == 1){
    input <- read_excel(x, col_names = TRUE)
  }else if(length(grep(pattern = ".xlsx", x = x, fixed = TRUE)) == 1){
    input <- read_excel(x, col_names = TRUE)
  } else{stop(say("STOP: Non supported format. Please try txt or excel tab separated with headers", by = "poop"))}
  input <- as.data.frame(input)
  if(opt$differential){ 
    if(opt$normalization == "none"){ input
    }else if(opt$normalization == "normalizeQuantiles"){
        input[,-1] <- normalizeQuantiles(as.matrix(input[,-1]), ties = TRUE)
        input
    }
  }else if(!opt$differential){ input}
})

if(length(inputs) == 1 & opt$intersect & !opt$differential){
  cat("\n Because only one file is detected and intersect is defined as TRUE, differential is forced as TRUE \n")
  opt$differential <- TRUE
}

##### 3. Diff expression / sva:OK #####

if(opt$differential){
cat("\n Performing differential expression ... \n")
Differential <- list()
for(l in 1:length(inputs)){
  
  # format, args and qc
  if(length(inputs) != length(conditions)){ stop(say("STOP: Same number of condition vectors and file inputs is required", by = "poop"))}
  groups <- as.numeric(strsplit(conditions[[l]], split = "-")[[1]])
  if(sum(groups) != ncol(inputs[[l]][,-1])){stop(say(paste0("STOP: Your condition vector dont match the number of replicates in", names(inputs[l])),by = "poop"))}
  comb.contrast <- combn(1:length(groups), 2)
  cat(paste0("\n" ,names(inputs[l]), sep = ": ", length(groups), " treatments and ", ncol(comb.contrast), " two-by-two contrasts \n"))
  if(ncol(inputs[[l]][,-1]) < 12){
    opt$sva <- FALSE
    cat("\n \t sva turned OFF because nsamples < 12 \n")
  }
  groups.list <- list()
  for(z in 1:length(groups)){
    if(z == 1){ groups.list[[z]] <- 1:groups[z]
    } else if(z != 1){groups.list[[z]] <- (max(groups.list[[z - 1]]) + 1):(max(groups.list[[z - 1]] + groups[z]))}}
  
  # differential
  if(opt$sva){
    # sva + limma: cat n.svas
    cat("\n \t sva: ON  \n")
    for(i in 1:ncol(comb.contrast)){
      # Build null-model and model
      samples <- c(groups.list[[comb.contrast[1,i]]], groups.list[[comb.contrast[2,i]]]) + 1
      PhenoData <- data.frame(samples = colnames(inputs[[l]][,samples]), Treatment = c(rep(paste0("Treatment", comb.contrast[1,i]),groups[comb.contrast[1,i]]), c(rep(paste0("Treatment", comb.contrast[2,i]),groups[comb.contrast[2,i]]))))
      mod = model.matrix(~0+as.factor(Treatment), data=PhenoData)
      colnames(mod) <- c(levels(PhenoData$Treatment))
      mod0 = model.matrix(~1,data=PhenoData)
      colnames(mod0) <- "samples"
      # Side_note01: In the future remove infinte and zero values if svas are giving problems
      # sva
      n.sv = num.sv(as.matrix(inputs[[l]][,samples]),mod,method="leek")
      cat(paste0("\n \t sva: ", n.sv, " unknown batch effects founded \n"))
      if(n.sv == 0){
        # Fit model
        fit <-lmFit(inputs[[l]][,samples],mod)
        # Select contrast
        contrastss <- c(paste0(colnames(mod)[2], "-", colnames(mod)[1]))
        contrast.matrix <- makeContrasts(contrastss, levels=mod)
        fit2 <- contrasts.fit(fit,contrast.matrix)
        fit2 <- eBayes(fit2)
        # Filter and export
        DE_filtered <- limma::topTable(fit2, number = Inf, coef = 1, p.value = 0.05, sort.by = "logFC")
        cat(paste0("\n \t ", contrastss, sep = " | ", nrow(inputs[[l]]), " total proteins , ", nrow(DE_filtered), " differential proteins (0.05) \n"))
        DE_filtered <- cbind(Accession = inputs[[l]][rownames(DE_filtered),1], DE_filtered)
        file.name <- paste0(names(inputs[l]), contrastss)
        write.table(DE_filtered, file = paste0(opt$out, file.name, ".txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
        Differential[[file.name]] <- DE_filtered[,1]
      }else if(n.sv > 0){
        # Model correction by sva
        svobj = sva(as.matrix(inputs[[l]][,samples]),mod,mod0,n.sv=n.sv)
        colnames(svobj$sv) <- paste(rep("col",ncol(svobj$sv)),c(1:ncol(svobj$sv)),sep="")
        modSv = cbind(mod,svobj$sv)
        # Fit model
        fit <-lmFit(inputs[[l]][,samples],modSv)
        # Select contrast
        contrastss <- c(paste0(colnames(modSv)[2], "-", colnames(modSv)[1]))
        contrast.matrix <- makeContrasts(contrastss, levels=modSv)
        fit2 <- contrasts.fit(fit,contrast.matrix)
        fit2 <- eBayes(fit2)
        # Filter and export
        DE_filtered <- limma::topTable(fit2, number = Inf, coef = 1, p.value = 0.05, sort.by = "logFC")
        cat(paste0("\n \t ", contrastss, sep = " | ", nrow(inputs[[l]]), " total proteins , ", nrow(DE_filtered), " differential proteins (0.05) \n"))
        DE_filtered <- cbind(Accession = inputs[[l]][rownames(DE_filtered),1], DE_filtered)
        file.name <- paste0(names(inputs[l]), contrastss)
        write.table(DE_filtered, file = paste0(opt$out, file.name, ".txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
        Differential[[file.name]] <- DE_filtered[,1]
      }

    }
    
  }else if(!opt$sva){
    # limma
    cat("\n \t sva: OFF \n")
    for(i in 1:ncol(comb.contrast)){ 
    # Build model
    samples <- c(groups.list[[comb.contrast[1,i]]], groups.list[[comb.contrast[2,i]]]) + 1
    PhenoData <- data.frame(samples = colnames(inputs[[l]][,samples]), Treatment = c(rep(paste0("Treatment", comb.contrast[1,i]),groups[comb.contrast[1,i]]), c(rep(paste0("Treatment", comb.contrast[2,i]),groups[comb.contrast[2,i]]))))
    mod = model.matrix(~0+as.factor(Treatment), data=PhenoData)
    colnames(mod) <- levels(PhenoData$Treatment)
    # Fit model
    fit <-lmFit(inputs[[l]][,samples],mod)
    # Select contrast
    contrastss <- c(paste0(colnames(mod)[2], "-", colnames(mod)[1]))
    contrast.matrix <- makeContrasts(contrastss, levels=mod)
    fit2 <- contrasts.fit(fit,contrast.matrix)
    fit2 <- eBayes(fit2)
    # Filter and export
    DE_filtered <- limma::topTable(fit2, number = Inf, coef = 1, p.value = 0.05, sort.by = "logFC")
    cat(paste0("\n \t ", contrastss, sep  = " | ", nrow(inputs[[l]]), " total proteins , ", nrow(DE_filtered), " differential proteins (0.05) \n"))
    DE_filtered <- cbind(Accession = inputs[[l]][rownames(DE_filtered),1], DE_filtered)
    file.name <- paste0(names(inputs[l]), contrastss)
    write.table(DE_filtered, file = paste0(opt$out, file.name, ".txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
    Differential[[file.name]] <- DE_filtered[,1]
    }
  }
}}

##### 4. Version ID liftover:OK #####

cat("\n Executing CreID liftover conversion to v5.5! \n")
DB <- read.table("./Data/DB/ChlamydomonasTranscriptNameConversionBetweenReleases.Mch12b.txt", header = F)
colnames(DB) <- c("5.5", "3.1", "Genbank","4", "4.3", "u5", "u9", "5.3.1")
inputs <- lapply(inputs, FUN = function(x){ x[,1]})

if(opt$differential){VictorgoestoBED <- list(inputs = inputs, differential = Differential)
}else if(!opt$differential){ VictorgoestoBED <- list(inputs = inputs)}

  
VictorgoestoBED <- lapply(VictorgoestoBED, FUN = function(y){
  lapply(y, FUN = function(x){  
  versions <- list()
  for(i in 1:ncol(DB)){
    versions[i] <- length(which((x %in% DB[,i]) == TRUE))
  }
  version <- colnames(DB)[which(unlist(versions,use.names = F) == max(unlist(versions,use.names = F)))]
  if(version == "5.5"){
    x <- x[x %in% DB[,1]]
    x <- x[!is.na(x)] 
    as.character(x)
  }else if(version != "5.5"){
    newIDs <- list()
    for(i in 1:length(x)){
      if(length(as.character(DB[which(as.character(DB[,version]) == as.character(x[i])),1])) != 0){
        newIDs[[i]] <- as.character(DB[which(as.character(DB[,version]) == as.character(x[i])),1])
        if(length(newIDs[[i]]) > 1){ newIDs[[i]] <- newIDs[[i]][[1]]}
      }else if(length(as.character(DB[which(as.character(DB[,version]) == as.character(x[i])),1])) == 0){
        NA
      }}
    newIDs <- unlist(newIDs, use.names = F)
    newIDs <- newIDs[!is.na(newIDs)]
    as.character(newIDs)
     }
  })
})

##### 5. Background selection:OK #####

cat("\n Universe background (whole Cre proteome/coding-transcriptome) is a good choice for your enrichments. You can find it at Chlamytina/Data/DB/Universe!")
# File background is VictorgoestoBED$inputs
cat("\n File background may be a good reference for your enrichments! \n")

if(length(inputs) > 1){
if(opt$intersect){ 
  Global_background <- unique(unlist(VictorgoestoBED$inputs, use.names = F))
  VictorgoestoBED$Global_background <- Global_background
  cat("\n For the args defined, Global background may be a good reference for your enrichments! \n")
  }
}

if(opt$intersect & opt$differential){
  Diff_background <- unique(unlist(VictorgoestoBED$differential, use.names = F))
  VictorgoestoBED$Diff_background <- Diff_background
  cat("\n For the args defined, Differential background may be interesting for your enrichments! \n")
}

##### 6. Intersection:OK #####

if(length(inputs) > 1){ 
if(opt$intersect){
  # intersect between files
  cat("\n Intersect is defined as TRUE \n")
  myV <- nVennR::plotVenn(VictorgoestoBED$inputs, nCycles = 7000, opacity = 0.2, borderWidth = 3, systemShow = T, fontScale = 2, outFile = paste0(opt$out, "intersection.svg"))
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
  # intersect between diffs
  cat("\n Differential is defined as TRUE \n")
  myV <- nVennR::plotVenn(VictorgoestoBED$differential, nCycles = 7000, opacity = 0.2, borderWidth = 3, systemShow = T, fontScale = 2, outFile = paste0(opt$out, "diff_intersection.svg"))
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
} else{cat("\n Intersect is defined as FALSE \n")}
}else if(length(inputs) == 1 & length(groups[1]) == 2){
  if(opt$differential & opt$intersect){ cat("Not enough conditions/files to do differential intersection/intersection")
    } else{cat("\n Intersect is defined as FALSE \n")}
}else if(length(inputs) == 1 & length(groups[1]) > 2){
  if(opt$differential & opt$intersect){
  myV <- nVennR::plotVenn(VictorgoestoBED$differential, nCycles = 7000, opacity = 0.2, borderWidth = 3, systemShow = T, fontScale = 2)
  myV2 <- nVennR::listVennRegions(myV)
  nonames <- names(myV2)
  nonames <- sapply(nonames, USE.NAMES = F,FUN = function(x){
    strsplit(x, split = "(", fixed = T)[[1]][2]
  })
  names(myV2) <- paste0("Diff_uniq_", nonames)
  names(myV2) <- gsub(")", "", names(myV2), fixed = T)
  names(myV2) <- gsub(", ", "", names(myV2), fixed = T)
  VictorgoestoBED$intersect_diff <- myV2
  } else{cat("\n Intersect is defined as FALSE \n")}  
}

##### 7. From ID table to bed #####

unzip(zipfile = "./Data/DB/Phytozome_download.zip", exdir = "./Data/DB/unzip")
if(!("Creinhardtii_281_v5.5.gene_exons.gff3" %in% list.files("./Data/DB/unzip/Phytozome/PhytozomeV12_unrestricted/Creinhardtii/annotation/"))){
  R.utils::gunzip(filename = "./Data/DB/unzip/Phytozome/PhytozomeV12_unrestricted/Creinhardtii/annotation/Creinhardtii_281_v5.5.gene_exons.gff3.gz", destname = "./Data/DB/unzip/Phytozome/PhytozomeV12_unrestricted/Creinhardtii/annotation/Creinhardtii_281_v5.5.gene_exons.gff3")
}
TxDb.Cre <- makeTxDbFromGFF(file = "./Data/DB/unzip/Phytozome/PhytozomeV12_unrestricted/Creinhardtii/annotation/Creinhardtii_281_v5.5.gene_exons.gff3")
Transcripts <- transcripts(TxDb.Cre)
cat("\n Computing CreIDs coherence ... \n")

QC <- lapply(VictorgoestoBED$inputs, function(x){
  x %in% Transcripts$tx_name})
if(length(which(unlist(QC,use.names = F) == FALSE)) > 0){stop(say("STOP: Something go wrong with CreIDs. Check", by = "poop"))}
if(length(which(unlist(QC,use.names = F) == FALSE)) == 0){cat("\n Passed \n")}

cat("\n Converting to .bed ... \n")
AllIDs <- unique(unlist(VictorgoestoBED$inputs, use.names = F))
chr <- list()
start <- list()
end <- list()
for(j in 1:length(AllIDs)){
    chr[j] <- as.character(Transcripts@seqnames[which(Transcripts$tx_name == AllIDs[j])])
    start[j] <- Transcripts@ranges@start[which(Transcripts$tx_name == AllIDs[j])]
    end[j] <- Transcripts@ranges@start[which(Transcripts$tx_name == AllIDs[j])] + Transcripts@ranges@width[which(Transcripts$tx_name == AllIDs[j])] - 1
  }
All.bed <- data.frame(row.names = AllIDs, unlist(chr, use.names = F), unlist(start, use.names = F), unlist(end, use.names = F))
All.bed[,1] <- gsub("chromosome_", "chr", All.bed[,1])
if(opt$chromosome){
  cat("\n Chromosome filtering is defined as TRUE \n")
  filter <- grep("scaffold", All.bed[,1])
  if(length(filter) == 0){
    All.bed <- All.bed
  }else if(length(filter) > 0){
    All.bed <- All.bed[-filter,]
  }
}
colnames(All.bed) <- c("chr", "start", "end")
levels <- paste0("chr", c(1:17))
All.bed$chr <- factor(All.bed$chr, levels, ordered = TRUE)
All.bed <- All.bed[order(All.bed$chr, All.bed$start),]

##### 8. Export beds #####

for(x in 1:length(VictorgoestoBED)){
  if(class(VictorgoestoBED[[x]]) == "list"){
    # all-beds
    for(y in 1:length(VictorgoestoBED[[x]])){
     VictorgoestoBED[[x]][[y]] <- unique(VictorgoestoBED[[x]][[y]])
     VictorgoestoBED[[x]][[y]] <- VictorgoestoBED[[x]][[y]][VictorgoestoBED[[x]][[y]] %in% row.names(All.bed)]
     VictorgoestoBED[[x]][[y]] <- VictorgoestoBED[[x]][[y]][order(VictorgoestoBED[[x]][[y]])]
     VictorgoestoBED[[x]][[y]] <- All.bed[VictorgoestoBED[[x]][[y]],]
     write.table(VictorgoestoBED[[x]][[y]], file = paste0(opt$out,names(VictorgoestoBED[[x]][y]), ".bed"), row.names = FALSE, col.names = F, sep = "\t", quote = FALSE)
    }
  }else if(class(VictorgoestoBED[[x]]) != "list"){
    # backgrounds except files
    VictorgoestoBED[[x]] <- unique(VictorgoestoBED[[x]])
    VictorgoestoBED[[x]] <- VictorgoestoBED[[x]][VictorgoestoBED[[x]] %in% row.names(All.bed)]
    VictorgoestoBED[[x]] <- VictorgoestoBED[[x]][order(VictorgoestoBED[[x]])]
    VictorgoestoBED[[x]] <- All.bed[VictorgoestoBED[[x]],]
    write.table(VictorgoestoBED[[x]], file = paste0(opt$out,names(VictorgoestoBED[x]), ".bed"), row.names = FALSE, col.names = F, sep = "\t", quote = FALSE)
  }
}

say(paste0("DataPrepare has finished ! ", Sys.time()), by = "rabbit", what_color = "white", by_color = "yellow")

