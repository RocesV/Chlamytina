##### Process #####

# Order here by chr and by start
# filter non-chromosome intervals

setwd("D:/AA_Project4_Chlamytina/Epigenetic marks/inputs/TFs/StressS/")

for(i in list.files()){
  file <- read.table(i, header = FALSE)
  file[,1] <- gsub("chr_", "chr", file[,1])
  file[,1] <- gsub("chromosome_", "chr", file[,1])
  filter <- grep("scaffold", file[,1])
  if(length(filter) == 0){
    file <- file
  }else if(length(filter) > 0){
    file <- file[-filter,]
  }
  colnames(file) <- c("chr", "start", "end")
  levels <- paste0("chr", c(1:17))
  levels <- c(levels, "chrC","chrM")
  file$chr <- factor(file$chr, levels, ordered = TRUE)
  file <- file[order(file$chr, file$start),]
  write.table(file, file = paste0("D:/AA_Project4_Chlamytina/Epigenetic marks/Process/TFs/StressS/", strsplit(i, split = ".", fixed = TRUE)[[1]][1] ,"_sortfilt",".bed") , row.names = FALSE, col.names = F, sep = "\t", quote = FALSE)
}

# MNase

nucleosome <- read.delim("clipboard", header = TRUE)
nucleosome <- nucleosome[,c(1,2,3,7,6)]
nucleosome[,1] <- gsub("chr_", "chr", nucleosome[,1])
nucleosome[,1] <- gsub("chromosome_", "chr", nucleosome[,1])
filter <- grep("scaffold", nucleosome[,1])
if(length(filter) == 0){
  nucleosome <- nucleosome
}else if(length(filter) > 0){
  nucleosome <- nucleosome[-filter,]
}
levels <- paste0("chr", c(1:17))
levels <- c(levels, "chrC","chrM")
nucleosome$chr <- factor(nucleosome$chr, levels, ordered = TRUE)
nucleosome <- nucleosome[order(nucleosome$chr, nucleosome$start),]
write.table(nucleosome, file = "D:/AA_Project4_Chlamytina/Epigenetic marks/Process/MNase/nucleosome_sortfilt.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

# 6mA

sixmA_light <- read.delim("clipboard", header = T)
sixmA_light <- sixmA_light[,c(1,2,3,10,8)]
sixmA_light[,1] <- gsub("chr_", "chr", sixmA_light[,1])
sixmA_light[,1] <- gsub("chromosome_", "chr", sixmA_light[,1])
filter <- grep("scaffold", sixmA_light[,1])
if(length(filter) == 0){
  sixmA_light <- sixmA_light
}else if(length(filter) > 0){
  sixmA_light <- sixmA_light[-filter,]
}
levels <- paste0("chr", c(1:17))
levels <- c(levels, "chrC","chrM")
sixmA_light$chr <- factor(sixmA_light$chr, levels, ordered = TRUE)
sixmA_light <- sixmA_light[order(sixmA_light$chr, sixmA_light$start),]
write.table(sixmA_light, file = "D:/AA_Project4_Chlamytina/Epigenetic marks/Process/6mA/6mA_dark_sortfilt.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

# 5mC: create regions, filter and sort


setwd("D:/AA_Project4_Chlamytina/Epigenetic marks/inputs/m5C/")

for(i in list.files()){
  file <- read.table(i, header = T)
  file <- file[,-8]
  file$chr <- gsub("chr_", "chr", file$chr)
  file$chr <- gsub("chromosome_", "chr", file$chr)
  file$chrBase <- gsub("chr_", "chr", file$chrBase)
  file$chrBase <- gsub("chromosome_", "chr", file$chrBase)
  filter <- grep("scaffold", file[,2])
  if(length(filter) == 0){
    file <- file
  }else if(length(filter) > 0){
    file <- file[-filter,]
  }
  levels <- paste0("chr", c(1:17))
  levels <- c(levels, "chrC","chrM")
  file$chr <- factor(file$chr, levels, ordered = TRUE)
  file <- file[order(file$chr, file$base),]
  write.table(file, paste0(i, "_formated.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  myobj=methRead(paste0(i, "_formated.txt"),
                 sample.id= i,assembly="chl", context = c("CpG", "CHG", "CHH"))
  res=methSeg(myobj,diagnostic.plot=TRUE,maxInt=100,minSeg=10,G=1:6, join.neighbours = TRUE)
  res <- res[which(res$seg.group == 6),]
  bed <- data.frame(chr = res@seqnames, start = res@ranges@start, end = (res@ranges@start + res@ranges@width - 1),
                    IDs = paste(res$ID, c(1:length(res$ID)), sep = "_"), score = res$seg.mean)
  
  write.table(bed, file = paste("D:/AA_Project4_Chlamytina/Epigenetic marks/Process/5mC/",i), row.names = F, col.names = F, sep = "\t", quote = FALSE)
}



