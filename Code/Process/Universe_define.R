##### Universe Chlamytina ###########
## Side_note01: For phytozome gff3 gene_exons match exactly with proteins.fa so we think that these files are filtered for only protein-coding transcripts.
## Whole-transcriptome and proteome Universe in this particular case are the same

setwd("D:/AA_GIT_myproject/Chlamytina/")

# 1 Whole-transcriptome/Proteome Universe (Phytozome-based)

unzip(zipfile = "./Data/DB/Phytozome_download.zip", exdir = "./Data/DB/unzip")
if(!("Creinhardtii_281_v5.5.gene_exons.gff3" %in% list.files("./Data/DB/unzip/Phytozome/PhytozomeV12_unrestricted/Creinhardtii/annotation/"))){
  R.utils::gunzip(filename = "./Data/DB/unzip/Phytozome/PhytozomeV12_unrestricted/Creinhardtii/annotation/Creinhardtii_281_v5.5.gene_exons.gff3.gz", destname = "./Data/DB/unzip/Phytozome/PhytozomeV12_unrestricted/Creinhardtii/annotation/Creinhardtii_281_v5.5.gene_exons.gff3")
}

TxDb.Cre <- makeTxDbFromGFF(file = "./Data/DB/unzip/Phytozome/PhytozomeV12_unrestricted/Creinhardtii/annotation/Creinhardtii_281_v5.5.gene_exons.gff3")
Transcripts <- transcripts(TxDb.Cre)

chr <- list()
start <- list()
end <- list()
tx_name <- list()

for(i in 1:length(Transcripts)){
  chr[i] <- as.character(Transcripts@seqnames[i])
  start[i] <- Transcripts@ranges@start[i]
  end[i] <- Transcripts@ranges@start[i] + Transcripts@ranges@width[i] - 1
  tx_name[i] <- Transcripts$tx_name[i]
  }

All.Transcripts <- data.frame(row.names = unlist(tx_name, use.names = F), chr = unlist(chr, use.names = F), start = unlist(start, use.names = F), 
                              end = unlist(end, use.names = F), tx_name = unlist(tx_name, use.names = F))
# chromosome filtering because scaffolds are highly variable between assemblies

All.Transcripts[,1] <- gsub("chromosome_", "chr", All.Transcripts[,1])
filter <- grep("scaffold", All.Transcripts[,1])
if(length(filter) == 0){
  All.Transcripts <- All.Transcripts
}else if(length(filter) > 0){
  All.Transcripts <- All.Transcripts[-filter,]
}
colnames(All.Transcripts) <- c("chr", "start", "end", "tx_name")
write.table(All.Transcripts, file =  "D:/AA_GIT_myproject/Chlamytina/Data/DB/Universe backgrounds/TranscriptUniverse.bed", row.names = FALSE, col.names = F, sep = "\t", quote = FALSE)
