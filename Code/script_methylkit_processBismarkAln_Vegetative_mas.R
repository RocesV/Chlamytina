#methylation call levels with MethylKit
library(methylKit)
library(tidyverse)

#Generación de los txt para la muestra para todos los contextos de metilación. processBismarkAln es en sí
#la fucnión que genera los niveles de metilación
setwd("/Users/ana/Desktop/Procesado_MethylKit")

file.list.bam.vegmas <- list(file="SRR2051057_rmdup_trimmed_bismark_bt2_sorted.bam")

vegmas <- processBismarkAln(location = file.list.bam.vegmas, sample.id = list("vegetative_mt+"),
                            assembly = "chl", save.folder = "methylation_calls_methylkit",
                            save.context = c("CpG","CHG","CHH"), nolap = FALSE, mincov = 10, 
                            minqual = 20, phred64 = FALSE, treatment = c(0))

#Workflow con la Vegetativa mt+ (SRR2051057)

#Leer las methylation calls, unir los tres contextos en una misma tabla y ordenarlas
setwd("/Users/ana/Desktop/Procesado_MethylKit/methylation_calls_methylkit")

vegmas_1 <- read.table(file = "vegetative_mt+_CHG.txt", sep="\t", header = T)
vegmas_1$con <- "CHG"
vegmas_2 <- read.table(file = "vegetative_mt+_CHH.txt", sep="\t", header = T)
vegmas_2$con <- "CHH"
vegmas_3 <- read.table(file = "vegetative_mt+_CpG.txt", sep="\t", header = T)
vegmas_3$con <- "CpG"

vegmas_all <- rbind(vegmas_1,vegmas_2,vegmas_3)
vegmas_all_sorted <- vegmas_all[order(vegmas_all[,2],vegmas_all[,3]),]

write.table(vegmas_all_sorted, file = "vegmas_all_sorted.txt",  sep = " ")