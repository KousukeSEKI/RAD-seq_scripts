library(Biostrings) 
library(data.table)

colnames <- c("sam1","sam2","sam3","sam4","sam5","sam6","sam7","sam8","sam9","sam10","sam11","sam12","sam13","sam14","sam15","sam16","sam17","sam18","sam19","sam20","sam21")
tagsP1 <- read.delim("path/RAD-seq/NGS_Data_Analysis_folder/parent/P1_unique_tag1over_mapped_aln.sam", header=F, fill=TRUE, na.strings="NA", col.names=colnames)
write.csv(tagsP1,"path/RAD-seq/NGS_Data_Analysis_folder/parent/P1_tags1.csv", row.names=F)
tagsP1 <- fread("path/RAD-seq/NGS_Data_Analysis_folder/parent/P1_tags1.csv", header=F, fill=TRUE, na.strings="NA", strip.white=TRUE, verbose=TRUE)
colnames(tagsP1) <- c("sam1","sam2","sam3","sam4","sam5","sam6","sam7","sam8","sam9","sam10","sam11","sam12","sam13","sam14","sam15","sam16","sam17","sam18","sam19","sam20","sam21")
tags1 <- tagsP1[sam2==16, sam10]
write.csv(tags1,"path/RAD-seq/NGS_Data_Analysis_folder/parent/P1_tags1.csv", row.names=F)
tags1 <- fread("path/RAD-seq/NGS_Data_Analysis_folder/parent/P1_tags1.csv", header=T)
file.remove("path/RAD-seq/NGS_Data_Analysis_folder//parent/P1_tags1.csv")
tags2 <- cbind(c(2*(1:nrow(tags1))),tags1[,1])
names(tags2) <- c("number","tags")
tags3 <- c(paste(">16_",1:nrow(tags1),sep=""))
write.csv(tags3,"path/RAD-seq/NGS_Data_Analysis_folder/parent/tags3.csv", row.names=F)
tags4 <- fread("path/RAD-seq/NGS_Data_Analysis_folder/parent/tags3.csv", header=T)
tags5 <- cbind(c(2*(1:nrow(tags1))-1),tags4)
names(tags5) <- c("number","tags")
tags6 <- rbind(tags2,tags5)
tags7 <- data.frame(tags6)
tags8 <- tags7[order(tags7$number),]
write(tags8$tags,"path/RAD-seq/NGS_Data_Analysis_folder/parent/P1_16.fasta")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/parent/tags3.csv")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/parent/P1_16.fasta"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/parent/P1_16_revcomp.fasta"
fasta <- readDNAStringSet(in_f, format="fasta")
fasta <- reverseComplement(fasta)
writeXStringSet(fasta, file=out_f, format="fasta", width=150)
system.time(fasta <- fread("path/RAD-seq/NGS_Data_Analysis_folder/parent/P1_16_revcomp.fasta", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotideP1 <- fasta[2+2*(0:((nrow(fasta)/2)-1))])
colnames(nucleotideP1) <- c("sam10")
tagsP1sam2_16f <- tagsP1[tagsP1$sam2=="16", sam1:9]
tagsP1sam2_16r <- tagsP1[tagsP1$sam2=="16", sam11:21]
tagsP1_16recomp <- cbind(tagsP1sam2_16f, nucleotideP1)
tagsP1_16recomp <- cbind(tagsP1_16recomp, tagsP1sam2_16r)
tagsP1sam2_0 <- tagsP1[tagsP1$sam2=="0", ]
tagsP1_0_16revcomp <- rbind(tagsP1sam2_0, tagsP1_16recomp)
write.csv(tagsP1_0_16revcomp,"path/RAD-seq/NGS_Data_Analysis_folder/parent/P1_unique_tag1over_mapped_0_16revcomp.csv", row.names=F)
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/parent/P1_16.fasta")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/parent/P1_16_revcomp.fasta")
tagsP1_0_16revcomp_37 <- tagsP1_0_16revcomp[tagsP1_0_16revcomp$sam5>=0 & substr(tagsP1_0_16revcomp$sam20, 1,5)=="", ]
tagsP1_0_16revcomp_37_2 <- tagsP1_0_16revcomp[tagsP1_0_16revcomp$sam5>=0 & substr(tagsP1_0_16revcomp$sam20, 1,5)=="MD:Z:", ]
tagsP1_0_16revcomp_37_3 <- rbind(tagsP1_0_16revcomp_37, tagsP1_0_16revcomp_37_2)
tagsP1_0_16revcomp_37 <- tagsP1_0_16revcomp_37_3[substr(tagsP1_0_16revcomp_37_3$sam21, 1,5)=="", ]
tagsP1_0_16revcomp_37$sam4 <- as.numeric(tagsP1_0_16revcomp_37$sam4)
tagsP1_0_16revcomp_37sort1 <- tagsP1_0_16revcomp_37[order(tagsP1_0_16revcomp_37$sam4),]
tagsP1_0_16revcomp_37sort2 <- tagsP1_0_16revcomp_37sort1[order(tagsP1_0_16revcomp_37sort1$sam3),]
write.csv(tagsP1_0_16revcomp_37sort2,"path/RAD-seq/NGS_Data_Analysis_folder/parent/P1_unique_tag1over_mapped_0_16revcomp_0-37.csv", row.names=F)

colnames <- c("sam1","sam2","sam3","sam4","sam5","sam6","sam7","sam8","sam9","sam10","sam11","sam12","sam13","sam14","sam15","sam16","sam17","sam18","sam19","sam20","sam21")
tagsP2 <- read.delim("path/RAD-seq/NGS_Data_Analysis_folder/parent/P2_unique_tag1over_mapped_aln.sam", header=F, fill=TRUE, na.strings="NA", col.names=colnames)
write.csv(tagsP2,"path/RAD-seq/NGS_Data_Analysis_folder/parent/P2_tags1.csv", row.names=F)
tagsP2 <- fread("path/RAD-seq/NGS_Data_Analysis_folder/parent/P2_tags1.csv", header=F, fill=TRUE, na.strings="NA", strip.white=TRUE, verbose=TRUE)
colnames(tagsP2) <- c("sam1","sam2","sam3","sam4","sam5","sam6","sam7","sam8","sam9","sam10","sam11","sam12","sam13","sam14","sam15","sam16","sam17","sam18","sam19","sam20","sam21")
tags1 <- tagsP2[sam2==16, sam10]
write.csv(tags1,"path/RAD-seq/NGS_Data_Analysis_folder/parent/P2_tags1.csv", row.names=F)
tags1 <- fread("path/RAD-seq/NGS_Data_Analysis_folder/parent/P2_tags1.csv", header=T)
file.remove("path/RAD-seq/NGS_Data_Analysis_folder//parent/P2_tags1.csv")
tags2 <- cbind(c(2*(1:nrow(tags1))),tags1[,1])
names(tags2) <- c("number","tags")
tags3 <- c(paste(">16_",1:nrow(tags1),sep=""))
write.csv(tags3,"path/RAD-seq/NGS_Data_Analysis_folder/parent/tags3.csv", row.names=F)
tags4 <- fread("path/RAD-seq/NGS_Data_Analysis_folder/parent/tags3.csv", header=T)
tags5 <- cbind(c(2*(1:nrow(tags1))-1),tags4)
names(tags5) <- c("number","tags")
tags6 <- rbind(tags2,tags5)
tags7 <- data.frame(tags6)
tags8 <- tags7[order(tags7$number),]
write(tags8$tags,"path/RAD-seq/NGS_Data_Analysis_folder/parent/P2_16.fasta")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/parent/tags3.csv")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/parent/P2_16.fasta"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/parent/P2_16_revcomp.fasta"
fasta <- readDNAStringSet(in_f, format="fasta")
fasta <- reverseComplement(fasta)
writeXStringSet(fasta, file=out_f, format="fasta", width=150)
system.time(fasta <- fread("path/RAD-seq/NGS_Data_Analysis_folder/parent/P2_16_revcomp.fasta", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotideP2 <- fasta[2+2*(0:((nrow(fasta)/2)-1))])
colnames(nucleotideP2) <- c("sam10")
tagsP2sam2_16f <- tagsP2[tagsP2$sam2=="16", sam1:9]
tagsP2sam2_16r <- tagsP2[tagsP2$sam2=="16", sam11:21]
tagsP2_16recomp <- cbind(tagsP2sam2_16f, nucleotideP2)
tagsP2_16recomp <- cbind(tagsP2_16recomp, tagsP2sam2_16r)
tagsP2sam2_0 <- tagsP2[tagsP2$sam2=="0", ]
tagsP2_0_16revcomp <- rbind(tagsP2sam2_0, tagsP2_16recomp)
write.csv(tagsP2_0_16revcomp,"path/RAD-seq/NGS_Data_Analysis_folder/parent/P2_unique_tag1over_mapped_0_16revcomp.csv", row.names=F)
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/parent/P2_16.fasta")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/parent/P2_16_revcomp.fasta")
tagsP2_0_16revcomp_37 <- tagsP2_0_16revcomp[tagsP2_0_16revcomp$sam5>=0 & substr(tagsP2_0_16revcomp$sam20, 1,5)=="", ]
tagsP2_0_16revcomp_37_2 <- tagsP2_0_16revcomp[tagsP2_0_16revcomp$sam5>=0 & substr(tagsP2_0_16revcomp$sam20, 1,5)=="MD:Z:", ]
tagsP2_0_16revcomp_37_3 <- rbind(tagsP2_0_16revcomp_37, tagsP2_0_16revcomp_37_2)
tagsP2_0_16revcomp_37 <- tagsP2_0_16revcomp_37_3[substr(tagsP2_0_16revcomp_37_3$sam21, 1,5)=="", ]
tagsP2_0_16revcomp_37$sam4 <- as.numeric(tagsP2_0_16revcomp_37$sam4)
tagsP2_0_16revcomp_37sort1 <- tagsP2_0_16revcomp_37[order(tagsP2_0_16revcomp_37$sam4),]
tagsP2_0_16revcomp_37sort2 <- tagsP2_0_16revcomp_37sort1[order(tagsP2_0_16revcomp_37sort1$sam3),]
write.csv(tagsP2_0_16revcomp_37sort2,"path/RAD-seq/NGS_Data_Analysis_folder/parent/P2_unique_tag1over_mapped_0_16revcomp_0-37.csv", row.names=F)

tagsP1P2_0_16revcomp <- rbind(tagsP1_0_16revcomp, tagsP2_0_16revcomp)
write.csv(tagsP1P2_0_16revcomp, "path/RAD-seq/NGS_Data_Analysis_folder/parent/P1P2_unique_tag1over_mapped_0_16revcomp.csv", row.names=FALSE)
tagsP1P2_0_16revcomp_37 <- rbind(tagsP1_0_16revcomp_37sort2, tagsP2_0_16revcomp_37sort2)
write.csv(tagsP1P2_0_16revcomp_37, "path/RAD-seq/NGS_Data_Analysis_folder/parent/P1P2_unique_tag1over_mapped_0_16revcomp_0-37.csv", row.names=FALSE)
mappinglist <- tagsP1P2_0_16revcomp_37[, sam1:5]
genotypelist <- tagsP1P2_0_16revcomp_37[, sam10]
mappinggenotypelist <- cbind(mappinglist, genotypelist)
colnames(mappinggenotypelist) <- c("Marker","Flag","Linkage_group","Position","Quality","tags")
write.csv(mappinggenotypelist, "path/RAD-seq/NGS_Data_Analysis_folder/parent/mappinggenotypelist.csv", row.names=FALSE)
