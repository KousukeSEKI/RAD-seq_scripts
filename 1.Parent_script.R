library(Biostrings)
library(ShortRead) 
library(QuasR)   
library(data.table)

dir.create("path/RAD-seq/NGS_Data_Analysis_folder/parent")
in_f <- c("path/RAD-seq/NGS_raw_data_folder/P1_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/P1_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/parent/P1.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/parent/P1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/parent/P_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(resP1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/parent/P_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/parent/P_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(resP2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/parent/P_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/parent/P_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(resP3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/parent/P_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/parent/AdapterTrim_P1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(resP4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res4 <- rbind(resP1,resP2,resP3,resP4)
write.csv(res4, "path/RAD-seq/NGS_Data_Analysis_folder/parent/P1.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/parent/P_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/parent/P_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/parent/P_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/parent/AdapterTrim_P1.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/parent/P1.tagcount", row.names = F))
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/parent/P1.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/P1_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/P1_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/parent/P2.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/parent/P2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/parent/P_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(resP1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/parent/P_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/parent/P_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(resP2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/parent/P_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/parent/P_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(resP3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/parent/P_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/parent/AdapterTrim_P2.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(resP4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res4 <- rbind(resP1,resP2,resP3,resP4)
write.csv(res4, "path/RAD-seq/NGS_Data_Analysis_folder/parent/P2.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/parent/P_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/parent/P_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/parent/P_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/parent/AdapterTrim_P2.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/parent/P2.tagcount", row.names = F))
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/parent/P2.tagcount1over", row.names = F))
data <- c("tags","freq")
P1csv <- read.csv("path/RAD-seq/NGS_Data_Analysis_folder/parent/P1.tagcount", header=F, skip=1, col.names=data)
P2csv <- read.csv("path/RAD-seq/NGS_Data_Analysis_folder/parent/P2.tagcount", header=F, skip=1, col.names=data)
P1P2tagcount <- rbind(P1csv,P2csv)
system.time(write.csv(sort(table(P1P2tagcount$tags), decreasing=F), file="path/RAD-seq/NGS_Data_Analysis_folder/parent/P1P2tagcount.table", row.names = F))
P1P2tagcount_table <- read.csv("path/RAD-seq/NGS_Data_Analysis_folder/parent/P1P2tagcount.table", header=F, skip=1, col.names=data)
system.time(write.csv(subset(P1P2tagcount_table, freq == 2, select=c(tags)), file="path/RAD-seq/NGS_Data_Analysis_folder/parent/P1P2common.tag", row.names = F))
system.time(write.csv(subset(P1P2tagcount_table, freq == 1, select=c(tags)), file="path/RAD-seq/NGS_Data_Analysis_folder/parent/P1P2unique.tag", row.names = F))
P1P2common <- read.csv("path/RAD-seq/NGS_Data_Analysis_folder/parent/P1P2common.tag", header=F, skip=1, col.names=data)
system.time(write.csv(P1csv[match(P1P2common$tags, P1csv$tags),1:2], file="path/RAD-seq/NGS_Data_Analysis_folder/parent/P1common_tag.csv", row.names = F))
system.time(write.csv(P2csv[match(P1P2common$tags, P2csv$tags),1:2], file="path/RAD-seq/NGS_Data_Analysis_folder/parent/P2common_tag.csv", row.names = F))
P1P2unique <- read.csv("path/RAD-seq/NGS_Data_Analysis_folder/parent/P1P2unique.tag", header=F, skip=1, col.names=data)
system.time(write.csv(P1csv[match(P1P2unique$tags, P1csv$tags),1:2], file="path/RAD-seq/NGS_Data_Analysis_folder/parent/P1unique_tag.csv", row.names = F))
system.time(write.csv(P2csv[match(P1P2unique$tags, P2csv$tags),1:2], file="path/RAD-seq/NGS_Data_Analysis_folder/parent/P2unique_tag.csv", row.names = F))
P1unique_tag <- read.csv("path/RAD-seq/NGS_Data_Analysis_folder/parent/P1unique_tag.csv", header=F, skip=1, col.names=data)
system.time(write.csv(subset(P1unique_tag, freq > 1, select=c(tags,freq)), file="path/RAD-seq/NGS_Data_Analysis_folder/parent/P1unique_tag1over.csv", row.names = F))
P2unique_tag <- read.csv("path/RAD-seq/NGS_Data_Analysis_folder/parent/P2unique_tag.csv", header=F, skip=1, col.names=data)
system.time(write.csv(subset(P2unique_tag, freq > 1, select=c(tags,freq)), file="path/RAD-seq/NGS_Data_Analysis_folder/parent/P2unique_tag1over.csv", row.names = F))
P1csv1over <- read.csv("path/RAD-seq/NGS_Data_Analysis_folder/parent/P1.tagcount1over", header=F, skip=1, col.names=data)
P2csv1over <- read.csv("path/RAD-seq/NGS_Data_Analysis_folder/parent/P2.tagcount1over", header=F, skip=1, col.names=data)
P1P2tagcount1over <- rbind(P1csv1over,P2csv1over)
system.time(write.csv(sort(table(P1P2tagcount1over$tags), decreasing=F), file="path/RAD-seq/NGS_Data_Analysis_folder/parent/P1P2tagcoun1over.table", row.names = F))
P1P2tagcount1over_table <- read.csv("path/RAD-seq/NGS_Data_Analysis_folder/parent/P1P2tagcoun1over.table", header=F, skip=1, col.names=data)
system.time(write.csv(subset(P1P2tagcount1over_table, freq == 2, select=c(tags)), file="path/RAD-seq/NGS_Data_Analysis_folder/parent/P1P2common1over.tag", row.names = F))
P1P2common1over <- read.csv("path/RAD-seq/NGS_Data_Analysis_folder/parent/P1P2common1over.tag", header=F, skip=1, col.names=data)
system.time(write.csv(P1csv[match(P1P2common1over$tags, P1csv1over$tags),1:2], file="path/RAD-seq/NGS_Data_Analysis_folder/parent/P1common_tag1over.csv", row.names = F))
system.time(write.csv(P2csv[match(P1P2common1over$tags, P2csv1over$tags),1:2], file="path/RAD-seq/NGS_Data_Analysis_folder/parent/P2common_tag1over.csv", row.names = F))
P1unique_tag1over <- read.csv("path/RAD-seq/NGS_Data_Analysis_folder/parent/P1unique_tag1over.csv", header=F, skip=1, col.names=data)
P2unique_tag1over <- read.csv("path/RAD-seq/NGS_Data_Analysis_folder/parent/P2unique_tag1over.csv", header=F, skip=1, col.names=data)
P1P2unique_tag1over <- rbind(P1unique_tag1over,P2unique_tag1over)
system.time(write.csv(subset(P1P2unique_tag1over, select=c(tags)), file="path/RAD-seq/NGS_Data_Analysis_folder/parent/All_genotypelist.csv", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/parent/P1P2tagcount.table")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/parent/P1P2common.tag")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/parent/P1P2unique.tag")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/parent/P1P2common1over.tag")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/parent/P1P2tagcoun1over.table")
tags <- fread("path/RAD-seq/NGS_Data_Analysis_folder/parent/P1unique_tag1over.csv", header=T)
tags2 <- cbind(c(2*(1:nrow(tags))),tags[,1])
names(tags2) <- c("number","tags")
tags3 <- c(paste(">P1unique_tags1over_",1:nrow(tags),sep=""))
write.csv(tags3,"path/RAD-seq/NGS_Data_Analysis_folder/tags3.csv", row.names=F)
tags4 <- fread("path/RAD-seq/NGS_Data_Analysis_folder/tags3.csv", header=T)
tags5 <- cbind(c(2*(1:nrow(tags))-1),tags4)
names(tags5) <- c("number","tags")
tags6 <- rbind(tags2,tags5)
tags7 <- data.frame(tags6)
tags8 <- tags7[order(tags7$number),]
write(tags8$tags,"path/RAD-seq/NGS_Data_Analysis_folder/parent/P1_unique_tag1over.fasta")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/tags3.csv")
tags <- fread("path/RAD-seq/NGS_Data_Analysis_folder/parent/P2unique_tag1over.csv", header=F,skip=1)
tags2 <- cbind(c(2*(1:nrow(tags))),tags[,1])
names(tags2) <- c("number","tags")
tags3 <- c(paste(">P2unique_tags1over_",1:nrow(tags),sep=""))
write.csv(tags3,"path/RAD-seq/NGS_Data_Analysis_folder/tags3.csv", row.names=F)
tags4 <- fread("path/RAD-seq/NGS_Data_Analysis_folder/tags3.csv", header=T)
tags5 <- cbind(c(2*(1:nrow(tags))-1),tags4)
names(tags5) <- c("number","tags")
tags6 <- rbind(tags2,tags5)
tags7 <- data.frame(tags6)
tags8 <- tags7[order(tags7$number),]
write(tags8$tags,"path/RAD-seq/NGS_Data_Analysis_folder/parent/P2_unique_tag1over.fasta")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/tags3.csv")
tags <- fread("path/RAD-seq/NGS_Data_Analysis_folder/parent/P1common_tag1over.csv", header=F,skip=1)
tags2 <- cbind(c(2*(1:nrow(tags))),tags[,1])
names(tags2) <- c("number","tags")
tags3 <- c(paste(">P1P2_common_tags1over_",1:nrow(tags),sep=""))
write.csv(tags3,"path/RAD-seq/NGS_Data_Analysis_folder/tags3.csv", row.names=F)
tags4 <- fread("path/RAD-seq/NGS_Data_Analysis_folder/tags3.csv", header=T)
tags5 <- cbind(c(2*(1:nrow(tags))-1),tags4)
names(tags5) <- c("number","tags")
tags6 <- rbind(tags2,tags5)
tags7 <- data.frame(tags6)
tags8 <- tags7[order(tags7$number),]
write(tags8$tags,"path/RAD-seq/NGS_Data_Analysis_folder/parent/P1P2_common_tag1over.fasta")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/tags3.csv")
