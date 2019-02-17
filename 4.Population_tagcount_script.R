library(Biostrings)
library(ShortRead)
library(QuasR)   
library(data.table)

in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_001_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_001_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_001.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_001.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_001.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_001.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_001.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_001.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_001.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_001.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_002_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_002_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_002.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_002.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_002.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_002.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_002.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_002.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_002.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_002.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_003_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_003_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_003.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_003.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_003.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_003.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_003.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_003.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_003.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_003.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_004_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_004_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_004.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_004.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_004.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_004.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_004.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_004.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_004.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_004.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_005_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_005_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_005.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_005.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_005.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_005.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_005.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_005.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_005.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_005.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_006_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_006_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_006.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_006.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_006.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_006.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_006.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_006.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_006.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_006.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_007_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_007_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_007.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_007.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_007.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_007.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_007.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_007.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_007.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_007.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_008_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_008_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_008.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_008.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_008.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_008.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_008.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_008.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_008.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_008.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_009_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_009_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_009.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_009.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_009.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_009.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_009.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_009.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_009.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_009.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_010_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_010_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_010.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_010.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_010.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_010.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_010.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_010.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_010.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_010.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_011_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_011_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_011.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_011.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_011.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_011.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_011.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_011.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_011.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_011.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_012_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_012_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_012.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_012.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_012.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_012.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_012.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_012.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_012.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_012.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_013_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_013_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_013.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_013.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_013.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_013.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_013.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_013.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_013.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_013.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_014_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_014_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_014.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_014.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_014.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_014.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_014.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_014.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_014.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_014.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_015_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_015_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_015.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_015.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_015.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_015.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_015.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_015.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_015.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_015.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_016_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_016_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_016.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_016.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_016.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_016.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_016.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_016.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_016.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_016.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_017_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_017_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_017.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_017.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_017.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_017.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_017.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_017.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_017.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_017.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_018_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_018_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_018.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_018.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_018.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_018.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_018.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_018.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_018.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_018.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_019_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_019_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_019.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_019.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_019.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_019.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_019.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_019.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_019.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_019.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_020_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_020_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_020.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_020.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_020.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_020.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_020.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_020.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_020.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_020.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_021_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_021_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_021.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_021.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_021.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_021.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_021.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_021.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_021.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_021.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_022_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_022_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_022.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_022.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_022.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_022.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_022.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_022.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_022.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_022.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_023_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_023_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_023.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_023.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_023.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_023.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_023.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_023.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_023.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_023.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_024_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_024_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_024.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_024.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_024.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_024.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_024.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_024.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_024.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_024.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_025_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_025_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_025.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_025.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_025.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_025.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_025.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_025.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_025.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_025.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_026_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_026_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_026.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_026.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_026.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_026.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_026.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_026.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_026.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_026.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_027_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_027_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_027.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_027.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_027.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_027.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_027.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_027.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_027.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_027.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_028_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_028_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_028.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_028.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_028.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_028.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_028.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_028.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_028.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_028.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_029_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_029_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_029.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_029.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_029.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_029.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_029.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_029.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_029.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_029.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_030_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_030_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_030.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_030.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_030.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_030.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_030.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_030.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_030.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_030.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_031_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_031_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_031.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_031.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_031.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_031.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_031.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_031.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_031.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_031.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_032_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_032_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_032.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_032.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_032.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_032.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_032.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_032.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_032.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_032.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_033_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_033_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_033.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_033.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_033.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_033.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_033.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_033.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_033.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_033.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_034_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_034_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_034.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_034.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_034.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_034.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_034.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_034.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_034.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_034.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_035_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_035_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_035.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_035.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_035.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_035.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_035.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_035.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_035.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_035.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_036_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_036_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_036.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_036.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_036.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_036.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_036.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_036.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_036.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_036.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_037_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_037_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_037.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_037.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_037.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_037.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_037.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_037.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_037.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_037.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_038_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_038_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_038.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_038.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_038.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_038.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_038.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_038.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_038.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_038.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_039_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_039_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_039.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_039.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_039.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_039.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_039.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_039.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_039.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_039.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_040_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_040_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_040.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_040.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_040.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_040.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_040.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_040.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_040.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_040.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_041_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_041_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_041.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_041.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_041.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_041.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_041.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_041.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_041.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_041.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_042_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_042_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_042.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_042.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_042.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_042.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_042.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_042.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_042.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_042.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_043_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_043_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_043.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_043.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_043.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_043.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_043.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_043.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_043.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_043.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_044_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_044_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_044.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_044.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_044.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_044.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_044.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_044.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_044.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_044.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_045_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_045_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_045.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_045.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_045.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_045.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_045.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_045.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_045.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_045.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_046_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_046_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_046.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_046.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_046.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_046.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_046.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_046.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_046.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_046.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_047_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_047_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_047.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_047.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_047.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_047.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_047.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_047.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_047.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_047.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_048_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_048_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_048.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_048.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_048.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_048.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_048.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_048.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_048.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_048.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_049_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_049_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_049.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_049.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_049.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_049.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_049.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_049.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_049.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_049.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_050_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_050_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_050.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_050.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_050.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_050.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_050.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_050.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_050.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_050.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_051_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_051_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_051.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_051.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_051.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_051.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_051.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_051.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_051.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_051.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_052_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_052_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_052.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_052.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_052.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_052.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_052.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_052.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_052.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_052.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_053_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_053_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_053.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_053.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_053.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_053.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_053.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_053.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_053.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_053.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_054_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_054_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_054.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_054.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_054.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_054.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_054.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_054.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_054.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_054.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_055_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_055_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_055.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_055.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_055.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_055.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_055.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_055.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_055.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_055.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_056_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_056_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_056.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_056.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_056.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_056.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_056.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_056.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_056.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_056.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_057_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_057_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_057.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_057.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_057.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_057.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_057.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_057.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_057.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_057.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_058_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_058_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_058.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_058.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_058.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_058.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_058.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_058.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_058.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_058.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_059_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_059_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_059.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_059.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_059.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_059.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_059.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_059.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_059.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_059.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_060_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_060_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_060.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_060.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_060.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_060.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_060.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_060.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_060.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_060.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_061_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_061_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_061.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_061.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_061.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_061.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_061.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_061.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_061.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_061.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_062_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_062_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_062.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_062.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_062.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_062.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_062.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_062.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_062.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_062.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_063_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_063_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_063.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_063.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_063.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_063.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_063.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_063.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_063.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_063.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_064_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_064_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_064.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_064.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_064.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_064.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_064.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_064.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_064.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_064.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_065_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_065_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_065.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_065.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_065.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_065.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_065.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_065.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_065.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_065.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_066_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_066_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_066.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_066.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_066.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_066.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_066.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_066.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_066.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_066.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_067_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_067_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_067.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_067.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_067.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_067.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_067.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_067.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_067.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_067.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_068_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_068_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_068.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_068.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_068.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_068.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_068.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_068.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_068.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_068.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_069_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_069_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_069.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_069.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_069.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_069.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_069.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_069.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_069.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_069.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_070_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_070_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_070.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_070.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_070.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_070.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_070.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_070.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_070.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_070.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_071_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_071_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_071.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_071.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_071.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_071.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_071.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_071.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_071.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_071.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_072_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_072_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_072.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_072.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_072.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_072.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_072.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_072.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_072.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_072.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_073_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_073_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_073.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_073.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_073.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_073.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_073.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_073.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_073.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_073.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_074_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_074_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_074.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_074.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_074.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_074.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_074.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_074.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_074.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_074.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_075_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_075_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_075.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_075.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_075.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_075.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_075.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_075.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_075.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_075.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_076_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_076_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_076.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_076.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_076.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_076.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_076.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_076.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_076.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_076.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_077_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_077_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_077.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_077.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_077.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_077.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_077.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_077.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_077.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_077.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_078_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_078_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_078.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_078.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_078.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_078.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_078.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_078.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_078.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_078.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_079_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_079_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_079.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_079.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_079.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_079.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_079.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_079.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_079.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_079.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_080_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_080_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_080.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_080.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_080.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_080.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_080.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_080.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_080.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_080.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_081_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_081_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_081.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_081.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_081.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_081.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_081.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_081.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_081.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_081.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_082_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_082_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_082.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_082.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_082.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_082.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_082.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_082.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_082.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_082.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_083_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_083_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_083.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_083.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_083.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_083.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_083.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_083.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_083.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_083.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_084_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_084_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_084.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_084.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_084.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_084.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_084.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_084.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_084.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_084.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_085_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_085_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_085.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_085.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_085.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_085.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_085.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_085.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_085.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_085.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_086_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_086_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_086.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_086.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_086.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_086.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_086.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_086.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_086.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_086.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_087_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_087_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_087.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_087.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_087.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_087.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_087.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_087.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_087.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_087.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_088_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_088_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_088.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_088.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_088.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_088.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_088.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_088.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_088.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_088.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_089_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_089_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_089.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_089.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_089.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_089.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_089.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_089.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_089.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_089.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_090_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_090_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_090.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_090.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_090.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_090.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_090.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_090.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_090.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_090.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_091_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_091_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_091.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_091.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_091.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_091.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_091.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_091.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_091.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_091.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_092_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_092_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_092.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_092.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_092.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_092.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_092.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_092.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_092.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_092.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_093_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_093_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_093.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_093.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_093.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_093.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_093.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_093.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_093.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_093.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_094_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_094_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_094.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_094.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_094.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_094.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_094.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_094.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_094.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_094.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_095_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_095_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_095.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_095.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_095.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_095.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_095.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_095.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_095.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_095.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_096_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_096_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_096.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_096.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_096.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_096.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_096.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_096.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_096.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_096.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_097_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_097_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_097.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_097.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_097.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_097.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_097.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_097.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_097.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_097.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_098_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_098_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_098.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_098.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_098.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_098.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_098.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_098.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_098.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_098.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_099_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_099_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_099.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_099.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_099.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_099.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_099.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_099.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_099.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_099.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_100_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_100_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_100.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_100.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_100.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_100.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_100.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_100.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_100.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_100.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_101_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_101_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_101.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_101.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_101.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_101.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_101.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_101.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_101.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_101.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_102_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_102_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_102.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_102.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_102.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_102.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_102.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_102.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_102.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_102.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_103_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_103_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_103.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_103.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_103.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_103.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_103.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_103.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_103.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_103.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_104_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_104_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_104.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_104.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_104.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_104.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_104.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_104.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_104.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_104.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_105_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_105_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_105.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_105.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_105.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_105.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_105.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_105.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_105.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_105.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_106_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_106_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_106.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_106.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_106.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_106.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_106.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_106.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_106.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_106.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_107_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_107_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_107.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_107.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_107.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_107.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_107.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_107.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_107.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_107.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_108_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_108_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_108.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_108.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_108.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_108.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_108.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_108.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_108.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_108.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_109_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_109_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_109.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_109.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_109.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_109.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_109.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_109.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_109.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_109.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_110_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_110_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_110.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_110.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_110.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_110.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_110.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_110.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_110.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_110.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_111_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_111_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_111.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_111.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_111.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_111.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_111.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_111.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_111.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_111.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_112_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_112_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_112.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_112.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_112.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_112.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_112.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_112.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_112.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_112.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_113_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_113_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_113.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_113.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_113.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_113.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_113.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_113.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_113.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_113.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_114_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_114_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_114.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_114.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_114.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_114.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_114.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_114.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_114.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_114.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_115_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_115_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_115.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_115.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_115.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_115.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_115.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_115.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_115.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_115.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_116_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_116_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_116.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_116.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_116.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_116.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_116.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_116.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_116.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_116.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_117_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_117_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_117.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_117.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_117.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_117.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_117.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_117.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_117.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_117.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_118_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_118_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_118.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_118.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_118.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_118.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_118.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_118.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_118.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_118.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_119_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_119_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_119.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_119.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_119.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_119.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_119.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_119.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_119.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_119.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_120_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_120_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_120.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_120.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_120.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_120.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_120.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_120.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_120.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_120.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_121_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_121_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_121.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_121.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_121.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_121.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_121.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_121.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_121.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_121.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_122_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_122_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_122.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_122.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_122.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_122.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_122.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_122.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_122.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_122.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_123_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_123_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_123.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_123.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_123.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_123.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_123.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_123.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_123.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_123.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_124_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_124_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_124.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_124.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_124.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_124.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_124.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_124.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_124.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_124.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_125_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_125_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_125.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_125.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_125.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_125.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_125.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_125.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_125.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_125.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_126_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_126_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_126.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_126.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_126.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_126.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_126.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_126.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_126.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_126.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_127_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_127_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_127.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_127.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_127.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_127.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_127.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_127.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_127.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_127.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_128_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_128_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_128.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_128.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_128.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_128.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_128.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_128.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_128.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_128.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_129_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_129_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_129.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_129.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_129.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_129.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_129.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_129.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_129.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_129.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_130_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_130_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_130.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_130.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_130.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_130.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_130.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_130.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_130.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_130.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_131_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_131_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_131.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_131.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_131.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_131.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_131.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_131.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_131.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_131.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_132_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_132_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_132.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_132.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_132.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_132.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_132.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_132.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_132.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_132.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_133_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_133_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_133.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_133.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_133.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_133.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_133.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_133.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_133.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_133.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_134_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_134_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_134.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_134.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_134.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_134.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_134.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_134.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_134.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_134.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_135_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_135_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_135.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_135.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_135.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_135.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_135.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_135.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_135.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_135.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_136_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_136_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_136.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_136.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_136.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_136.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_136.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_136.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_136.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_136.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_137_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_137_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_137.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_137.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_137.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_137.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_137.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_137.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_137.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_137.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_138_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_138_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_138.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_138.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_138.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_138.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_138.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_138.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_138.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_138.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_139_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_139_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_139.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_139.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_139.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_139.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_139.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_139.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_139.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_139.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_140_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_140_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_140.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_140.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_140.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_140.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_140.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_140.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_140.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_140.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_141_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_141_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_141.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_141.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_141.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_141.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_141.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_141.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_141.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_141.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_142_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_142_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_142.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_142.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_142.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_142.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_142.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_142.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_142.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_142.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_143_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_143_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_143.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_143.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_143.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_143.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_143.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_143.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_143.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_143.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_144_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_144_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_144.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_144.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_144.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_144.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_144.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_144.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_144.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_144.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_145_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_145_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_145.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_145.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_145.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_145.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_145.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_145.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_145.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_145.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_146_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_146_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_146.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_146.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_146.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_146.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_146.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_146.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_146.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_146.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_147_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_147_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_147.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_147.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_147.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_147.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_147.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_147.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_147.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_147.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_148_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_148_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_148.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_148.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_148.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_148.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_148.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_148.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_148.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_148.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_149_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_149_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_149.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_149.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_149.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_149.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_149.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_149.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_149.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_149.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_150_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_150_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_150.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_150.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_150.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_150.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_150.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_150.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_150.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_150.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_151_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_151_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_151.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_151.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_151.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_151.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_151.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_151.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_151.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_151.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_152_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_152_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_152.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_152.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_152.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_152.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_152.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_152.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_152.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_152.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_153_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_153_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_153.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_153.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_153.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_153.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_153.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_153.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_153.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_153.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_154_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_154_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_154.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_154.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_154.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_154.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_154.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_154.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_154.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_154.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_155_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_155_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_155.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_155.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_155.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_155.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_155.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_155.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_155.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_155.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_156_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_156_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_156.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_156.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_156.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_156.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_156.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_156.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_156.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_156.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_157_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_157_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_157.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_157.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_157.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_157.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_157.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_157.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_157.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_157.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_158_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_158_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_158.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_158.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_158.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_158.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_158.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_158.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_158.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_158.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_159_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_159_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_159.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_159.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_159.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_159.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_159.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_159.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_159.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_159.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_160_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_160_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_160.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_160.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_160.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_160.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_160.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_160.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_160.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_160.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_161_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_161_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_161.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_161.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_161.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_161.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_161.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_161.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_161.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_161.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_162_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_162_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_162.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_162.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_162.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_162.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_162.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_162.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_162.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_162.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_163_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_163_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_163.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_163.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_163.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_163.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_163.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_163.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_163.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_163.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_164_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_164_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_164.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_164.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_164.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_164.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_164.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_164.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_164.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_164.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_165_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_165_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_165.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_165.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_165.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_165.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_165.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_165.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_165.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_165.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_166_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_166_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_166.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_166.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_166.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_166.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_166.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_166.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_166.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_166.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_167_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_167_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_167.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_167.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_167.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_167.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_167.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_167.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_167.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_167.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_168_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_168_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_168.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_168.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_168.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_168.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_168.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_168.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_168.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_168.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_169_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_169_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_169.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_169.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_169.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_169.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_169.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_169.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_169.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_169.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_170_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_170_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_170.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_170.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_170.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_170.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_170.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_170.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_170.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_170.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_171_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_171_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_171.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_171.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_171.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_171.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_171.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_171.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_171.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_171.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_172_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_172_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_172.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_172.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_172.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_172.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_172.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_172.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_172.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_172.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_173_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_173_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_173.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_173.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_173.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_173.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_173.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_173.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_173.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_173.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_174_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_174_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_174.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_174.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_174.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_174.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_174.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_174.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_174.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_174.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_175_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_175_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_175.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_175.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_175.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_175.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_175.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_175.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_175.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_175.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_176_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_176_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_176.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_176.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_176.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_176.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_176.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_176.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_176.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_176.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_177_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_177_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_177.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_177.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_177.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_177.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_177.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_177.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_177.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_177.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_178_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_178_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_178.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_178.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_178.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_178.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_178.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_178.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_178.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_178.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_179_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_179_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_179.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_179.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_179.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_179.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_179.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_179.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_179.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_179.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_180_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_180_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_180.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_180.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_180.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_180.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_180.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_180.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_180.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_180.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_181_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_181_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_181.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_181.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_181.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_181.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_181.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_181.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_181.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_181.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_182_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_182_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_182.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_182.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_182.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_182.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_182.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_182.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_182.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_182.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_183_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_183_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_183.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_183.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_183.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_183.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_183.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_183.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_183.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_183.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_184_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_184_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_184.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_184.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_184.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_184.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_184.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_184.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_184.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_184.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_185_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_185_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_185.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_185.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_185.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_185.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_185.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_185.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_185.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_185.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_186_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_186_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_186.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_186.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_186.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_186.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_186.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_186.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_186.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_186.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_187_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_187_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_187.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_187.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_187.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_187.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_187.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_187.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_187.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_187.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_188_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_188_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_188.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_188.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_188.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_188.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_188.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_188.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_188.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_188.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_189_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_189_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_189.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_189.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_189.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_189.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_189.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_189.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_189.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_189.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_190_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_190_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_190.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_190.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_190.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_190.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_190.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_190.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_190.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_190.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_191_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_191_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_191.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_191.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_191.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_191.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_191.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_191.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_191.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_191.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_192_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_192_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_192.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_192.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_192.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_192.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_192.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_192.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_192.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_192.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_193_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_193_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_193.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_193.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_193.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_193.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_193.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_193.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_193.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_193.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_194_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_194_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_194.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_194.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_194.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_194.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_194.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_194.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_194.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_194.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_195_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_195_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_195.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_195.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_195.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_195.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_195.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_195.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_195.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_195.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_196_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_196_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_196.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_196.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_196.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_196.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_196.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_196.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_196.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_196.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_197_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_197_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_197.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_197.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_197.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_197.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_197.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_197.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_197.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_197.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_198_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_198_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_198.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_198.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_198.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_198.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_198.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_198.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_198.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_198.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_199_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_199_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_199.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_199.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_199.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_199.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_199.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_199.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_199.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_199.tagcount1over", row.names = F))
in_f <- c("path/RAD-seq/NGS_raw_data_folder/F2_200_1.fq.gz", "path/RAD-seq/NGS_raw_data_folder/F2_200_2.fq.gz")
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_200.fastq"
param_trim <- 0
fastq <- readFastq(in_f)              
hoge <- width(sread(fastq)) - param_trim
hoge[hoge < 1] <- 1                   
hoge1 <- DNAStringSet(sread(fastq), start=1, end=hoge)
hoge2 <- BStringSet(quality(quality(fastq)), start=1, end=hoge)
fastq <- ShortReadQ(hoge1, hoge2, id(fastq))
writeFastq(fastq, out_f, compress=F)   
fastq <- readDNAStringSet(out_f, format="fastq")
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_200.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
param_adapter <- "CATGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 1
system.time(res1 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
param_adapter <- "TAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 2
system.time(res2 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
param_adapter <- "TAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  #Adapter sequence 3
system.time(res3 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
in_f <- "path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq"
out_f <- "path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_200.fastq"
param_adapter <- "CATGAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"  #Adapter sequence 4
system.time(res4 <- preprocessReads(filename=in_f,  outputFilename=out_f,Rpattern=param_adapter, max.Rmismatch=rep(0, nchar(param_adapter)), nBases=0, minLength=20, nrec=1000000))      
res <- rbind(res1,res2,res3,res4)
write.csv(res, "path/RAD-seq/NGS_Data_Analysis_folder/F2_200.fastq.csv")
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-1.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-2.fastq","path/RAD-seq/NGS_Data_Analysis_folder/F2_Adap-3.fastq")
system.time(fastq <- fread("path/RAD-seq/NGS_Data_Analysis_folder/AdapterTrim_F2_200.fastq", header=F, sep="\t", na.strings="NA", strip.white=TRUE, verbose=TRUE))
system.time(nucleotide <- fastq[2+4*(0:((nrow(fastq)/4)-1))])
system.time(RADtag <- sort(table(nucleotide), decreasing=TRUE))
system.time(write.csv(RADtag, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_200.tagcount", row.names = F))
file.remove("path/RAD-seq/NGS_Data_Analysis_folder/F2_200.fastq")
RAD_1over <- RADtag[RADtag>1]
system.time(write.csv(RAD_1over, file="path/RAD-seq/NGS_Data_Analysis_folder/F2_200.tagcount1over", row.names = F))
