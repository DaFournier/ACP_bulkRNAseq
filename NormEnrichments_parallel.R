#!/usr/bin/env Rscript

# Author: David Fournier - December 2017
# A script to count the reads from bam files falling into peaks 
# Two inputs: bam file list and peak list. 

###########################################################################################
#
# syntax:
#
# R CMD BATCH '--args bamlist peakfile ncores outfile' NormEnrichments_parallel.R /dev/null &
#
# syntax with output file ReadCounts.Rout:
#
# R CMD BATCH '--args bamlist peakfile ncores outfile' NormEnrichments_parallel.R &

## instructions

#
# - bamlist: a text file containing the list of bam files to count reads from. 
#   - header should contain words "Filename" and "SampleName" separated by a tabulation
#   - other lines should contain the name of each bam file and the name that will
#     appear in the final enrichment table, these two fiels being separated by a tabulation. 
#   - one line should contain the ChIP-seq associated to the input of this series. 
#     In this case, the name of the sample has to be "Input" with a upper-case I. 
#   - example: 
#
#FileName	SampleName 
#chip_1.bam	Sample1    
#chip_2.bam	Sample2   
#input_filename.bam	Input


#  note regarding replicates: replicates have to be merged prior to launch this script. One bam file per feature.

# - peakfile: put a file name containing the peak list in a bed format. 

# - ncores: number of cores used to generate the count tables. 1 means no parallelization. 
 
#############################################################################################



library(QuasR)
library(foreach)
library(doParallel)

args = commandArgs(trailingOnly=TRUE)
#args = c("bamFiles.txt","master.bed",7,"Peak_counts2.txt") # for debug

cores=detectCores(); cl=makeCluster(as.numeric(args[3]))  # ncores threads 
registerDoParallel(cl)

# Input: bam files
bams=read.table(args[1], header=T, sep="\t")
ii=which(bams$SampleName=="Input")
is=which(bams$SampleName!="Input")

if(args[3]>1){
  bamfl=list()
  for(i in is){
   write.table(bams[c(i,ii),], "file.tmp", sep="\t",quote=F, row.names=F)
   bamf=qAlign("file.tmp", "BSgenome.Hsapiens.UCSC.hg19", checkOnly=T, paired="no") # genome does not matter
   bamfl=c(bamfl, bamf)
  }
  unlink("file.tmp")
}
bamf=qAlign(args[1], "BSgenome.Hsapiens.UCSC.hg19", checkOnly=T, paired="no") # genome does not matter
stat<-alignmentStats(bamf)
snames=bamf@alignments$SampleName

# Input: peaks list
peaks=read.table(args[2], h=F, sep="\t")

# Read counts 
 
grpk=GRanges(seqnames=as.factor(peaks$V1), range=IRanges(peaks$V2, peaks$V3), peakID=peaks$V4)
if(args[3]==1){ # no parallelization
 counts=as.data.frame(qCount(bamf, grpk)) # counts reads from bam files in the different peaks grpk
}else{ # parallelization
 counts=foreach(i=1:length(bamfl), .combine='cbind', .packages='QuasR') %dopar% qCount(bamfl[[i]], grpk)
 counts=as.data.frame(counts[,c(1,3*c(1:length(bamfl))-1,3)]) #gets rid of redundant cols
} 
counts[counts==0]=1 # pseudocount

# Normalization-Enrichment calculation

write.table(counts, "master_tmp.tab", sep="\t",quote=F, row.names=F)


li=stat[ii,2]

for(i in is){ 
    ls=stat[i,2]
    sname=snames[i]
    eval(parse(text=paste0("counts","$",sname,"_norm=(counts$",sname,"/",ls,")","/(counts$Input/",li,")")))#debug
}	


# Output file

counts=cbind(peaks[,1:3],counts)
write.table(counts, args[4], sep="\t",quote=F, row.names=F)

