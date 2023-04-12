#!/usr/bin/env Rscript

# Author: David Fournier - January 2018
# generation of wig tracks

###################################################
#
# syntax:
#
# Rscript wig.R $bamlistfile $Nb_cores >wig.Rout 2>wig.Rerr
#
################################################### 

library(QuasR)
library(foreach)
library(doParallel)

args = commandArgs(trailingOnly=TRUE)

cores=detectCores(); cl=makeCluster(as.numeric(args[2]))  # ncores threads 
registerDoParallel(cl)

# Input: bam files
bams=read.table(args[1], header=T, sep="\t")
ii=which(bams$SampleName=="Input")
is=which(bams$SampleName!="Input")

#bamf=qAlign(args[1], "BSgenome.Hsapiens.UCSC.hg19", checkOnly=T, paired="no") # genome does not matter
#print(args)
if(args[2]>1){
  bamfl=list()
  for(i in is){
   write.table(bams[c(i,ii),], "file.tmp", sep="\t",quote=F, row.names=F)
   bamf=qAlign("file.tmp", "BSgenome.Hsapiens.UCSC.hg19", checkOnly=T, paired="no") # genome does not matter
   bamfl=c(bamfl, bamf)
  }
  unlink("file.tmp")
}
bamf=qAlign(args[1], "BSgenome.Hsapiens.UCSC.hg19", checkOnly=T, paired="no") # genome does not matter

if(args[2]==1){ # no parallelization
	for (i in 1:length(bamfl)){
		qExportWig(bamfl[[i]],  binsize=100L, scaling=TRUE, collapseBySample=TRUE)
	}
}else{ # parallelization
 foreach(i=1:length(bamf), .packages='QuasR') %dopar% qExportWig(bamfl[[i]],  binsize=100L, scaling=TRUE, collapseBySample=TRUE)
}

#

