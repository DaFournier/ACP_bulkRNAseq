#!/bin/bash

######################################################
#
# Pipeline NextGenBender for ChIP-seq processing
# Version without scheduler
# By David Fournier - January 2018
#
# programs to install: check file installation.sh
#
######################################################


#   Requirements
#   ============
#
# - Need enough space to content all the files plus 17G
#   for the macs2 step that generates temp files in large amounts -- up to 17G
# - Do not use sample names with dashes (symbol "-") in them. bowtie & R QuasR
#   package do not process files with such names correctly.


#   Input file format
#   =================
#
# Data to be downloaded or imported should be placed
# in the file fastq_files.csv with the following , delimited format:
#
# SampleName,SRAid,URL
#
# - The URL is optional as fastq-dump can import data with the SRA id only.
# - SRAid with the same SampleName will be considered as replicates to merge
#   whether they are biological or technical. If you want to separate biological
#   replicates, give them different values of the SampleName field.
# - Input/control file should be referred to as "Input" with uppercase I.
# - Example:
#
# H3K27me3,,ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR001/SRR001363/SRR001363.sra
# H3K4me3,SRR001431,
# H3K4me3,SRR001432,
# H3K9me3_1,SRR001424,
# H3K9me3_2,SRR001425,
# Input,SRR001426,
#
# will generate one final bam file and one final peak file for both H3K27me3
# and H3K4me3, a technical replicate, while H3K9me3_1 and H3K9me3_2, two
# biological replicates whose separate information is crucial to the analysis,
# will be kept separated.
# H3K27me3 will be downloaded from the URL mentioned above, the other samples
# will be downloaded based on SRA id using the fastq-dump tool. 
# SRR001426 is a control experiment that will be treated as such in
# the enrichment calculation step at the end of the pipeline, providing
# that parameter "enrich" is set to 1. 

#   Output files
#   ============
#
# The pipeline generates bam and bed files for each separate sample,
# id est which is not merged. On top of that, two files summarize the
# data and can be used for subsequent analysis:
#  - master.bed contains the peaks from all experiments, the overlapping
#   ones by more than 50% being agglomerated;
#  - master.tab is a large table containing the enrichments for all the
#   features (but the input) at the different peaks. This table can be
#   used immediately to perform genome analysis, such as correlation b.
#   marks, chromatin states definition or quantification at genomic regions. 

# Parameters

download=1 # switch step on/off
fastqc=0
align=0
alignqc=0
wig=0
peaks=0
master=0
enrich=0
corr=0

datafile="fastq_files.csv"
bamfiles="bamfiles.txt"
#datadir="../bam_files_epig_ESC/epigenetic_code"
datadir="."
bamdir="."
genome="hs" # hs: Homo sapiens; mm: Mus musculus
index=/path/to/index/folder/hg19/hg19 # index location-- points at the basename of index files (e.g.here: hg19)
PARALLEL_NPROCS=8

# 1 -- Download data

if [ "$download" -eq 1 ]; then
    
    snames=$(awk -F ',' '{print $1}' $datafile)
    SRA=$(awk -F ',' '{print match($2, /[^ ]/) ? $2 : "empty"}' $datafile)
    URL=$(awk -F ',' '{print match($3, /[^ ]/) ? $3 : "empty"}' $datafile)
    i=1
    for i in $snames
    do
	if [ "${SRA[$i]}" == "empty" ] && [ "${URL[$i]}" != "empty" ]; then
	    wget ${URL[$i]} > ${snames[$i]}.sra
	    fastq-dump ${snames[$i]}.sra
	    rm ${snames[$i]}.sra
	fi
	if [ "${SRA[$i]}" != "empty" ]; then
	    #wget ${URL[$i]} > ${snames[$i]}.sra
	    fastq-dump ${snames[$i]}.sra
	    rm ${snames[$i]}.sra
	fi
   
    done
        
fi

# 2 -- Raw files QC

if [ "$fastqc" -eq 1 ]; then

    PARALLEL_NPROCS=20
    
    if [ -e "*.fastq" ];then 
	seq $PARALLEL_NPROCS | parallel -n0 fastqc *.fastq > /dev/null
    fi
    if [ -e "*.fastq.gz" ];then
	seq $PARALLEL_NPROCS | parallel -n0 fastqc <(gunzip -c '*.fastq.gz') > /dev/null
    fi

    text='Sample\tBasic_stats\tPer_base_seq_quality\tPer_tile_seq\tPer_seq_quality_scores\tPer_base_seq_content\tPer_seq_GC_content\tPer_base_N_content\tSeq_Length_Distri\tSeq_Duplication_Levels\tOverrepresented_seqs\tAdapter_Content\tKmer_Content\tStatus'
    echo -e $text > QCreport.txt

    for f in *.fastq 
    do
	filename=$(basename "$f" .fastq); fn="${filename%.fastq}"
	#echo $f
	tests="$(unzip -cq ${fn}_fastqc.zip ${fn}_fastqc/summary.txt | cat | awk '{print $1}')"
	set -- $tests
        t='\t' ; echo -en ${fn}$t$tests >> QCreport.txt
    done
    
fi

# 3 -- Aligning files

if [ "$align" -eq 1 ]; then
    
    PARALLEL_NPROCS=20 # for GNU parallel
    $nbt=$PARALLEL_NPROCS/10 # GNU parallel bowtie instances (10 threads/instance)

    if [ -e "*.fastq" ];then
	(parallel -t -j $nbt bowtie2 -p 10 -x $index <(find . -name '*.fastq') | samtools view -bS -o *.bam - ) 2> align.err
    fi
    if [ -e "*.fastq.gz" ];then 
	(parallel -t -j $nbt bowtie2 -p 10 -x $index <(gunzip -c '*.fastq.gz')  | samtools view -bS -o *.bam - ) 2>> align.err
    fi
     
    #    find . -name '*.fastq' | parallel -t -j $nbt "bowtie -p 10 -x $index {/}" | samtools view -Sb  {}  {}.bam  #old
fi

# 4 -- Alignment QC

if [ "$alignqc" -eq 1 ]; then
    for i in *.bam
    do
	echo -en ${i} >> QCreport.txt
	samtools stats $i | grep "\%" >> QCreport.txt   
	echo -en  >> QCreport.txt
	
    done
fi
    
if [ "$align" -eq 1 ]; then
    
    if [ -e "$datafile" ];then  # merging for identical SampleNames
    
	snames=$(awk -F ',' '{print $1}' $datafile)
	snames2=$(awk -F ',' '{print $1}' $datafile | uniq)
	
	for i in $snames2
	do
	    samtools merge ${names2}_merged.bam <(ls *.bam | grep ${names2})      
	done
	
    fi
    
fi

# 5 -- Generation of wig files

if [ "$wig" -eq 1 ] || [ "$peaks" -eq 1 ]; then  #indexing bam files
    find . -name '*.bam' | parallel -t -j $PARALLEL_NPROCS  "samtools index {/} " 
fi

if [ "$wig" -eq 1 ]; then
    
    PARALLEL_NPROCS=40 # for GNU parallel
    Rscript wig.R $bamfiles $PARALLEL_NPROCS >wig.Rout 2>wig.Rerr
fi

# 6 -- Peak calling

if [ "$peaks" -eq 1 ]; then
    
    PARALLEL_NPROCS=40 # for GNU parallel
    ls -1 *bam | sort | sed -r 's/.bam//g' | sort | uniq > tmp.txt
    cat tmp.txt | parallel --max-procs=$PARALLEL_NPROCS 'macs2 callpeak -t {}.bam -g hs -n {}_peaks -q 0.01 --tempdir . 2> {}.stderr'
fi

# 7 -- generation of master.bed file

if [ "$master" -eq 1 ]; then
    
    > $datadir/master.bed
    for i in $datadir/*.narrowPeak
    do
	if [ -e $datadir/master.bed ]; then
	    
	    filename=$(basename "$i" _peaks.narrowPeak);	
	    cut -f 1-5 $datadir/${filename}_peaks.narrowPeak > $datadir/tmp_peaksfile.bed
	    intersectBed -a $datadir/tmp_peaksfile.bed -b $datadir/master.bed -f 0.5 -r -wa > $datadir/overlapping_peaks.bed # peaks w. 50% overlap w. master.bed
	    subtractBed -a $datadir/tmp_peaksfile.bed -b $datadir/overlapping_peaks.bed -f 1 -r > $datadir/new_peaks.bed # peaks from tmp_peaksfile.bed not in overlapping_peaks.bed -- only new peaks are present
	    cat $datadir/new_peaks.bed >> $datadir/master.bed
	    sortBed -i $datadir/master.bed > $datadir/out.log
	    rm $datadir/new_peaks.bed $datadir/overlapping_peaks.bed
	else
	    cat $datadir/${filename}_peaks.narrowPeak > $datadir/master.bed
	fi
    done
    
fi
	
# 8 -- Normalization and enrichment above input

if [ "$enrich" -eq 1 ]; then

    echo -e 'FileName\tSampleName' > $datadir/bam.tmp    
    awk -F ',' '{print $2 "\t" $1}' $datafile >> $datadir/bam.tmp
    bamfile=$datadir/bam.tmp
    peakfile=$datadir/master.bed
    outfile=$datadir/master.tab    
    Rscript NormEnrichments_parallel.R $bamfile $peakfile $PARALLEL_NPROCS $outfile >NormEnrichments_parallel.Rout 2>NormEnrichments_parallel.Rerr
    rm bam.tmp
fi
