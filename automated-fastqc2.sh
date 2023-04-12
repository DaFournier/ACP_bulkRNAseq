#!/bin/sh


# variables to set at beginning of the pipeline

PARALLEL_NPROCS=6 # for GNU parallel
QUALITY_SCORE=10 # for trim_galore
ADAPTER_SEQ="" # specific sequence adapter to trim
TRIM=1 # 1: trimming if low quality or adapter bias 0: no trimming
MASK=1 # 1: mask if tile bias only (no quality bias) 0: no masking

# more to does on the pipeline: 

if the diagnosis turns to be zero, the pipeline has to be interrupted with 
a message: "please check the low-quality samples reported in report.txt and
eventually do not include them when you launch the pipeline again. "

# variables of this script

diagnosis=1


# creates study report

text='Sample\tBasic_stats\tPer_base_seq_quality\tPer_tile_seq\tPer_seq_quality_scores\tPer_base_seq_content\tPer_seq_GC_content\tPer_base_N_content\tSeq_Length_Distri\tSeq_Duplication_Levels\tOverrepresented_seqs\tAdapter_Content\tKmer_Content\tStatus'
echo -e $text > report.txt


# fastqc

seq $PARALLEL_NPROCS | parallel -n0 fastqc *.fastq > /dev/null


if [ $TRIM == 1 ]; then

  for f in *.fastq 
  do
      filename=$(basename "$f" .fastq); fn="${filename%.fastq}"
      #echo $f
      tests="$(unzip -cq ${fn}_fastqc.zip ${fn}_fastqc/summary.txt | cat | awk '{print $1}')"
      set -- $tests
      if [ $2 == "FAIL" ] && [ $3 == "PASS" ]; then #combined adapter and quality trimming
         t='\t' ; echo -en ${fn}_old$t$tests >> report.txt
         if [ -z "$ADAPTER_SEQ" ]; then
          seq $PARALLEL_NPROCS | parallel -n0 trim_galore -q $QUALITY_SCORE --length 15 --phred33 $f > /dev/null  
        else
          seq $PARALLEL_NPROCS | parallel -n0 trim_galore -q $QUALITY_SCORE -a $ADAPTER_SEQ --length 15 --phred33 $f > /dev/null  #custom seq
        fi
        fastqc ${fn}_trimmed.fq > /dev/null

        # report the result
        tests="$(unzip -cq ${fn}_trimmed_fastqc.zip ${fn}_trimmed_fastqc/summary.txt | cat | awk '{print $1}')"
        set -- $tests
        tmp=$2
#SRR001449_H3K4me2_trimmed_fastqc.zip
        seqs="$(unzip -cq ${fn}_trimmed_fastqc.zip ${fn}_trimmed_fastqc/fastqc_data.txt | cat | grep Total\ Sequences)"
        set -- $seqs
        if [ $tmp != "FAIL" ] && [ $3 > 1000000 ]; then
     #   if [ $tmp == "FAIL" ] && [ $3 > 1000000 ]; then
          printf "\t%s" "Efficient trimming" >> report.txt
          t='\t' ; echo -en ${fn}$t$tests >> report.txt  # tests for trimmed seq
          mv ${fn}.fastq ${fn}_old.fastq
          mv ${fn}_trimmed.fq ${fn}.fastq # trimmed files become the new ones       
        else
          rm ${fn}_trimmed.fq
          diagnosis=0  # will stop the pipeline 
        fi
      elif [ $3 == "FAIL" ]; then
        # do the masking
        #fastqc *.fastq > /dev/null
        # report the result
      elif [ $2 == "PASS" ] && [ $3 == "PASS" ]; then 
        printf "\t%s" "QC_PASSED" >> report.txt
      fi
  done
  #seq $PARALLEL_NPROCS | parallel -n0 fastqc *.fastq > /dev/null
  rm *_fastqc.html _fastqc.zip *_trimming_report.txt

else # if no trimming, simple report with advise
  for f in *.fastq # debug
  do
    filename=$(basename "$f"); fn="${filename%.fastq}"
    tests="$(unzip -cq ${fn}_fastqc.zip ${fn}_fastqc/summary.txt | cat | awk '{print $1}')"
    t='\t' ; n='\n'; echo -en $n$fn.fastq$t$tests >> report.txt
    set -- $tests
    if [ $2 == "FAIL" ] && [ $3 == "PASS" ]; then
      printf "\t%s" "Trimming advised" >> report.txt
    elif [ $3 == "FAIL" ]; then
      printf "\t%s" "Masking advised" >> report.txt
    fi
  done
fi

# for recollection: the different tests done by fastqc
#Sample	Basic_stats 1	Per_base_seq_quality 2	Per_tile_seq 3	Per_seq_quality_scores 4	Per_base_seq_content 5	Per_seq_GC_content 6	Per_base_N_content 7	Seq_Length_Distri 8	Seq_Duplication_Levels 9	Overrepresented_seqs 10	Adapter_Content 11	Kmer_Content 12


