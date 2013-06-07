# ===============================================================
# Pipeline for calculating binned density of histone modification (and any other features in BigWig format)
# cluster manner
# Author: Xianjun.Dong @ umassmed.edu
# Date: Thu Mar 10
# Usage:
# qsub bigWigAverageOverBed_81bins.sh ../data/controlSet3/controlSet3.alluniq.mRNA.81bins 81
# qsub bigWigAverageOverBed_81bins.sh ../data/controlSet3/controlSet3.alluniq.genome.81bins 81
# Requirement:
# Step1: cat ../data/controlSet2/*.bed* | sort -u | awk '{OFS="\t"; $4=$4"__"$5; $5=0;$9=0; print;}' | awk -v option=genome -f bigWigAverageOverBed_generate_81bins.awk > ../data/controlSet2/controlSet2.alluniq.genome.81bins
# Step2: prepare the bigwig file (e.g. Download BIGWIG FILES: encode-box-01@fasp.encode.ebi.ac.uk:byDataType/signal/jan2011/bigwig/ or ftp://ftp-private.ebi.ac.uk/byDataType/signal/jan2011/bigwig Login:encode-box-01 Passwd: enc*deDOWN
# or scp from Jessie's folder at: hpcc03:/isilon/weng/Jie/bw_from_Anshul_JanFreeze/)

## CHANGE: use cut -f1,5 to get mean0 (instead of mean) from bigWigAverageOverBed- average over bases with non-covered bases counting as zeroes
# ===============================================================


#!/bin/sh
#$ -V
#$ -cwd
#$ -pe openmpi 2
#$-o $HOME/sge_jobs_output/sge_job.$JOB_ID.out -j y
#$ -S /bin/bash
#$ -l mem_free=2G
#$ -t 1-3
# decided by echo ${#array[@]}

# for mouse transcriptome data
# for Phil's grant
array=(`ls /home/dongx/nearline/Xin/smallRNA/Jia/smallRNA-term.TAP.GCAC.raw.uniqmap.xkxh.norm.bed.gz.*.bw /home/dongx/nearline/Xin/smallRNA/Jia/Phil.SRA.wt.ox.6wk.testes.raw.uniqmap.xkxh.norm.bed.gz.*.bw`)
# for the paper Fig1
#array=(`ls /home/dongx/nearline/Xin/smallRNA/Jia/Phil.SRA.wt.ox.6wk.testes.raw.uniqmap.xkxh.norm.bed.gz.*.bw /home/dongx/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand_R1*.bw /home/dongx/scratch/ucsc/mouse_adult_wt_CAGE_*unique*.bw /home/dongx/scratch/ucsc/mouse_adult_wt_PAS*unique*.bw`)
# for mouse chipseq data
# array=(/home/dongx/scratch/ucsc/mouse_adult_wt_H3K4me3_50nt.accepted_hits.bw /home/dongx/scratch/ucsc/mouse_adult_wt_Pol2_50nt.accepted_hits.bw /home/wangj2/scratch/bill/mouse/PolIII_testes_1900_mmu0SQ121.bw)  # not used anymore
array=(/home/dongx/scratch/mouse_adult_wt_InputIP_50nt/input.mouse.by.chris/H3K4-mouse_R1.gz.m1.sam.sorted.bam.bigwig /home/dongx/scratch/mouse_adult_wt_InputIP_50nt/input.mouse.by.chris/pol2.chip.v3.fastq.sam.sorted.bam.bigWig /home/wangj2/scratch/bill/mouse/PolIII_testes_1900_mmu0SQ121.bw)  # update to unique
# for BIP
#array=(`ls /home/dongx/scratch/ucsc/mouse_adult_wt_smallRNAseq_76nt_strand_ox*.bw /home/dongx/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand_R1*.bw /home/dongx/scratch/ucsc/mouse_adult_wt_CAGE_*unique*.bw /home/dongx/scratch/ucsc/mouse_adult_wt_PAS*unique*.bw /home/dongx/scratch/ucsc/mouse_adult_wt_Degradome*unique*.bw` /home/wangj2/scratch/bill/mouse/PolIII_testes_1900_mmu0SQ121.bw /home/dongx/nearline/Xin/others/mmENCODE/wgEncodeCshlLongRnaSeqTestisAdult8wks*Rep1.bigWig /home/dongx/scratch/ucsc/mouse_adult_wt_ChIP_Amyb_50nt.accepted_hits.bw /home/dongx/scratch/ucsc/mouse_adult_wt_H3K4me3_50nt.accepted_hits.bw /home/dongx/scratch/ucsc/mouse_adult_wt_Pol2_50nt.accepted_hits.bw)

inputfile=$1
option=$2 # bin number: 80 or 41  (80bins: 100bp bins for both TSS and TTS; 41bins: 100bp bins for TSS-flanking, [2k, TTS] as one big bin;)

filename="${array[$SGE_TASK_ID-1]}"  # array is 0-indexed

export JOBOUTPUT_DIR=/home/dongx/projects/piRNA/data/for_aggregationplot
export JOBOUTPUT=$JOBOUTPUT_DIR/output_`basename $inputfile`_of_`basename $filename`
export JOBOUTPUT_HEADER=$HOME/scratch/output_header_`basename $inputfile`_of_`basename $filename`

# Make the directory for the job ID you are running if it does not exist
[ -d $JOBOUTPUT_DIR ] || mkdir -p $JOBOUTPUT_DIR

# 'header' column  # sort by TranID
#awk '{if($1!="chr"){print $6}}' $GENOME/mm9/Annotation/Genes/NCBIM37.biomart67.transcripts.cpg.tab | sort > $JOBOUTPUT_HEADER
#awk '{if($1!="chr"){print $5}}' ../data/piRNA.clusters.coordinates.cpg.bed | sort > $JOBOUTPUT_HEADER
cut -f4 -d' ' $inputfile | sed 's/.bin.*//g' | sort -u > $JOBOUTPUT_HEADER

if [[ "$option" == '80' ]]; then
    bigWigAverageOverBed $filename $inputfile stdout | cut -f1,5 | sed 's/.bin/\t/g' | sort -k1,1 -k2,2n | awk '{printf $3" "; if(NR%80==0) printf $1"\n";}' | sort -k81,81 | cut -f1-80 -d' ' | paste $JOBOUTPUT_HEADER - > $JOBOUTPUT
fi
if [[ "$option" == '41' ]]; then
    bigWigAverageOverBed $filename $inputfile stdout | cut -f1,5 | sed 's/.bin/\t/g' | sort -k1,1 -k2,2n | awk '{printf $3" "; if(NR%41==0) printf $1"\n";}' | sort -k42,42 | cut -f1-41 -d' ' | paste $JOBOUTPUT_HEADER - > $JOBOUTPUT
fi
if [[ "$option" == '81' ]]; then
    bigWigAverageOverBed $filename $inputfile stdout | cut -f1,5 | sed 's/.bin/\t/g' | sort -k1,1 -k2,2n | awk '{printf $3" "; if(NR%81==0) printf $1"\n";}' | sort -k82,82 | cut -f1-81 -d' ' | paste $JOBOUTPUT_HEADER - > $JOBOUTPUT
fi
