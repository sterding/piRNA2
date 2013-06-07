# ===============================================================
# Pipeline for calculating binned density of histone modification (and any other features in BigWig format)
# cluster manner
# Author: Xianjun.Dong @ umassmed.edu
# Date: Thu Mar 10
# Usage:
# ./bigWig_81bins.sh ../data/controlSet2/controlSet2.alluniq.mRNA.81bins
# ./bigWig_81bins.sh ../data/controlSet2/controlSet2.alluniq.genome.81bins 
# Requirement:
# Step1: cat ../data/controlSet2/*.bed* | sort -u | awk '{OFS="\t"; $4=$4"__"$5; $5=0;$9=0; print;}' | awk -v option=genome -f bigWigAverageOverBed_generate_81bins.awk > ../data/controlSet2/controlSet2.alluniq.genome.81bins
# Step2: prepare the bigwig file (e.g. Download BIGWIG FILES: encode-box-01@fasp.encode.ebi.ac.uk:byDataType/signal/jan2011/bigwig/ or ftp://ftp-private.ebi.ac.uk/byDataType/signal/jan2011/bigwig Login:encode-box-01 Passwd: enc*deDOWN
# or scp from Jessie's folder at: hpcc03:/isilon/weng/Jie/bw_from_Anshul_JanFreeze/)

## CHANGE: use cut -f1,5 to get mean0 (instead of mean) from bigWigAverageOverBed- average over bases with non-covered bases counting as zeroes
# ===============================================================


#!/bin/sh

# for mouse transcriptome data
#array=(`ls /home/dongx/scratch/ucsc/mouse_adult_wt_smallRNAseq_76nt_strand_ox*.bw /home/dongx/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand_R1*.bw /home/dongx/scratch/ucsc/mouse_adult_wt_CAGE_*unique*.bw /home/dongx/scratch/ucsc/mouse_adult_wt_PAS*unique*.bw /home/dongx/nearline/Xin/others/mmENCODE/wgEncodeCshlLongRnaSeqTestisAdult8wks*Rep1.bigWig`)
# for mouse chipseq data
#array=(/home/wangj2/scratch/bill/mouse/PolIII_testes_1900_mmu0SQ121.bw /home/dongx/scratch/ucsc/mouse_adult_wt_ChIP_Amyb_50nt.accepted_hits.bw /home/dongx/scratch/ucsc/mouse_adult_wt_H3K4me3_50nt.accepted_hits.bw /home/dongx/scratch/ucsc/mouse_adult_wt_Pol2_50nt.accepted_hits.bw)
array=(/home/wangj2/scratch/bill/mouse/PolIII_testes_1900_mmu0SQ121.bw)

inputfile=$1 #/home/dongx/projects/piRNA/data/controlSet2/controlSet2.alluniq.genome.81bins
cut -f4 -d' ' $inputfile | sed 's/.bin.*//g' | sort -u > /tmp/header_`basename $inputfile`

JOBOUTPUT_DIR=/home/dongx/projects/piRNA/data/for_aggregationplot
[ -d $JOBOUTPUT_DIR ] || mkdir -p $JOBOUTPUT_DIR

for filename in "${array[@]}"
do
    JOBOUTPUT=$JOBOUTPUT_DIR/output_`basename $inputfile`_of_`basename $filename`
    echo $filename
    bigWigAverageOverBed $filename $inputfile stdout | cut -f1,5 | sed 's/.bin/\t/g' | sort -k1,1 -k2,2n | awk '{printf $3" "; if(NR%81==0) printf $1"\n";}' | sort -k82,82 | cut -f1-81 -d' ' | paste /tmp/header_`basename $inputfile` - > $JOBOUTPUT &
done
