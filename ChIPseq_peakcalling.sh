## interactive-version of pipeline for running ChIP-seq peak calling etc.
## Author: Xianjun Dong (xianjun.dong@umassmed.edu)
## Date: 2012-July-10
## Version: 0.1
## Usage:
#qsub ChIPseq_peakcalling.sh ~/scratch/chicken_adult_wt_ChIP_Amyb_50nt/accepted_hits.bam ~/scratch/chicken_adult_wt_InputIP_50nt/accepted_hits.bam gg chicken_adult_wt_ChIP_Amyb_50nt TF
#qsub ChIPseq_peakcalling.sh ~/scratch/rat_adult_wt_ChIP_Amyb_50nt/accepted_hits.bam ~/scratch/rat_adult_wt_InputIP_50nt/accepted_hits.bam rn rat_adult_wt_ChIP_Amyb_50nt TF
#qsub ChIPseq_peakcalling.sh ~/scratch/frog_adult_wt_ChIP_Amyb_50nt/accepted_hits.bam ~/scratch/frog_adult_wt_InputIP_50nt/accepted_hits.bam xt frog_adult_wt_ChIP_Amyb_50nt TF
#qsub ChIPseq_peakcalling.sh ~/scratch/mouse_adult_wt_H3K4me3_50nt/accepted_hits.bam ~/scratch/mouse_adult_wt_InputIP_50nt/accepted_hits.bam mm mouse_adult_wt_H3K4me3_50nt HM
#qsub ChIPseq_peakcalling.sh ~/scratch/chicken_adult_wt_H3K4me3_50nt/accepted_hits.bam ~/scratch/chicken_adult_wt_InputIP_50nt/accepted_hits.bam gg chicken_adult_wt_H3K4me3_50nt HM
#qsub ChIPseq_peakcalling.sh ~/scratch/frog_adult_wt_H3K4me3_50nt/accepted_hits.bam ~/scratch/frog_adult_wt_InputIP_50nt/accepted_hits.bam xt frog_adult_wt_H3K4me3_50nt HM
#qsub ChIPseq_peakcalling.sh ~/scratch/rat_adult_wt_H3K4me3_50nt/accepted_hits.bam ~/scratch/rat_adult_wt_InputIP_50nt/accepted_hits.bam rn rat_adult_wt_H3K4me3_50nt HM

#!/bin/sh
#$ -V
#$ -pe single 4
#$ -cwd
#$ -o $HOME/sge_jobs_output/sge_job.$JOB_ID.out -j y
#$ -S /bin/bash
#$ -M $EMAIL
#$ -S /bin/bash
#$ -m e
#$ -l mem_free=2G

###########################################
############## 1. Configuring
###########################################

treatment=$1
input=$2

# genome size data from: http://useast.ensembl.org/Xenopus_tropicalis/Info/StatsTable?db=core
genomeSize='hs'
case $3 in
        hs) genomeSize='hs' ;;
	mm) genomeSize='mm' ;;
	gg) genomeSize='1.11e9' ;;
	rn) genomeSize='2.72e9' ;;
	xt) genomeSize='1.52e9' ;;
	*) echo "INVALID NUMBER!"; exit ;;
esac

sampleName=$4
sampleType=$5 #'TF' or 'HM'

###########################################
###############  2. call peaks
###########################################

if [ $sampleType == 'TF' ]; then
    macs14 -t $treatment -c $input -g $genomeSize -w -n $sampleName --space=10
fi
if [ $sampleType == 'HM' ]; then
    macs14 -t $treatment -c $input -g $genomeSize -w -n $sampleName --space=10 --nomodel
fi

R --vanilla < ${sampleName}_model.r

###########################################
###############  3. UCSC display
###########################################

echo "track type=narrowPeak name=$sampleName description=$samepleName useScore=1" > $HOME/scratch/ucsc/${sampleName}_peaks.encodePeak
cat ${sampleName}_peaks.encodePeak >> $HOME/scratch/ucsc/${sampleName}_peaks.encodePeak
cp ${sampleName}_peaks.xls $HOME/scratch/ucsc/
