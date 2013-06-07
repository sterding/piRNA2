## sge script to merge adult RNAseq bam files and do assembly
## Author: Xianjun Dong (xianjun.dong@umassmed.edu)
## Date: 2012-Oct-2
## Version: 0.2
## Usage: qsub mergeRNAseq.sh

#!/bin/sh
#$ -V
#$ -pe single 8
#$ -cwd
#$ -o $HOME/sge_jobs_output/sge_job.$JOB_ID.out -j y
#$ -S /bin/bash
#$ -l mem_free=5G

export BOWTIE_INDEXES=$GENOME/mm9/Sequence/BowtieIndex/
export ANNOTATION=$GENOME/mm9/Annotation/Genes

cpu=8

[ -d ~/scratch/mouse_adult_wt_RNAseq_PE50nt_strand_merged ] || mkdir ~/scratch/mouse_adult_wt_RNAseq_PE50nt_strand_merged;
cd ~/scratch/mouse_adult_wt_RNAseq_PE50nt_strand_merged

ln -fs $HOME/sge_jobs_output/sge_job.$JOB_ID.out sge.log
#
##  2. merge BAM
#ln -s /home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/accepted_hits.bam r1.bam
#ln -s /home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand_R2/merge_reads/accepted_hits.bam r2.bam
#ln -s /home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand_R3/merge_reads/accepted_hits.bam r3.bam
#samtools merge -f accepted_hits.bam -r -h $GENOME/mm9/Sequence/BowtieIndex/sorted.header.sam r1.bam r2.bam r3.bam
### an equivelent but faster way: samtools cat
##samtools cat -o accepted_hits.bam r1.bam r2.bam r3.bam
##
## 3. generate tracks for UCSC browsing
#samtools sort accepted_hits.bam accepted_hits.sorted
#mv accepted_hits.sorted.bam accepted_hits.bam
#samtools index accepted_hits.bam
bam2bw accepted_hits.bam

cp ~/scratch/mouse_adult_wt_RNAseq_PE50nt_strand_merged/*.bam* ~/scratch/mouse_adult_wt_RNAseq_PE50nt_strand_merged/*s.bw ~/scratch/ucsc

###########################################
echo "################# assembly"
###########################################

## run cufflinks to get FPKM
cufflinks -v --no-update-check --frag-len-mean 400 --frag-len-std-dev 100 --library-type fr-firststrand -o cufflinks_o250 -p $cpu -g $ANNOTATION/genes.gtf -M $ANNOTATION/chrM.rRNA.tRNA.gtf accepted_hits.bam -u -j 0.2 --min-frags-per-transfrag 40 --overlap-radius 250 2&> cufflinks_o250.run.log

# gtf of assembly
echo "track name=mouse_adult_wt_RNAseq_PE50nt_strand_merged_o250 description=mouse_adult_wt_RNAseq_PE50nt_strand_merged_o250 visibility=pack colorByStrand='200,100,0 0,100,200'" > $HOME/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand_merged.cufflinks_o250.transcripts.gtf
cat cufflinks_o250/transcripts.gtf >>  $HOME/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand_merged.cufflinks_o250.transcripts.gtf
gzip -f  $HOME/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand_merged.cufflinks_o250.transcripts.gtf

## run cufflinks to get FPKM
cufflinks -v --no-update-check --frag-len-mean 400 --frag-len-std-dev 100 --library-type fr-firststrand -o cufflinks_o100 -p $cpu -g $ANNOTATION/genes.gtf -M $ANNOTATION/chrM.rRNA.tRNA.gtf accepted_hits.bam -u -j 0.2 --min-frags-per-transfrag 40 --overlap-radius 100 2&> cufflinks_o100.run.log

# gtf of assembly
echo "track name=mouse_adult_wt_RNAseq_PE50nt_strand_merged_o100 description=mouse_adult_wt_RNAseq_PE50nt_strand_merged_o100 visibility=pack colorByStrand='200,100,0 0,100,200'" > $HOME/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand_merged.cufflinks_o100.transcripts.gtf
cat cufflinks_o100/transcripts.gtf >>  $HOME/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand_merged.cufflinks_o100.transcripts.gtf
gzip -f  $HOME/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand_merged.cufflinks_o100.transcripts.gtf
