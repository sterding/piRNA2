## sge script to run NorahDesk to re-construct ncRNA
## Author: Xianjun Dong (xianjun.dong@umassmed.edu)
## Date: 2012-July-11
## Version: 0.1
## Usage: qsub run_NorahDesk.sh

#!/bin/sh
#$ -V
#$ -pe single 4
#$ -cwd
#$ -o $HOME/sge_jobs_output/sge_job.$JOB_ID.out -j y
#$ -S /bin/bash
#$ -M sterding.hpcc@gmail.com
#$ -S /bin/bash
#$ -m e
#$ -l mem_free=4G

export CLASSPATH=$HOME/bin/NorahDesk/bin
export ANNOTATION=$GENOME/mm9/Annotation/Genes

samplename="mouse_adult_wt_smallRNAseq_merge"

mkdir ~/scratch/$samplename; cd ~/scratch/$samplename

# 1. merge BAM from smallRNA-seq libs
samtools cat -o $samplename.bam `ls ~/scratch/ucsc/mouse_adult_wt_smallRNAseq_76nt_strand_*.bam`

# 2. generate tracks for UCSC browsing
# index (sorted) bam
samtools sort $samplename.bam $samplename.sorted
mv $samplename.sorted.bam $samplename.bam
samtools index $samplename.bam
#bam -> bigwig
#bamToBed -i mouse_adult_wt_RNAseq_merge.bam -split | sort -k1,1 | bedItemOverlapCount mm9 -chromSize=$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n > mouse_adult_wt_RNAseq_merge.bedGraph
bedtools genomecov -ibam -strand + -split -bg -dz -i $samplename.bam -g $ANNOTATION/ChromInfo.txt
bedGraphToBigWig $samplename.bedGraph $ANNOTATION/ChromInfo.txt $samplename.bw

echo " 3. Assembly using cufflinks ........";
## for assmelby
cufflinks --overlap-radius 150 --no-update-check -o ./cufflinks -p 8 -g $ANNOTATION/genes.gtf -M $ANNOTATION/chrM.rRNA.tRNA.gtf $samplename.bam

# gtf of assembly
echo "track name=$samplename.cufflinks description=$samplename.cufflinks visibility=pack colorByStrand='200,100,0 0,100,200'" > $HOME/scratch/ucsc/$samplename.cufflinks.transcripts.gtf
cat ./cufflinks/transcripts.gtf >> $HOME/scratch/ucsc/$samplename.cufflinks.transcripts.gtf
gzip -f $HOME/scratch/ucsc/$samplename.cufflinks.transcripts.gtf

echo " 4. Assembly using NorahDesk ........"

java ncRNA_prediction.PredictRNA NUMCHR=19 GENOMESEQ=$GENOME/mm9/Sequence/Chromosomes INFILE=$samplename.bam
##java ncRNA_annotation.AnnotateRNA NCRNA=<file name> ENSEMBLE=<file name> INFILE=<file name>


cp ~/scratch/mouse_adult_wt_RNAseq_merge/*.bam* ~/scratch/mouse_adult_wt_RNAseq_merge/*.bw ~/scratch/ucsc
