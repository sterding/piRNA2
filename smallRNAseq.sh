## interactive-version of pipeline for running RNAseq mapping/assmelby on cluster
## Author: Xianjun Dong (xianjun.dong@umassmed.edu)
## Usage: ./smallRNAseq.sh ~/nearline/Xin/smallRNA/Jia/Phil.SRA.wt.ox.6wk.testes.raw.xkxh.norm.all0.bam Phil.SRA.wt.ox.6wk.testes.raw.xkxh.norm.all0
# ./smallRNAseq.sh ~/nearline/Xin/smallRNA/Jia/Phil.SRA.wt.unox.6wk.testes.raw.xkxh.norm.all0.bam Phil.SRA.wt.unox.6wk.testes.raw.xkxh.norm.all0
## Date: 2012-March-13
## Version: 0.1

#!/bin/sh

#1.Prepare input/parameters for the script

mappingbam=$1
samplename=$2
junctionReads=$HOME/scratch/mouse_adult_wt_RNAseq_merge/mouse_adult_wt_RNAseq_merge.junction.bam
junciton2intron_bed=

dir=$HOME/scratch
index=mm9

[ -d $dir/$samplename ] || mkdir $dir/$samplename;

echo "
###########################################
## bash script for running RNAseq pipeline on cluster
###########################################
#!/bin/sh
#$ -V
#$ -pe single 8
#$ -cwd
#$ -o \$HOME/sge_jobs_output/sge_job.\$JOB_ID.out -j y
#$ -S /bin/bash
#$ -m e
#$ -l mem_free=4G

echo \$0;
echo \$@;
###########################################
echo \"############## 1. Configuring\";
###########################################

cd $dir/$samplename

export BOWTIE_INDEXES=\$GENOME/$index/Sequence/BowtieIndex
export BOWTIE2_INDEXES=\$GENOME/$index/Sequence/Bowtie2Index
export ANNOTATION=\$GENOME/$index/Annotation/Genes

############################################
echo \"############# 2. remap RNAseq (PE50, stranded) with provided junction\";
############################################

~/bin/tophat-2.0.4.Linux_x86_64/tophat -o $dir/$samplename --no-convert-bam -p $cpu --read-mismatches $mm $tophat $PE --library-type fr-firststrand --min-anchor-length 8 --min-intron-length 30 --max-intron-length 50000 --splice-mismatches 1 --max-multihits 100 --no-coverage-search genome_offrate3 $readsfile
if [ -f $junctionReads ]; then
    echo \"Junction file found\";
    echo \"Oringal accepted_hits.bam from Tophat is renamed as accepted_hits.beforeMergeWithJunction.bam!\";
    samtools reheader \$BOWTIE_INDEXES/sorted.header.sam $mappingbam > accepted_hits.beforeMergeWithJunction.bam
    echo \"Merging with the junction reads: $junctionReads!\";
    #samtools cat -o accepted_hits.bam accepted_hits.beforeMergeWithJunction.bam $junctionReads
    samtools merge -f accepted_hits.bam accepted_hits.beforeMergeWithJunction.bam $junctionReads
else
    echo \"Junction file not found, just do reheader\";
    samtools reheader \$BOWTIE_INDEXES/sorted.header.sam $mappingbam > accepted_hits.bam
fi
echo \"Sorting bam file ...\";
samtools sort accepted_hits.bam accepted_hits.sorted
mv accepted_hits.sorted.bam accepted_hits.bam
samtools index accepted_hits.bam

############################################
echo \"############# 6. prepare for tracks files to display on UCSC\";
############################################

[ -d $dir/ucsc ] || mkdir $dir/ucsc

# make index for the (sorted) BAM
cp accepted_hits.bam $dir/ucsc/$samplename.accepted_hits.bam
cp accepted_hits.bam.bai $dir/ucsc/$samplename.accepted_hits.bam.bai

# bigwig
bamToBed -i $dir/ucsc/$samplename.accepted_hits.bam -split > $dir/ucsc/$samplename.accepted_hits.bed
awk '{if(\$6==\"+\") print}' $dir/ucsc/$samplename.accepted_hits.bed | sort -k1,1 | bedItemOverlapCount mm9 -chromSize=\$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n > $dir/ucsc/$samplename.accepted_hits.+.bedGraph
awk '{if(\$6==\"-\") print}' $dir/ucsc/$samplename.accepted_hits.bed | sort -k1,1 | bedItemOverlapCount mm9 -chromSize=\$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n | awk '{OFS=\"\t\"; print \$1,\$2,\$3,\"-\"\$4}' > $dir/ucsc/$samplename.accepted_hits.-.bedGraph
bedGraphToBigWig $dir/ucsc/$samplename.accepted_hits.+.bedGraph \$ANNOTATION/ChromInfo.txt $dir/ucsc/$samplename.accepted_hits.+.bw
bedGraphToBigWig $dir/ucsc/$samplename.accepted_hits.-.bedGraph \$ANNOTATION/ChromInfo.txt $dir/ucsc/$samplename.accepted_hits.-.bw
rm $dir/ucsc/$samplename.accepted_hits.bed $dir/ucsc/$samplename.accepted_hits*.bedGraph

###########################################
echo \"################# 7. assembly\";
###########################################

## run cufflinks to get FPKM
## cufflinks2 --overlap-radius 150 --no-update-check -o $dir/$samplename/cufflinks -p 8 -g \$ANNOTATION/genes.gtf -M \$ANNOTATION/chrM.rRNA.tRNA.gtf -b \$BOWTIE_INDEXES/genome.fa --multi-read-correct accepted_hits.bam
## for assmelby
#cufflinksOut=$dir/$samplename/cufflinks_maskRef
#cufflinks --overlap-radius 250 --max-bundle-frags 100000000 --no-update-check -o \$cufflinksOut -p 8 -M \$ANNOTATION/genes.gtf accepted_hits.bam
#cufflinksOut=$dir/$samplename/cufflinks_noRef
#cufflinks --overlap-radius 250 --max-bundle-frags 100000000 --no-update-check -o \$cufflinksOut -p 8 accepted_hits.bam
#cufflinksOut=$dir/$samplename/cufflinks_norm
#cufflinks --overlap-radius 250 --max-bundle-frags 100000000 --no-update-check -o \$cufflinksOut -p 8 -g \$ANNOTATION/genes.gtf -M \$ANNOTATION/chrM.rRNA.tRNA.gtf --frag-len-mean 70 --frag-len-std-dev 20 --upper-quartile-norm --no-effective-length-correction accepted_hits.bam
#cufflinksOut=$dir/$samplename/cufflinks_multi
#cufflinks --overlap-radius 250 --max-bundle-frags 100000000 --no-update-check -o \$cufflinksOut -p 8 -g \$ANNOTATION/genes.gtf -M \$ANNOTATION/chrM.rRNA.tRNA.gtf --multi-read-correct accepted_hits.bam
cufflinksOut=$dir/$samplename/cufflinks_multinorm
cufflinks --overlap-radius 200 --max-bundle-frags 100000000 --no-update-check -o \$cufflinksOut -p 8 -g \$ANNOTATION/genes.gtf -M \$ANNOTATION/chrM.rRNA.tRNA.gtf --upper-quartile-norm --no-effective-length-correction --multi-read-correct accepted_hits.bam

# gtf of assembly
echo \"track name=$samplename description=$samplename visibility=pack colorByStrand='200,100,0 0,100,200'\" > $dir/ucsc/$samplename.transcripts.gtf
cat \$cufflinksOut/transcripts.gtf >> $dir/ucsc/$samplename.transcripts.gtf
gzip -f $dir/ucsc/$samplename.transcripts.gtf


echo \"JOBDONE!\"

">$dir/$samplename/$samplename.sge

cat $dir/$samplename/$samplename.sge
qsub $dir/$samplename/$samplename.sge

# when the above is done
# rsync -azv $dir/ucsc/ dongx@zlab.umassmed.edu:~dongx/public_html/tracks/RNAseq
