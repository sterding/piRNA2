#!/bin/sh
#$ -V
#$ -pe single 8
#$ -cwd
#$ -o $HOME/sge_jobs_output/sge_job.$JOB_ID.out -j y
#$ -S /bin/bash
#$ -l mem_free=2G

# usage: qsub ~/projects/piRNA/src/junction_analysis.sh

# include
source $HOME/projects/mylib.sh

# arguments
junctions2intron_bed=$1  # bed6 format
RNAseq_R1=$2 # $HOME/nearline/Xin/RNAseq/mouse_adult_wt_RNAseq_PE50nt_strand/Adult_testis_NoIndex_L004_R1.fq
RNAseq_R2=$3 # $HOME/nearline/Xin/RNAseq/mouse_adult_wt_RNAseq_PE50nt_strand/Adult_testis_NoIndex_L004_R2.fq
smallRNA_R1=$4 # $HOME/nearline/Xin/smallRNA/Jia/Phil.SRA.wt.ox.6wk.testes.raw.uniq.reads.fa.gz

# configure
index=mm9
FLANKING=200

# log
ln -fs $HOME/sge_jobs_output/sge_job.$JOB_ID.out sge.log


echo "# ======================================"
echo "# get donor-acceptor signal"
echo "# ======================================"

awk -v fk=2 '{OFS="\t";left=($6=="+")?($2-2):($3-2); right=($6=="+")?($3-2):($2-2); print $1,left,left+4,$4,$5,$6;print $1,right,right+4,$4,$5,$6;}' $junctions2intron_bed | fastaFromBed -tab -s -name -fi $GENOME/$index/Sequence/BowtieIndex/genome.fa -bed stdin -fo stdout | groupBy -g 1 -c 2 -o collapse | sort -k1,1 | awk '{OFS="\t";split($2,a,"");print $1,tolower(a[3])tolower(a[4])"-"tolower(a[6])tolower(a[7]), toupper(a[1])toupper(a[2])"("tolower(a[3])tolower(a[4])"-"tolower(a[6])tolower(a[7])")"toupper(a[8])toupper(a[9]);}' > ${junctions2intron_bed/bed/intronsignal.tab}

echo "# ======================================"
echo "##### map smallRNA/RNAseq reads to the splicing sites"
echo "# ======================================"

# left + right splicing site
awk -v fk=$FLANKING '{OFS="\t"; flanking=($2>fk)?fk:$2; print $1,$2-flanking,$2+flanking,$4"_1", $5,$6; flanking=($3>fk)?fk:$3; print $1,$3-flanking,$3+flanking,$4"_2", $5,$6;}' $junctions2intron_bed > ${junctions2intron_bed/bed/splicing.bed}

################################
echo "## map RNAseq reads to splicing sites ..."
################################
mapping_junction_reads mapped2splicing_rnaseq -q $RNAseq_R1 1 ${junctions2intron_bed/bed/splicing.bed}
awk 'BEGIN{FS="\t";OFS="\t";}{print $1,$3,$2,$5,$4}' mapped2splicing_rnaseq/reads.species.tab > ${junctions2intron_bed/bed/mapped2splicing_rnaseq_1.bed}

mapping_junction_reads mapped2splicing_rnaseq -q $RNAseq_R2 1 ${junctions2intron_bed/bed/splicing.bed}
cp mapped2splicing_rnaseq/reads.species.tab ${junctions2intron_bed/bed/mapped2splicing_rnaseq_2.bed}
paste ${junctions2intron_bed/bed/mapped2splicing_rnaseq_1.bed} ${junctions2intron_bed/bed/mapped2splicing_rnaseq_2.bed} | awk 'BEGIN{FS="\t";OFS="\t";}{print $1,$2+$7,$3+$8,$4+$9,$5+$10}' | sed 's/_/\t/g' | sort -k1,1 | groupBy -g 1 -c 3,4,5,6 -o sum,sum,sum,sum > ${junctions2intron_bed/bed/mapped2splicing_rnaseq.bed}

################################
echo "## map piRNA reads to splicing site ..."
################################
mapping_junction_reads mapped2splicing_pirna -f $smallRNA_R1 1 ${junctions2intron_bed/bed/splicing.bed}
sed 's/_/\t/g' mapped2splicing_pirna/reads.species.tab | sort -k1,1 | groupBy -g 1 -c 3,4,5,6 -o sum,sum,sum,sum > ${junctions2intron_bed/bed/mapped2splicing_pirna.bed}

echo "# ======================================"
echo "##### map smallRNA/RNAseq reads to the junctions"
echo "# ======================================"

echo "# extend it to bed12"
awk -v fk=$FLANKING '{OFS="\t"; flanking=($2>fk)?fk:$2; print $1,$2-flanking,$3+flanking,$4, $5,$6,$2-flanking,$3+flanking,"255,0,0",2,flanking","flanking,"0,"(flanking+$3-$2);}' $junctions2intron_bed > ${junctions2intron_bed/bed/extended.bed}

################################
echo "# map RNAseq reads to junction ..."
################################

# bowtie2 does not support gz for -U, so do this first
# zcat Adult_testis_NoIndex_L004_R1*.gz > Adult_testis_NoIndex_L004_R1.fq;  zcat Adult_testis_NoIndex_L004_R2*.gz > Adult_testis_NoIndex_L004_R2.fq;
# or
# zcat ~/nearline/Xin/RNAseq/mouse_adult_wt_RNAseq_PE50nt_strand/Adult_testis_NoIndex_L004_*.gz | sed 's/ /\//g' > ~/nearline/Xin/RNAseq/mouse_adult_wt_RNAseq_PE50nt_strand/Adult_testis_NoIndex_L004.fq &

# since PE50 is from fr-firststrand, so /1 is from antisense-strand, and /2 is from sense strand.
mapping_junction_reads mapped2junction_rnaseq -q $RNAseq_R1 1 ${junctions2intron_bed/bed/extended.bed} mapped2genome_rnaseq_accepted_hits_R1.sam
awk 'BEGIN{FS="\t";OFS="\t";}{print $1,$3,$2,$5,$4}' mapped2junction_rnaseq/reads.species.tab > ${junctions2intron_bed/bed/mapped2junction_rnaseq_1.bed}

mapping_junction_reads mapped2junction_rnaseq -q $RNAseq_R2 1 ${junctions2intron_bed/bed/extended.bed} mapped2genome_rnaseq_accepted_hits_R2.sam
cp mapped2junction_rnaseq/reads.species.tab ${junctions2intron_bed/bed/mapped2junction_rnaseq_2.bed}
paste ${junctions2intron_bed/bed/mapped2junction_rnaseq_1.bed} ${junctions2intron_bed/bed/mapped2junction_rnaseq_2.bed} | awk 'BEGIN{FS="\t";OFS="\t";}{print $1,$3+$7,$2+$8,$5+$9,$4+$10}' > ${junctions2intron_bed/bed/mapped2junction_rnaseq.bed}

################################
echo "# map smallRNA reads to junction ..."
################################

mapping_junction_reads mapped2junction_pirna -f $smallRNA_R1 1 ${junctions2intron_bed/bed/extended.bed} mapped2genome_pirna_accepted_hits.sam
paste mapped2junction_pirna/reads.species.tab > ${junctions2intron_bed/bed/mapped2junction_pirna.bed}

echo "## merge together"

sort -k4,4 $junctions2intron_bed | paste - ${junctions2intron_bed/bed/intronsignal.tab}  ${junctions2intron_bed/bed/mapped2junction_rnaseq.bed} ${junctions2intron_bed/bed/mapped2junction_pirna.bed} ${junctions2intron_bed/bed/mapped2splicing_rnaseq.bed} ${junctions2intron_bed/bed/mapped2splicing_pirna.bed} | sort -k1,1 -k2,2n -k3,3n | cut -f1-6,8-9,11-14,16-19,21-24,26-29 > ${junctions2intron_bed/bed/final.tab}

echo " ======= DONE =========="
