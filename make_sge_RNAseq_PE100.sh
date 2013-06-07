###########################################
## bash script for running RNAseq pipeline on cluster
###########################################

#!/bin/sh
#$ -V
#$ -pe openmpi 8
#$ -cwd
#$ -o $HOME/sge_jobs_output/sge_job.$JOB_ID.out -j y
#$ -S /bin/bash
#$ -l mem_free=4G

# usage:
# qsub ~/projects/piRNA/src/make_sge_RNAseq_PE100.sh mouse_04dpp_wt_RNAseq_PE100nt_strand ~/nearline/Xin/RNAseq/01-mm/wild.type.04.rep1.rnaseq.r1.paired.tworun.fastq.gz ~/nearline/Xin/RNAseq/01-mm/wild.type.04.rep1.rnaseq.r2.paired.tworun.fastq.gz
# qsub ~/projects/piRNA/src/make_sge_RNAseq_PE100.sh mouse_07dpp_wt_RNAseq_PE100nt_strand ~/nearline/Xin/RNAseq/01-mm/wild.type.07.rep1.rnaseq.r1.paired.tworun.fastq.gz ~/nearline/Xin/RNAseq/01-mm/wild.type.07.rep1.rnaseq.r2.paired.tworun.fastq.gz

###########################################
echo "############## 1. Configuring"
###########################################

samplename=$1 #mouse_adult_wt_RNAseq_PE50nt_strand_R2
read1=$2 # R1.fastq
read2=$3 # R2.fastq
smallRNAreads=$4 #$HOME/nearline/Xin/smallRNA/Jia/Phil.SRA.wt.ox.6wk.testes.raw.uniq.reads.fa.gz

trinity_output=$HOME/nearline/Xin/RNAseq/$samplename/trinity.introns.bed #/home/wangj2/wnearline/xin_rnaseq/junc_pirna_cluster/Phil.Rnaseq.mouse_adult_testis.npa.trinity.blat.maplen95p.mm1p.allmap.introns.bed # intron in bed6

export BOWTIE2_INDEXES=$GENOME/mm9/Sequence/Bowtie2Index/

cpu=8

export ANNOTATION=$GENOME/mm9/Annotation/Genes
[ -d $HOME/scratch/ucsc ] || mkdir $HOME/scratch/ucsc
[ -d $HOME/scratch/$samplename ] || mkdir -p $HOME/scratch/$samplename

ln -fs $HOME/sge_jobs_output/sge_job.$JOB_ID.out $HOME/scratch/$samplename/sge2.log


phred=`getphred $read1`; tophat_scoreoption=""; # default
[ "$phred" == "Phred+64" ] && tophat_scoreoption="--solexa1.3-quals";
[ "$phred" == "Solexa+64" ] && tophat_scoreoption="--solexa-quals";
echo "tophat_scoreoption: $tophat_scoreoption"


############################################
#echo "################  2. quality filter: adaptor removal/clip"
############################################
#
## QC
#[ -d /home/dongx/scratch/$samplename/preQC ] || mkdir -p /home/dongx/scratch/$samplename/preQC
#fastqc --outdir /home/dongx/scratch/$samplename/preQC --extract --threads 4 -f fastq Adult_testis_NoIndex_L004_R1.fq Adult_testis_NoIndex_L004_R2.fq
#rm /home/dongx/scratch/$samplename/preQC/*fastqc.zip
#
####### adaptor removal
#[ -d filtered ] || mkdir filtered
#fastq-mcf -o filtered/Adult_testis_NoIndex_L004_R1.fq -o filtered/Adult_testis_NoIndex_L004_R2.fq -l 16 -q 15 -w 4 -x 10 -u -P 33 adaptor.fa Adult_testis_NoIndex_L004_R1.fq Adult_testis_NoIndex_L004_R2.fq
#
##############################################
#echo "################# 3. QC"
#############################################
#
#cd filtered
#[ -d /home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/postQC ] || mkdir -p /home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/postQC
#fastqc --outdir /home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/postQC --extract --threads 4 -f fastq Adult_testis_NoIndex_L004_R1.fq Adult_testis_NoIndex_L004_R2.fq
#rm /home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/postQC/*fastqc.zip
#cp /home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/*_fastqc /home/dongx/scratch/ucsc

#############################################
echo "################# 4. mapping"
#############################################
## mapped to genome using Tophat (incl. junctions)
# run tophat here...
cd $HOME/scratch/$samplename
export BOWTIE2_INDEXES=$GENOME/mm9/Sequence/Bowtie2Index/
tophat -o $HOME/scratch/$samplename $tophat_scoreoption --no-convert-bam -p $cpu --read-mismatches 2 --mate-inner-dist 400 --mate-std-dev 100 --library-type fr-firststrand --min-anchor-length 8 --min-intron-length 30 --max-intron-length 50000 --splice-mismatches 1 --max-multihits 100 --no-coverage-search --keep-tmp genome_offrate3 $read1 $read2

## check & fix the XS:A tag (require: _fixXStag.awk under the same folder)
_fixXStag -v libtype='fr-firststrand' -v save_discrepancy_to_file=discrepant_reads.sam $HOME/scratch/$samplename/accepted_hits.sam > $HOME/scratch/$samplename/_accepted_hits.sam

##SAM->BAM->index
samtools view -Sbut $BOWTIE2_INDEXES/genome.fai _accepted_hits.sam | samtools sort - accepted_hits.sorted
mv accepted_hits.sorted.bam accepted_hits.bam
samtools index accepted_hits.bam

MAP2GENOME_SAM=$HOME/scratch/$samplename/_accepted_hits.sam

exit;

#############################################
#echo "############  4.1 run junction analysis here ..."
#############################################
#cd $HOME/scratch/$samplename
#[ -d junction_analysis_result ] || mkdir junction_analysis_result;
#cd $HOME/scratch/$samplename/junction_analysis_result
## get junction table (200bp_exon + intron + 200bp_exon)
#tophat_output=$HOME/scratch/$samplename/junctions.bed  # exon+intron+exon format
#grep -v track $tophat_output | awk '{OFS="\t"; if($10==2) {split($11,sz,","); split($12,st, ","); print $1,$2+sz[1],$2+st[2],$4,0,$6;}}' > junctions2intron.tophat.bed
#
## only for mouse adult
##awk '{OFS="\t"; if($10==2) {split($11,sz,","); split($12,st, ","); print $1,$2+sz[1],$2+st[2],".",0,$6;}}' junctions.tophat.bed12 | sort -u | awk '{OFS="\t"; printf("%s\t%s\t%s\tJUNC%.8d\t%s\t%s\n", $1,$2,$3,NR,$5,$6);}' > junctions2intron.tophat.bed
#
#if [ -e $trinity_output ];
#then
#    grep -v track $trinity_output | awk '{OFS="\t"; print $1, $2, $3, "junc"NR, 0, $6;}' > junctions2intron.trinity.bed
#    cat junctions2intron.tophat.bed junctions2intron.trinity.bed | sort -k1,1 -k2,2n -k3,3n -k6,6 | groupBy -g 1,2,3,6 -c 4, -o collapse | awk '{OFS="\t"; s=$5;sub(/junc[0-9]+/,"Tri",s); sub(/JUNC[0-9]+/,"Top",s); printf("%s\t%s\t%s\tJUNC%.10d\t%d\t%s\t%s\t%s\n", $1,$2,$3,NR,0,$4,s,$5);}' > junctions2intron.bed
#else
#    cat junctions2intron.tophat.bed | sort -k1,1 -k2,2n -k3,3n -k6,6 | groupBy -g 1,2,3,6 -c 4, -o collapse | awk '{OFS="\t"; s=$5;sub(/junc[0-9]+/,"Tri",s); sub(/JUNC[0-9]+/,"Top",s); printf("%s\t%s\t%s\tJUNC%.10d\t%d\t%s\t%s\t%s\n", $1,$2,$3,NR,0,$4,s,$5);}' > junctions2intron.bed
#fi
#
#bedtools slop -b 2000 -i $HOME/projects/piRNA/result/10-02-12.clusters.with.assigned.name.bed -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a junctions2intron.bed -b stdin -s -wo | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$12,$7,$8;}' | sort -u > junctions2intron.pi.bed9
#
### require the smallRNAseq reads file
#cut -f1-6 junctions2intron.pi.bed9 | sort -u > junctions2intron.pi.bed
#jobid=`qsub ~/projects/piRNA/src/junction_analysis.sh junctions2intron.pi.bed $HOME/nearline/Xin/RNAseq/$samplename/$read1 $HOME/nearline/Xin/RNAseq/$samplename/$read2 $smallRNAreads | cut -f3 -d' '`
#echo "####### junction_analysis job is submiited with jobid $jobid"
#ln -fs $HOME/sge_jobs_output/sge_job.$jobid.out $HOME/scratch/$samplename/junction_analysis_result/sge.log

#############################################
echo "################# 4.2 get reads mapped to the Trinity junctions"
#############################################
cd $HOME/scratch/$samplename

echo "# stat ..."
sam2stat _accepted_hits.sam &

[ -d merge_reads ] || mkdir merge_reads

if [ -e $trinity_output ];
then
    echo "# make gtf file for Trinity intron..."
    cat junction_analysis_result/junctions2intron.trinity.bed | awk -v fk=500 '{OFS="\t"; flanking=($2>fk)?fk:$2; print $1,$2-flanking,$3+flanking,$4, $5, $6, $2-flanking, $3+flanking, "255,0,0", 2, flanking","flanking,"0,"(flanking+$3-$2);}' | bedToGenePred stdin stdout | genePredToGtf file stdin junction_analysis_result/junctions2intron.trinity.gtf

    JUNCTION2INTRON_gtf=$HOME/scratch/$samplename/junction_analysis_result/junctions2intron.trinity.gtf

    echo "# mapping reads to trinity junctions..."
    export BOWTIE2_INDEXES=$GENOME/mm9/Sequence/Bowtie2Index/
    $HOME/bin/tophat-2.0.4.Linux_x86_64/tophat -o ./map2trinityjunctions --no-convert-bam -p $cpu --read-mismatches 1  --mate-inner-dist 300 $tophat_scoreoption --library-type fr-firststrand --min-anchor-length 8 --min-intron-length 30 --max-intron-length 50000 --splice-mismatches 1 --max-multihits 100 --no-coverage-search -G $JUNCTION2INTRON_gtf --transcriptome-index=$HOME/scratch/$samplename/map2junctions/junctions --no-novel-juncs --transcriptome-only genome $HOME/nearline/Xin/RNAseq/$samplename/$read1 $HOME/nearline/Xin/RNAseq/$samplename/$read2

    echo "# intersect with unmapped ..."
    samtools view unmapped.bam | cut -f1 | sed 's/\/[12]$//g' | sort -u > unmapped.ID
    # paired reads mapped to junctions (incl. mates)
    awk '{if($6~/^[0-9]+M[0-9]+N[0-9]+M$/) print $1}' ./map2trinityjunctions/accepted_hits.sam | sort -u | sort -m unmapped.ID - | uniq -d | fgrep -w -f - ./map2trinityjunctions/accepted_hits.sam | _fixXStag -v libtype=fr-firststrand | sort > ./map2trinityjunctions/accepted_hits_junctions.sam
    echo "## merge reads with genome mapper"
    cat $MAP2GENOME_SAM map2trinityjunctions/accepted_hits_junctions.sam > merge_reads/accepted_hits.sam
else
    cp $MAP2GENOME_SAM  merge_reads/accepted_hits.sam
fi

exit

###########################################
echo "############### 5. post-processing, format converting"
###########################################
[ -d merge_reads ] || mkdir merge_reads
cd $HOME/scratch/$samplename/merge_reads;

### sam -> bam -> sorted -> index
#samtools view -Sbut $BOWTIE2_INDEXES/genome.fai accepted_hits.sam | samtools sort - accepted_hits.sorted
#mv accepted_hits.sorted.bam accepted_hits.bam
#samtools index accepted_hits.bam
#
## make index for the (sorted) BAM
#cp accepted_hits.bam  $HOME/scratch/ucsc/$samplename.accepted_hits.bam
#cp accepted_hits.bam.bai  $HOME/scratch/ucsc/$samplename.accepted_hits.bam.bai

## for unique mappers
grep -P "NH:i:1\t" accepted_hits.sam > accepted_hits_unique.sam
samtools view -Sbut $BOWTIE2_INDEXES/genome.fai accepted_hits_unique.sam | samtools sort - accepted_hits_unique.sorted
mv accepted_hits_unique.sorted.bam accepted_hits_unique.bam
samtools index accepted_hits_unique.bam

# make index for the (sorted) BAM
cp accepted_hits_unique.bam  $HOME/scratch/ucsc/$samplename.accepted_hits_unique.bam
cp accepted_hits_unique.bam.bai  $HOME/scratch/ucsc/$samplename.accepted_hits_unique.bam.bai

###########################################
echo "################# 6. assembly"
###########################################

cd $HOME/scratch/$samplename/merge_reads;



### run cufflinks to get FPKM
#cufflinks -v --no-update-check --frag-len-mean 400 --frag-len-std-dev 100 --library-type fr-firststrand -o $HOME/scratch/$samplename/cufflinks2_o250 -p $cpu -g $ANNOTATION/genes.gtf -M $ANNOTATION/chrM.rRNA.tRNA.gtf accepted_hits.bam -u -j 0.2 --min-frags-per-transfrag 40 --overlap-radius 250 2&> $HOME/scratch/$samplename/cufflinks2_o250.run.log
#
## gtf of assembly
#echo "track name=${samplename}_o250 description=${samplename}_o250 visibility=pack colorByStrand='200,100,0 0,100,200'" > $HOME/scratch/ucsc/${samplename}_o250.transcripts.gtf
#cat $HOME/scratch/$samplename/cufflinks2_o250/transcripts.gtf >>  $HOME/scratch/ucsc/${samplename}_o250.transcripts.gtf
#gzip -f  $HOME/scratch/ucsc/${samplename}_o250.transcripts.gtf

### run cufflinks to get FPKM
#cufflinks -v --no-update-check --frag-len-mean 400 --frag-len-std-dev 100 --library-type fr-firststrand -o $HOME/scratch/$samplename/cufflinks2_o100 -p $cpu -g $ANNOTATION/genes.gtf -M $ANNOTATION/chrM.rRNA.tRNA.gtf accepted_hits.bam -u -j 0.2 --min-frags-per-transfrag 40 --overlap-radius 100 2&> $HOME/scratch/$samplename/cufflinks2_o100.run.log
#
## gtf of assembly
#echo "track name=${samplename}_o100 description=${samplename}_o100 visibility=pack colorByStrand='200,100,0 0,100,200'" > $HOME/scratch/ucsc/${samplename}_o100.transcripts.gtf
#cat $HOME/scratch/$samplename/cufflinks2_o100/transcripts.gtf >>  $HOME/scratch/ucsc/${samplename}_o100.transcripts.gtf
#gzip -f  $HOME/scratch/ucsc/${samplename}_o100.transcripts.gtf

###########  unique mapper
## run cufflinks to get FPKM
cufflinks -v --no-update-check --frag-len-mean 400 --frag-len-std-dev 100 --library-type fr-firststrand -o $HOME/scratch/$samplename/cufflinks_o100_accepted_hits_unique -p $cpu -g $ANNOTATION/genes.gtf -M $ANNOTATION/chrM.rRNA.tRNA.gtf accepted_hits_unique.bam -u -j 0.2 --min-frags-per-transfrag 40 --overlap-radius 100 2&> $HOME/scratch/$samplename/cufflinks_o100_accepted_hits_unique.run.log

# gtf of assembly
echo "track name=${samplename}_o100_accepted_hits_unique description=${samplename}_o100_accepted_hits_unique visibility=pack colorByStrand='200,100,0 0,100,200'" > $HOME/scratch/ucsc/${samplename}_o100_accepted_hits_unique.transcripts.gtf
cat $HOME/scratch/$samplename/cufflinks_o100_accepted_hits_unique/transcripts.gtf >>  $HOME/scratch/ucsc/${samplename}_o100_accepted_hits_unique.transcripts.gtf
gzip -f  $HOME/scratch/ucsc/${samplename}_o100_accepted_hits_unique.transcripts.gtf



###########################################
echo "############## 7. prepare for tracks files to display on UCSC"
###########################################

# bam -> bigwig
bam2bw  $HOME/scratch/ucsc/$samplename.accepted_hits_unique.bam

echo "JOBDONE!"; exit;
