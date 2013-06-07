################################################################
# The wrapper script for the RNAseq analysis pipeline
# Usage: TOADD
# version: 1.0
# date: 2013-02-28
# Author: Xianjun Dong
################################################################


#!/bin/sh
#$ -V
#$ -pe single 8
#$ -cwd
#$ -o $HOME/sge_jobs_output/sge_job.$JOB_ID.out -j y
#$ -S /bin/bash
#$ -l mem_free=4G

# include
source $HOME/projects/mylib.sh

###########################################
echo "############## 0. Configuration"
###########################################

samplename=$1 #mouse_adult_wt_RNAseq_PE50nt_strand_R2
read1=$2 # R1.fastq
read2=$3 # R2.fastq
smallRNAreads=$4 #$HOME/nearline/Xin/smallRNA/Jia/Phil.SRA.wt.ox.6wk.testes.raw.uniq.reads.fa.gz

trinity_output=$HOME/nearline/Xin/RNAseq/$samplename/trinity.introns.bed
#link to /home/wangj2/wnearline/xin_rnaseq/junc_pirna_cluster/Phil.Rnaseq.mouse_adult_testis.npa.trinity.blat.maplen95p.mm1p.allmap.introns.bed # intron in bed6

cpu=8

export BOWTIE2_INDEXES=$GENOME/mm9/Sequence/Bowtie2Index/
export ANNOTATION=$GENOME/mm9/Annotation/Genes
[ -d $HOME/scratch/ucsc ] || mkdir $HOME/scratch/ucsc
[ -d $HOME/scratch/$samplename ] || mkdir -p $HOME/scratch/$samplename

ln -fs $HOME/sge_jobs_output/sge_job.$JOB_ID.out $HOME/scratch/$samplename/sge2.log

cd $HOME/nearline/Xin/RNAseq/$samplename

phred=`getphred $read1`; tophat_scoreoption=""; # default
[ "$phred" == "Phred+64" ] && tophat_scoreoption="--solexa1.3-quals";
[ "$phred" == "Solexa+64" ] && tophat_scoreoption="--solexa-quals";
echo "tophat_scoreoption: $tophat_scoreoption"

###########################################
echo "############## 1. QC"
###########################################
run_QC R1.fastq R2.fastq

###########################################
echo "############## 2. Mapping"
###########################################
run_mapping R1.fastq R2.fastq

# mapping stat
get_mapping_stat mapping.sam

###########################################
echo "############## 3. Assembly"
###########################################
run_assembly mapping.bam

###########################################
echo "############## 3.2 Refine with CAGE/PAS"
###########################################
run_5refinement annotation.gtf CAGE_narrowpeak
run_3refinement annotation.gtf PAS_narrowpeak

###########################################
echo "############## 4. quantify mRNA"
###########################################
run_quanfiy_mRNA transcripts.gtf mapping.bam

###########################################
echo "############## 5. quantify smallRNA"
###########################################
./run_quantify_smallRNA transcripts.gtf $HOME/scratch/Gan2011_mouse_TypeASpermatogonia_8dpp_wt_smallRNA/map2piRNA.bed
./run_quantify_smallRNA transcripts.gtf $HOME/scratch/Gan2011_mouse_TypeASpermatogonia_8dpp_wt_smallRNA/map2piRNA_u.bed

./run_quantify_smallRNA ~/scratch/mouse_07dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/transcripts.gtf $HOME/scratch/Gan2011_mouse_TypeASpermatogonia_8dpp_wt_smallRNA/map2piRNA.uniqmap.bed exon > ~/scratch/mouse_07dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/transcripts.uniqmap.ppm.rpkm &
./run_quantify_smallRNA ~/scratch/mouse_07dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/transcripts.gtf $HOME/scratch/Gan2011_mouse_TypeASpermatogonia_8dpp_wt_smallRNA/map2piRNA.allmap.bed exon > ~/scratch/mouse_07dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/transcripts.allmap.ppm.rpkm &
./run_quantify_smallRNA ~/scratch/mouse_07dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/transcripts.gtf $HOME/scratch/Gan2011_mouse_TypeASpermatogonia_8dpp_wt_smallRNA/map2piRNA.uniqmap.bed gene > ~/scratch/mouse_07dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/genes.uniqmap.ppm.rpkm &
./run_quantify_smallRNA ~/scratch/mouse_07dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/transcripts.gtf $HOME/scratch/Gan2011_mouse_TypeASpermatogonia_8dpp_wt_smallRNA/map2piRNA.allmap.bed gene > ~/scratch/mouse_07dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/genes.allmap.ppm.rpkm &
./run_quantify_smallRNA ~/scratch/mouse_07dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/transcripts.gtf $HOME/scratch/Gan2011_mouse_TypeASpermatogonia_8dpp_wt_smallRNA/map2piRNA.uniqmap.bed gene yes > ~/scratch/mouse_07dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/genes.intronyes.uniqmap.ppm.rpkm &
./run_quantify_smallRNA ~/scratch/mouse_07dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/transcripts.gtf $HOME/scratch/Gan2011_mouse_TypeASpermatogonia_8dpp_wt_smallRNA/map2piRNA.allmap.bed gene yes > ~/scratch/mouse_07dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/genes.intronyes.allmap.ppm.rpkm &


cd $HOME/scratch/mouse_10dpp_wt_smallRNA_ox_rep1/
$HOME/projects/piRNA/src/run_quantify_smallRNA ~/scratch/mouse_10.5dpp_wt_RNAseq_PE50nt_strand/cufflinks2_o100/transcripts.gtf map2piRNA.uniqmap.bed gene no uniq > genes.uniqmap.nointron.ppm.rpkm &
$HOME/projects/piRNA/src/run_quantify_smallRNA ~/scratch/mouse_10.5dpp_wt_RNAseq_PE50nt_strand/cufflinks2_o100/transcripts.gtf map2piRNA.allmap.bed gene no all > genes.allmap.nointron.ppm.rpkm &
sum_reads_all__sum_species_all__sum_reads_uniq__sum_species_uniq=`awk '{split($4,a,"_");if($5==1) uniq[$4]=a[2]; all[$4]=a[2];}END{n=0;s=0;n2=0;s2=0;for(i in all) {n++;s=s+all[i]}; for(i in uniq) {n2++;s2=s2+uniq[i];} print s,n, s2, n2;}' map2piRNA.allmap.bed`
Rscript $HOME/projects/piRNA/src/plot_quantify_smallRNA genes.allmap.nointron.ppm.rpkm genes.uniqmap.nointron.ppm.rpkm mouse_10dpp_wt_smallRNA_ox_rep1.gene.nointron.piRNA $sum_reads_all__sum_species_all__sum_reads_uniq__sum_species_uniq
# screenshot for top loci
gtf2bed ~/scratch/mouse_10.5dpp_wt_RNAseq_PE50nt_strand/cufflinks2_o100/transcripts.gtf | grep -f <(awk '{if($7>100 || $8>100)print}' genes.allmap.nointron.ppm.rpkm | awk '{print $1"__"}') > selected.bed

cd $HOME/scratch/mouse_10dpp_wt_smallRNA_ox_rep1/
$HOME/projects/piRNA/src/run_quantify_smallRNA ~/scratch/mouse_10.5dpp_wt_RNAseq_PE50nt_strand/cufflinks2_o100/transcripts.gtf map2piRNA.uniqmap.bed gene yes uniq > genes.uniqmap.yesintron.ppm.rpkm &
$HOME/projects/piRNA/src/run_quantify_smallRNA ~/scratch/mouse_10.5dpp_wt_RNAseq_PE50nt_strand/cufflinks2_o100/transcripts.gtf map2piRNA.allmap.bed gene yes all > genes.allmap.yesintron.ppm.rpkm &
sum_reads_all__sum_species_all__sum_reads_uniq__sum_species_uniq=`awk '{split($4,a,"_");if($5==1) uniq[$4]=a[2]; all[$4]=a[2];}END{n=0;s=0;n2=0;s2=0;for(i in all) {n++;s=s+all[i]}; for(i in uniq) {n2++;s2=s2+uniq[i];} print s,n, s2, n2;}' map2piRNA.allmap.bed`
Rscript $HOME/projects/piRNA/src/plot_quantify_smallRNA genes.allmap.yesintron.ppm.rpkm genes.uniqmap.yesintron.ppm.rpkm mouse_10dpp_wt_smallRNA_ox_rep1.gene.yesintron.piRNA $sum_reads_all__sum_species_all__sum_reads_uniq__sum_species_uniq

cd $HOME/scratch/mouse_42dpp_wt_smallRNA_ox_rep1/
$HOME/projects/piRNA/src/run_quantify_smallRNA ~/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks_o100/transcripts.gtf map2piRNA.uniqmap.bed gene no uniq > genes.uniqmap.nointron.ppm.rpkm &
$HOME/projects/piRNA/src/run_quantify_smallRNA ~/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks_o100/transcripts.gtf map2piRNA.allmap.bed gene no all > genes.allmap.nointron.ppm.rpkm &
sum_reads_all__sum_species_all__sum_reads_uniq__sum_species_uniq=`awk '{split($4,a,"_");if($5==1) uniq[$4]=a[2]; all[$4]=a[2];}END{n=0;s=0;n2=0;s2=0;for(i in all) {n++;s=s+all[i]}; for(i in uniq) {n2++;s2=s2+uniq[i];} print s,n, s2, n2;}' map2piRNA.allmap.bed`
Rscript $HOME/projects/piRNA/src/plot_quantify_smallRNA genes.allmap.nointron.ppm.rpkm genes.uniqmap.nointron.ppm.rpkm mouse_adult_wt_RNAseq_PE50nt_strand.gene.nointron.piRNA $sum_reads_all__sum_species_all__sum_reads_uniq__sum_species_uniq
