################################################################
# The wrapper script for the smallRNAseq analysis pipeline
# Usage:

#qsub $HOME/projects/piRNA/src/smallRNAseq.pipeline.sh /home/dongx/nearline/Xin/smallRNA/mouse_00dpp_wt_smallRNA_ox_rep1.fastq mouse_00dpp_wt_smallRNA_ox_rep1_mm1 1 /home/dongx/scratch/mouse_04dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/transcripts.gtf /home/wangj2/wnearline/xin_rnaseq/blat/Phil.Rnaseq.nonPolyA.B6.4dpp.testis.trinity.blat.maplen95p.mm1p.allmap.bed #9065593 (Tue May 28 17:33:38 EDT 2013)

#qsub $HOME/projects/piRNA/src/smallRNAseq.pipeline.sh /home/dongx/nearline/Xin/smallRNA/mouse_00dpp_wt_smallRNA_ox_rep1.fastq mouse_00dpp_wt_smallRNA_ox_rep1 0 /home/dongx/scratch/mouse_04dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/transcripts.gtf /home/wangj2/wnearline/xin_rnaseq/blat/Phil.Rnaseq.nonPolyA.B6.4dpp.testis.trinity.blat.maplen95p.mm1p.allmap.bed #9065276
#qsub $HOME/projects/piRNA/src/smallRNAseq.pipeline.sh /home/dongx/nearline/Xin/smallRNA/mouse_00dpp_wt_smallRNA_unox_rep1.fastq mouse_00dpp_wt_smallRNA_unox_rep1 0 /home/dongx/scratch/mouse_04dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/transcripts.gtf /home/wangj2/wnearline/xin_rnaseq/blat/Phil.Rnaseq.nonPolyA.B6.4dpp.testis.trinity.blat.maplen95p.mm1p.allmap.bed #
#
#qsub $HOME/projects/piRNA/src/smallRNAseq.pipeline.sh /home/dongx/nearline/Xin/smallRNA/mouse_02dpp_wt_smallRNA_ox_rep1.fastq mouse_02dpp_wt_smallRNA_ox_rep1 0 /home/dongx/scratch/mouse_04dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/transcripts.gtf /home/wangj2/wnearline/xin_rnaseq/blat/Phil.Rnaseq.nonPolyA.B6.4dpp.testis.trinity.blat.maplen95p.mm1p.allmap.bed #9065277
#qsub $HOME/projects/piRNA/src/smallRNAseq.pipeline.sh /home/dongx/nearline/Xin/smallRNA/mouse_02dpp_wt_smallRNA_unox_rep1.fastq mouse_02dpp_wt_smallRNA_unox_rep1 0 /home/dongx/scratch/mouse_04dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/transcripts.gtf /home/wangj2/wnearline/xin_rnaseq/blat/Phil.Rnaseq.nonPolyA.B6.4dpp.testis.trinity.blat.maplen95p.mm1p.allmap.bed #9064594
#
#qsub $HOME/projects/piRNA/src/smallRNAseq.pipeline.sh /home/dongx/nearline/Xin/smallRNA/mouse_04dpp_wt_smallRNA_ox_rep1.fastq mouse_04dpp_wt_smallRNA_ox_rep1 0 /home/dongx/scratch/mouse_04dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/transcripts.gtf /home/wangj2/wnearline/xin_rnaseq/blat/Phil.Rnaseq.nonPolyA.B6.4dpp.testis.trinity.blat.maplen95p.mm1p.allmap.bed #9065278
#qsub $HOME/projects/piRNA/src/smallRNAseq.pipeline.sh /home/dongx/nearline/Xin/smallRNA/mouse_04dpp_wt_smallRNA_unox_rep1.fastq mouse_04dpp_wt_smallRNA_unox_rep1 0 /home/dongx/scratch/mouse_04dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/transcripts.gtf /home/wangj2/wnearline/xin_rnaseq/blat/Phil.Rnaseq.nonPolyA.B6.4dpp.testis.trinity.blat.maplen95p.mm1p.allmap.bed #NO DATA
#
#qsub $HOME/projects/piRNA/src/smallRNAseq.pipeline.sh /home/dongx/nearline/Xin/smallRNA/mouse_07dpp_wt_smallRNA_ox_rep1.fastq mouse_07dpp_wt_smallRNA_ox_rep1 0 /home/dongx/scratch/mouse_07dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/transcripts.gtf /home/wangj2/wnearline/xin_rnaseq/blat/Phil.Rnaseq.nonPolyA.B6.7dpp.testis.trinity.blat.maplen95p.mm1p.allmap.bed #9065279
#qsub $HOME/projects/piRNA/src/smallRNAseq.pipeline.sh /home/dongx/nearline/Xin/smallRNA/mouse_07dpp_wt_smallRNA_unox_rep1.fastq mouse_07dpp_wt_smallRNA_unox_rep1 0 /home/dongx/scratch/mouse_07dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/transcripts.gtf /home/wangj2/wnearline/xin_rnaseq/blat/Phil.Rnaseq.nonPolyA.B6.7dpp.testis.trinity.blat.maplen95p.mm1p.allmap.bed  #9064586

# qsub $HOME/projects/piRNA/src/smallRNAseq.pipeline.sh /home/dongx/nearline/Xin/smallRNA/mouse_10dpp_wt_smallRNA_ox_rep1.fastq mouse_10dpp_wt_smallRNA_ox_rep1 0 /home/dongx/scratch/mouse_10.5dpp_wt_RNAseq_PE50nt_strand/cufflinks2_o100/transcripts.gtf /home/wangj2/wnearline/xin_rnaseq/blat/Phil.Rnaseq.nonPolyA.B6.10.5dpp.testis.trinity.blat.maplen90p.mm5p.allmap.bed #9065280
# qsub $HOME/projects/piRNA/src/smallRNAseq.pipeline.sh /home/dongx/nearline/Xin/smallRNA/mouse_10dpp_wt_smallRNA_unox_rep1.fastq mouse_10dpp_wt_smallRNA_unox_rep1 0 /home/dongx/scratch/mouse_10.5dpp_wt_RNAseq_PE50nt_strand/cufflinks2_o100/transcripts.gtf /home/wangj2/wnearline/xin_rnaseq/blat/Phil.Rnaseq.nonPolyA.B6.10.5dpp.testis.trinity.blat.maplen90p.mm5p.allmap.bed #9050829

# qsub $HOME/projects/piRNA/src/smallRNAseq.pipeline.sh ~/nearline/Xin/smallRNA/028.mm.wild.type.42.rep1.sraOX.fastq mouse_42dpp_wt_smallRNA_ox_rep1 0 /home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks_o100/transcripts.gtf /home/wangj2/wnearline/xin_rnaseq/blat/Phil.Rnaseq.nonPolyA.B6.adult.testis.pooled.trinity.blat.maplen90p.mm5p.allmap.bed  #9065281
# qsub $HOME/projects/piRNA/src/smallRNAseq.pipeline.sh ~/nearline/Xin/smallRNA/006.mm.wild.type.42.rep1.sraUNox.fastq mouse_42dpp_wt_smallRNA_unox_rep1 0 /home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks_o100/transcripts.gtf /home/wangj2/wnearline/xin_rnaseq/blat/Phil.Rnaseq.nonPolyA.B6.adult.testis.pooled.trinity.blat.maplen90p.mm5p.allmap.bed  #9064077


# map smallRNA to early stage of annotation to see if more covered
# qsub $HOME/projects/piRNA/src/smallRNAseq.pipeline.sh /home/dongx/nearline/Xin/smallRNA/mouse_10dpp_wt_smallRNA_ox_rep1.fastq mouse_10dpp_wt_smallRNA_ox_rep1_to4dpp 0 /home/dongx/scratch/mouse_04dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/transcripts.gtf /home/wangj2/wnearline/xin_rnaseq/blat/Phil.Rnaseq.nonPolyA.B6.4dpp.testis.trinity.blat.maplen95p.mm1p.allmap.bed #9050292
# qsub $HOME/projects/piRNA/src/smallRNAseq.pipeline.sh ~/nearline/Xin/smallRNA/006.mm.wild.type.42.rep1.sraUNox.fastq mouse_42dpp_wt_smallRNA_unox_rep1_to4dpp 0 /home/dongx/scratch/mouse_04dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/transcripts.gtf /home/wangj2/wnearline/xin_rnaseq/blat/Phil.Rnaseq.nonPolyA.B6.4dpp.testis.trinity.blat.maplen95p.mm1p.allmap.bed #9050293
# qsub $HOME/projects/piRNA/src/smallRNAseq.pipeline.sh ~/nearline/Xin/smallRNA/028.mm.wild.type.42.rep1.sraOX.fastq mouse_42dpp_wt_smallRNA_ox_rep1_to4dpp 0 /home/dongx/scratch/mouse_04dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/transcripts.gtf /home/wangj2/wnearline/xin_rnaseq/blat/Phil.Rnaseq.nonPolyA.B6.4dpp.testis.trinity.blat.maplen95p.mm1p.allmap.bed #9050830
# qsub $HOME/projects/piRNA/src/smallRNAseq.pipeline.sh /home/dongx/nearline/Xin/fromPublication/Gan2011_mouse_TypeASpermatogonia_8dpp_wt_smallRNA.fastq Gan2011_mouse_TypeASpermatogonia_8dpp_wt_smallRNA 0 /home/dongx/scratch/mouse_07dpp_wt_RNAseq_PE100nt_strand/cufflinks2_o100/transcripts.gtf /home/wangj2/wnearline/xin_rnaseq/blat/Phil.Rnaseq.nonPolyA.B6.7dpp.testis.trinity.blat.maplen95p.mm1p.allmap.bed   #9050294


# version: 1.0
# date: 2013-02-28
# Author: Xianjun Dong
# Limit: this only works for single-end smallRNA fastq input
# LOG:
# 1. add rev-complementary to adaptor.fa (May25-2013)
################################################################

#!/bin/sh
#$ -V
#$ -pe openmpi 24
#$ -cwd
#$ -o $HOME/sge_jobs_output/sge_job.$JOB_ID.out -j y
#$ -S /bin/bash
#$ -l mem_free=8G

export BOWTIE_INDEXES=$GENOME/mm9/Sequence/BowtieIndex/
export ANNOTATION=$GENOME/mm9/Annotation
#ncRNA=$ANNOTATION/ncRNA/ncRNA.wt_miRNA_piRNA.bed12  # ncRNA from fRNAdb (http://www.ncrna.org/frnadb/download); it seems fRNAdb has mixed tRNA with miRNA...
ncRNA=$ANNOTATION/ncRNA/ncRNA.bed  # mergen tRNA from tRNAscan-SE and rRNA/snRNA/snoRNA from biomart
miRNA=$ANNOTATION/SmallRNA/mm10.liftover.to.mm9.bed  # precursor from miRBase
repeats=$ANNOTATION/Variation/repeatmasker.mm9.ucsc.sorted.bed  # repeat from repeatmasker
genes=$ANNOTATION/Genes/genes.gtf  # mRNA
pseudogenes=$ANNOTATION/Genes/pseudoYale60.bed # pseudogenes from Yale sets
exons=$ANNOTATION/Genes/genes.exons.bed12  # mRNA
cds=$ANNOTATION/Genes/genes.cds.bed12  # mRNA
introns=$ANNOTATION/Genes/genes.introns.bed12
three_prime_utr=$ANNOTATION/Genes/genes.3utr.bed12
five_prime_utr=$ANNOTATION/Genes/genes.5utr.bed12

cpu=24

###########################################
echo "############## 0. Configuration"
###########################################
# include
source $HOME/projects/mylib.sh
export PATH=$PATH:$HOME/projects/piRNA/src

fullfilepath=$1
samplename=$2
mm=$3
annotation_from_cufflinks=$4 # gtf format of cufflinks
annotation_from_trinity=$5 # bed12 format of trinity

annotation="$(mktemp ~/scratch/cufflinkstrinity.XXXXXXXXXX)"
echo "# making temp file in bed12 format of trinity+cufflinks (only used for mapping reads to junctions)"
echo $annotation
[[ "$annotation_from_trinity" == "" ]] && ln -fs $annotation_from_cufflinks $annotation
[[ "$annotation_from_trinity" != "" ]] && gtf2bed $annotation_from_cufflinks | cat - $annotation_from_trinity | sort-bed - > $annotation

## adaptor sequence
#PCRId1: 0.5dppox
#5Õ-CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA-3Õ
[[ "$samplename" =~ "mouse_00dpp_wt_smallRNA_ox" ]] && adaptersequence="CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
#PCRId2: 0.5dppunox
#5Õ-CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA-3Õ
[[ "$samplename" =~ "mouse_00dpp_wt_smallRNA_unox" ]] && adaptersequence="CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
#PCRId3: 2.5dppox
#5Õ-CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA-3Õ
[[ "$samplename" =~ "mouse_02dpp_wt_smallRNA_ox" ]] && adaptersequence="CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
#PCRId4: 2.5dppunox
#5Õ-CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA-3Õ
[[ "$samplename" =~ "mouse_02dpp_wt_smallRNA_unox" ]] && adaptersequence="CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
#PCRId5: 4.5dppox
#5Õ-CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA-3Õ
[[ "$samplename" =~ "mouse_04dpp_wt_smallRNA_ox" ]] && adaptersequence="CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
#PCRId6: 4.5dppunox
#5Õ-CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA-3Õ
[[ "$samplename" =~ "mouse_04dpp_wt_smallRNA_unox" ]] && adaptersequence="CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
#PCRId7: 7.5dppox
#5Õ-CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA-3Õ
[[ "$samplename" =~ "mouse_07dpp_wt_smallRNA_ox" ]] && adaptersequence="CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
#PCRId8: 7.5dppunox
#5Õ-CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA-3Õ
[[ "$samplename" =~ "mouse_07dpp_wt_smallRNA_unox" ]] && adaptersequence="CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA";
# 10dpp/adult: TCGTATGCCGTCTTCTGCTTG
[[ "$samplename" =~ "mouse_10dpp_wt_smallRNA|mouse_42dpp_wt_smallRNA" ]] && adaptersequence="TCGTATGCCGTCTTCTGCTTG";
echo "adaptor sequence: " $adaptersequence

ext="${fullfilepath##*.}"
[[ "$samplename" == "" ]] && { samplename=$(basename "$fullfilepath"); samplename="${samplename%.*}";}
[[ "$mm" == "" ]] && mm=0;

# check the file format by file extension
[[ "$ext" =~ "fq|fastq|FASTQ" ]] || die "The input is not in fastq format. You have to convert it into fastq in order to run the pipeline.";

Min=16; Max=40

[ -d $HOME/scratch/ucsc ] || mkdir $HOME/scratch/ucsc
[ -d $HOME/scratch/$samplename ] || mkdir -p $HOME/scratch/$samplename

cd $HOME/scratch/$samplename
echo -e ">adapter1\n$adaptersequence" > adaptor.fa; faRc -keepCase adaptor.fa stdout | fasta_formatter >> adaptor.fa;
ln -fs $fullfilepath $samplename.fq
ln -fs $HOME/sge_jobs_output/sge_job.$JOB_ID.out run.log

############################################
#echo "############## 1. QC"
############################################
run_QC_smallRNA $samplename.fq

###########################################
echo "############## 1.2 Collapse reads into uniq species [optional for smallRNA only]"
###########################################
_insert2uniqreads filtered/$samplename.fq $Min $Max > $samplename.fa

##########################################
echo "############## 2. Mapping"
##########################################
## map to genome
[ "$mm" == "0" ] && bowtie -v 0 -a --un unmapped.fa -p 8 -f genome $samplename.fa 2>mapping.log | awk '{OFS="\t";print $3,$4,$4+length($5),$1,$7,$2,$4,$4+length($5),0,1,length($5),0;}' > map2genome.bed
[ "$mm" != "0" ] && bowtie -v $mm -a --best --un unmapped.fa -p 8 -f genome $samplename.fa 2>mapping.log | awk '{OFS="\t";print $3,$4,$4+length($5),$1,$7,$2,$4,$4+length($5),0,1,length($5),0;}' > map2genome.bed

echo "### map to transcriptome"
# 0 mismatch, and at least 3nt overhang
mapping_junction_reads2 map2junctions -f $samplename.fa 0 $annotation
## merge reads mapped to genome and transcriptome (using sort-merge)
##sort -m -k4,4 <(sort -k4,4 map2junctions/junctions.bed) map2genome.bed > allmap.bed
LC_ALL=C sort -m --parallel=$cpu --buffer-size=2G --temporary-directory=$HOME/scratch -k4,4 <(sort -k4,4 map2junctions/junctions.bed) map2genome.bed > allmap.bed

echo "# add mapping count (like NH tag)"
join -1 4 -2 1 -a 1 -a 2 allmap.bed <(groupBy -g 4 -c 4 -o count -i allmap.bed) | awk 'BEGIN{OFS="\t";}{print $2,$3,$4,$1,$13,$6,$7,$8,$9,$10,$11,$12}' > allmap.bed2
#LC_ALL=C sort --parallel=$cpu --buffer-size=2G -k1,1 -k2,2n allmap.bed2 > allmap.bed  # parallel version
sort-bed --max-mem 8G allmap.bed2 > allmap.bed
rm allmap.bed2

# unique mapper only
awk '{if($5==1) print}' allmap.bed > uniqmap.bed

###########################################
echo "############## 3. distribution of reads mapped to different categories of annotation (e.g. miRNA(precursor) ncRNA(tRNA/rRNA/snRNA/snoRNA), piRNA (genic/TE/pseudogene/intergenic etc.))"
###########################################
rm -f *.pdf allmap.bed_2_piRNA.unaccountable*
echo "### handle uniq-mappers"
bash reads.map2where.sh uniqmap.bed
echo "### handle all-mappers (which takes time)"
bash reads.map2where.sh allmap.bed

###########################################
echo "############## 4. quantify the assembled loci with piRNA reads"
###########################################
f=`bash run_quantify_smallRNA $annotation_from_cufflinks uniqmap.bed_2_piRNA gene no uniq`
sort -k2,2nr $f > genes.uniqmap.nointron.ppm.rpkm
f=`bash run_quantify_smallRNA $annotation_from_cufflinks allmap.bed_2_piRNA gene no all`
sort -k2,2nr $f > genes.allmap.nointron.ppm.rpkm
sum__readsall_speciesall_readsuniq_speciesuniq=`awk '{split($4,a,"_");if($5==1) uniq[$4]=a[2]; all[$4]=a[2];}END{n=0;s=0;n2=0;s2=0;for(i in all) {n++;s=s+all[i]}; for(i in uniq) {n2++;s2=s2+uniq[i];} print s,n, s2, n2;}' allmap.bed_2_piRNA`
Rscript $HOME/projects/piRNA/src/plot_quantify_smallRNA genes.allmap.nointron.ppm.rpkm genes.uniqmap.nointron.ppm.rpkm $samplename.gene.nointron.piRNA $sum__readsall_speciesall_readsuniq_speciesuniq

echo "### handle allmap.bed_2_piRNA.unaccountable"
bash reads.map2where.sh allmap.bed_2_piRNA.unaccountable piRNA

for i in *.mapping.statistics*.pdf; do mv $i $samplename.$i; done

# screenshot for top loci
#gtf2bed $annotation_from_cufflinks | grep -f <(awk '{if($7>100 || $8>100)print}' genes.allmap.nointron.ppm.rpkm | awk '{print $1"__"}') > selected.bed
