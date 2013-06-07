################################################################
# The wrapper script to see where the mapping reads mapped to
# Usage:
# qsub -l mem_free=5G -pe openmpi 12 $HOME/projects/piRNA/src/reads.map2where.sh ~/scratch/mouse_10dpp_wt_smallRNA_ox_rep1/map2piRNA.allmap.bed.unaccountable #9031518
# bash $HOME/projects/piRNA/src/reads.map2where.sh ~/scratch/mouse_10dpp_wt_smallRNA_ox_rep1/allmap.bed #9032073
# bash $HOME/projects/piRNA/src/reads.map2where.sh ~/scratch/mouse_10dpp_wt_smallRNA_ox_rep1/uniqmap.bed #9032074
# version: 1.0
# date: 2013-05-16
# Author: Xianjun Dong
## Requirement: the input bed files are already sorted (e.g. by sort-bed -)
################################################################

#!/bin/sh

##$ -V
##$ -pe openmpi 12
##$ -cwd
##$ -o $HOME/sge_jobs_output/sge_job.$JOB_ID.out -j y
##$ -S /bin/bash
##$ -l mem_free=5G

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

# include
source $HOME/projects/mylib.sh
export PATH=$PATH:$HOME/projects/piRNA/src

BEDfile_fullpath=$1
type=$2 # piRNA

## ============ divide the reads into different categories ==============

if [ "$type" != "piRNA" ]
then
    # 1. mapped to miRNA precusor
    intersectBed -a $BEDfile_fullpath -b $miRNA -s -wo -sorted -f 0.5 > $BEDfile_fullpath"_2_miRNA"

    # 2. mapped to the other ncRNA (tRNA/rRNA/snRNA/snoRNA)
    intersectBed -a $BEDfile_fullpath -b $ncRNA -s -wo -sorted -f 0.5 > $BEDfile_fullpath"_2_ncRNA"

    # 3: divide the reads unmapped to miRNA/ncRNA into piRNA and non-piRNA according to their lenght;
    # Note: if a read maps to miRNA/ncRNA annotation, it won't go to the next step for calculating piRNA coveragex
    ## Alternative: do the size-selection later
    > $BEDfile_fullpath"_2_piRNA"; > $BEDfile_fullpath"_2_non24_40nt"
    fgrep -f <(cat $BEDfile_fullpath"_2_ncRNA" $BEDfile_fullpath"_2_miRNA" | awk '{a[$4]=1;}END{for(i in a) print i;}') -vw $BEDfile_fullpath | awk -v path=$BEDfile_fullpath '{if(($3-$2)>=24 && ($3-$2)<=40) print >> path"_2_piRNA"; else print >> path"_2_non24_40nt";}'
fi
[[ "$type" == "piRNA" ]] && ln -sf $BEDfile_fullpath $BEDfile_fullpath"_2_piRNA"

## futher divide the putative piRNA reads to repeats/exons/UTR/intron/intergenic
# 4. reads2repeats in individual class
rm -f $BEDfile_fullpath"_2_piRNA_"*
## divide to the 10 classes used in UCSC (http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=337043749&c=chr4&g=rmsk)
# intersectBed -a $BEDfile_fullpath"_2_piRNA" -b $repeats -s -wo -sorted -split -f 0.5| awk -v path=$BEDfile_fullpath"_2_piRNA" '{OFS="\t"; split($16, n, "__"); print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12>> path"_repeat_"n[3];}'
## divide to the 7 classes used in RepBase (http://www.girinst.org/repbase/update/browse.php): TE(Class1: (LTR, LINE/SINE), Class2(DNA), Simple_repeat, pseudogene, others)
#intersectBed -a $BEDfile_fullpath"_2_piRNA" -b $repeats -s -wo -sorted -split -f 0.5| awk -v path=$BEDfile_fullpath"_2_piRNA" '{OFS="\t"; split($16, n, "__"); class=n[3]; if(class~/Satellite|Low_complexity/) class="Simple_repeat"; if(class~/RNA/) class="Pseudogene"; if(class~/Unknown|RC/) class="Other"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12>> path"_repeat_"class;}'
# add L1 and IAP (May28-2013)
intersectBed -a $BEDfile_fullpath"_2_piRNA" -b $repeats -s -wo -sorted -split -f 0.5| awk -v path=$BEDfile_fullpath"_2_piRNA" '{OFS="\t"; split($16, n, "__"); class=n[3]; if($16~/IAP.*__ERVK__LTR/) class="LTR_IAP"; if($16~/L1__LINE/) class="LINE_L1"; if(class~/Satellite|Low_complexity/) class="Simple_repeat"; if(class~/RNA/) class="Pseudogene"; if(class~/Unknown|RC/) class="Other"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12>> path"_repeat_"class;}'

intersectBed -a $BEDfile_fullpath"_2_piRNA" -b $repeats -s -wa -sorted -split -v -f 0.5 | sort-bed - > $BEDfile_fullpath"-2_piRNA_nonrepeats"
# 5. reads2genic
intersectBed -a $BEDfile_fullpath"-2_piRNA_nonrepeats" -b <(sort-bed $five_prime_utr) -sorted -s -wa -split -u -f 0.5> $BEDfile_fullpath"_2_piRNA_genic_5UTR"
intersectBed -a $BEDfile_fullpath"-2_piRNA_nonrepeats" -b <(sort-bed $cds) -sorted -s -wa -split -u -f 0.5> $BEDfile_fullpath"_2_piRNA_genic_CDS"
intersectBed -a $BEDfile_fullpath"-2_piRNA_nonrepeats" -b <(sort-bed $introns) -sorted -s -wa -split -u -f 0.5> $BEDfile_fullpath"_2_piRNA_genic_intron"
intersectBed -a $BEDfile_fullpath"-2_piRNA_nonrepeats" -b <(sort-bed $three_prime_utr) -sorted -s -wa -split -u -f 0.5> $BEDfile_fullpath"_2_piRNA_genic_3UTR"
# 6. reads2pseudogene
intersectBed -a $BEDfile_fullpath"-2_piRNA_nonrepeats" -b <(sort-bed $pseudogenes) -sorted -s -f 0.5> $BEDfile_fullpath"_2_piRNA_pseudogene"
# 7. read2intergenic
intersectBed -a $BEDfile_fullpath"-2_piRNA_nonrepeats" -b <(gtf2bed $genes | sort-bed - $pseudogenes) -sorted -s -v -f 0.5> $BEDfile_fullpath"_2_piRNA_intergenic"

## ============ get stat of species/normalized_reads on each category ==============
> $BEDfile_fullpath".mapping.statistics.txt"
[[ "$type" == "piRNA" ]] || awk -v type=`basename $BEDfile_fullpath` '{OFS="\t"; split($4,a,"_"); species[$4]=1; reads=reads+(a[2]/$5);}END{for(i in species) n_species++; if(n_species) print type, n_species, reads;}' $BEDfile_fullpath >> $BEDfile_fullpath".mapping.statistics.txt"
for i in $BEDfile_fullpath"_2_"*; do awk -v type=`basename $i` '{OFS="\t"; split($4,a,"_"); species[$4]=1; reads=reads+(a[2]/$5);}END{for(i in species) n_species++; if(n_species) print type, n_species, reads;}' $i >> $BEDfile_fullpath".mapping.statistics.txt"; done

# divide to the 7 classes used in RepBase (http://www.girinst.org/repbase/update/browse.php): TE(Class1: (LTR, LINE/SINE), Class2(DNA), Simple_repeat, pseudogene, others)
#sed 's/Satellite\|Low_complexity/Simple_repeat/;s/_tRNA\|_srpRNA\|_snRNA\|_scRNA\|_rRNA\|_RNA/_Pseudogene/;s/_Unknown\|_RC/_Other/;' $BEDfile_fullpath".mapping.statistics.txt" | sort | groupBy -g 1 -c 2,3 -o sum,sum > $BEDfile_fullpath".mapping.statistics.2.txt"

## ============ plot the statistics ==============

#
echo "
args <- commandArgs(TRUE)
filename=args[1]
title=sub('.txt', '', filename)
df=read.table(filename, header=F)
rownames(df)=df[,1]; df=df[,-1]; colnames(df)=c('species_count', 'partitioned_reads_count')

# total number of reads/species
total_s=df[grep('_2_', rownames(df), invert =T), 1]
total_r=df[grep('_2_', rownames(df), invert =T), 2]

df=df[grep('_2_', rownames(df)),]
rownames(df)=sub('.*_2_', '', rownames(df))
d=c('miRNA','ncRNA', rownames(df)[c(grep('repeat_LTR', rownames(df)), grep('SINE|LINE', rownames(df)), grep('repeat_DNA', rownames(df)), grep('Simple_repeat|Low_complexity|Satellite', rownames(df)), grep('repeat.*RNA|repeat_Pseudogene', rownames(df)), grep('repeat_Unknown|repeat_Other|repeat_RC', rownames(df)), grep('genic_5UTR', rownames(df)), grep('genic_CDS', rownames(df)), grep('genic_intron', rownames(df)), grep('genic_3UTR', rownames(df)), grep('piRNA_pseudogene', rownames(df)), grep('intergenic', rownames(df)))])

df=df[match(d,rownames(df), nomatch=0),]

pdf(paste(title, 'pdf', sep='.'), width=8, height=8)
#pie(mapping)
par(mar=c(8,5,2,0.5), mgp=c(4,0.1,0), tck = 0.03, las=2, mfrow=c(2,1));
barplot(df[,1], cex.name=0.8, names.arg=rownames(df), main=title, ylab=colnames(df)[1])
barplot(df[,2], names.arg=rownames(df),  cex.name=0.8, main='', ylab=colnames(df)[2])
dev.off()
" > /tmp/stat.R
Rscript /tmp/stat.R $BEDfile_fullpath".mapping.statistics.txt"
#Rscript /tmp/stat.R $BEDfile_fullpath".mapping.statistics.2.txt"

echo "Done: reads.map2where.sh $1 $2";
