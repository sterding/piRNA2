#!/bin/sh

FLANKING=50

zcat ~/nearline/Xin/smallRNA/Jia/smallRNA-term.TAP.GCAC.raw.uniqmap.xkxh.norm.bed.gz | awk '{OFS="\t";if($1~/^chr/) {tss=($4=="+")?$2:($3-1); print $1, tss, tss+1, $5, $6, $4;}}' > ~/nearline/Xin/smallRNA/Jia/smallRNA-term.TAP.GCAC.raw.uniqmap.xkxh.norm.bed

for i in ../data/controlSet3/piRNA.prepachytene.bed ../data/controlSet3/piRNA.hybrid.bed ../data/controlSet3/piRNA.pachytene.bi.bed ../data/controlSet3/piRNA.pachytene.uni.bed ../data/controlSet3/most.abundant.x100912.transcript.bed.NM.bi ../data/controlSet3/most.abundant.x100912.transcript.bed.NM.uni ../data/controlSet3/most.abundant.x100912.transcript.bed.NR;
do
    echo $i;
    #grep -v track $i | awk '{OFS="\t"; tss=($6=="+")?$2:$3; print $1,tss,tss,$4,$5,$6;}' | slopBed -l $FLANKING -r $FLANKING -s -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b ~/nearline/Xin/smallRNA/Jia/smallRNA-term.TAP.GCAC.raw.uniqmap.xkxh.norm.bed -s -wao | groupBy -g 1,2,3,4,5,6 -c 11 -o sum | sort > /tmp/`basename $i`
    #grep -v track $i | cut -f1-6 | slopBed -b 1000 -s -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b ~/nearline/Xin/smallRNA/Jia/smallRNA-term.TAP.GCAC.raw.uniqmap.xkxh.norm.bed -s -wao | groupBy -g 1,2,3,4,5,6 -c 11 -o sum | sort | cut -f7 | paste /tmp/`basename $i` - > $i.TSS$FLANKING.gene1000.count.tab

    paste <(sort -k4,4 -k5,5 $i.TSS$FLANKING.count.tab) <(sort -k4,4 -k5,5 $i.gene1000.count.tab | cut -f7) >  $i.TSS$FLANKING.gene1000.count.tab
done

## CAGE
awk '{OFS="\t";$4=".";$5=1;print}' /home/dongx/scratch/mouse_adult_wt_CAGE_PE100nt_strand/mouse_adult_wt_CAGE_PE100nt_strand.unique.bed | groupBy -g 1,2,3,4,6 -c 5 -o sum | awk '{OFS="\t"; print $1,$2,$3,$4,$6,$5}'> /tmp/mouse_adult_wt_CAGE_PE100nt_strand.unique.collapsed.bed
sort -k1,1 -k2,2n -k6,6 /tmp/mouse_adult_wt_CAGE_PE100nt_strand.unique.collapsed.bed | groupBy -g 1,2,3,4,6 -c 5 -o sum | awk '{OFS="\t"; print $1,$2,$3,$4,$6,$5}' > /home/dongx/scratch/mouse_adult_wt_CAGE_PE100nt_strand/mouse_adult_wt_CAGE_PE100nt_strand.unique.collapsed.bed

for i in ../data/controlSet3/piRNA.prepachytene.bed ../data/controlSet3/piRNA.hybrid.bed ../data/controlSet3/piRNA.pachytene.bi.bed ../data/controlSet3/piRNA.pachytene.uni.bed ../data/controlSet3/most.abundant.x100912.transcript.bed.NM.bi ../data/controlSet3/most.abundant.x100912.transcript.bed.NM.uni ../data/controlSet3/most.abundant.x100912.transcript.bed.NR;
do
    echo $i;
    grep -v track $i | awk '{OFS="\t"; tss=($6=="+")?$2:$3; print $1,tss,tss,$4,$5,$6;}' | slopBed -l $FLANKING -r $FLANKING -s -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b ~/scratch/mouse_adult_wt_CAGE_PE100nt_strand/mouse_adult_wt_CAGE_PE100nt_strand.unique.collapsed.bed -s -wao | groupBy -g 1,2,3,4,5,6 -c 11 -o sum | sort -k4,4 -k5,5 > /tmp/`basename $i`
    grep -v track $i | cut -f1-6 | slopBed -b 1000 -s -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b ~/scratch/mouse_adult_wt_CAGE_PE100nt_strand/mouse_adult_wt_CAGE_PE100nt_strand.unique.collapsed.bed -s -wao | groupBy -g 1,2,3,4,5,6 -c 11 -o sum | sort -k4,4 -k5,5 | cut -f7 | paste /tmp/`basename $i` - > $i.CAGE.TSS$FLANKING.gene1000.count.tab
    #paste <(sort -k4,4 -k5,5 $i.TSS$FLANKING.count.tab) <(sort -k4,4 -k5,5 $i.gene1000.count.tab | cut -f7) >  $i.TSS$FLANKING.gene1000.count.tab
done

exit;

## R
pdf('smallRNA.TAP.CAGE.density.TSS50.gene1000.count.pdf', width=21,height=8, colormodel='cmyk')
par(mfrow=c(2,7),mar=c(8,4,4,2), oma = c(0, 0, 2, 0))
#TAP
for(i in paste(c('../data/controlSet3/piRNA.prepachytene.bed', '../data/controlSet3/piRNA.hybrid.bed', '../data/controlSet3/piRNA.pachytene.bi.bed', '../data/controlSet3/piRNA.pachytene.uni.bed', '../data/controlSet3/most.abundant.x100912.transcript.bed.NM.bi', '../data/controlSet3/most.abundant.x100912.transcript.bed.NM.uni', '../data/controlSet3/most.abundant.x100912.transcript.bed.NR'),'.TSS50.gene1000.count.tab', sep=''))
{
    df=read.table(i, header=F)
    df=df[,c(7,8)]
    colnames(df)=c('count_TSS50','count_gene1000')
    df$count_TSS50[df$count_TSS50==-1]=0
    df$count_gene1000[df$count_gene1000==-1]=1

    ii=sub('.*controlSet3/(.*).TSS.*','\\1', i)
    ii=sub('most.abundant.x100912.transcript.bed.','',ii)
    ii=sub('.bed','',ii)
    ii=sub('.bi','\nbidirectional',ii)
    ii=sub('.uni','\nunidirectional',ii)
    print(ii)
    boxplot(100*df$count_TSS50/df$count_gene1000, 100-100*df$count_TSS50/df$count_gene1000, names=c("5' in [-50,50]", "5' not in [-50,50]"), ylab="Percent(%)")
    title(main=ii, sub=paste("Wilcox test p-value:", format(wilcox.test(100*df$count_TSS50/df$count_gene1000, 100-100*df$count_TSS50/df$count_gene1000, paired=T)$p.value,3,3)))
}
mtext("PARS", outer = TRUE, cex = 1.5)
# CAGE
for(i in paste(c('../data/controlSet3/piRNA.prepachytene.bed', '../data/controlSet3/piRNA.hybrid.bed', '../data/controlSet3/piRNA.pachytene.bi.bed', '../data/controlSet3/piRNA.pachytene.uni.bed', '../data/controlSet3/most.abundant.x100912.transcript.bed.NM.bi', '../data/controlSet3/most.abundant.x100912.transcript.bed.NM.uni', '../data/controlSet3/most.abundant.x100912.transcript.bed.NR'),'.CAGE.TSS50.gene1000.count.tab', sep=''))
{
    df=read.table(i, header=F)
    df=df[,c(7,8)]
    colnames(df)=c('count_TSS50','count_gene1000')
    df$count_TSS50[df$count_TSS50==-1]=0
    df$count_gene1000[df$count_gene1000==-1]=1

    ii=sub('.*controlSet3/(.*).TSS.*','\\1', i)
    ii=sub('most.abundant.x100912.transcript.bed.','',ii)
    ii=sub('.bed','',ii)
    ii=sub('.bi','\nbidirectional',ii)
    ii=sub('.uni','\nunidirectional',ii)
    ii=sub('.CAGE','',ii)
    print(ii)
    boxplot(100*df$count_TSS50/df$count_gene1000, 100-100*df$count_TSS50/df$count_gene1000, names=c("5' in [-50,50]", "5' not in [-50,50]"), ylab="Percent(%)")
    title(main=ii, sub=paste("Wilcox test p-value:", format(wilcox.test(100*df$count_TSS50/df$count_gene1000, 100-100*df$count_TSS50/df$count_gene1000, paired=T)$p.value,3,3)))
}
mtext("CAGE", outer = TRUE, line=-30, cex = 1.5)

dev.off()



