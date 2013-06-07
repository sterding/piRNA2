#!/bin/sh

# include
source $HOME/projects/mylib.sh

FLANKINGup=500
FLANKINGdown=50
CENTER=150
D=20

inputfile1=/home/lix/nearline/new.assembly/x100910.rnaseq.transcripts.bed
inputfile2=$HOME/projects/piRNA/result/100312.transcript.with.transcripts.removed.new.names.bed
CAGE=/home/dongx/scratch/mouse_adult_wt_CAGE_PE100nt_strand/mouse_adult_wt_CAGE_PE100nt_strand.unique
#CAGE=/home/dongx/scratch/mouse_adult_wt_CAGE_PE100nt_strand/mouse_adult_wt_CAGE_PE100nt_strand  # all mappers
#PAS=/home/dongx/scratch/mouse_adult_wt_PAS200_100nt/mouse_adult_wt_PAS200_100nt.unique
PAS=/home/dongx/scratch/mouse_adult_wt_PAS200_100nt/mouse_adult_wt_PAS200_100nt  # all mappers
annotationFile_all=../data/x100910.rnaseq.transcripts.ALL.bed
annotationFile=../data/100312.transcript.with.transcripts.removed.new.names.bed
annotationFileORF=../data/100312.transcript.with.transcripts.removed.new.names.orf.bed
TSSregions=$annotationFile_all.TSSregion.bed
TTSregions=$annotationFile_all.TTSregion.bed
NM=../data/controlSet2/x100910.rnaseq.transcripts.NM.all.final.bed.c
NR=../data/controlSet2/x100910.rnaseq.transcripts.NR.all.final.bed.c
genic=../data/controlSet2/piRNA.genic.bed
intergenic=../data/controlSet2/piRNA.intergenic.bed
# previous annotation
lau=../data/lau.kingston.2006.pachytene.clusters.mm9.bed.stranded.bed
girard=../data/girard.hannon.2006.pachytene.clusters.mm9.bed
aravin=../data/aravin.hannon.2007.prepachytene.mm9.simpleBed.bed
#H3K4me3=/home/dongx/scratch/mouse_adult_wt_H3K4me3_50nt/mouse_adult_wt_H3K4me3_50nt_peaks.encodePeak  # OLD one
H3K4me3=/home/dongx/scratch/mouse_adult_wt_InputIP_50nt/input.mouse.by.chris/H3K4_Mouse.250.summit.macs.bed  # Update one from Chris
Amyb=~/scratch/mouse_adult_wt_ChIP_Amyb_50nt/A-Myb-mouse_R2.gz.k_summits.bed  # from Chris: /isilon/moorelab/isilon/RoyC/Cell.Submission/01-mm/chip/amyb.mouse

##################################################
### CAGE_Peak <-> 5' and 3' end <->  PAS_PEAK distances from previous annotations (Lau or Girard).
## e.g. CAGE/PAS peak around [-3k, +3k] of the 5'/3' end
##################################################
dos2unix $girard; dos2unix $lau; dos2unix $aravin
cat $lau | awk '{OFS="\t"; if($6=="+") print $1,$2,$2,$4,$5,$6; if($6=="-") print $1,$3,$3,$4,$5,$6;}' | slopBed -b 3000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b $H3K4me3 -wao | sort -k4,4 -k11,11nr | awk 'BEGIN{id="";}{OFS="\t";if($4!=id) {tss=$2+3000;peak=$8+$16;dis=peak-tss;if(dis<0) dis=-dis; if(peak<0 || dis>3000) dis=4000; print "H3K4me3", "lau",$4, dis;id=$4}}' > ../data/distanceToCAGEPASpeak.tab
cat $girard | awk '{OFS="\t"; if($6=="+") print $1,$2,$2,$4,$5,$6; if($6=="-") print $1,$3,$3,$4,$5,$6;}' | slopBed -b 3000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b $H3K4me3 -wao | sort -k4,4 -k11,11nr | awk 'BEGIN{id="";}{OFS="\t";if($4!=id) {tss=$2+3000;peak=$8+$16;dis=peak-tss;if(dis<0) dis=-dis; if(peak<0 || dis>3000) dis=4000; print "H3K4me3","girard",$4, dis ;id=$4}}' >> ../data/distanceToCAGEPASpeak.tab
cat $aravin | awk '{OFS="\t"; if($6=="+") print $1,$2,$2,$4,$5,$6; if($6=="-") print $1,$3,$3,$4,$5,$6;}' | slopBed -b 3000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b $H3K4me3 -wao | sort -k4,4 -k11,11nr | awk 'BEGIN{id="";}{OFS="\t";if($4!=id) {tss=$2+3000;peak=$8+$16;dis=peak-tss;if(dis<0) dis=-dis; if(peak<0 || dis>3000) dis=4000; print "H3K4me3","aravin",peak, $4, dis ;id=$4}}' >> ../data/distanceToCAGEPASpeak.tab
grep -P "\.1\t" $annotationFileORF | awk '{OFS="\t"; if($6=="+") print $1,$2,$2,$4,$5,$6; if($6=="-") print $1,$3,$3,$4,$5,$6;}' | slopBed -b 3000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b $H3K4me3 -wao | sort -k4,4 -k11,11nr | awk 'BEGIN{id="";}{OFS="\t";if($4!=id) {tss=$2+3000;peak=$8+$16;dis=peak-tss;if(dis<0) dis=-dis; if(peak<0 || dis>3000) dis=4000; print "H3K4me3","li",$4, dis ;id=$4}}' >> ../data/distanceToCAGEPASpeak.tab

# update for Chris
cat $lau | awk '{OFS="\t"; if($6=="+") print $1,$2,$2,$4,$5,$6; if($6=="-") print $1,$3,$3,$4,$5,$6;}' | slopBed -b 3000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b $H3K4me3 -wao | sort -k4,4 -k11,11nr | awk 'BEGIN{id="";}{OFS="\t";if($4!=id) {tss=$2+3000;peak=$8;dis=peak-tss;if(dis<0) dis=-dis; if(peak<0 || dis>3000) dis=4000; print "H3K4me3", "lau",$4, dis;id=$4}}' > ../data/distanceToH3K4peak.tab
cat $girard | awk '{OFS="\t"; if($6=="+") print $1,$2,$2,$4,$5,$6; if($6=="-") print $1,$3,$3,$4,$5,$6;}' | slopBed -b 3000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b $H3K4me3 -wao | sort -k4,4 -k11,11nr | awk 'BEGIN{id="";}{OFS="\t";if($4!=id) {tss=$2+3000;peak=$8;dis=peak-tss;if(dis<0) dis=-dis; if(peak<0 || dis>3000) dis=4000; print "H3K4me3","girard",$4, dis ;id=$4}}' >> ../data/distanceToH3K4peak.tab
cat $aravin | awk '{OFS="\t"; if($6=="+") print $1,$2,$2,$4,$5,$6; if($6=="-") print $1,$3,$3,$4,$5,$6;}' | slopBed -b 3000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b $H3K4me3 -wao | sort -k4,4 -k11,11nr | awk 'BEGIN{id="";}{OFS="\t";if($4!=id) {tss=$2+3000;peak=$8;dis=peak-tss;if(dis<0) dis=-dis; if(peak<0 || dis>3000) dis=4000; print "H3K4me3","aravin",$4, dis ;id=$4}}' >> ../data/distanceToH3K4peak.tab
grep -P "\.1\t" $annotationFileORF | awk '{OFS="\t"; if($6=="+") print $1,$2,$2,$4,$5,$6; if($6=="-") print $1,$3,$3,$4,$5,$6;}' | slopBed -b 3000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b $H3K4me3 -wao | sort -k4,4 -k11,11nr | awk 'BEGIN{id="";}{OFS="\t";if($4!=id) {tss=$2+3000;peak=$8;dis=peak-tss;if(dis<0) dis=-dis; if(peak<0 || dis>3000) dis=4000; print "H3K4me3","li",$4, dis ;id=$4}}' >> ../data/distanceToH3K4peak.tab


cat $lau | awk '{OFS="\t"; if($6=="+") print $1,$2,$2,$4,$5,$6; if($6=="-") print $1,$3,$3,$4,$5,$6;}' | slopBed -b 3000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b $Amyb -wao | sort -k4,4 -k11,11nr | awk 'BEGIN{id="";}{OFS="\t";if($4!=id) {tss=$2+3000;peak=$8;dis=peak-tss;if(dis<0) dis=-dis;if(peak<0 || dis>3000) dis=4000; print "A-Myb", "lau",$4, dis;id=$4}}' >> ../data/distanceToCAGEPASpeak.tab
cat $girard | awk '{OFS="\t"; if($6=="+") print $1,$2,$2,$4,$5,$6; if($6=="-") print $1,$3,$3,$4,$5,$6;}' | slopBed -b 3000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b $Amyb -wao | sort -k4,4 -k11,11nr | awk 'BEGIN{id="";}{OFS="\t";if($4!=id) {tss=$2+3000;peak=$8;dis=peak-tss;if(dis<0) dis=-dis;if(peak<0 || dis>3000) dis=4000; print "A-Myb", "girard",$4, dis;id=$4}}' >> ../data/distanceToCAGEPASpeak.tab
cat $aravin | awk '{OFS="\t"; if($6=="+") print $1,$2,$2,$4,$5,$6; if($6=="-") print $1,$3,$3,$4,$5,$6;}' | slopBed -b 3000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b $Amyb -wao | sort -k4,4 -k11,11nr | awk 'BEGIN{id="";}{OFS="\t";if($4!=id) {tss=$2+3000;peak=$8;dis=peak-tss;if(dis<0) dis=-dis;if(peak<0 || dis>3000) dis=4000; print "A-Myb", "aravin",$4, dis;id=$4}}' >> ../data/distanceToCAGEPASpeak.tab
grep -P "\.1\t" $annotationFileORF | awk '{OFS="\t"; if($6=="+") print $1,$2,$2,$4,$5,$6; if($6=="-") print $1,$3,$3,$4,$5,$6;}' | slopBed -b 3000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b $Amyb -wao | sort -k4,4 -k11,11nr | awk 'BEGIN{id="";}{OFS="\t";if($4!=id) {tss=$2+3000;peak=$8;dis=peak-tss;if(dis<0) dis=-dis;if(peak<0 || dis>3000) dis=4000; print "A-Myb", "li",$4, dis;id=$4}}' >> ../data/distanceToCAGEPASpeak.tab

cat $lau | awk '{OFS="\t"; if($6=="+") print $1,$2,$2,$4,$5,$6; if($6=="-") print $1,$3,$3,$4,$5,$6;}' | slopBed -b 3000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b $CAGE.narrowpeak -s -wao | sort -k4,4 -k11,11nr | awk 'BEGIN{id="";}{OFS="\t";if($4!=id) {tss=$2+3000;peak=$8+$16;dis=peak-tss;if(dis<0) dis=-dis; if(peak<0 || dis>3000) dis=4000; print "CAGE", "lau",$4, dis;id=$4}}' >> ../data/distanceToCAGEPASpeak.tab
cat $girard | awk '{OFS="\t"; if($6=="+") print $1,$2,$2,$4,$5,$6; if($6=="-") print $1,$3,$3,$4,$5,$6;}' | slopBed -b 3000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b $CAGE.narrowpeak -s -wao | sort -k4,4 -k11,11nr | awk 'BEGIN{id="";}{OFS="\t";if($4!=id) {tss=$2+3000;peak=$8+$16;dis=peak-tss;if(dis<0) dis=-dis; if(peak<0 || dis>3000) dis=4000; print "CAGE","girard",$4, dis ;id=$4}}' >> ../data/distanceToCAGEPASpeak.tab
cat $aravin | awk '{OFS="\t"; if($6=="+") print $1,$2,$2,$4,$5,$6; if($6=="-") print $1,$3,$3,$4,$5,$6;}' | slopBed -b 3000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b $CAGE.narrowpeak -s -wao | sort -k4,4 -k11,11nr | awk 'BEGIN{id="";}{OFS="\t";if($4!=id) {tss=$2+3000;peak=$8+$16;dis=peak-tss;if(dis<0) dis=-dis; if(peak<0 || dis>3000) dis=4000; print "CAGE","aravin",$4, dis ;id=$4}}' >> ../data/distanceToCAGEPASpeak.tab
grep -P "\.1\t" $annotationFileORF | awk '{OFS="\t"; if($6=="+") print $1,$2,$2,$4,$5,$6; if($6=="-") print $1,$3,$3,$4,$5,$6;}' | slopBed -b 3000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b $CAGE.narrowpeak -s -wao | sort -k4,4 -k11,11nr | awk 'BEGIN{id="";}{OFS="\t";if($4!=id) {tss=$2+3000;peak=$8+$16;dis=peak-tss;if(dis<0) dis=-dis;if(peak<0 || dis>3000) dis=4000; print "CAGE","li",$4, dis ;id=$4}}' >> ../data/distanceToCAGEPASpeak.tab

cat $lau | awk '{OFS="\t"; if($6=="-" || $6==".") print $1,$2,$2,$4,$5,"-"; if($6=="+" || $6==".") print $1,$3,$3,$4,$5,"+"; }' | slopBed -b 3000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b $PAS.narrowpeak -s -wao | sort -k4,4 -k11,11nr | awk 'BEGIN{id="";}{OFS="\t";if($4!=id) {tss=$2+3000;peak=$8+$16;dis=peak-tss;if(dis<0) dis=-dis; if(peak<0 || dis>3000) dis=4000;print "PAS","lau",$4, dis ;id=$4}}' >> ../data/distanceToCAGEPASpeak.tab
cat $girard | awk '{OFS="\t"; if($6=="-" || $6==".") print $1,$2,$2,$4,$5,"-"; if($6=="+" || $6==".") print $1,$3,$3,$4,$5,"+"; }' | slopBed -b 3000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b $PAS.narrowpeak -s -wao | sort -k4,4 -k11,11nr | awk 'BEGIN{id="";}{OFS="\t";if($4!=id) {tss=$2+3000;peak=$8+$16;dis=peak-tss;if(dis<0) dis=-dis; if(peak<0 || dis>3000) dis=4000;print "PAS","girard",$4, dis ;id=$4}}' >> ../data/distanceToCAGEPASpeak.tab
cat $aravin | awk '{OFS="\t"; if($6=="-" || $6==".") print $1,$2,$2,$4,$5,"-"; if($6=="+" || $6==".") print $1,$3,$3,$4,$5,"+"; }' | slopBed -b 3000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b $PAS.narrowpeak -s -wao | sort -k4,4 -k11,11nr | awk 'BEGIN{id="";}{OFS="\t";if($4!=id) {tss=$2+3000;peak=$8+$16;dis=peak-tss;if(dis<0) dis=-dis; if(peak<0 || dis>3000) dis=4000;print "PAS","aravin",$4, dis ;id=$4}}' >> ../data/distanceToCAGEPASpeak.tab
grep -P "\.1\t" $annotationFileORF | awk '{OFS="\t"; if($6=="-" || $6==".") print $1,$2,$2,$4,$5,"-"; if($6=="+" || $6==".") print $1,$3,$3,$4,$5,"+"; }' | slopBed -b 3000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b $PAS.narrowpeak -s -wao | sort -k4,4 -k11,11nr | awk 'BEGIN{id="";}{OFS="\t";if($4!=id) {tss=$2+3000;peak=$8+$16;dis=peak-tss;if(dis<0) dis=-dis; if(peak<0 || dis>3000) dis=4000;print "PAS","li",$4, dis ;id=$4}}' >> ../data/distanceToCAGEPASpeak.tab


echo "
df=read.table('../data/distanceToCAGEPASpeak.tab', header=F)
colnames(df)=c('type','source','name','distance')
df$source=as.character(df$source)
df$type=as.character(df$type)
df$source[df$source=='li']='This Study'
df$source[df$source=='lau']='Lau et al. 2006'
df$source[df$source=='girard']='Girard et al. 2006'
df$source[df$source=='aravin']='Aravin et al. 2007'
df$source=factor(df$source, levels=c('This Study','Lau et al. 2006','Girard et al. 2006', 'Aravin et al. 2007'))
df$type=factor(df$type, levels=c('H3K4me3','A-Myb', 'CAGE','PAS'))

pdf('Fig2B.distanceToCAGEPASpeak.pdf',width=9, height=4, colormodel='cmyk')
par(lwd=2/3, mfrow=c(2,1))
library(ggplot2)
ggplot(df, aes(x=distance, fill=source)) + xlab('Distance to the peak summit (bp)') + geom_bar(position='dodge', binwidth=100) + facet_wrap(~ type, ncol=2) + aes(y = ..density..)
ggplot(df, aes(x=distance, fill=source)) + xlab('Distance to the peak summit (bp)') + geom_bar(position='dodge', binwidth=100) + facet_wrap(~ type, ncol=2)
dev.off()

pdf('Fig2S.distanceToCAGEPASpeak.separated.pdf',width=12, height=12, colormodel='cmyk')
par(lwd=2/3, mfrow=c(4,4))
hist(df$distance[grepl('This', df$source) & df$type=='H3K4me3'], breaks=30, main='', xlab='Distance to H3K4me3 peaks (this study)');
hist(df$distance[grepl('This', df$source) & df$type=='A-Myb'], breaks=30, main='', xlab='Distance to A-Myb peaks (this study)');
hist(df$distance[grepl('This', df$source) & df$type=='CAGE'], breaks=30, main='', xlab='Distance to CAGE peaks (this study)');
hist(df$distance[grepl('This', df$source) & df$type=='PAS'], breaks=30, main='', xlab='Distance to PAS peaks (this study)');

hist(df$distance[grepl('Lau', df$source) & df$type=='H3K4me3'], breaks=30, main='', xlab='Distance to H3K4me3 peaks (Lau et al. 2006)');
hist(df$distance[grepl('Lau', df$source) & df$type=='A-Myb'], breaks=30, main='', xlab='Distance to A-Myb peaks (Lau et al. 2006)');
hist(df$distance[grepl('Lau', df$source) & df$type=='CAGE'], breaks=30, main='', xlab='Distance to CAGE peaks (Lau et al. 2006)');
hist(df$distance[grepl('Lau', df$source) & df$type=='PAS'], breaks=30, main='', xlab='Distance to PAS peaks (Lau et al. 2006)');

hist(df$distance[grepl('Girard', df$source) & df$type=='H3K4me3'], breaks=30, main='', xlab='Distance to H3K4me3 peaks (Girard et al. 2006)');
hist(df$distance[grepl('Girard', df$source) & df$type=='A-Myb'], breaks=30, main='', xlab='Distance to A-Myb peaks (Girard et al. 2006)');
hist(df$distance[grepl('Girard', df$source) & df$type=='CAGE'], breaks=30, main='', xlab='Distance to CAGE peaks (Girard et al. 2006)');
hist(df$distance[grepl('Girard', df$source) & df$type=='PAS'], breaks=30, main='', xlab='Distance to PAS peaks (Girard et al. 2006)');

hist(df$distance[grepl('Aravin', df$source) & df$type=='H3K4me3'], breaks=30, main='', xlab='Distance to H3K4me3 peaks (Aravin et al. 2007)');
hist(df$distance[grepl('Aravin', df$source) & df$type=='A-Myb'], breaks=30, main='', xlab='Distance to A-Myb peaks (Aravin et al. 2007)');
hist(df$distance[grepl('Aravin', df$source) & df$type=='CAGE'], breaks=30, main='', xlab='Distance to CAGE peaks (Aravin et al. 2007)');
hist(df$distance[grepl('Aravin', df$source) & df$type=='PAS'], breaks=30, main='', xlab='Distance to PAS peaks (Aravin et al. 2007)');

dev.off()

" > /tmp/Fig2B.R
Rscript /tmp/Fig2B.R

##################################################
# CAGE/PAS tag clustering
##################################################
cat $CAGE.plus.bedGraph  | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); sum=0; for(i=1;i<=length(b);i++) sum+=b[i]; cum=0; for(i=1;i<=length(b);i++) {cum+=b[i]; split(a[i],c,":|-"); if(cum>=(0.1*sum) && cum<(0.9*sum)) print c[1], c[2],c[3], a[i], b[i]; if(cum>=(0.9*sum)) {print c[1], c[2],c[3], a[i], b[i]; break;}}}' | mergeBed -d $D -scores collapse -nms | awk -v strand="plus" -f cage2narrowPeak > $CAGE.narrowpeak
cat $CAGE.minus.bedGraph  | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, -$4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); sum=0; for(i=1;i<=length(b);i++) sum+=b[i]; cum=0; for(i=1;i<=length(b);i++) {cum+=b[i]; split(a[i],c,":|-"); if(cum>=(0.1*sum) && cum<(0.9*sum)) print c[1], c[2],c[3], a[i], b[i]; if(cum>=(0.9*sum)) {print c[1], c[2],c[3], a[i], b[i]; break;}}}' | mergeBed -d $D -scores collapse -nms | awk -v strand="minus" -f cage2narrowPeak >> $CAGE.narrowpeak

#cat $CAGE.plus.bedGraph | intersectBed -a stdin -b <(awk '{OFS="\t"; tss=($6=="+")?$2:$3; if($1!~/random/) print $1,tss,tss,"pi_"NR,1000,$6;}' $inputfile2 | slopBed -b 2000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt) -wa -u | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); sum=0; for(i=1;i<=length(b);i++) sum+=b[i]; cum=0; for(i=1;i<=length(b);i++) {cum+=b[i]; split(a[i],c,":|-"); if(cum>=(0.1*sum) && cum<(0.9*sum)) print c[1], c[2],c[3], a[i], b[i]; if(cum>=(0.9*sum)) {print c[1], c[2],c[3], a[i], b[i]; break;}}}' | mergeBed -d $D -scores collapse -nms | awk -v strand="plus" -f cage2narrowPeak > $CAGE.piRNA.narrowpeak
#cat $CAGE.minus.bedGraph | intersectBed -a stdin -b <(awk '{OFS="\t"; tss=($6=="+")?$2:$3; if($1!~/random/) print $1,tss,tss,"pi_"NR,1000,$6;}' $inputfile2 | slopBed -b 2000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt) -wa -u | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, -$4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); sum=0; for(i=1;i<=length(b);i++) sum+=b[i]; cum=0; for(i=1;i<=length(b);i++) {cum+=b[i]; split(a[i],c,":|-"); if(cum>=(0.1*sum) && cum<(0.9*sum)) print c[1], c[2],c[3], a[i], b[i]; if(cum>=(0.9*sum)) {print c[1], c[2],c[3], a[i], b[i]; break;}}}' | mergeBed -d $D -scores collapse -nms | awk -v strand="minus" -f cage2narrowPeak >> $CAGE.piRNA.narrowpeak

# PAS peak calling
cat $PAS.plus.bedGraph  | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); sum=0; for(i=1;i<=length(b);i++) sum+=b[i]; cum=0; for(i=1;i<=length(b);i++) {cum+=b[i]; split(a[i],c,":|-"); if(cum>=(0.1*sum) && cum<(0.9*sum)) print c[1], c[2],c[3], a[i], b[i]; if(cum>=(0.9*sum)) {print c[1], c[2],c[3], a[i], b[i]; break;}}}' | mergeBed -d $D -scores collapse -nms | awk -v strand="plus" -f cage2narrowPeak > $PAS.narrowpeak
cat $PAS.minus.bedGraph  | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, -$4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); sum=0; for(i=1;i<=length(b);i++) sum+=b[i]; cum=0; for(i=1;i<=length(b);i++) {cum+=b[i]; split(a[i],c,":|-"); if(cum>=(0.1*sum) && cum<(0.9*sum)) print c[1], c[2],c[3], a[i], b[i]; if(cum>=(0.9*sum)) {print c[1], c[2],c[3], a[i], b[i]; break;}}}' | mergeBed -d $D -scores collapse -nms | awk -v strand="minus" -f cage2narrowPeak >> $PAS.narrowpeak

#cat $PASallmapper.plus.bedGraph  | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); sum=0; for(i=1;i<=length(b);i++) sum+=b[i]; cum=0; for(i=1;i<=length(b);i++) {cum+=b[i]; split(a[i],c,":|-"); if(cum>=(0.1*sum) && cum<(0.9*sum)) print c[1], c[2],c[3], a[i], b[i]; if(cum>=(0.9*sum)) {print c[1], c[2],c[3], a[i], b[i]; break;}}}' | mergeBed -d $D -scores collapse -nms | awk -v strand="plus" -f cage2narrowPeak > $PASallmapper.narrowpeak
#cat $PASallmapper.minus.bedGraph  | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, -$4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); sum=0; for(i=1;i<=length(b);i++) sum+=b[i]; cum=0; for(i=1;i<=length(b);i++) {cum+=b[i]; split(a[i],c,":|-"); if(cum>=(0.1*sum) && cum<(0.9*sum)) print c[1], c[2],c[3], a[i], b[i]; if(cum>=(0.9*sum)) {print c[1], c[2],c[3], a[i], b[i]; break;}}}' | mergeBed -d $D -scores collapse -nms | awk -v strand="minus" -f cage2narrowPeak >> $PASallmapper.narrowpeak

##################################################
# CAGE/PAS tag correlation
##################################################
# correlation of PAS vs. CAGE for mRNA and piRNA
rm ../data/controlSet2/*summit
for i in ../data/controlSet2/*bed*;
do
    cat $i | awk '{OFS="\t"; s=($6=="+")?$2:$3; print $1,s-1000,s+1000,$4"___"$5"___"$9,1000, $6;}'  | intersectBed -a stdin -b $CAGE.narrowpeak -s -wao | sort -k4,4 -k11,11nr | groupBy -g 4 -c 11 -o max > /tmp/CAGE.summit
    echo $i CAGE
    cat $i | awk '{OFS="\t"; s=($6=="+")?$3:$2; print $1,s-1000,s+1000,$4"___"$5"___"$9,1000, $6;}'  | intersectBed -a stdin -b $PAS.narrowpeak -s -wao | sort -k4,4 -k11,11nr | groupBy -g 4 -c 11 -o max | paste /tmp/CAGE.summit - | cut -f1,2,4 > $i.CAGEPAS.summit
    echo $i PAS
done
echo "
pdf('/tmp/CAGE.PAS.correlation.pdf', width=5, height=25)
format.value <- function(R=R, where='slide'){
             if(where=='slide') return(ifelse(R<2.2e-26, '2.2e-16', ifelse(R<0.001, format(R, digits=2, scienticif=TRUE), format(R, digits=2))));
             if(where=='table') return(ifelse(R<2.2e-26, '2.2e-16', ifelse(R<0.001, format(R, digits=2, scienticif=TRUE), format(R, digits=2))));
}
par(mfrow=c(5,1))
# piRNA-all
df=rbind(read.table('../data/controlSet2/piRNA.genic.bed.CAGEPAS.summit'), read.table('../data/controlSet2/piRNA.intergenic.bed.CAGEPAS.summit'))
rownames(df)=df[,1]; df=df[,-1]; colnames(df)=c('CAGE','PAS'); df=df+2;
plot(df, col='#0000ff22', log='xy', asp=1, pch=16, main='piRNA-all')
legend('topright',c(paste('Spearman\'s rho =', round(cor(df[,1],df[,2], method='spearman'),3)), paste('Pearson\'s r =', round(cor(df[,1],df[,2], method='pearson'),3)), paste('Wilcox test p-value = ', format.value(wilcox.test(df[,1],df[,2], paired=T)\$p.value))))
# piRNA-genic
df=read.table('../data/controlSet2/piRNA.genic.bed.CAGEPAS.summit')
rownames(df)=df[,1]; df=df[,-1]; colnames(df)=c('CAGE','PAS'); df=df+2
plot(df, col='#0000ff22', log='xy', asp=1, pch=16, main='piRNA-genic')
legend('topright',c(paste('Spearman\'s rho =', round(cor(df[,1],df[,2], method='spearman'),3)), paste('Pearson\'s r =', round(cor(df[,1],df[,2], method='pearson'),3)), paste('Wilcox test p-value = ', format.value(wilcox.test(df[,1],df[,2], paired=T)\$p.value))))
# piRNA-intergenic
df=read.table('../data/controlSet2/piRNA.intergenic.bed.CAGEPAS.summit')
rownames(df)=df[,1]; df=df[,-1]; colnames(df)=c('CAGE','PAS'); df=df+2
plot(df, col='#0000ff22', log='xy', asp=1, pch=16, main='piRNA-intergenic')
legend('topright',c(paste('Spearman\'s rho =', round(cor(df[,1],df[,2], method='spearman'),3)), paste('Pearson\'s r =', round(cor(df[,1],df[,2], method='pearson'),3)), paste('Wilcox test p-value = ', format.value(wilcox.test(df[,1],df[,2], paired=T)\$p.value))))
# NM
df=read.table('../data/controlSet2/x100910.rnaseq.transcripts.NM.all.final.bed.c.CAGEPAS.summit')
rownames(df)=df[,1]; df=df[,-1]; colnames(df)=c('CAGE','PAS'); df=df+2
plot(df, col='#0000ff22', log='xy', asp=1, pch=16, main='mRNA')
legend('topright',c(paste('Spearman\'s rho =', round(cor(df[,1],df[,2], method='spearman'),3)), paste('Pearson\'s r =', round(cor(df[,1],df[,2], method='pearson'),3)), paste('Wilcox test p-value = ', format.value(wilcox.test(df[,1],df[,2], paired=T)\$p.value))))
# NR
df=read.table('../data/controlSet2/x100910.rnaseq.transcripts.NR.all.final.bed.c.CAGEPAS.summit')
rownames(df)=df[,1]; df=df[,-1]; colnames(df)=c('CAGE','PAS'); df=df+2
plot(df, col='#0000ff22', log='xy', asp=1, pch=16, main='ncRNA')
legend('topright',c(paste('Spearman\'s rho =', round(cor(df[,1],df[,2], method='spearman'),3)), paste('Pearson\'s r =', round(cor(df[,1],df[,2], method='pearson'),3)), paste('Wilcox test p-value = ', format.value(wilcox.test(df[,1],df[,2], paired=T)\$p.value))))
dev.off()
" > /tmp/CAGE.PAS.correlation.R
Rscript /tmp/CAGE.PAS.correlation.R; scp /tmp/CAGE.PAS.correlation.pdf zlab:~/public_html/temp



awk '{OFS="\t"; if($1!~/random/) print $1,$2,$3,"nonpi_"NR,1000,$6, $0}' $inputfile1 > ../data/x100910.rnaseq.transcripts.bed
awk '{OFS="\t"; if($1!~/random/) print $1,$2,$3,"pi_"NR,1000,$6, $0}' $inputfile2 > ../data/100312.transcript.with.transcripts.removed.new.names.bed

intersectBed -a ../data/x100910.rnaseq.transcripts.bed -b ../data/100312.transcript.with.transcripts.removed.new.names.bed -s -v | cut -f1-6 | cat - <(cut -f1-6 ../data/100312.transcript.with.transcripts.removed.new.names.bed) | sort -k1,1 -k2,2n > ../data/x100910.rnaseq.transcripts.ALL.bed


# 1. TSS-flanking region [-500, +50], and TTS-flanking region [-50,+500]
grep -v track $annotationFile_all | awk '{OFS="\t"; tss=($6=="+")?$2:$3; print $1,tss,tss,$4,$5,$6;}' | slopBed -l $FLANKINGup -r $FLANKINGdown -s -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt > $TSSregions
grep -v track $annotationFile_all | awk '{OFS="\t"; tts=($6=="+")?$3:$2; print $1,tts,tts,$4,$5,$6;}' | slopBed -l $FLANKINGdown -r $FLANKINGup -s -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt > $TTSregions

# 2. BIP : two Tx with different orientation, and distance of TSS is between [-100,+1000]
# only for piRNA-BIP
grep -v track $annotationFile | awk '{OFS="\t"; tss=($6=="+")?$2:$3; print $1,tss,tss,$4,$5,$6;}' | slopBed -l 1000 -r 100 -s -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt

intersectBed -a <(awk '{OFS="\t"; if($6=="-")print}' $TSSregions) -b <(awk '{OFS="\t"; if($6=="+")print}' $TSSregions) -S -wo | awk -v up=$FLANKINGup -v down=$FLANKINGdown '{OFS="\t"; s=$2+down; e=$8+up; print $4,$1":"((s<e)?s:e)"-"((s<e)?e:s), e-s, $10, (($4~/^pi/)?"pi":"pc")":"(($10~/^pi/)?"pi":"pc"); print $10,$1":"((s<e)?s:e)"-"((s<e)?e:s), e-s, $4, (($10~/^pi/)?"pi":"pc")":"(($4~/^pi/)?"pi":"pc");}' | sort -k1,1 -k3,3n | awk 'BEGIN{ID=""}{if($1!=ID) {print; ID=$1;}}' > $annotationFile_all.BIP

intersectBed -a <(awk '{OFS="\t"; if($6=="-")print}' $TSSregions) -b <(awk '{OFS="\t"; if($6=="+")print}' $TSSregions) -S -wo | awk -v up=$FLANKINGup -v down=$FLANKINGdown '{OFS="\t"; s=$2+down; e=$8+up; if($4~/^pi/) print $4,$1":"((s<e)?s:e)"-"((s<e)?e:s), e-s,$10, "pi:"(($10~/^pi/)?"pi":"pc"); if($10~/^pi/) print $10,$1":"((s<e)?s:e)"-"((s<e)?e:s),e-s, $4, "pi:"(($4~/^pi/)?"pi":"pc");}' | sort -k1,1 -k3,3n | awk 'BEGIN{ID=""}{if($1!=ID) {OFS="\t"; print; ID=$1;}}' > $annotationFile.BIP

# BIP length distribution
textHistogram -col=2 -maxBinCount=1200 -binSize=20 <(cut -f2,3 $annotationFile.BIP | sort -u) > $annotationFile.BIP.lengthstat
textHistogram -col=2 -maxBinCount=1200 -binSize=20 <(awk '{if($1~/^nonpi/ && $4~/^nonpi/) print $2,$3}' $annotationFile_all.BIP | sort -u) > $annotationFile_all.BIP.lengthstat

# take uniq (sort -u before awk)
cut -f6,13,19-20 ../data/TSS_table_piRNA_adult_mm9.tab | grep "chr" | awk '{OFS="\t";if($4=="pi:pi") $1="+"; print}' | sort -u | awk '{OFS="\t"; split($3,a,":|-"); print a[1], a[2], a[3], "piBIP_"NR, a[3]-a[2], $1, $2, $4}' > ../data/TSS_table_piRNA_adult_mm9.tab.BIP
cut -f6,13,19-20 ../data/TSS_table_piRNA_adult_mm9.tab | grep "chr" | awk '{OFS="\t";if($4=="pi:pi") $1="+"; print}' | sort -u | awk '{OFS="\t"; split($3,a,":|-"); print a[1], int((a[3]+a[2])/2), int((a[3]+a[2])/2), "piBIP_"NR, 1000,$1}' | slopBed -b 1000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo stdout | perl -ne '$h=$_; $h=~s/^>//; chomp($h); $s=<>; @n1 = ($s=~/CG/gi); $n1=@n1; $n2 = length($s); @n3 = (($s=~/C/gi), ($s=~/G/gi)); $n3=@n3; printf("%s\t%.3f\n", $h, 4*$n2*$n1/($n3**2+1));' | paste ../data/TSS_table_piRNA_adult_mm9.tab.BIP - | cut -f1-8,10 > ../data/TSS_table_piRNA_adult_mm9.tab.BIP.cpg.bed

### pc:pc
cut -f2,3,5 $annotationFile_all.BIP | sort -u | awk '{OFS="\t"; if($3=="pc:pc") {split($1,a,":|-"); print a[1], a[2],a[3],"pcBIP_"NR, a[3]-a[2],"+", "all", "pc:pc";}}' > $annotationFile_all.pcpc.BIP
cut -f2,3,5 $annotationFile_all.BIP | sort -u | awk '{OFS="\t"; if($3=="pc:pc") {split($1,a,":|-"); print a[1], int((a[3]+a[2])/2), int((a[3]+a[2])/2), "pcBIP_"NR, 1000,"+";}}'  | slopBed -b 1000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo stdout | perl -ne '$h=$_; $h=~s/^>//; chomp($h); $s=<>; @n1 = ($s=~/CG/gi); $n1=@n1; $n2 = length($s); @n3 = (($s=~/C/gi), ($s=~/G/gi)); $n3=@n3; printf("%s\t%.3f\n", $h, 4*$n2*$n1/($n3**2+1));' | paste $annotationFile_all.pcpc.BIP - | cut -f1-8,10 > $annotationFile_all.pcpc.BIP.cpg.bed
# divide into NM, NR
cat ../data/controlSet2/x100910.rnaseq.transcripts.NM.all.final.bed | awk '{OFS="\t"; tss=($6=="+")?$2:$3; if($1!~/random/) print $1,tss,tss,$4,$5,$6;}' | slopBed -l $FLANKINGup -r $FLANKINGdown -s -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a $annotationFile_all.pcpc.BIP.cpg.bed -b stdin -wa -u | awk '{OFS="\t"; $7="NM"; print }' >  /tmp/a
cat ../data/controlSet2/x100910.rnaseq.transcripts.NR.final.bed | awk '{OFS="\t"; tss=($6=="+")?$2:$3; if($1!~/random/) print $1,tss,tss,$4,$5,$6;}' | slopBed -l $FLANKINGup -r $FLANKINGdown -s -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a $annotationFile_all.pcpc.BIP.cpg.bed -b stdin -wa -u | awk '{OFS="\t"; $7="NR"; print }' >>  /tmp/a
sort -k4,4 /tmp/a | groupBy -g 1,2,3,4,5,6,8,9 -c 7 -o collapse | awk '{OFS="\t"; print $1, $2, $3,$4,$5,$6,$9,$7,$8;}' > $annotationFile_all.pcpc.BIP.cpg.bed.NMNR

## pc + piRNA
cat ../data/TSS_table_piRNA_adult_mm9.tab.BIP.cpg.bed $annotationFile_all.pcpc.BIP.cpg.bed.NMNR > ../data/BIP.all.piNMNR.cpg.bed
# get binned signal for aggregation
awk '{OFS="\t"; print $1,$2,$3,$6,$4}' ../data/BIP.all.piNMNR.cpg.bed | awk -f bigWigAverageOverBed_generate_81bins.awk > ../data/BIP.all.piNMNR.cpg.bed.81bins
qsub bigWigAverageOverBed_81bins.sh ../data/BIP.all.piNMNR.cpg.bed.81bins 81
Rscript draw_aggregation_plot.R BIP

## nucleotide frequency of 150nt region of BIP
cut -f13,19-20 ../data/TSS_table_piRNA_adult_mm9.tab | grep "chr" | sort -u | awk '{OFS="\t"; split($2,a,":|-"); print a[1], int((a[3]+a[2])/2), int((a[3]+a[2])/2), "piBIP_"NR, 1000,"+";}' | slopBed -b 150 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo ../data/TSS_table_piRNA_adult_mm9.tab.BIP.150nt.fa
cut -f13,19-20 ../data/TSS_table_piRNA_adult_mm9.tab | grep "chr" | grep "pi:pc" | sort -u | awk '{OFS="\t"; split($2,a,":|-"); print a[1], int((a[3]+a[2])/2), int((a[3]+a[2])/2), "piBIP_"NR, 1000,"+";}' | slopBed -b 150 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo ../data/TSS_table_piRNA_adult_mm9.tab.BIP.pipc.150nt.fa
cut -f13,19-20 ../data/TSS_table_piRNA_adult_mm9.tab | grep "chr" | grep "pi:pi" | sort -u | awk '{OFS="\t"; split($2,a,":|-"); print a[1], int((a[3]+a[2])/2), int((a[3]+a[2])/2), "piBIP_"NR, 1000,"+";}' | slopBed -b 150 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo ../data/TSS_table_piRNA_adult_mm9.tab.BIP.pipi.150nt.fa
cut -f13,19-20 ../data/TSS_table_piRNA_adult_mm9.tab | grep "chr" | grep "pi:nc" | sort -u | awk '{OFS="\t"; split($2,a,":|-"); print a[1], int((a[3]+a[2])/2), int((a[3]+a[2])/2), "piBIP_"NR, 1000,"+";}' | slopBed -b 150 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo ../data/TSS_table_piRNA_adult_mm9.tab.BIP.pinc.150nt.fa
cut -f2,3,5 $annotationFile_all.BIP | sort -u | awk '{OFS="\t"; if($3=="pc:pc") {split($1,a,":|-"); print a[1], int((a[3]+a[2])/2), int((a[3]+a[2])/2), "pcBIP_"NR, 1000,"+";}}'  | slopBed -b 150 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo $annotationFile_all.BIP.150nt.fa

# motif analysis of BIP region
echo "
#!/bin/sh
#$ -V
#$ -cwd
#$ -pe openmpi 4
#$ -o \$HOME/sge_jobs_output/sge_job.\$JOB_ID.out -j y
#$ -S /bin/bash
#$ -l mem_free=2G

maxsize=\`wc -c \$1 | cut -f1 -d' '\`

id=\`basename \$1\`
meme_p -oc meme_\$id -mod zoops -nmotifs 3 -minw 6 -maxw 20 -revcomp -nostatus -dna -p 4 -maxsize \$maxsize \$1
tomtom -oc tomtom_\$id -min-overlap 5 -verbosity 1 meme_\$id/meme.txt motifdb/JASPAR_CORE_2009.meme motifdb/transfac.meme
fimo --text --verbosity 1 --output-pthresh 0.005 meme_\$id/meme.txt \$1 > \$1.fimo
" > /tmp/meme.sge
qsub /tmp/meme.sge ../data/TSS_table_piRNA_adult_mm9.tab.BIP.pipc.150nt.fa
qsub /tmp/meme.sge ../data/TSS_table_piRNA_adult_mm9.tab.BIP.pipi.150nt.fa
qsub /tmp/meme.sge ../data/TSS_table_piRNA_adult_mm9.tab.BIP.pinc.150nt.fa
qsub /tmp/meme.sge ../data/TSS_table_piRNA_adult_mm9.tab.BIP.150nt.fa
qsub /tmp/meme.sge $annotationFile_all.BIP.150nt.fa

echo "
args <- commandArgs(TRUE)
pdf('BIP.motif.density.pdf', width=20, height=6)
df=read.table(args[1])
plot(df[,1],df[,2], type='h')
abline(v=150, col='gray')
dev.off()
" > /tmp/BIPmotifdraw.R

fimo --text --verbosity 1 --output-pthresh 0.005 --motif 1 meme_TSS_table_piRNA_adult_mm9.tab.BIP.150nt.fa/meme.txt ../data/TSS_table_piRNA_adult_mm9.tab.BIP.150nt.fa | cut -f3-6 | grep -v "Start" | sort -k1,1n -k2,2n -k3,3 | groupBy -g 1,2,3 -c 4 -o sum | awk '{OFS="\t"; if($3=="+") print $1,$4; if($3=="-") print $2,-$4;}' | Rscript /tmp/BIPmotifdraw.R stdin
scp BIP.motif.density.pdf zlab:~/public_html/temp

# 3. CpG score of [-1000,+1000] region of TSS
grep -v track $annotationFile | awk '{OFS="\t"; tss=($6=="+")?$2:$3; print $1, tss, tss,$4,$5,$6;}' | slopBed -b 1000 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo $annotationFile.cpg.fa
cat $annotationFile.cpg.fa | perl -ne '$h=$_; $h=~s/^>//; chomp($h); $s=<>; @n1 = ($s=~/CG/gi); $n1=@n1; $n2 = length($s); @n3 = (($s=~/C/gi), ($s=~/G/gi)); $n3=@n3; printf("%s\t%.3f\n", $h, 4*$n2*$n1/($n3**2+1));' > $annotationFile.cpg

# 4. TATA : scan the region -35bp to -22bp (expected position for the TATA-box) upstream of TSS detected by CAGE.
awk -v flk=1000 '{if(NR%2==1) print; if(NR%2==0) print substr($1,flk-40,20);}' $annotationFile.cpg.fa | fimo --text --verbosity 1 --output-pthresh 0.005 ../data/motif/TATA.meme - | awk  '{OFS="\t"; if($5=="+") print $2,($3-40)"."$8"."($4-40), $6}' | sort -k1,1 -k3,3gr | groupBy -g 1 -c 2 -o collapse > $annotationFile.TATA

# 5. T2T : two Tx with different orientation, and distance of TTS is between [-100,+1000]
intersectBed -a <(awk '{OFS="\t"; if($6=="+")print}' $TTSregions) -b <(awk '{OFS="\t"; if($6=="-")print}' $TTSregions) -S -wo | awk -v up=$FLANKINGup -v down=$FLANKINGdown '{OFS="\t"; s=$2+down; e=$8+up; print $4,$1":"((s<e)?s:e)"-"((s<e)?e:s), e-s; print $10,$1":"((s<e)?s:e)"-"((s<e)?e:s), e-s;}' | sort -k1,1 -k3,3n | awk 'BEGIN{ID=""}{if($1!=ID) {print; ID=$1;}}' > $annotationFile.T2T

# T2T length distribution
textHistogram -col=3 -maxBinCount=1200 -binSize=20 $annotationFile.T2T > $annotationFile.T2T.lengthstat

##################################################
### PAS
##################################################
# 6. motif : PAS signal
# get sequence of [-50,0] of 3'
grep -v track $annotationFile | awk '{OFS="\t"; tts=($6=="-")?$2:$3; print $1, tts, tts,$4,$5,$6;}' | slopBed -b 50 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo $annotationFile.pas.fa
grep -v track $annotationFile | awk '{OFS="\t"; tts=($6=="-")?$2:$3; if($10~/^pi/) print $1, tts, tts,$4,$5,$6;}' | slopBed -b 50 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo $annotationFile.pas.genic.fa
grep -v track $annotationFile | awk '{OFS="\t"; tts=($6=="-")?$2:$3; if($10!~/^pi/) print $1, tts, tts,$4,$5,$6;}' | slopBed -b 50 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo $annotationFile.pas.intergenic.fa

# pas for pc genes

# get sequence of [-50,50] of top 1000 PAS peaks
# only NR with overlaps with any PAS peaks
cat $PAS.narrowpeak | awk '{OFS="\t"; print $1,$2+$10,$2+$10,$4,$5,$6; }' | slopBed -b 50 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b <(awk '{OFS="\t"; tts=($6=="+")?$3:$2; print $1,tts,tts,$9,0,$6;}' $NR | slopBed -b 50 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt) -s -wo | sort -k10,10 -k13,13nr | awk 'BEGIN{id="";}{if($10!=id) {print;id=$10;}}' | cut -f1-6 | sort -u| fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo ../data/mouse_adult_wt_PAS200_100nt.unique.all.50PAS50.inNR.fa
# only genes overlap with the top1000 PAS peaks
cat $PAS.narrowpeak | sort -k5,5nr -k7,7nr | head -n1000 | awk '{OFS="\t"; print $1,$2+$10,$2+$10,$4,$5,$6; }' | slopBed -b 50 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b <(awk '{OFS="\t"; tts=($6=="+")?$3:$2; print $1,tts,tts,$9,0,$6;}' $NM | slopBed -b 50 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt) -s -wo | sort -k10,10 -k13,13nr | awk 'BEGIN{id="";}{if($10!=id) {print;id=$10;}}' | cut -f1-6 | sort -u| fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo ../data/mouse_adult_wt_PAS200_100nt.unique.top1000.50PAS50.inNM.fa
# all NM
cat $PAS.narrowpeak | awk '{OFS="\t"; print $1,$2+$10,$2+$10,$4,$5,$6; }' | slopBed -b 50 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b <(awk '{OFS="\t"; tts=($6=="+")?$3:$2; print $1,tts,tts,$9,0,$6;}' $NM | slopBed -b 50 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt) -s -wo | sort -k10,10 -k13,13nr | awk 'BEGIN{id="";}{if($10!=id) {print;id=$10;}}' | cut -f1-6 | sort -u| fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo ../data/mouse_adult_wt_PAS200_100nt.unique.all.50PAS50.inNM.fa &

#cat $PAS.narrowpeak | sort -k5,5nr -k7,7nr | awk '{OFS="\t"; print $1,$2+$10,$2+$10,$4,$5,$6; }' | slopBed -b 50 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b $NR -s -wa -u | fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo ../data/mouse_adult_wt_PAS200_100nt.unique.all.50PAS50.inNR.fa
#cat $PAS.narrowpeak | sort -k5,5nr -k7,7nr | awk '{OFS="\t"; print $1,$2+$10,$2+$10,$4,$5,$6; }' | slopBed -b 50 -g $GENOME/mm9/Annotation/Genes/ChromInfo.txt | intersectBed -a stdin -b $NM -s -wa -u | fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo ../data/mouse_adult_wt_PAS200_100nt.unique.all.50PAS50.inNM.fa

memepas() {
    fa=$1;
    maxsize=`wc -c $fa | cut -f1 -d' '`
    id=`basename $fa`
    meme -oc meme_$id -w 6 -mod zoops -nmotifs 1 -dna -maxsize $maxsize $fa
    fimo --text --verbosity 1 --output-pthresh 0.005 --motif 1 meme_$id/meme.txt $fa | awk '{if($5=="+") print;}' > $fa.PAS
}
memepas $annotationFile.pas.fa &
memepas $annotationFile.pas.intergenic.fa &
memepas $annotationFile.pas.genic.fa &
memepas ../data/mouse_adult_wt_PAS200_100nt.unique.top1000.50PAS50.inNR.fa &
memepas ../data/mouse_adult_wt_PAS200_100nt.unique.top1000.50PAS50.inNM.fa &
memepas ../data/mouse_adult_wt_PAS200_100nt.unique.all.50PAS50.inNM.fa &

echo "
args <- commandArgs(TRUE);
pdf('PAS.position.pdf',width=8, height=5, colormodel='cmyk');
df=read.table(args[1]);
N=as.numeric(args[2]);
motifs=as.character(unique(df\$V1));
cls=colors()[c(552,81,78, 40,90,117,11,259,10,86,12,84,94,122,200,51,139,100)];
mtfs=c('AAUAAA','AUUAAA','AGUAAA','AAGAAA','AAUACA','AGAAAA','ACUAAA','AAAAAA','CAUAAA','AACAAA','UAUAAA','AGGAAA','AUGAAA','AUAAAA','AAUAGA','GAUAAA','AAUAAG','AUCAAA')
plot(df[df\$V1==motifs[1],2],df[df\$V1==motifs[1],3],type='h',xlim=c(1,50),ylim=range(0, df[,3]), xaxt='n',las=1, tck=0.01, lwd=2,col=cls[1],xlab='Position relative to 3\' end of RNA', ylab='Motif occurrence');
n=sum(df[df\$V1==motifs[1],3]);
cl=cls[1];
j=1;
mtf=mtfs[1];
for(i in 2:length(mtfs)){
    if(mtfs[i] %in% motifs) {
        lines(df[df\$V1==mtfs[i],2]+0.2*j,df[df\$V1==mtfs[i],3], type='h',col=cls[i],lwd=2);
        n=c(n,sum(df[df\$V1==mtfs[i],3]))
        cl=c(cl,cls[i]);
        mtf=c(mtf,mtfs[i])
        j=j+1;
    }
}
n0=c(n,N-sum(n))
cl=c(cl,'black')
axis(1,at=seq(0,50,1),tcl=0.2, tck=0.01,labels=F);
axis(1,at=seq(0,50,5),tcl=0.5, tck=0.01,labels=seq(-50,0,5));
#library(Hmisc)
#minor.tick(ny=n, tick.ratio=n)
legend('topleft', paste(c(as.character(mtf),'none'), ' (n=', n0, ', ',round(100*n0/N,1) , '%)',sep=''), text.col=cl, bty='n');
dev.off();" > /tmp/memepas.R

## only motif1 and its variants
# piRNA
#cat $annotationFile.pas.fa.PAS | awk '{OFS="\t"; if($5=="+") print}' | cut -f1,3 | sort -k1,1n -k2,2n | groupBy -g 1,2 -c 2 -o count | Rscript /tmp/memepas.R stdin; cp PAS.position.pdf ../result/PAS.piRNA.position.pdf
N=`wc -l $annotationFile.pas.intergenic.fa | awk '{print $1/2}'`
sort -k2,2 -k7,7g $annotationFile.pas.intergenic.fa.PAS | awk 'BEGIN{id=""}{if($3<50 && $2!=id) {print;id=$2;}}' | sort -k7,7g -k3,3 | groupBy -g 8,3 -c 3 -o count | sed 's/T/U/g' | Rscript /tmp/memepas.R stdin $N; cp PAS.position.pdf ../result/PAS.piRNA.intergenic.position.pdf
N=`wc -l $annotationFile.pas.genic.fa | awk '{print $1/2}'`
sort -k2,2 -k7,7g $annotationFile.pas.genic.fa.PAS | awk 'BEGIN{id=""}{if($3<50 && $2!=id) {print;id=$2;}}' | sort -k7,7g -k3,3 | groupBy -g 8,3 -c 3 -o count | sed 's/T/U/g' | Rscript /tmp/memepas.R stdin $N; cp PAS.position.pdf ../result/PAS.piRNA.genic.position.pdf

# NM/NR with top 1000 PAS peaks
N=`wc -l ../data/mouse_adult_wt_PAS200_100nt.unique.top1000.50PAS50.inNM.fa | awk '{print $1/2}'`
sort -k2,2 -k7,7g ../data/mouse_adult_wt_PAS200_100nt.unique.top1000.50PAS50.inNM.fa.PAS | awk 'BEGIN{id=""}{if($3<50 && $2!=id) {print;id=$2;}}' | sort -k7,7g -k3,3 | groupBy -g 8,3 -c 3 -o count | sed 's/T/U/g' | Rscript /tmp/memepas.R stdin $N; cp PAS.position.pdf ../result/PAS.NM.position.pdf
N=`wc -l ../data/mouse_adult_wt_PAS200_100nt.unique.top1000.50PAS50.inNR.fa | awk '{print $1/2}'`
sort -k2,2 -k7,7g ../data/mouse_adult_wt_PAS200_100nt.unique.top1000.50PAS50.inNR.fa.PAS | awk 'BEGIN{id=""}{if($3<50 && $2!=id) {print;id=$2;}}' | sort -k7,7g -k3,3 | groupBy -g 8,3 -c 3 -o count | sed 's/T/U/g' | Rscript /tmp/memepas.R stdin $N; cp PAS.position.pdf ../result/PAS.NR.position.pdf


#### with GU box
#echo "
#args <- commandArgs(TRUE);
#pdf('PAS.position.pdf',width=8, height=5, colormodel='cmyk');
#df=read.table(args[1]);
#motifs=unique(df\$V1);
#cls=colors()[c(552,81,78, 40,90,117,259,10,86)];
#plot(df[df\$V1==motifs[1],2],df[df\$V1==motifs[1],3],type='h',xlim=c(1,100),xaxt='n',las=1, tck=0.01, lwd=2,col=cls[1],xlab='Position relative to 3\' end of RNA', ylab='Motif occurrence');
#n=sum(df[df\$V1==motifs[1],3]);
#for(i in 3:length(motifs)){
#    # different axis for the last one (GT content)
#    if(i==length(motifs)) par(usr=c(par('usr')[1:2], range(df[df\$V1==motifs[i],3], 10)))
#    lines(df[df\$V1==motifs[i],2]+0.3*(i-1),df[df\$V1==motifs[i],3], type='h',col=cls[i],lwd=2);
#    n=c(n,sum(df[df\$V1==motifs[i],3]))
#    if(i==length(motifs)) axis(4, col.axis=cls[i])
#}
#axis(1,at=seq(0,100,10),col='gray',labels=seq(-50,50,10));
#legend('topleft', paste(c('AAUAAA', 'GU'), '(n=',n,')',sep=""), col=cls[c(1,3)], bty='n', lty=1, lwd=2);
#
#dev.off();" > /tmp/memepas.R
## piRNA
#cat $annotationFile.pas.fa.PAS <(perl -slne '$h=$_; $h=~s/^>//; chomp($h); $s=<>; while ($s=~/GT/gi) {printf("3\t%s\t%d\t%d\t\+\t1\t0\tGT\n", $h, pos($s)-3, pos($s));}' $annotationFile.pas.fa) | awk '{OFS="\t"; if($5=="+") print}' | cut -f1,3 | sort -k1,1n -k2,2n | groupBy -g 1,2 -c 2 -o count | Rscript /tmp/memepas.R stdin; cp PAS.position.pdf ../result/PAS.piRNA.position.pdf
#cat $annotationFile.pas.intergenic.fa.PAS <(perl -slne '$h=$_; $h=~s/^>//; chomp($h); $s=<>; while ($s=~/GT/gi) {printf("3\t%s\t%d\t%d\t\+\t1\t0\tGT\n", $h, pos($s)-3, pos($s));}' $annotationFile.pas.intergenic.fa) | awk '{OFS="\t"; if($5=="+") print}' | cut -f1,3 | sort -k1,1n -k2,2n | groupBy -g 1,2 -c 2 -o count | Rscript /tmp/memepas.R stdin; cp PAS.position.pdf ../result/PAS.piRNA.intergenic.position.pdf
#cat $annotationFile.pas.genic.fa.PAS <(perl -slne '$h=$_; $h=~s/^>//; chomp($h); $s=<>; while ($s=~/GT/gi) {printf("3\t%s\t%d\t%d\t\+\t1\t0\tGT\n", $h, pos($s)-3, pos($s));}' $annotationFile.pas.genic.fa) | awk '{OFS="\t"; if($5=="+") print}' | cut -f1,3 | sort -k1,1n -k2,2n | groupBy -g 1,2 -c 2 -o count | Rscript /tmp/memepas.R stdin; cp PAS.position.pdf ../result/PAS.piRNA.genic.position.pdf
## NM/NR with top 1000 PAS peaks
#cat ../data/mouse_adult_wt_PAS200_100nt.unique.top1000.50PAS50.inNM.fa.PAS <(perl -slne '$h=$_; $h=~s/^>//; chomp($h); $s=<>; while ($s=~/GT/gi) {printf("3\t%s\t%d\t%d\t\+\t1\t0\tGT\n", $h, pos($s)-3, pos($s));}' ../data/mouse_adult_wt_PAS200_100nt.unique.top1000.50PAS50.inNM.fa) | awk '{OFS="\t"; if($5=="+") print}' | cut -f1,3 | sort -k1,1n -k2,2n | groupBy -g 1,2 -c 2 -o count | Rscript /tmp/memepas.R stdin; cp PAS.position.pdf ../result/PAS.NM.top1000.position.pdf
#cat ../data/mouse_adult_wt_PAS200_100nt.unique.top1000.50PAS50.inNR.fa.PAS <(perl -slne '$h=$_; $h=~s/^>//; chomp($h); $s=<>; while ($s=~/GT/gi) {printf("3\t%s\t%d\t%d\t\+\t1\t0\tGT\n", $h, pos($s)-1, pos($s));}' ../data/mouse_adult_wt_PAS200_100nt.unique.top1000.50PAS50.inNR.fa) | awk '{OFS="\t"; if($5=="+") print}' | cut -f1,3 | sort -k1,1n -k2,2n | groupBy -g 1,2 -c 2 -o count | Rscript /tmp/memepas.R stdin; cp PAS.position.pdf ../result/PAS.NR.top1000.position.pdf

# copy to Dropbox
#scp hpcc:/home/dongx/projects/piRNA/result/PAS.*.pdf


# join
join -a 1 -1 1 -2 1 -e "NA" -o '0,1.2,2.2' <(sort -k1,1 $annotationFile.cpg) <(sort -k1,1 $annotationFile.TATA) | sed 's/ /\t/g' | sort -k1,1 | join -a 1 -1 1 -2 1 -e "NA" -o '0,1.2,1.3,2.2,2.5' - <(sort -k1,1 $annotationFile.BIP) | sed 's/ /\t/g;s/pi_//g' | sort -k1,1n | paste $annotationFile - | cut -f7-18,20- > $annotationFile.tab


##################################################
## relative position of EEJ on meta-mRNA
##################################################
echo "args <- commandArgs(TRUE);pdf('EEJ.position.pdf',width=8, height=5); df=read.table(args[1]); plot(df\$V1,df\$V2,type='h',xlim=c(0,100),xaxt='n',lwd=2,xlab='EEJ position relative to the full-length of mRNA (%)', ylab='EEJ count');axis(1,at=seq(0,100,10),col='gray',labels=seq(0,100,10));dev.off()" > /tmp/EEJ.R
# for piRNA
cut -f7- $annotationFile | awk '{split($11,a,","); printf("%s\t", $0); if($10>1) {sumLen=0; for(i=1;i<=$10;i++) {sumLen+=a[i];sL[i]=sumLen;} for(i=1;i<$10;i++) printf("%d,",($6=="+")?(100*sL[i]/sumLen+0.5):(100-100*sL[i]/sumLen+0.5))} print "";}' | cut -f13 | sed 's/,/\n/g;' | sed '/^$/d' | sort -k1,1n | groupBy -g 1 -c 1 -o count | Rscript /tmp/EEJ.R stdin; scp EEJ.position.pdf zlab:~/public_html/temp/EEJ.position.piRNA.pdf
cut -f7- $annotationFile | awk '{split($11,a,","); printf("%s\t", $0); if($10>1 && $4~/^pi/) {sumLen=0; for(i=1;i<=$10;i++) {sumLen+=a[i];sL[i]=sumLen;} for(i=1;i<$10;i++) printf("%d,",($6=="+")?(100*sL[i]/sumLen+0.5):(100-100*sL[i]/sumLen+0.5))} print "";}' | cut -f13 | sed 's/,/\n/g;' | sed '/^$/d' | sort -k1,1n | groupBy -g 1 -c 1 -o count | Rscript /tmp/EEJ.R stdin; scp EEJ.position.pdf zlab:~/public_html/temp/EEJ.position.piRNA.genic.pdf
cut -f7- $annotationFile | awk '{split($11,a,","); printf("%s\t", $0); if($10>1 && $4!~/^pi/) {sumLen=0; for(i=1;i<=$10;i++) {sumLen+=a[i];sL[i]=sumLen;} for(i=1;i<$10;i++) printf("%d,",($6=="+")?(100*sL[i]/sumLen+0.5):(100-100*sL[i]/sumLen+0.5))} print "";}' | cut -f13 | sed 's/,/\n/g;' | sed '/^$/d' | sort -k1,1n | groupBy -g 1 -c 1 -o count | Rscript /tmp/EEJ.R stdin; scp EEJ.position.pdf zlab:~/public_html/temp/EEJ.position.piRNA.intergenic.pdf
# for NM/NR
cat ../data/controlSet2/x100910.rnaseq.transcripts.NM.all.final.bed | awk '{split($11,a,","); printf("%s\t", $0); if($10>1) {sumLen=0; for(i=1;i<=$10;i++) {sumLen+=a[i];sL[i]=sumLen;} for(i=1;i<$10;i++) printf("%d,",($6=="+")?(100*sL[i]/sumLen+0.5):(100-100*sL[i]/sumLen+0.5))} print "";}' | cut -f13 | sed 's/,/\n/g;' | sed '/^$/d' | sort -k1,1n | groupBy -g 1 -c 1 -o count | Rscript /tmp/EEJ.R stdin; scp EEJ.position.pdf zlab:~/public_html/temp/EEJ.position.NM.pdf
cat ../data/controlSet2/x100910.rnaseq.transcripts.NM.all.final.bed | awk '{split($11,a,","); printf("%s\t", $0); if($10==2) {sumLen=0; for(i=1;i<=$10;i++) {sumLen+=a[i];sL[i]=sumLen;} for(i=1;i<$10;i++) printf("%d,",($6=="+")?(100*sL[i]/sumLen+0.5):(100-100*sL[i]/sumLen+0.5))} print "";}' | cut -f13 | sed 's/,/\n/g;' | sed '/^$/d' | sort -k1,1n | groupBy -g 1 -c 1 -o count | Rscript /tmp/EEJ.R stdin; scp EEJ.position.pdf zlab:~/public_html/temp/EEJ.position.NM.1intron.pdf
cat ../data/controlSet2/x100910.rnaseq.transcripts.NR.final.bed | awk '{split($11,a,","); printf("%s\t", $0); if($10>1) {sumLen=0; for(i=1;i<=$10;i++) {sumLen+=a[i];sL[i]=sumLen;} for(i=1;i<$10;i++) printf("%d,",($6=="+")?(100*sL[i]/sumLen+0.5):(100-100*sL[i]/sumLen+0.5))} print "";}' | cut -f13 | sed 's/,/\n/g;' | sed '/^$/d' | sort -k1,1n | groupBy -g 1 -c 1 -o count | Rscript /tmp/EEJ.R stdin; scp EEJ.position.pdf zlab:~/public_html/temp/EEJ.position.NR.pdf

cut -f7- $annotationFile | awk '{split($11,a,","); printf("%s\t", $0); if($10>1) {sumLen=0; for(i=1;i<=$10;i++) {sumLen+=a[i];sL[i]=sumLen;} for(i=1;i<$10;i++) printf("%d,",($6=="+")?(100*sL[i]/sumLen+0.5):(100-100*sL[i]/sumLen+0.5))} print "";}' | cut -f13 | sed 's/,/\n/g;' | sed '/^$/d' | awk '{OFS="\t"; print "pi.all", $1;}' > ../data/EEJ.pct.cvs
cut -f7- $annotationFile | awk '{split($11,a,","); printf("%s\t", $0); if($10>1 && $4~/^pi/) {sumLen=0; for(i=1;i<=$10;i++) {sumLen+=a[i];sL[i]=sumLen;} for(i=1;i<$10;i++) printf("%d,",($6=="+")?(100*sL[i]/sumLen+0.5):(100-100*sL[i]/sumLen+0.5))} print "";}' | cut -f13 | sed 's/,/\n/g;' | sed '/^$/d'  | awk '{OFS="\t"; print "pi.genic", $1;}' >> ../data/EEJ.pct.cvs
cut -f7- $annotationFile | awk '{split($11,a,","); printf("%s\t", $0); if($10>1 && $4!~/^pi/) {sumLen=0; for(i=1;i<=$10;i++) {sumLen+=a[i];sL[i]=sumLen;} for(i=1;i<$10;i++) printf("%d,",($6=="+")?(100*sL[i]/sumLen+0.5):(100-100*sL[i]/sumLen+0.5))} print "";}' | cut -f13 | sed 's/,/\n/g;' | sed '/^$/d'  | awk '{OFS="\t"; print "pi.intergenic", $1;}' >> ../data/EEJ.pct.cvs
cat ../data/controlSet2/x100910.rnaseq.transcripts.NM.all.final.bed | awk '{split($11,a,","); printf("%s\t", $0); if($10>1) {sumLen=0; for(i=1;i<=$10;i++) {sumLen+=a[i];sL[i]=sumLen;} for(i=1;i<$10;i++) printf("%d,",($6=="+")?(100*sL[i]/sumLen+0.5):(100-100*sL[i]/sumLen+0.5))} print "";}' | cut -f13 | sed 's/,/\n/g;' | sed '/^$/d'  | awk '{OFS="\t"; print "NM", $1;}' >> ../data/EEJ.pct.cvs
cat ../data/controlSet2/x100910.rnaseq.transcripts.NR.final.bed | awk '{split($11,a,","); printf("%s\t", $0); if($10>1) {sumLen=0; for(i=1;i<=$10;i++) {sumLen+=a[i];sL[i]=sumLen;} for(i=1;i<$10;i++) printf("%d,",($6=="+")?(100*sL[i]/sumLen+0.5):(100-100*sL[i]/sumLen+0.5))} print "";}' | cut -f13 | sed 's/,/\n/g;' | sed '/^$/d'  | awk '{OFS="\t"; print "NR", $1;}' >> ../data/EEJ.pct.cvs

# test in R
df=read.table('../data/EEJ.pct.cvs');
wilcox.test(df[df$V1=='pi.all',2],df[df$V1=='NM',2])
wilcox.test(df[df$V1=='NR',2],df[df$V1=='NM',2])
wilcox.test(df[df$V1=='pi.genic',2],df[df$V1=='pi.intergenic',2])


##################################################
## ORF scanning
##################################################
echo "
args <- commandArgs(TRUE);
df=read.table(args[1], fill=T, header=F);
pdf('ORF.pdf',width=10, height=round(nrow(df)/7));
par(mar=c(4,9,2,9))
plot(0:100,ylim=c(0,nrow(df)+1),xlab='relative position to 5\' end of mature transcript(%)', ylab='', type='n', frame.plot=F, xaxs='r', yaxs='i', axes=F)
for(i in c(1:nrow(df))){
    len=df[i,2]
    EEJ=as.numeric(unlist(strsplit(as.character(df[i,3]),',')))
    ORF=as.numeric(unlist(strsplit(as.character(df[i,4]),':')))
    uORFs=unlist(strsplit(as.character(df[i,5]),';'))
    uATGs=as.numeric(unlist(strsplit(as.character(df[i,6]),',')))
    # draw meta gene
    lines(c(0,100),c(i,i), col='gray')
    # draw ORF
    if(ORF[1]!=-1) lines(c(100*ORF[1]/len,100*sum(ORF)/len),c(i,i), col='darkgray', lwd=6)
    # draw EEJ
    if(EEJ[1]!=-1) points(100*EEJ/len,rep(i,length(EEJ)), col='blue', pch='|')
    # draw uORF
    ji=rep(i,length(uORFs))+runif(length(uORFs),0,0.3)
    k=1
    for(uORF in uORFs){
        orf=as.numeric(unlist(strsplit(uORF,':')))
        if(orf[1]!=-1 && orf[2]>=30) lines(c(100*orf[1]/len,100*sum(orf)/len),c(ji[k],ji[k]), col='red', lwd=2)
        k=k+1
    }
    # draw uATG
    if(uATGs[1]!=-1) points(100*uATGs[uATGs<ORF[1]]/len, rep(i,sum(uATGs<ORF[1])), col='green', pch='.')
}
axis(1, at=c(0:10)*10, labels=c(0:10)*10)

df[df[,7]<0,7]=0
df[df[,8]<0,8]=0
df[df[,9]<0,9]=0

axis(2, at=c(1:nrow(df))[grepl('^pi', df[,1])], labels=df[grepl('^pi', df[,1]),1], las=2, cex.axis=.7, col.axis='#D8A93B')
axis(2, at=c(1:nrow(df))[!grepl('^pi', df[,1])], labels=df[!grepl('^pi', df[,1]),1], las=2, cex.axis=.7, col.axis='#643A82')
axis(4, at=c(1:nrow(df)), labels=paste(df[,2],' (',df[,7], ',', df[,8],',', df[,9],')',sep=''), las=2, tick = F, line=-2, cex.axis=.7)
dev.off()
" > /tmp/orf.R

mkdir ../data/ORF
# for piRNA
#cat $annotationFileORF > /tmp/pi.annotation.bed
cut -f7- $annotationFile > /tmp/pi.annotation.bed
ORFanalysis /tmp/pi.annotation.bed | awk '{OFS="\t"; split($4,a,":"); split($6,uAUG,","); uaugs=""; for(i=1;i<length(uATG);i++) if(uAUG[i]<a[1]) uaugs=uaugs""uAUG[i]","; if(uaugs=="") uaugs="-1,"; $6=uaugs; print $0,a[1],a[2],(a[1]==-1)?-1:($2-a[1]-a[2])}' > ../data/ORF/pi.EEJ.ORF.uORF.uATG.tab
cat ../data/ORF/pi.EEJ.ORF.uORF.uATG.tab | grep -P "^pi" | sort -k8,8nr > /tmp/orf.tab
cat ../data/ORF/pi.EEJ.ORF.uORF.uATG.tab | grep -v -P "^pi" | sort -k8,8nr >> /tmp/orf.tab
cat /tmp/orf.tab | Rscript /tmp/orf.R stdin; cp ORF.pdf ../data/ORF/ORF.piRNA.type_length_sorted.pdf
cat ../data/ORF/pi.EEJ.ORF.uORF.uATG.tab | sort -k8,8nr | Rscript /tmp/orf.R stdin; cp ORF.pdf ../data/ORF/ORF.piRNA.length_sorted.pdf
cat ../data/ORF/pi.EEJ.ORF.uORF.uATG.tab | Rscript /tmp/orf.R stdin; cp ORF.pdf ../data/ORF/ORF.piRNA.unsorted.pdf

# for NM
intersectBed -a ../data/controlSet2/x100910.rnaseq.transcripts.NM.all.final.bed.c -b /tmp/pi.annotation.bed -s -v | awk '{OFS="\t"; name=$4"__"$5;$5=1000;$4=name;$9=0; print}'> /tmp/nm.annotation.bed
ORFanalysis /tmp/nm.annotation.bed | awk '{OFS="\t"; split($4,a,":"); split($6,uAUG,","); uaugs=""; for(i=1;i<length(uATG);i++) if(uAUG[i]<a[1]) uaugs=uaugs""uAUG[i]","; if(uaugs=="") uaugs="-1,";$6=uaugs; print $0,a[1],a[2],(a[1]==-1)?-1:($2-a[1]-a[2])}'> ../data/ORF/NM.EEJ.ORF.uORF.uATG.tab
cat ../data/ORF/NM.EEJ.ORF.uORF.uATG.tab | sort -k8,8nr | Rscript /tmp/orf.R stdin; cp ORF.pdf ../data/ORF/ORF.NM.length_sorted.pdf

# for NR
intersectBed -a ../data/controlSet2/x100910.rnaseq.transcripts.NR.all.final.bed.c -b /tmp/pi.annotation.bed -s -v | awk '{OFS="\t"; name=$4"__"$5;$5=1000;$4=name;$9=0; print}' > /tmp/nr.annotation.bed
ORFanalysis /tmp/nr.annotation.bed | awk '{OFS="\t"; split($4,a,":"); split($6,uAUG,","); uaugs=""; for(i=1;i<length(uATG);i++) if(uAUG[i]<a[1]) uaugs=uaugs""uAUG[i]","; if(uaugs=="") uaugs="-1,";$6=uaugs; print $0,a[1],a[2],(a[1]==-1)?-1:($2-a[1]-a[2])}'> ../data/ORF/NR.EEJ.ORF.uORF.uATG.tab
cat ../data/ORF/NR.EEJ.ORF.uORF.uATG.tab | sort -k8,8nr | Rscript /tmp/orf.R stdin; cp ORF.pdf ../data/ORF/ORF.NR.length_sorted.pdf

# random intergenic region
awk '{OFS="\t"; if($4!~/^pi/) {$4="RIR_"$4; print}}' /tmp/pi.annotation.bed | bedtools shuffle -excl <(cut -f1-3 $GENOME/mm9/Annotation/Genes/refGene.bed | cat - $GENOME/mm9/Annotation/Genes/mm9.gap.bed) -g /home/dongx/nearline/genomes/mm9/Annotation/Genes/ChromInfo.txt| awk '{OFS="\t"; $7=$2;$8=$3;$10=1;$11=($3-$2)","; $12="0,"; print}'> /tmp/random.intergenic.annotation.bed
ORFanalysis /tmp/random.intergenic.annotation.bed | awk '{OFS="\t"; split($4,a,":"); split($6,uAUG,","); uaugs=""; for(i=1;i<length(uATG);i++) if(uAUG[i]<a[1]) uaugs=uaugs""uAUG[i]","; if(uaugs=="") uaugs="-1,";$6=uaugs; print $0,a[1],a[2],(a[1]==-1)?-1:($2-a[1]-a[2])}' > ../data/ORF/random.intergenic.EEJ.ORF.uORF.uATG.tab
cat ../data/ORF/random.intergenic.EEJ.ORF.uORF.uATG.tab | sort -k8,8nr | Rscript /tmp/orf.R stdin; cp ORF.pdf ../data/ORF/ORF.RIR.length_sorted.pdf

#novel transcript
grep "#" ../data/dataforBIP_fromXin2/x100912.filtered.bed | awk '{OFS="\t"; $7=$2;$8=$3;$9=0;$10=1;$11=($3-$2); $12="0,"; print}' | sed 's/___#N\/A//g' > /tmp/novel.NA.annotation.bed
ORFanalysis /tmp/novel.NA.annotation.bed | awk '{OFS="\t"; split($4,a,":"); split($6,uAUG,","); uaugs=""; for(i=1;i<length(uATG);i++) if(uAUG[i]<a[1]) uaugs=uaugs""uAUG[i]","; if(uaugs=="") uaugs="-1,";$6=uaugs; print $0,a[1],a[2],(a[1]==-1)?-1:($2-a[1]-a[2])}' > ../data/ORF/novel.NA.EEJ.ORF.uORF.uATG.tab
cat ../data/ORF/novel.NA.EEJ.ORF.uORF.uATG.tab | sort -k8,8nr | Rscript /tmp/orf.R stdin; cp ORF.pdf ../data/ORF/ORF.novel.NA.length_sorted.pdf


## germline genes
#intersectBed -a ../data/controlSet/germline.specifc.nonAMyb.bed -b /tmp/pi.annotation.bed -s -v | awk '{OFS="\t"; name=$4"__"$5;$5=1000;$4=name;$9=0; print}' > /tmp/germline.annotation.bed
#ORFanalysis /tmp/germline.annotation.bed | awk '{OFS="\t"; split($4,a,":"); split($6,uAUG,","); uaugs=""; for(i=1;i<length(uATG);i++) if(uAUG[i]<a[1]) uaugs=uaugs""uAUG[i]","; if(uaugs=="") uaugs="-1,";$6=uaugs; print $0,a[1],a[2],(a[1]==-1)?-1:($2-a[1]-a[2])}'> ../data/ORF/germline.EEJ.ORF.uORF.uATG.tab
#cat ../data/ORF/germline.EEJ.ORF.uORF.uATG.tab | sort -k8,8nr | Rscript /tmp/orf.R stdin; cp ORF.pdf ../data/ORF/ORF.germline.length_sorted.pdf
#
## housekeeping
#intersectBed -a ../data/controlSet/housekeeping.bed -b /tmp/pi.annotation.bed -s -v | awk '{OFS="\t"; name=$4"__"$5;$5=1000;$4=name;$9=0; print}' > /tmp/housekeeping.annotation.bed
#ORFanalysis /tmp/housekeeping.annotation.bed | awk '{OFS="\t"; split($4,a,":"); split($6,uAUG,","); uaugs=""; for(i=1;i<length(uATG);i++) if(uAUG[i]<a[1]) uaugs=uaugs""uAUG[i]","; if(uaugs=="") uaugs="-1,";$6=uaugs; print $0,a[1],a[2],(a[1]==-1)?-1:($2-a[1]-a[2])}'> ../data/ORF/housekeeping.EEJ.ORF.uORF.uATG.tab
#cat ../data/ORF/housekeeping.EEJ.ORF.uORF.uATG.tab | sort -k8,8nr | Rscript /tmp/orf.R stdin; cp ORF.pdf ../data/ORF/ORF.housekeeping.length_sorted.pdf
#
## housekeeping
#intersectBed -a ../data/controlSet/germline.silenced.bed -b /tmp/pi.annotation.bed -s -v | awk '{OFS="\t"; name=$4"__"$5;$5=1000;$4=name;$9=0; print}' > /tmp/germline.silenced.annotation.bed
#ORFanalysis /tmp/germline.silenced.annotation.bed | awk '{OFS="\t"; split($4,a,":"); split($6,uAUG,","); uaugs=""; for(i=1;i<length(uATG);i++) if(uAUG[i]<a[1]) uaugs=uaugs""uAUG[i]","; if(uaugs=="") uaugs="-1,";$6=uaugs; print $0,a[1],a[2],(a[1]==-1)?-1:($2-a[1]-a[2])}'> ../data/ORF/germline.silenced.EEJ.ORF.uORF.uATG.tab
#cat ../data/ORF/germline.silenced.EEJ.ORF.uORF.uATG.tab | sort -k8,8nr | Rscript /tmp/orf.R stdin; cp ORF.pdf ../data/ORF/ORF.germline.silenced.length_sorted.pdf

scp ../data/ORF/*.pdf zlab:~/public_html/temp

# compare annotated ORF with predicted ORF for NM
awk '{OFS="\t"; $4=$4"__"NR; print}' /home/dongx/nearline/genomes/mm9/Annotation/Genes/refGene.bed > ../data/ORF/refGene.bed
ORFanalysis ../data/ORF/refGene.bed | awk '{OFS="\t"; split($4,a,":"); split($6,uAUG,","); uaugs=""; for(i=1;i<length(uATG);i++) if(uAUG[i]<a[1]) uaugs=uaugs""uAUG[i]","; if(uaugs=="") uaugs="-1,";$6=uaugs; print $0,a[1],a[2],(a[1]==-1)?-1:($2-a[1]-a[2])}' > ../data/ORF/refGene.EEJ.ORF.uORF.uATG.tab &
awk '{OFS="\t"; split($11,l,","); split($12,s,","); cdsS=$7; cdsE=$8; sum=0; utr5=0; utr3=0; for(i=1;i<=$10;i++) {sum=sum+l[i]; if(($2+s[i]+l[i])<=cdsS) utr5+=l[i]; if(($2+s[i])<cdsS && ($2+s[i]+l[i])>cdsS) utr5+=(cdsS-$2-s[i]); if(($2+s[i])<cdsE && ($2+s[i]+l[i])>cdsE) utr3+=($2+s[i]+l[i]-cdsE); if(($2+s[i])>cdsE) utr3+=l[i];} print $4, ($6=="+")?utr5:utr3, sum-utr5-utr3, ($6=="+")?utr3:utr5;}' ../data/ORF/refGene.bed > ../data/ORF/refGene.5UTR.ORF.3UTR.tab
join -1 1 -2 1 -a 1 -a 2 <(sort -k1,1 ../data/ORF/refGene.EEJ.ORF.uORF.uATG.tab) <(sort -k1,1 ../data/ORF/refGene.5UTR.ORF.3UTR.tab) | sed 's/ /\t/g' > ../data/ORF/refGene.EEJ.ORF.uORF.uATG.tab2

echo "
args <- commandArgs(TRUE);
#df=read.table(args[1], header=F, col.names=c('id','mRNA_length','EEJ','ORF','uORF','uAUG','UTR5_length','ORF_length','UTR3_length','Refseq_5UTR_length','Refseq_ORF_length','Refseq_3UTR_length'));
df=read.table('../data/ORF/refGene.EEJ.ORF.uORF.uATG.tab2', header=F, col.names=c('id','mRNA_length','EEJ','ORF','uORF','uAUG','UTR5_length','ORF_length','UTR3_length','Refseq_5UTR_length','Refseq_ORF_length','Refseq_3UTR_length'));
row.names(df)=df[,1]
pdf('orf.compare.NM.pdf', width=12,height=4)
df=df[grep('^NM',df\$id),c(7:12)]
df=df[df\$Refseq_ORF_length %% 3 ==0, ]
df[df==-1]=0
par(mfrow=c(1,3))
plot(df\$UTR5_length, df\$Refseq_5UTR_length, col='#0000ff22', pch=16, asp=1, xlim=range(df\$UTR5_length, df\$RefSeq_5UTR_length), ylim=range(df\$UTR5_length, df\$RefSeq_5UTR_length))
legend('topright',paste(round(sum(df\$UTR5_length==df\$Refseq_5UTR_length)*100/length(df\$Refseq_5UTR_length),2), '% of ',length(df\$Refseq_5UTR_length),' are same',sep=''))
plot(df\$ORF_length, df\$Refseq_ORF_length, col='#0000ff22', pch=16, asp=1, xlim=range(df\$ORF_length, df\$Refseq_ORF_length), ylim=range(df\$ORF_length, df\$Refseq_ORF_length))
legend('topright',paste(round(sum(df\$ORF_length==df\$Refseq_ORF_length)*100/length(df\$Refseq_ORF_length),2), '% of ',length(df\$Refseq_ORF_length),' are same',sep=''))
plot(df\$UTR3_length, df\$Refseq_3UTR_length, col='#0000ff22', pch=16, asp=1, xlim=range(df\$UTR3_length, df\$Refseq_3UTR_length), ylim=range(df\$UTR3_length, df\$Refseq_3UTR_length))
legend('topright',paste(round(sum(df\$UTR3_length==df\$Refseq_3UTR_length)*100/length(df\$Refseq_3UTR_length),2), '% of ', length(df\$Refseq_3UTR_length),' are same',sep=''))
dev.off()
" > /tmp/orfcompare.R
Rscript /tmp/orfcompare.R ../data/ORF/refGene.EEJ.ORF.uORF.uATG.tab2; scp orf.compare.NM.pdf zlab:~/public_html/temp

# statistics  (TODO: check uAUG count, make sure not including all AUG)
> ../data/ORF/merge.tab;
for i in ../data/ORF/NM.EEJ.ORF.uORF.uATG.tab ../data/ORF/NR.EEJ.ORF.uORF.uATG.tab ../data/ORF/pi.EEJ.ORF.uORF.uATG.tab ../data/ORF/random.intergenic.EEJ.ORF.uORF.uATG.tab ../data/ORF/novel.NA.EEJ.ORF.uORF.uATG.tab; do j=`basename $i`; k=${j/.EEJ.ORF.uORF.uATG.tab/}; echo $k; awk -v k=$k '{OFS="\t"; j=(k=="pi")?(($1~/^pi/)?"pi.genic":"pi.intergenic"):k; print $0,j;}' $i >> ../data/ORF/merge.tab; done
awk '{OFS="\t"; split($3,a,","); split($5,b,";"); split($6,c,","); l_sum=0; n=0; for(i=1;i<=length(b);i++) {split(b[i],d,":");if(d[2]>0) {n++; l_sum=l_sum+d[2];}} print $0, length(a)-1,n, (n==0)?0:l_sum/n, length(c)-1;}' ../data/ORF/merge.tab > ../data/ORF/merge2.tab
echo "
args <- commandArgs(TRUE);
df=read.table(args[1], fill=T, header=F, col.names=c('id','mRNA_length','EEJ','ORF','uORF','uAUG','UTR5_length','ORF_length','UTR3_length','Category','EEJ_count','uORF_count','uORF_length','uAUG_count'));
#df=read.table('../data/ORF/merge2.tab', fill=T, header=F, col.names=c('id','mRNA_length','EEJ','ORF','uORF','uAUG','UTR5_length','ORF_length','UTR3_length','Category','EEJ_count','uORF_count','uORF_length','uAUG_count'));
df\$Category=as.character(df\$Category)
df=df[df\$Category %in% c('NM','NR','pi.intergenic','pi.genic','random.intergenic','novel.NA'), ]
df\$Category=factor(df\$Category, levels=c('pi.genic','NM','pi.intergenic','NR','novel.NA','random.intergenic'))
pdf('ORFboxplot.pdf',width=12, height=9);
par(mfcol=c(3,2), mar=c(4,3,2,1))
# orf length comparison
boxplot(UTR5_length ~ Category, df, xaxt='n', outline=F, main='5\' UTR length')
boxplot(ORF_length ~ Category, df, xaxt='n', outline=F, main='ORF length')
boxplot(UTR3_length ~ Category, df, xaxt='n', outline=F, main='3\' UTR length')
axis(1, at=c(1:length(levels(df\$Category))), sapply(levels(df\$Category), function(x) paste(x, paste('(n=',sum(df\$Category==x),')',sep=""),sep='\n')), padj=0.5)
boxplot(uAUG_count ~ Category, df, xaxt='n', outline=F, main='uAUG count')
boxplot(uORF_count ~ Category, df, xaxt='n', outline=F, main='uORF count')
boxplot(uORF_length ~ Category, df, xaxt='n', outline=F, main='uORF length')
axis(1, at=c(1:length(levels(df\$Category))), sapply(levels(df\$Category), function(x) paste(x, paste('(n=',sum(df\$Category==x),')',sep=""),sep='\n')), padj=0.5)
dev.off();
"> /tmp/orfboxplot.R
Rscript /tmp/orfboxplot.R ../data/ORF/merge2.tab; scp ORFboxplot.pdf zlab:~/public_html/temp

##################################################
### dinucleotide frequency
##################################################

awk '{OFS="\t"; if($4!~/^pi/) {print}}' /tmp/pi.annotation.bed  | bedtools getfasta -s -tab -name -fi /home/dongx/nearline/genomes/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo stdout | awk '{s=toupper($2); if(s~/N/) next;  len=length(s); split("", di); for(i=1;i<length(s);i++) di[substr(s,i,2)]++; i=1; for(j in di) print $i, j, int(0.5+100*di[j]/(len-1)), "piRNA-intergenic";}' > /tmp/di.tab
awk '{OFS="\t"; if($4~/^pi/) {print}}' /tmp/pi.annotation.bed  | bedtools getfasta -s -tab -name -fi /home/dongx/nearline/genomes/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo stdout | awk '{s=toupper($2); if(s~/N/) next;  len=length(s); split("", di); for(i=1;i<length(s);i++) di[substr(s,i,2)]++; i=1; for(j in di) print $i, j, int(0.5+100*di[j]/(len-1)), "piRNA-genic";}' >> /tmp/di.tab
cat /tmp/nm.annotation.bed  | bedtools getfasta -s -tab -name -fi /home/dongx/nearline/genomes/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo stdout | awk '{s=toupper($2); if(s~/N/) next; len=length(s); split("", di); for(i=1;i<length(s);i++) di[substr(s,i,2)]++; i=1; for(j in di) print $i, j, int(0.5+100*di[j]/(len-1)), "mRNA";}' >> /tmp/di.tab
cat /tmp/nr.annotation.bed  | bedtools getfasta -s -tab -name -fi /home/dongx/nearline/genomes/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo stdout | awk '{s=toupper($2); if(s~/N/) next;  len=length(s); split("", di); for(i=1;i<length(s);i++) di[substr(s,i,2)]++; i=1; for(j in di) print $i, j, int(0.5+100*di[j]/(len-1)), "ncRNA";}' >> /tmp/di.tab
cat /tmp/random.intergenic.annotation.bed | bedtools getfasta -s -tab -name -fi /home/dongx/nearline/genomes/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo stdout | awk '{s=toupper($2); if(s~/N/) next;  len=length(s); split("", di); for(i=1;i<length(s);i++) di[substr(s,i,2)]++; i=1; for(j in di) print $1, j, int(0.5+100*di[j]/(len-1)), "ramdom.intergenic.regions";}' >> /tmp/di.tab
awk '{OFS="\t"; if($4!~/^pi/) {s=($6=="+")?($2-1000):($3-500); print $1,s,s+1500,$4,$5,$6}}' /tmp/pi.annotation.bed  | bedtools getfasta -s -tab -name -fi /home/dongx/nearline/genomes/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo stdout | awk '{s=toupper($2); if(s~/N/) next;  len=length(s); split("", di); for(i=1;i<length(s);i++) di[substr(s,i,2)]++; i=1; for(j in di) print $i, j, int(0.5+100*di[j]/(len-1)), "piRNA-intergenic.promoter";}' >> /tmp/di.tab
awk '{OFS="\t"; if($4~/^pi/) {s=($6=="+")?($2-1000):($3-500); print $1,s,s+1500,$4,$5,$6}}' /tmp/pi.annotation.bed  | bedtools getfasta -s -tab -name -fi /home/dongx/nearline/genomes/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo stdout | awk '{s=toupper($2); if(s~/N/) next;  len=length(s); split("", di); for(i=1;i<length(s);i++) di[substr(s,i,2)]++; i=1; for(j in di) print $i, j, int(0.5+100*di[j]/(len-1)), "piRNA-genic.promoter";}' >> /tmp/di.tab


#df1=read.table('/tmp/pi.annotation.bed.di.tab')
#colnames(df1)=c('id', 'di','percentage')
## AA', 'CA', 'AC', 'CC', 'TA', 'AT', 'AG', 'GA', 'CG', 'GC', 'TC', 'CT', 'TG', 'TT', 'GT', 'GG
#df2=read.table('/tmp/random.intergenic.annotation.bed.di.tab')
#colnames(df2)=c('id', 'di','percentage')
#df3=read.table('/tmp/pi.intergenic.promoter.bed.di.tab')
#colnames(df3)=c('id', 'di','percentage')
#pdf('/tmp/di.pdf', colormodel='cmyk', width=10,height=12)
#par(mfrow=c(3,1))
#boxplot(percentage ~ di, df1, ylim=range(df1\$percentage, df2\$percentage, df3\$percentage), ylab='percentage', main='intergenic piRNA loci')
#boxplot(percentage ~ di, df2, ylim=range(df1\$percentage, df2\$percentage, df3\$percentage), ylab='percentage', main='random intergenic regions')
#boxplot(percentage ~ di, df3, ylim=range(df1\$percentage, df2\$percentage, df3\$percentage), ylab='percentage', main='promoters of intergenic piRNA loci')
#dev.off()
echo "
df1=read.table('/tmp/di.tab')
colnames(df1)=c('id', 'di','percentage','type')
pdf('/tmp/di.pdf', colormodel='cmyk', width=20,height=5)
library(ggplot2)
df1\$type <- factor(df1\$type, levels = c('ramdom.intergenic.regions', 'piRNA-intergenic','ncRNA', 'piRNA-genic', 'mRNA', 'piRNA-intergenic.promoter','piRNA-genic.promoter'))
ggplot(df1, aes(factor(di), percentage)) + geom_boxplot(aes(fill = factor(type)), outlier.colour = 'darkgray', outlier.size = 1)
dev.off()
" > /tmp/di.R
Rscript /tmp/di.R; scp /tmp/di.pdf zlab:~/public_html/temp

##################################################
### BIP
##################################################
### define BIP using dataset from Xin: /home/lix/nearline/new.assembly/alternative/promoter/for.xiangjun/final.list
# pi:pi, pi:union, pi:NM, pi:NR, pi:none
# union:union, union:NM, union:NR, union:none,
# NM:NM, NM:NR, NM:none
# NR:NR, NR:none

pi=../data/dataforBIP_fromXin2/piRNA.true.promoter.bed
NM=../data/dataforBIP_fromXin2/NM.promoter.2.bed
NR=../data/dataforBIP_fromXin2/NR.promoter.2.bed
union=../data/dataforBIP_fromXin2/x100912.filtered.bed
none=../data/dataforBIP_fromXin2/noneunion.cage.narrowpeak
intersectBed -a <(awk '{OFS="\t"; if($5>5 || $7>10) {left=($6=="+")?($2+$10-500):($2+$10-50); if(left<0) left=0; print $1,left,left+550,$4,$5,$6;}}' $CAGE.narrowpeak) -b <(awk '{OFS="\t"; left=($6=="+")?($2-500):($3-50);if(left<0) left=0; print $1,left,left+550,$4,$5,$6; }' $union) -s -v > ../data/dataforBIP_fromXin2/noneunion.cage.narrowpeak

getBIP $pi $pi | awk '{OFS="\t"; print $0, "pi:pi";}' > ../data/dataforBIP_fromXin2/BIP.pi.pi.bed &
getBIP $pi $NM | awk '{OFS="\t"; print $0, "pi:NM";}'  > ../data/dataforBIP_fromXin2/BIP.pi.NM.bed &
getBIP $pi $NR | awk '{OFS="\t"; print $0, "pi:NR";}'  > ../data/dataforBIP_fromXin2/BIP.pi.NR.bed &
getBIP $pi $union | awk '{OFS="\t"; print $0, "pi:union";}'  > ../data/dataforBIP_fromXin2/BIP.pi.union.bed &
getBIP $pi $none | awk '{OFS="\t"; print $0, "pi:none";}'  > ../data/dataforBIP_fromXin2/BIP.pi.none.bed &

getBIP $union $union | awk '{OFS="\t"; print $0, "union:union";}'  > ../data/dataforBIP_fromXin2/BIP.union.union.bed &
getBIP $NM $union | awk '{OFS="\t"; print $0, "NM:union";}'  > ../data/dataforBIP_fromXin2/BIP.NM.union.bed &
getBIP $NR $union | awk '{OFS="\t"; print $0, "NR:union";}'  > ../data/dataforBIP_fromXin2/BIP.NR.union.bed &
getBIP $none $union | awk '{OFS="\t"; print $0, "none:union";}'  > ../data/dataforBIP_fromXin2/BIP.none.union.bed &

#getBIP $NM $NM | awk '{OFS="\t"; print $0, "NM:NM";}'  > ../data/dataforBIP_fromXin2/BIP.NM.NM.bed &
#getBIP $NM $NR | awk '{OFS="\t"; print $0, "NM:NR";}'  > ../data/dataforBIP_fromXin2/BIP.NM.NR.bed &
#getBIP $NM $none | awk '{OFS="\t"; print $0, "NM:none";}'  > ../data/dataforBIP_fromXin2/BIP.NM.none.bed &
#
#getBIP $NR $NR | awk '{OFS="\t"; print $0, "NR:NR";}'  > ../data/dataforBIP_fromXin2/BIP.NR.NR.bed &
#getBIP $NR $none | awk '{OFS="\t"; print $0, "NR:none";}'  > ../data/dataforBIP_fromXin2/BIP.NR.none.bed &

for i in ../data/dataforBIP_fromXin2/BIP.*.bed;
do
    sort -k2,2 $i | paste - <(getCpGscore <(awk '{OFS="\t";s=($3<$6)?$3:$6; e=($3<$6)?$6:$3; print $1,s,e,$2,0,"+"}' $i) | sort -k1,1 | cut -f2) | join -1 2 -2 1 -a 1 -a 2 -e'-' -o '1.1,0,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,2.2,2.3' - <(scanTATA <(awk '{OFS="\t";s=($4=="+")?($3-50):$3; print $1,s,s+50,$2,0,$4}' $i)) | sed 's/ /\t/g' | awk '{if($0~/^chr/) print}' | join -1 2 -2 1 -a 1 -a 2 -e'-' -o '1.1,0,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,2.2,2.3' - <(scanTATA <(awk '{OFS="\t"; s=($7=="+")?($6-50):$6; print $1,s,s+50,$2,0,$7}' $i)) | sed 's/ /\t/g' | awk '{if($0~/^chr/) print}'  | paste - <(scanCpGisland <(awk '{OFS="\t";s=($3<$6)?$3:$6; e=($3<$6)?$6:$3; print $1,s,e,$2,0,"+"}' $i) | sort -k1,1 | cut -f2) > $i.TATA.CpG.CGI.tab &
done


##################################################
### Conservation
##################################################
# compare conservation between NM/NR/pi-inger/pi-genic
for i in $GENOME/mm9/Conservation/phyloP30way/placentalMammals/*.wigFix.gz; do zcat $i | wigToBigWig stdin /home/dongx/nearline/genomes/mm9/Annotation/Genes/ChromInfo.txt $i.bw & done
# convert phyloP bigwig file
cat $annotationFileORF > /tmp/pi.annotation.bed
intersectBed -a ../data/ORF/refGene.bed -b /tmp/nm.annotation.bed -wa -u -s | grep NM_ > /tmp/refseq.annotation.bed
for j in pi refseq nr random.intergenic novel.NA;
do
    > /tmp/$j.annotation.bed.exons.phyloP
    > /tmp/$j.annotation.bed.introns.phyloP
    > /tmp/$j.annotation.bed.promoters.phyloP
    > /tmp/$j.annotation.bed.5UTR.phyloP
    > /tmp/$j.annotation.bed.CDS.phyloP
    > /tmp/$j.annotation.bed.3UTR.phyloP
    for i in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrM chrX chrY;
    do
        bw=$GENOME/mm9/Conservation/phyloP30way/placentalMammals/$i.phyloP30way.placental.wigFix.gz.bw
        echo $j, $i, 'exon';
        grep -w $i /tmp/$j.annotation.bed | bigWigAverageOverBed $bw stdin stdout | cut -f1,6 >> /tmp/$j.annotation.bed.exons.phyloP
        echo $j, $i, 'intron';
        grep -w $i /tmp/$j.annotation.bed | awk '{OFS="\t";split($11,a,","); split($12,b,","); A=""; B=""; for(i=1;i<length(a)-1;i++) {A=A""(b[i+1]-b[i]-a[i])",";B=B""(b[i]+a[i]-(b[1]+a[1]))",";} if($10>1) print $1,$2+a[1], $3-a[length(a)-1], $4,$5,$6,$2+a[1], $3-a[length(a)-1],$9,$10-1,A,B;}' | bigWigAverageOverBed $bw stdin stdout | cut -f1,6 >> /tmp/$j.annotation.bed.introns.phyloP
        echo $j, $i, '5UTR';
        grep -w $i /tmp/$j.annotation.bed | awk '{OFS="\t";split($11,a,","); split($12,b,","); A=""; B="";  if($7==$8) next; if($6=="+" && $2<$7) {for(i=1;i<length(a);i++) if(($2+b[i]+a[i])<=$7) {A=A""a[i]",";B=B""b[i]",";} else {A=A""($7-$2-b[i])",";B=B""b[i]","; break; } print $1,$2,$7,$4,$5,$6,$2,$7,$9,i,A,B;} if($6=="-" && $8<$3) {for(i=length(a)-1;i>0;i--) if(($2+b[i])>=$8) {A=a[i]","A;B=($2+b[i]-$8)","B;} else {A=($2+b[i]+a[i]-$8)","A;B=0","B; break; } print $1,$8,$3,$4,$5,$6,$8,$3,$9,length(a)-i,A,B;}}' | bigWigAverageOverBed $bw stdin stdout | cut -f1,6 >> /tmp/$j.annotation.bed.5UTR.phyloP
        echo $j, $i, 'CDS';
        grep -w $i /tmp/$j.annotation.bed | awk '{OFS="\t";split($11,a,","); split($12,b,","); A=""; B=""; if($7==$8) next; j=0; for(i=1;i<length(a);i++) if(($2+b[i]+a[i])>$7 && ($2+b[i])<$8) {j++; start=$2+b[i]-$7; size=a[i]; if(($2+b[i])<=$7) {start=0;size=size-($7-($2+b[i]));} if(($2+a[i]+b[i])>=$8) {size=size-($2+a[i]+b[i]-$8);} A=A""size",";B=B""start",";} print $1,$7,$8,$4,$5,$6,$7,$8,$9,j,A,B;}' | bigWigAverageOverBed $bw stdin stdout | cut -f1,6 >> /tmp/$j.annotation.bed.CDS.phyloP
        echo $j, $i, '3UTR';
        grep -w $i /tmp/$j.annotation.bed | awk '{OFS="\t";split($11,a,","); split($12,b,","); A=""; B=""; if($7==$8) next; if($6=="-" && $8<$3) {for(i=1;i<length(a);i++) if(($2+b[i]+a[i])<=$7) {A=A""a[i]",";B=B""b[i]",";} else {A=A""($7-$2-b[i])",";B=B""b[i]","; break; } print $1,$2,$7,$4,$5,$6,$2,$7,$9,i,A,B;} if($6=="+" && $2<$7) {for(i=length(a)-1;i>0;i--) if(($2+b[i])>=$8) {A=a[i]","A;B=($2+b[i]-$8)","B;} else {A=($2+b[i]+a[i]-$8)","A;B=0","B; break; } print $1,$8,$3,$4,$5,$6,$8,$3,$9,length(a)-i,A,B;}}' | bigWigAverageOverBed $bw stdin stdout | cut -f1,6 >> /tmp/$j.annotation.bed.3UTR.phyloP
        echo $j, $i, 'promoter';
        grep -w $i /tmp/$j.annotation.bed | awk '{OFS="\t";s=($6=="+")?($2-500):$3; print $1,s,s+500,$4,$5,$6}' | bigWigAverageOverBed $bw stdin stdout | cut -f1,6 >> /tmp/$j.annotation.bed.promoters.phyloP
    done
done

# merge
>../data/conservation.exon.intron.utr.promoter.tab
for i in /tmp/*.annotation.bed.*.phyloP; do j=`basename $i | sed 's/.annotation.bed//;s/.phyloP//'`; awk -v j=$j '{OFS="\t"; split(j,a,"."); type=(a[1]=="pi")?(($1~/^pi/)?"piRNA-genic":"piRNA-intergenic"):a[1]; if(type=="refseq") type="mRNA"; if(type=="nr") type="ncRNA"; if(type=="random") type="random.intergenic.region"; if(type=="novel") type="novel.transcript"; print type,a[length(a)],$0}' $i >> ../data/conservation.exon.intron.utr.promoter.tab; done
# R script
echo "
df=read.table('../data/conservation.exon.intron.utr.promoter.tab', header=F)
colnames(df)=c('type', 'regions','id','phyloP')
pdf('/tmp/conservation.pdf', colormodel='cmyk', width=20,height=5)
library(ggplot2)
df\$type <- factor(df\$type, levels = c('random.intergenic.region', 'piRNA-intergenic','ncRNA', 'novel.transcript','piRNA-genic', 'mRNA'))
df\$regions <- factor(df\$regions, levels = c('exons','introns', 'promoters','5UTR','CDS', '3UTR'))
ggplot(df, aes(factor(regions), phyloP)) + geom_boxplot(aes(fill = factor(type)), outlier.colour = 'darkgray', outlier.size = 1) #+stat_summary(fun.y=mean, geom='point')
#ggplot(df, aes(factor(regions), phyloP)) + geom_ribbon(aes(fill = factor(type), ymax = ..density.., ymin = -..density..), stat = 'density')

dev.off()
" > /tmp/conservation.R
Rscript /tmp/conservation.R; scp /tmp/conservation.pdf zlab:~/public_html/temp



##################################################
### junction table (In Supplemental Table S1, tab (H))
##################################################
junctions=/home/dongx/projects/piRNA/data/for_junctiontable/mouse_adult_wt_RNAseq_PE50nt_strand.junctions2intron.all.final.merged.tab
NM=../data/controlSet2/x100910.rnaseq.transcripts.NM.all.final.bed.c
NR=../data/controlSet2/x100910.rnaseq.transcripts.NR.all.final.bed.c
pi_genic=../data/controlSet2/piRNA.genic.bed
pi_intergenic=../data/controlSet2/piRNA.intergenic.bed
juntions_only_in_4setoftranscript=../data/for_junctiontable/mouse_adult_wt_RNAseq_PE50nt_strand.junctions2intron.OnlyUsed.in.assembly.4types

# only junctions used in the above transctipt sets
awk '{OFS="\t"; $4=$1"_"$2"_"$3"_"$6; print}' $junctions > /tmp/junctions.bed  # no redundant
cat $pi_genic | awk '{OFS="\t"; split($11,a,","); split($12,b,","); for(i=1;i<$10;i++) print $1"_"($2+b[i]+a[i])"_"($2+b[i+1])"_"$6;}' | fgrep -f - /tmp/junctions.bed | cut -f4,7,9 | awk '{OFS="\t"; print $0, "piRNA-genic";}' > $juntions_only_in_4setoftranscript
cat $pi_intergenic | awk '{OFS="\t"; split($11,a,","); split($12,b,","); for(i=1;i<$10;i++) print $1"_"($2+b[i]+a[i])"_"($2+b[i+1])"_"$6;}' | fgrep -f - /tmp/junctions.bed | cut -f4,7,9 | awk '{OFS="\t"; print $0, "piRNA-intergenic";}' >> $juntions_only_in_4setoftranscript
cat $NM | awk '{OFS="\t"; split($11,a,","); split($12,b,","); for(i=1;i<$10;i++) print $1"_"($2+b[i]+a[i])"_"($2+b[i+1])"_"$6;}' | fgrep -f - /tmp/junctions.bed | cut -f4,7,9 | awk '{OFS="\t"; print $0, "mRNA";}' >> $juntions_only_in_4setoftranscript
cat $NR | awk '{OFS="\t"; split($11,a,","); split($12,b,","); for(i=1;i<$10;i++) print $1"_"($2+b[i]+a[i])"_"($2+b[i+1])"_"$6;}' | fgrep -f - /tmp/junctions.bed | cut -f4,7,9 | awk '{OFS="\t"; print $0, "ncRNA";}' >> $juntions_only_in_4setoftranscript

## all junctions
#intersectBed -a $junctions -b $NM -s -wo | cut -f1-7,9,11,15,24,25,29 | awk '{OFS="\t"; if($9>0) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11"___"$12"___"$13;}'|sort -k4,4 | groupBy -g 1,2,3,4,5,6,7,8,9,10 -c 11 -o collapse | awk '{OFS="\t"; print $0, "NM";}' > ../data/for_junctiontable/mouse_adult_wt_RNAseq_PE50nt_strand.junctions2intron.all.final.merged.tab.4types
#intersectBed -a $junctions -b $NR -s -wo | cut -f1-7,9,11,15,24,25,29 | awk '{OFS="\t"; if($9>0) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11"___"$12"___"$13;}'|sort -k4,4 | groupBy -g 1,2,3,4,5,6,7,8,9,10 -c 11 -o collapse | awk '{OFS="\t"; print $0, "NR";}' >> ../data/for_junctiontable/mouse_adult_wt_RNAseq_PE50nt_strand.junctions2intron.all.final.merged.tab.4types
#intersectBed -a $junctions -b $pi_genic -s -wo | cut -f1-7,9,11,15,24,25 | awk '{OFS="\t"; if($9>0) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11"___"$12;}'|sort -k4,4 | groupBy -g 1,2,3,4,5,6,7,8,9,10 -c 11 -o collapse | awk '{OFS="\t"; print $0, "pi_genic";}' >> ../data/for_junctiontable/mouse_adult_wt_RNAseq_PE50nt_strand.junctions2intron.all.final.merged.tab.4types
#intersectBed -a $junctions -b $pi_intergenic -s -wo | cut -f1-7,9,11,15,24,25 | awk '{OFS="\t"; if($9>0) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11"___"$12;}'|sort -k4,4 | groupBy -g 1,2,3,4,5,6,7,8,9,10 -c 11 -o collapse | awk '{OFS="\t"; print $0, "pi_intergenic";}' >> ../data/for_junctiontable/mouse_adult_wt_RNAseq_PE50nt_strand.junctions2intron.all.final.merged.tab.4types

cut -f2-4 $juntions_only_in_4setoftranscript | sort -k3,3 -k1,1 -k2,2 | groupBy -g 3,1,2 -c 2 -o count | sort -k1,1 -k2,2 -k4,4nr > $juntions_only_in_4setoftranscript.stat

##################################################
### Intron signal (asked by Melissa)
##################################################
# For each intron, I want the following:  10 nts upstream (i.e. in the upstream exon) and 20 nts downstream (i.e. in the intron) of the 5' SS, and 50 nts upstream (i.e. in the intron) and 10 nts downstream (i.e. in the downstream exon).
cat $annotationFileORF | awk '{OFS="\t"; split($11,a,","); split($12,b,","); for(i=1;i<$10;i++) {if($6=="+") {print $1, $2+b[i]+a[i]-10, $2+b[i]+a[i], $4"_intron"i, i,$6; print $1, $2+b[i]+a[i], $2+b[i]+a[i]+50, $4"_intron"i, i,$6; print $1, $2+b[i+1]-50, $2+b[i+1], $4"_intron"i, i,$6;print $1, $2+b[i+1], $2+b[i+1]+10, $4"_intron"i, i,$6;} if($6=="-") {print $1, $2+b[i+1], $2+b[i+1]+10, $4"_intron"i, i,$6; print $1, $2+b[i+1]-50, $2+b[i+1], $4"_intron"i, i,$6; print $1, $2+b[i]+a[i], $2+b[i]+a[i]+50, $4"_intron"i, i,$6; print $1, $2+b[i]+a[i]-10, $2+b[i]+a[i], $4"_intron"i, i,$6; }}}' | fastaFromBed -tab -s -name -fi $GENOME/mm9/Sequence/BowtieIndex/genome.fa -bed stdin -fo stdout | groupBy -g 1 -c 2 -o collapse | awk '{OFS="\t";split($2,a,","); print $1,toupper(a[1])"|"tolower(a[2]), tolower(a[3])"|"toupper(a[4]);}' > $annotationFileORF.intronsSignal.tab

##################################################
### PAS (based on scanned result from Bo)
# copied from Xin's ~/nearline/transcript.structure/PAS
##################################################
echo "
args <- commandArgs(TRUE);
df=read.table(args[1], header=T)
#print(args[1])
rownames(df)=df[,1]; df=df[,-1]
df[is.na(df)]=-1
df0=df[apply(df>0,1,sum)>0,]
df0=cbind(df0, cloestMotif=apply(df0,1,function(x) colnames(df0)[which.max(x)]))
write.table(df0, paste('../data/for_PAS/fromXin/', basename(args[1]), '.2', sep=''), col.names = NA, quote =F, sep ='\t')
" > /tmp/fixpas.R
for i in ../data/for_PAS/fromXin/*.txt;
do
    echo $i;
    Rscript /tmp/fixpas.R $i;
    grep -v -P "^\t" $i".2" | cut -f15 | sort | uniq -c | sort -k2,2 | sed 's/^[ \t]*//g;s/ /\t/g'
done


























FLANKING=500
CENTER=150
D=20

annotationFile=$HOME/projects/piRNA/result/100312.transcript.with.transcripts.removed.new.names.bed  # bed12 format
#CAGE=/home/dongx/scratch/mouse_10dpp_wt_CAGE_100nt_strand/mouse_10dpp_wt_CAGE_100nt_strand.unique
CAGE=/home/dongx/scratch/mouse_adult_wt_CAGE_PE100nt_strand/mouse_adult_wt_CAGE_PE100nt_strand.unique
PAS=/home/dongx/scratch/mouse_adult_wt_PAS200_100nt/mouse_adult_wt_PAS200_100nt.unique
CAGE_bw_plus=$CAGE.plus.bw; CAGE_bw_minus=$CAGE.minus.bw
TSSregions=../data/annotation.TSSregion.bed
TTSregions=../data/annotation.TTSregion.bed

# clustering CTSS to define CT with IQR (0.1-0.9 range)
# if multiple summits with equal height, take the most 5' end one
# calculate kurtosis/skewness/peakheight/width/sum_of_reads for each TC
cat $CAGE.plus.bedGraph  | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); sum=0; for(i=1;i<=length(b);i++) sum+=b[i]; cum=0; for(i=1;i<=length(b);i++) {cum+=b[i]; split(a[i],c,":|-"); if(cum>=(0.1*sum) && cum<(0.9*sum)) print c[1], c[2],c[3], a[i], b[i]; if(cum>=(0.9*sum)) {print c[1], c[2],c[3], a[i], b[i]; break;}}}' | mergeBed -d $D -scores collapse -nms | awk -v strand="plus" -f cage2narrowPeak > $CAGE.+.narrowpeak
cat $CAGE.minus.bedGraph  | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, -$4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); sum=0; for(i=1;i<=length(b);i++) sum+=b[i]; cum=0; for(i=1;i<=length(b);i++) {cum+=b[i]; split(a[i],c,":|-"); if(cum>=(0.1*sum) && cum<(0.9*sum)) print c[1], c[2],c[3], a[i], b[i]; if(cum>=(0.9*sum)) {print c[1], c[2],c[3], a[i], b[i]; break;}}}' | mergeBed -d $D -scores collapse -nms | awk -v strand="minus" -f cage2narrowPeak > $CAGE.-.narrowpeak

# get TSS/TTS flanking region from annotation
grep -v track $annotationFile | awk -v flanking=$FLANKING '{OFS="\t"; tss=($6=="+")?$2:$3; tts=($6=="+")?$3:$2; print $1, tss-flanking,tss+flanking,$4,$5, $6 > "../data/annotation.TSSregion.bed"; print $1, tts-flanking,tts+flanking,$4,$5, $6 > "../data/annotation.TTSregion.bed";}'
# optional: take unique TSS/TTS only


## get BIP regions: two CAGE TCs with 1) peak distance <1000nt 2) on different orinetation (i.e. <--|-->) and 3) both peak height>=5)
# get 500nt-flanking regions of summit (h>=5)
awk '{OFS="\t";if($6=="+") print $1, $2+$10,$2+$10+1,$4"."$2,$5,$6}' $CAGE.narrowpeak | slopBed -b 500 -g /home/dongx/nearline/genomes/mm9/Annotation/Genes/ChromInfo.txt > ../data/cage.summit.500nt.+.bed
awk '{OFS="\t";if($6=="-") print $1, $2+$10,$2+$10+1,$4"."$3,$5,$6}' $CAGE.narrowpeak | slopBed -b 500 -g /home/dongx/nearline/genomes/mm9/Annotation/Genes/ChromInfo.txt > ../data/cage.summit.500nt.-.bed
# total # of bidirectional TC (minus.summit to plus.summit, midpoint of two TCs)
intersectBed -a ../data/cage.summit.500nt.-.bed -b ../data/cage.summit.500nt.+.bed -S -wo | awk '{OFS="\t"; split($4,a,"."); split($10,b,"."); if($2<=$8) print $1,$2+500,$9-500,a[1]":"b[1],$13,".",$13,$5,$11,int((a[2]+b[2])/2)-$2-500;}' > ../data/BIP.all.narrowPeak
# ucsc track
echo "track name=BIP_mouse_adult type=narrowPeak description=bidirectional_promoter_defined_by_adult_CAGE visibility=2 " > $CAGE.BIP.all.narrowPeak
cat ../data/BIP.all.narrowPeak >> $CAGE.BIP.all.narrowPeak
scp $CAGE.BIP.all.narrowPeak zlab:~/public_html/tracks/CAGE

# BIP in piRNA cluster
cat ../data/cage.summit.500nt.+.bed ../data/cage.summit.500nt.-.bed | intersectBed -a $TSSregions -b stdin -wo -S | awk '{OFS="\t"; split($10,a,"."); if(($6=="+" && $8<=($2+20)) || ($6=="-" && $2<=($8+20))) print}' | sort -k4,4 -k11,11nr -k13,13nr | awk 'BEGIN{id=""}{if($4!=id) {id=$4; if($11>=5) print}}' > $CAGE.BIP.inpirna.dis1000nt.off20nt.peak5.tab

##manual add
#chr12	99657758	99658953	12-qE-23911.1,12-qE-23911.2;12-qE-7089.2,12-qE-7089.1	2
cat ../data/cage.summit.500nt.+.bed ../data/cage.summit.500nt.-.bed | intersectBed -a $TSSregions -b stdin -wo -S | awk '{OFS="\t"; split($10,a,"."); if(($6=="+" && $8<=($2+20)) || ($6=="-" && $2<=($8+20))) print}' | sort -k4,4 -k11,11nr -k13,13nr | grep 12-qE-23911 | awk 'BEGIN{id=""}{if($4!=id) {id=$4; if($11>=1) print}}' >> $CAGE.BIP.inpirna.dis1000nt.off20nt.peak5.tab
cat ../data/cage.summit.500nt.+.bed ../data/cage.summit.500nt.-.bed | intersectBed -a $TSSregions -b stdin -wo -S | awk '{OFS="\t"; split($10,a,"."); if(($6=="+" && $8<=($2+20)) || ($6=="-" && $2<=($8+20))) print}' | sort -k4,4 -k11,11nr -k13,13nr | grep TC_p26389.20032169 >> $CAGE.BIP.inpirna.dis1000nt.off20nt.peak5.tab
cat ../data/cage.summit.500nt.+.bed ../data/cage.summit.500nt.-.bed | intersectBed -a $TSSregions -b stdin -wo -S | awk '{OFS="\t"; split($10,a,"."); if(($6=="+" && $8<=($2+20)) || ($6=="-" && $2<=($8+20))) print}' | sort -k4,4 -k11,11nr -k13,13nr | grep TC_p26393.20032727 >> $CAGE.BIP.inpirna.dis1000nt.off20nt.peak5.tab


# total number of piRNA genes with BIP n=95, Tx.n=214
cat $CAGE.BIP.inpirna.dis1000nt.off20nt.peak5.tab | cut -f4 | awk '{gsub(/\.[0-9a-d]+$/, "", $0); print}' | sort -u | wc -l

awk '{OFS="\t"; s=($2<$8)?($2+500):($8+500+1); e=($2<$8)?($8+500):($2+500); m=int((s+e)/2); print $4, $1":"s"-"e;}' $CAGE.BIP.inpirna.dis1000nt.off20nt.peak5.tab > ${TSSregions/bed/BIP}

# 150nt region flanking around the midpoint of BIP
awk '{OFS="\t"; s=$2+500; e=$8+500; m=int((s+e)/2); print $1, m-150, m+150;}' $CAGE.BIP.inpirna.dis1000nt.off20nt.peak5.tab | sort -u | awk '{OFS="\t"; print $0, "BIP_"NR;}' > ${TSSregions/bed/BIP.bed}
fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed ${TSSregions/bed/BIP.bed} -fo ${TSSregions/bed/BIP.fa}





## pi-Map3k9.1/2 doesnot have CAGE signal
## to uscs
#echo "track name=BIP_inpiRNA_mouse_adult type=narrowPeak description=BIP_inpiRNA_mouse_adult visibility=2 color=255,0,0" > $CAGE.BIP.inpirna.dis1000nt.off20nt.peak5.narrowPeak
#cat ../data/cage.summit.500nt.+.bed ../data/cage.summit.500nt.-.bed | intersectBed -a $CAGE.BIP.inpirna.dis1000nt.off20nt.peak5.tab -b stdin -wo -s | awk '{OFS="\t"; if($18>=2) print }' |sort -k4,4 -k20,20nr |  awk 'BEGIN{id=""}{OFS="\t"; if($4!=id) {id=$4; split($10,a,"."); split($17,b,"."); BIP=($6=="+")?(a[1]"."b[1]):(b[1]"."a[1]); s=($6=="+")?($8+500):($15+500); e=($6=="+")?($15+500):($8+500); m=int((b[2]+a[2])/2); print $1,$2,$3,$4,$5,$6,$7,s,e,BIP,1000,".", 1000,($6=="+")?$11:$18, ($6=="+")?$18:$11, m-s;}}' | cut -f7- | sort -u | awk '{if($2<$3) print }'>> $CAGE.BIP.inpirna.dis1000nt.off20nt.peak5.narrowPeak
#
### manual add
##chr7	77053969	77055149	7-qD1-9417.1;7-qD1-16444.1,7-qD1-16444.2	2
##chr12	99657758	99658953	12-qE-23911.1,12-qE-23911.2;12-qE-7089.2,12-qE-7089.1	2
#cat ../data/cage.summit.500nt.+.bed ../data/cage.summit.500nt.-.bed | intersectBed -a $CAGE.BIP.inpirna.dis1000nt.off20nt.peak5.tab -b stdin -wo -s | grep 7-qD1-9417.1 | sort -k4,4 -k20,20nr |  awk 'BEGIN{id=""}{OFS="\t"; if($4!=id) {id=$4; split($10,a,"."); split($17,b,"."); BIP=($6=="+")?(a[1]"."b[1]):(b[1]"."a[1]); s=($6=="+")?($8+500):($15+500); e=($6=="+")?($15+500):($8+500); m=int((b[2]+a[2])/2); print $1,$2,$3,$4,$5,$6,$7,s,e,BIP,1000,".", 1000,($6=="+")?$11:$18, ($6=="+")?$18:$11, m-s;}}' | cut -f7- | sort -u | awk '{if($2<$3) print }'>> $CAGE.BIP.inpirna.dis1000nt.off20nt.peak5.narrowPeak
#cat ../data/cage.summit.500nt.+.bed ../data/cage.summit.500nt.-.bed | intersectBed -a $CAGE.BIP.inpirna.dis1000nt.off20nt.peak5.tab -b stdin -wo -s | grep 12-qE-23911 | sort -k4,4 -k20,20nr |  awk 'BEGIN{id=""}{OFS="\t"; if($4!=id) {id=$4; split($10,a,"."); split($17,b,"."); BIP=($6=="+")?(a[1]"."b[1]):(b[1]"."a[1]); s=($6=="+")?($8+500):($15+500); e=($6=="+")?($15+500):($8+500); m=int((b[2]+a[2])/2); print $1,$2,$3,$4,$5,$6,$7,s,e,BIP,1000,".", 1000,($6=="+")?$11:$18, ($6=="+")?$18:$11, m-s;}}' | cut -f7- | sort -u | awk '{if($2<$3) print }'>> $CAGE.BIP.inpirna.dis1000nt.off20nt.peak5.narrowPeak
#
#scp $CAGE.BIP.inpirna.dis1000nt.off20nt.peak5.narrowPeak zlab:~/public_html/tracks/CAGE
#
#
#intersectBed -a $TSSregions -b ../data/BIP.all.narrowPeak -wo | awk '{OFS="\t"; if(($6=="+" && ($8+$16)<=($2+500)) || ($6=="-" && ($8+$16)>=($2+500))) print $0, $14+$15}' | sort -k4,4 -k18,18nr -k11,11nr | awk 'BEGIN{id=""}{if($4!=id) {id=$4; print}}' | awk -v fk=$CENTER '{OFS="\t"; summit=int($8+$16); print $1,summit-fk, summit+fk, $10}' | sort -u > ${TSSregions/bed/BIP.bed}
#
#
##intersectBed -a $TSSregions -b ../data/mouse_adult_wt_CAGE_PE100nt_strand.BIP.all.narrowPeak -wo | awk '{OFS="\t"; a=($14>$15)?$14:$15;b=($14>$15)?$15:$14; if(a/b>20) print $0, $14+$15}' | sort -k4,4 -k18,18nr -k11,11nr | awk 'BEGIN{id=""}{if($4!=id) {id=$4; print}}' | awk -v fk=$CENTER '{OFS="\t"; summit=int($8+$16); print $1,summit-fk, summit+fk, $10}' | sort -u > ${TSSregions/bed/BIP.bed}
#fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed ${TSSregions/bed/BIP.bed} -fo ${TSSregions/bed/BIP.fa}

## CAGE density for BIP region and its nucleotide content at each position of [-150, +150] region
# like the figure: http://www.plosgenetics.org/article/info:doi/10.1371/journal.pgen.0020047?imageURI=info:doi/10.1371/journal.pgen.0020047.g006


# promoter sequence
fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed $TSSregions -fo ${TSSregions/bed/fa}
# cpg score
cat ${TSSregions/bed/fa} | perl -ne '$h=$_; $h=~s/^>//; chomp($h); $s=<>; @n1 = ($s=~/CG/gi); $n1=@n1; $n2 = length($s); @n3 = (($s=~/C/gi), ($s=~/G/gi)); $n3=@n3; printf("%s\t%.3f\n", $h, 4*$n2*$n1/($n3**2+1));' > ${TSSregions/bed/cpg}
# mouse BIP from NIH
echo -e "chr\tstart\tend\tname\tscore\tstrand\tnormalizedCpG" > ../data/PC.BIP.cpg.bed
grep -v "^#" ../data/wgEncodeNhgriBip.mm9.bed8 | awk '{OFS="\t"; print $1, $7,$8,$4, 800, "+";}' | sort -u > /tmp/PC.BIP.cpg.bed
awk '{OFS="\t"; m=int(($3+$2)/2); print $1, m-1000,m+1000,$4,$5,$6}' /tmp/PC.BIP.cpg.bed | fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo stdout | perl -ne '$h=$_; $h=~s/^>//; chomp($h); $s=<>; @n1 = ($s=~/CG/gi); $n1=@n1; $n2 = length($s); @n3 = (($s=~/C/gi), ($s=~/G/gi)); $n3=@n3; printf("%s\t%.3f\n", $h, 4*$n2*$n1/($n3**2+1));' | paste /tmp/PC.BIP.cpg.bed - | cut -f1-6,8 >> ../data/PC.BIP.cpg.bed
echo -e "chr\tstart\tend\tname\tscore\tstrand\tintergenic_genic\tnormalizedCpG" > ../data/piRNA.BIP.cpg.bed
grep -v "^chrom" ../data/10092012.annotation.TSSregion.table | cut -f13,19 | grep -v "NA" | sort -u | sed 's/[:|-]/\t/g' | awk '{OFS="\t"; ; print $2, $3,$4,"BIP_"NR, 800, "+", $1}' > /tmp/piRNA.BIP.cpg.bed
awk '{OFS="\t"; m=int(($3+$2)/2); print $1, m-1000,m+1000,$4,$5,$6,$7}' /tmp/piRNA.BIP.cpg.bed | fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo stdout | perl -ne '$h=$_; $h=~s/^>//; chomp($h); $s=<>; @n1 = ($s=~/CG/gi); $n1=@n1; $n2 = length($s); @n3 = (($s=~/C/gi), ($s=~/G/gi)); $n3=@n3; printf("%s\t%.3f\n", $h, 4*$n2*$n1/($n3**2+1));' | paste /tmp/piRNA.BIP.cpg.bed - | cut -f1-7,9 >> ../data/piRNA.BIP.cpg.bed

echo "df1=read.table('../data/piRNA.BIP.cpg.bed', header=T); df2=read.table('../data/PC.BIP.cpg.bed', header=T); pdf('../result/BIP.piRNA_vs_BIP.pc.pdf', width=10, height=5); par(mfrow=c(1,2)); boxplot(df1[,8], df2[,7], names=c('BIP.piRNA', 'BIP.protein-coding'), ylab='normalized CpG score'); hist(df1[,3]-df1[,2], freq=F, breaks=50, xlim=c(0,800),col='red', main='', xlab='BIP length(bp)'); lines(density(df1[,3]-df1[,2], bw=25), col='green'); hist(df2[,3]-df2[,2], breaks=50,freq=F, add=T, col=rgb(.5, .5, .5, 0.5)); lines(density(df2[,3]-df2[,2], bw=25), col='black'); legend('topright', c('piRNA BIP', 'protein-coding BIP'), col=c('red','gray'), pch=15); dev.off()" > /tmp/cpg.R
Rscript /tmp/cpg.R; scp ../result/BIP.piRNA_vs_BIP.pc.pdf zlab:~/public_html/temp

# 150nt region flanking around the midpoint of BIP
awk '{OFS="\t"; m=int(($3+$2)/2); print $1, m-150,m+150,$4,$5,$6}' /tmp/PC.BIP.cpg.bed | sort -u | fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo ../data/PC.BIP.150nt.fa
awk '{OFS="\t"; m=int(($3+$2)/2); print $1, m-150,m+150,$4,$5,$6}' /tmp/piRNA.BIP.cpg.bed | sort -u | fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo ../data/piRNA.BIP.150nt.all.fa
awk '{OFS="\t"; m=int(($3+$2)/2); if($7=="intergenic") print $1, m-150,m+150,$4,$5,$6}' /tmp/piRNA.BIP.cpg.bed | sort -u | fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo ../data/piRNA.BIP.150nt.intergenic.fa
awk '{OFS="\t"; m=int(($3+$2)/2); if($7=="genic") print $1, m-150,m+150,$4,$5,$6}' /tmp/piRNA.BIP.cpg.bed | sort -u | fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo ../data/piRNA.BIP.150nt.genic.fa

# create frequency logo for the alignment (on weblogo)
# get CAGE density at each position
let CENTER2=CENTER*2
> ${TSSregions/bed/BIP.plus}; > ${TSSregions/bed/BIP.minus}
while read -r -a arr;
do
    region="${arr[3]}"
    a=`bigWigSummary $CAGE_bw_plus ${arr[0]} ${arr[1]} ${arr[2]} $CENTER2 | sed 's/n\/a/0/g'`
    [[ "$a" != "" ]] && echo $region $a >> ${TSSregions/bed/BIP.plus}
    [[ "$a" == "" ]] && awk -v r=$region -v n=$CENTER2 'BEGIN{printf("%s\t",r); for(i=1;i<n;i++) printf("%s\t", 0); printf("%s\n", 0);}' >> ${TSSregions/bed/BIP.plus}
    a=`bigWigSummary $CAGE_bw_minus ${arr[0]} ${arr[1]} ${arr[2]} $CENTER2 | sed 's/n\/a/0/g'`
    [[ "$a" != "" ]] && echo $region $a >> ${TSSregions/bed/BIP.minus}
    [[ "$a" == "" ]] && awk -v r=$region -v n=$CENTER2 'BEGIN{printf("%s\t",r); for(i=1;i<n;i++) printf("%s\t", 0); printf("%s\n", 0);}' >> ${TSSregions/bed/BIP.minus}
done <${TSSregions/bed/BIP.bed}


echo "
args <- commandArgs(TRUE)

pdf('BIP.density.pdf', width=20, height=6)
par(mfrow=c(2,1), mar=c(0.25,4,0.25,0.5))
df=read.table(args[1])
df=data.frame(df[,-1], row.names=df[,1])
df1=apply(df, 2, sum)
df=read.table(args[2])
df=data.frame(df[,-1], row.names=df[,1])
df2=apply(df, 2, sum)
#plot(df1, ylim=range(df1,-df2), type='h', xaxt='n', ylab='Plus strand')
plot(df1, type='h', xaxt='n', ylab='Plus strand')
abline(v=$CENTER, col='gray')
#plot(df2, ylim=range(-df1,df2),  type='h', xaxt='n', ylab='Minus strand')
plot(df2,  type='h', xaxt='n', ylab='Minus strand')
abline(v=$CENTER, col='gray')
dev.off()
" > /tmp/draw.R

Rscript /tmp/draw.R ${TSSregions/bed/BIP.plus} ${TSSregions/bed/BIP.minus}
scp BIP.density.pdf zlab:~/public_html/temp



## TATA
#cat ${TSSregions/bed/fa} | perl -slne '$h=$_; $h=~s/^>//; chomp($h); $s=<>; while ($s=~/TATAAA/gi) {printf("TATAAA\t%s\t%d\t%d\t\+\n", $h, pos($s)-6, pos($s));}' | awk -v flk=$FLANKING '{OFS="\t"; if($5=="+" && $3>=(flk-40) && $4<=(flk-20)) print $2,($3-flk)".TATAAA."($4-flk)}' > ${TSSregions/bed/TATA}
# use FIMO
# scan the region -35bp to -22bp (expected position for the TATA-box) upstream of TSS detected by CAGE.
mkdir ../data/motif
echo "
A  [16 352   3 354 268 360 222 ]
C  [46   0  10   0   0   3   2 ]
G  [18   2   2   5   0  20  44 ]
T  [309  35 374  30 121   6 121 ]
" | sed '/^$/d;s/  \[/| /;s/ \]//;s/ +/ /g'> ../data/motif/TATA.cm  # copied from JASPAR
jaspar2meme -cm ../data/motif/ > ../data/motif/TATA.meme
fimo --text --verbosity 1 --output-pthresh 0.005 ../data/motif/TATA.meme ${TSSregions/bed/fa}| awk -v flk=$FLANKING '{OFS="\t"; if($5=="+" && $3>=(flk-40) && $4<=(flk-20)) print $2,($3-flk)"."$8"."($4-flk)}' > ${TSSregions/bed/TATA.fimo}

# use MEME
# We selected a 14bp long sequence spanning the region from -35bp to -22bp (expected position for the TATA-box) upstream of TSS detected by CAGE and performed a de novo motif discovery using MEME
# see W-box?
awk '{if(NR%2==1) print; if(NR%2==0) print substr($1,460,20);}' ${TSSregions/bed/fa} | meme -oc meme_tata -o meme_tata -w 7 -nmotifs 10 -text -dna -nostatus stdin &
cat ${TSSregions/bed/BIP.fa} | meme -oc meme_bip -w 7 -nmotifs 10 -dna -nostatus stdin &

cat ../data/PC.BIP.150nt.fa | meme -oc meme_PC_BIP -w 8 -nmotifs 10 -dna -nostatus stdin &
cat ../data/piRNA.BIP.150nt.all.fa | meme -oc meme_pi_BIP_all -w 8 -nmotifs 10 -dna -nostatus stdin &
cat ../data/piRNA.BIP.150nt.intergenic.fa | meme -oc meme_pi_BIP_intergenic -w 8 -nmotifs 10 -dna -nostatus stdin &
cat ../data/piRNA.BIP.150nt.genic.fa | meme -oc meme_pi_BIP_genic -w 8 -nmotifs 10 -dna -nostatus stdin &

echo "
#!/bin/sh
#$ -V
#$ -cwd
#$ -pe single 8
#$ -o \$HOME/sge_jobs_output/sge_job.\$JOB_ID.out -j y
#$ -S /bin/bash
#$ -l mem_free=2G

id=\`basename \$1\`
echo \$id
meme_p -oc meme_\$id -minw 7 -maxw 16 -nmotifs 10 -dna -p 8 -maxsize 200000 \$1
tomtom -oc tomtom_\$id -min-overlap 5 -verbosity 1 meme_\$id/meme.txt motifdb/JASPAR_CORE_2009.meme motifdb/transfac.meme
" > /tmp/meme.sge
qsub /tmp/meme.sge ../data/PC.BIP.150nt.fa
qsub /tmp/meme.sge ../data/piRNA.BIP.150nt.all.fa
qsub /tmp/meme.sge ../data/piRNA.BIP.150nt.intergenic.fa
qsub /tmp/meme.sge ../data/piRNA.BIP.150nt.genic.fa

# join table
join -a 1 -1 1 -2 1 -e "NA" -o '0,1.2,2.2' <(sort ${TSSregions/bed/cpg}) <(sort ${TSSregions/bed/TATA.fimo}) | sed 's/ /\t/g' | sort | join -a 1 -1 1 -2 1 -e "NA" -o '0,1.2,1.3,2.2' - <(sort ${TSSregions/bed/BIP}) | sed 's/ /\t/g' | sort | paste <(sort -k4,4 $annotationFile) - | awk '{OFS="\t"; $13=($4~/^[0-9]+/)?"intergenic":"genic"; print }' > ${TSSregions/bed/tab}

## aggregation plot of RNAseq/smallRNA/CAGE/PAS/H3k4/Pol2/Pol3/amyb... at TSS/TTS flanking region
Rscript draw_aggregation_plot.R piRNA

## aggregation plot of RNAseq/smallRNA/CAGE/PAS/H3k4/Pol2/Pol3/amyb... at the flanking region of BIP midpoint
grep -v "^chrom" ../data/10092012.annotation.TSSregion.table | cut -f19 | grep -v "NA" | sort -u | sed 's/[:|-]/\t/g' | awk '{OFS="\t"; print $0, "+", "BIP_"NR;}' | awk -f bigWigAverageOverBed_generate_81bins.awk > ../data/10092012.annotation.TSSregion.table.BIP.81bins
qsub bigWigAverageOverBed_81bins.sh ../data/10092012.annotation.TSSregion.table.BIP.81bins 81
Rscript draw_aggregation_plot.BIP.R piRNA

awk '{OFS="\t"; print $1,$7,$8,$6, $4;}' ../data/wgEncodeNhgriBip.mm9.cpg.bed8 >   | awk -f bigWigAverageOverBed_generate_81bins.awk > ../data/wgEncodeNhgriBip.mm9.cpg.bed8.81bins
qsub bigWigAverageOverBed_81bins.sh ../data/wgEncodeNhgriBip.mm9.cpg.bed8.81bins 81
Rscript draw_aggregation_plot.BIP.R protein-coding

scp ../result/aggregation.*.pdf zlab:~/public_html/temp

## aggregation plot of RNAseq/smallRNA/CAGE/PAS/H3k4/Pol2/Pol3/amyb... at the flanking region of T2T midpoint


## define TC, and promoter type (using -20 and -100 as cutoff)


## Do the same for PAS region
# fimo
# meme

# distance to stopcodon



## special for BIP genes

## and the tail-to-tail genes

## intergenic vs. genic

## pachytene vs. pre-pachytene

## a table, each row is a piRNA precursor, columns are the feature (e.g. intergenic, pachytene, BIP, CpG, sharp/broad promoter)




#binned regions (into 100 windows, 20bp per window)
N=100
bedtools makewindows -b $TSSregions -n $N -i srcwinnum | fastaFromBed -name -s -fi $GENOME/mm9/Sequence/WholeGenomeFasta/genome.fa -bed stdin -fo stdout | perl -ne '$h=$_; $h=~s/^>//; chomp($h); $s=<>; @n1 = ($s=~/CG/g); $n1=@n1; $n2 = length($s); @n3 = (($s=~/C/g), ($s=~/G/g)); $n3=@n3; printf("%s\t%.3f\n", $h, 4*$n2*$n1/($n3**2)));'

## is piRNA precusor the result of re-purpose of mRNA? is the ORF broken? proteimics evidence?
