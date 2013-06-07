## script for merging piRNA cluster annotation using different data sources
# Author: Xianjun Dong
# Date: 2012.07.05
# Usage:
#!/bin/bash


##INPUTBED=../data/cufflinksOutput.in.piRNA_cluster.bed  # bed12 format
#
## get bed format of transcript annotation in piRNA cluster
#cufflinksOutput=$HOME/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand.transcripts.gtf.gz
## convert gtf to bed (Ref: https://lists.soe.ucsc.edu/pipermail/genome/2011-April/025696.html)
#zcat $cufflinksOutput | grep -v track | perl gtf2bed.pl - > ${cufflinksOutput/gz/bed}
##more +2 ../data/piRNA.clusters.coordinates.cpg.bed | intersectBed -wb -a stdin -b ${cufflinksOutput/gz/bed} -s | cut -f7- | sort -u | awk '{if(($3-$2)>100) print}' > $INPUTBED
#
## for all
##INPUTBED=../data/cufflinksOutput.all.bed  # bed12 format
#awk '{if(($3-$2)>100) print}' ${cufflinksOutput/gz/bed} > $INPUTBED

# for o100
INPUTBED=../data/cufflinksOutput.all.o100.bed  # bed12 format
cat ~/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks_o100/transcripts.gtf | grep -v track | perl gtf2bed.pl - | awk '{if(($3-$2)>100) print}' > $INPUTBED

samplename="mouse_adult_wt_RNAseq_PE50_strand"

# UCSC track
echo "track name=cufflinks_o100 description=\"$samplename (cufflinks_o100)\" visibility=pack colorByStrand='200,100,0 0,100,200'" > ${INPUTBED/bed/cufflinks+CAGE+PAS_final.bed}
cat $INPUTBED  >>  ${INPUTBED/bed/cufflinks+CAGE+PAS_final.bed}

D=20
CAGE=$HOME/scratch/ucsc/mouse_adult_wt_CAGE_PE100nt_strand.peak.Dis$D.bed
PAS=$HOME/scratch/ucsc/mouse_adult_wt_PAS200_100nt.peak.Dis$D.bed

# Step 1 (Optional): clustering the CAGE TSS and PAS site
# get CAGE peak (clustering CAGE TSS with 50nt distance thereshold, output all if multiple summits)
cat ~/scratch/mouse_adult_wt_CAGE_PE100nt_strand/mouse_adult_wt_CAGE_PE100nt_strand.unique.plus.bedGraph  | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; for(i=1;i<=length(b);i++) {if(b[i]>max) {max=b[i]; ii=1; delete imax; imax[ii++]=a[i];} else if(b[i]==max) {imax[ii++]=a[i];}} for(j=1;j<ii;j++) print $1,$2,$3,imax[j],max,"+";}' > $CAGE
cat ~/scratch/mouse_adult_wt_CAGE_PE100nt_strand/mouse_adult_wt_CAGE_PE100nt_strand.unique.minus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; for(i=1;i<=length(b);i++) {if(b[i]<max) {max=b[i]; ii=1; delete imax; imax[ii++]=a[i];} else if(b[i]==max) {imax[ii++]=a[i];}} for(j=1;j<ii;j++) print $1,$2,$3,imax[j],-max,"-";}' >> $CAGE
# get PAS peak (clustering PAS site with 50nt distance thereshold, output all if multiple summits)
cat ~/scratch/mouse_adult_wt_PAS200_100nt/mouse_adult_wt_PAS200_100nt.unique.plus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; for(i=1;i<=length(b);i++) {if(b[i]>max) {max=b[i]; ii=1; delete imax; imax[ii++]=a[i];} else if(b[i]==max) {imax[ii++]=a[i];}} for(j=1;j<ii;j++) print $1,$2,$3,imax[j],max,"+";}' > $PAS
cat ~/scratch/mouse_adult_wt_PAS200_100nt/mouse_adult_wt_PAS200_100nt.unique.minus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; for(i=1;i<=length(b);i++) {if(b[i]<max) {max=b[i]; ii=1; delete imax; imax[ii++]=a[i];} else if(b[i]==max) {imax[ii++]=a[i];}} for(j=1;j<ii;j++) print $1,$2,$3,imax[j],-max,"-";}' >> $PAS

## OLD VERSION: output the most 5' (for CAGE) and most 3' (for PAS) one if multiple
#cat ~/scratch/mouse_adult_wt_CAGE_PE100nt_strand/mouse_adult_wt_CAGE_PE100nt_strand.plus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; imax=1; for(i=1;i<=length(b);i++) if(b[i]>=max) {max=b[i];imax=i;} print $1,$2,$3,a[imax],max,"+";}' > $CAGE
#cat ~/scratch/mouse_adult_wt_CAGE_PE100nt_strand/mouse_adult_wt_CAGE_PE100nt_strand.minus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; imax=1; for(i=1;i<=length(b);i++) if(b[i]<max) {max=b[i];imax=i;} print $1,$2,$3,a[imax],-max,"-";}' >> $CAGE
#cat ~/scratch/mouse_adult_wt_PAS200_100nt/mouse_adult_wt_PAS200_100nt.plus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; imax=1; for(i=1;i<=length(b);i++) if(b[i]>=max) {max=b[i];imax=i;} print $1,$2,$3,a[imax],max,"+";}' > $PAS
#cat ~/scratch/mouse_adult_wt_PAS200_100nt/mouse_adult_wt_PAS200_100nt.minus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; imax=1; for(i=1;i<=length(b);i++) if(b[i]<max) {max=b[i];imax=i;} print $1,$2,$3,a[imax],-max,"-";}' >> $PAS


UP=${UPs[i]}; DOWN=${DOWNs[i]};
echo "D:$D; UP:$UP; DOWN:$DOWN";
###########################################
################# refine the TSS/TTS position using peak submit of CAGE/H3K4me3/PAS
###########################################

##use h/(abs(d)+1)/1000 to measure the proximal peaks, choose the one with best score
echo "correcting 5' end"
cat $INPUTBED | awk -v UP=1000 -v DOWN=500 '{OFS="\t"; s=($6=="+")?($2-UP):((($3-DOWN)>$2)?($3-DOWN):$2); e=($6=="+")?((($2+DOWN)<$3)?($2+DOWN):$3):($3+UP); print $1,(s>0)?s:0,e,$4,$5,$6,$0;}' | intersectBed -a stdin -b $CAGE -wao -s | awk '{OFS="\t"; if($19==".") {$19=$1;$20=$2;$21=$3;$22=($6=="+")?($7":"$8"-"($8+1)):($7":"($9-1)"-"$9);$23=1;$24=$6;} split($22,a,":|-"); $20=a[2]; $21=a[3]; d=($6=="+")?($8-$20):($21-$9); $25=d;if(d<0)d=-d; print $0, $23/(d+1)/1000;}' | sort -k4,4 -k26,26gr -k25,25n | awk 'BEGIN{id=""}{if($10!=id) {print;id=$10;} }' | awk '{OFS="\t"; for(i=19;i<25;i++)printf "%s\t", $i; for(i=1;i<19;i++) printf "%s\t", $i; print $25, $26;}' | ./changeTSS | sort -k1,1 -k2,2n > ${INPUTBED/bed/merge_CAGE.bed}

echo "correcting 3' end"
cut -f1-12 ${INPUTBED/bed/merge_CAGE.bed} | awk -v UP=1000 -v DOWN=2000 '{OFS="\t"; s=($6=="+")?((($3-UP)>$2)?($3-UP):$2):($2-DOWN); e=($6=="+")?($3+DOWN):((($2+UP)<$3)?($2+UP):$3); print $1,(s>0)?s:0,e,$4,$5,$6,$0;}' | intersectBed -a stdin -b $PAS -wao -s | awk '{OFS="\t"; if($19==".") {$19=$1;$20=$2;$21=$3;$22=($6=="-")?($7":"$8"-"($8+1)):($7":"($9-1)"-"$9);$23=1;$24=$6;} split($22,a,":|-"); $20=a[2]; $21=a[3]; d=($6=="-")?($8-$20):($21-$9); $25=d; if(d<0)d=-d;print $0, $23/(d+1)/1000;}' | sort -k4,4 -k26,26gr -k25,25n | awk 'BEGIN{id=""}{if($10!=id) {print;id=$10;} }' | awk '{for(i=19;i<25;i++)printf "%s\t", $i; for(i=1;i<19;i++) printf "%s\t", $i; print $25, $26;}' | ./changeTTS | sort -k1,1 -k2,2n > ${INPUTBED/bed/merge_CAGE_PAS.bed}

echo "merge the result"
paste ${INPUTBED/bed/merge_CAGE_PAS.bed} ${INPUTBED/bed/merge_CAGE.bed} $INPUTBED | awk '{for(i=27;i<29;i++) printf "%s\t", $i; for(i=13;i<15;i++)printf "%s\t", $i; for(i=1;i<13;i++)printf "%s\t", $i; for(i=15;i<27;i++)printf "%s\t", $i; for(i=29;i<40;i++)printf "%s\t", $i; print $40;}' > ${INPUTBED/bed/cufflinks+CAGE+PAS_final.txt}

echo "track name=cufflinks_o100_CAGE_PAS description=\"$samplename (cufflinks_o100 + CAGE + PAS correction)\" visibility=pack colorByStrand='200,100,0 0,100,200'" >> ${INPUTBED/bed/cufflinks+CAGE+PAS_final.bed}
cut -f1-12 ${INPUTBED/bed/merge_CAGE_PAS.bed} | grep -v -P "\-[0-9]" >> ${INPUTBED/bed/cufflinks+CAGE+PAS_final.bed}

echo "JOB $D + 1000.CAGE.500 + 1000.PAS.2000 DONE";

scp ${INPUTBED/bed/cufflinks+CAGE+PAS_final.bed} zlab:~/public_html/tracks/RNAseq

# for UCSC2PDF.R
sort -k1,1 -k2,2n ../data/cufflinksOutput.in.piRNA_cluster.merge_CAGE_PAS.bed | bedtools merge -s -nms -i stdin | sed 's/;/|/g' > ../data/piRNA_cluster.merge_CAGE_PAS.mergeBed.bed


echo "JOB DONE!"; exit


###########

UPs=(1000 2000 3000 4000)
DOWNs=(500 1000 2000 3000)
for D in 20
do
    CAGE=$HOME/scratch/ucsc/mouse_adult_wt_CAGE_PE100nt_strand.peak.Dis$D.bed
    PAS=$HOME/scratch/ucsc/mouse_adult_wt_PAS200_100nt.peak.Dis$D.bed

    # Step 1 (Optional): clustering the CAGE TSS and PAS site
    # get CAGE peak (clustering CAGE TSS with 50nt distance thereshold, output all if multiple summits)
    cat ~/scratch/mouse_adult_wt_CAGE_PE100nt_strand/mouse_adult_wt_CAGE_PE100nt_strand.unique.plus.bedGraph  | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; for(i=1;i<=length(b);i++) {if(b[i]>max) {max=b[i]; ii=1; delete imax; imax[ii++]=a[i];} else if(b[i]==max) {imax[ii++]=a[i];}} for(j=1;j<ii;j++) print $1,$2,$3,imax[j],max,"+";}' > $CAGE
    cat ~/scratch/mouse_adult_wt_CAGE_PE100nt_strand/mouse_adult_wt_CAGE_PE100nt_strand.unique.minus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; for(i=1;i<=length(b);i++) {if(b[i]<max) {max=b[i]; ii=1; delete imax; imax[ii++]=a[i];} else if(b[i]==max) {imax[ii++]=a[i];}} for(j=1;j<ii;j++) print $1,$2,$3,imax[j],-max,"-";}' >> $CAGE
    # get PAS peak (clustering PAS site with 50nt distance thereshold, output all if multiple summits)
    cat ~/scratch/mouse_adult_wt_PAS200_100nt/mouse_adult_wt_PAS200_100nt.unique.plus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; for(i=1;i<=length(b);i++) {if(b[i]>max) {max=b[i]; ii=1; delete imax; imax[ii++]=a[i];} else if(b[i]==max) {imax[ii++]=a[i];}} for(j=1;j<ii;j++) print $1,$2,$3,imax[j],max,"+";}' > $PAS
    cat ~/scratch/mouse_adult_wt_PAS200_100nt/mouse_adult_wt_PAS200_100nt.unique.minus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; for(i=1;i<=length(b);i++) {if(b[i]<max) {max=b[i]; ii=1; delete imax; imax[ii++]=a[i];} else if(b[i]==max) {imax[ii++]=a[i];}} for(j=1;j<ii;j++) print $1,$2,$3,imax[j],-max,"-";}' >> $PAS

    ## OLD VERSION: output the most 5' (for CAGE) and most 3' (for PAS) one if multiple
    #cat ~/scratch/mouse_adult_wt_CAGE_PE100nt_strand/mouse_adult_wt_CAGE_PE100nt_strand.plus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; imax=1; for(i=1;i<=length(b);i++) if(b[i]>=max) {max=b[i];imax=i;} print $1,$2,$3,a[imax],max,"+";}' > $CAGE
    #cat ~/scratch/mouse_adult_wt_CAGE_PE100nt_strand/mouse_adult_wt_CAGE_PE100nt_strand.minus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; imax=1; for(i=1;i<=length(b);i++) if(b[i]<max) {max=b[i];imax=i;} print $1,$2,$3,a[imax],-max,"-";}' >> $CAGE
    #cat ~/scratch/mouse_adult_wt_PAS200_100nt/mouse_adult_wt_PAS200_100nt.plus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; imax=1; for(i=1;i<=length(b);i++) if(b[i]>=max) {max=b[i];imax=i;} print $1,$2,$3,a[imax],max,"+";}' > $PAS
    #cat ~/scratch/mouse_adult_wt_PAS200_100nt/mouse_adult_wt_PAS200_100nt.minus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; imax=1; for(i=1;i<=length(b);i++) if(b[i]<max) {max=b[i];imax=i;} print $1,$2,$3,a[imax],-max,"-";}' >> $PAS

    for i in 0 1 2 3
    do
        UP=${UPs[i]}; DOWN=${DOWNs[i]};
        echo "D:$D; UP:$UP; DOWN:$DOWN";
        ###########################################
        ################# refine the TSS/TTS position using peak submit of CAGE/H3K4me3/PAS
        ###########################################

        ##use h/(abs(d)+1)/1000 to measure the proximal peaks, choose the one with best score
        cat $INPUTBED | awk -v UP=$UP -v DOWN=$DOWN '{OFS="\t"; s=($6=="+")?($2-UP):((($3-DOWN)>$2)?($3-DOWN):$2); e=($6=="+")?((($2+DOWN)<$3)?($2+DOWN):$3):($3+UP); print $1,s,e,$4,$5,$6,$0;}' | intersectBed -a stdin -b $CAGE -wao -s | awk '{OFS="\t"; if($19==".") {$19=$1;$20=$2;$21=$3;$22=($6=="+")?($7":"$8"-"($8+1)):($7":"($9-1)"-"$9);$23=1;$24=$6;} split($22,a,":|-"); $20=a[2]; $21=a[3]; d=($6=="+")?($8-$20):($21-$9); $25=d;if(d<0)d=-d; print $0, $23/(d+1)/1000;}' | sort -k4,4 -k26,26gr -k25,25n | awk 'BEGIN{id=""}{if($10!=id) {print;id=$10;} }' | awk '{OFS="\t"; for(i=19;i<25;i++)printf "%s\t", $i; for(i=1;i<19;i++) printf "%s\t", $i; print $25, $26;}' | ./changeTSS | sort -k1,1 -k2,2n > ${INPUTBED/bed/merge_CAGE.bed}

        cut -f1-12 ${INPUTBED/bed/merge_CAGE.bed} | awk -v UP=$DOWN -v DOWN=$UP '{OFS="\t"; s=($6=="+")?((($3-UP)>$2)?($3-UP):$2):($2-DOWN); e=($6=="+")?($3+DOWN):((($2+UP)<$3)?($2+UP):$3); print $1,s,e,$4,$5,$6,$0;}' | intersectBed -a stdin -b $PAS -wao -s | awk '{OFS="\t"; if($19==".") {$19=$1;$20=$2;$21=$3;$22=($6=="-")?($7":"$8"-"($8+1)):($7":"($9-1)"-"$9);$23=1;$24=$6;} split($22,a,":|-"); $20=a[2]; $21=a[3]; d=($6=="-")?($8-$20):($21-$9); $25=d; if(d<0)d=-d;print $0, $23/(d+1)/1000;}' | sort -k4,4 -k26,26gr -k25,25n | awk 'BEGIN{id=""}{if($10!=id) {print;id=$10;} }' | awk '{for(i=19;i<25;i++)printf "%s\t", $i; for(i=1;i<19;i++) printf "%s\t", $i; print $25, $26;}' | ./changeTTS | sort -k1,1 -k2,2n > ${INPUTBED/bed/merge_CAGE_PAS.bed}

        paste ${INPUTBED/bed/merge_CAGE_PAS.bed} ${INPUTBED/bed/merge_CAGE.bed} $INPUTBED | awk '{for(i=27;i<29;i++) printf "%s\t", $i; for(i=13;i<15;i++)printf "%s\t", $i; for(i=1;i<13;i++)printf "%s\t", $i; for(i=15;i<27;i++)printf "%s\t", $i; for(i=29;i<40;i++)printf "%s\t", $i; print $40;}' > ${INPUTBED/bed/cufflinks+dis$D+$UP.CAGE.$DOWN+$UP.PAS.$DOWN.txt}

        echo "track name=cufflinks.D$D.${UP}CAGE$DOWN.${UP}PAS$DOWN description=\"$samplename (CUFFLINKS + Dis$D + $UP.CAGE.$DOWN + $UP.PAS.$DOWN)\" visibility=pack colorByStrand='200,100,0 0,100,200'" >> ${INPUTBED/bed/cufflinks+CAGE+PAS.bed}
        cut -f1-12 ${INPUTBED/bed/merge_CAGE_PAS.bed} | grep -v -P "\-[0-9]" >> ${INPUTBED/bed/cufflinks+CAGE+PAS.bed}

        echo "JOB $D + $UP.CAGE.$DOWN + $UP.PAS.$DOWN DONE";
    done
done

scp ${INPUTBED/bed/cufflinks+CAGE+PAS.bed} zlab:~/public_html/tracks/RNAseq

echo "JOB DONE!"; exit















## OLD VERSION

###########################################
################# 0. configuring
###########################################

mergeDistance=200

repeatFile=$HOME/projects/piRNA/data/mm9.rmsk.LINE_LTR.tab
cufflinksOutput1=$HOME/scratch/ucsc/mouse_adult_wt_smallRNAseq_76nt_strand_ox.transcripts.gtf.gz
cufflinksOutput2=$HOME/scratch/ucsc/mouse_adult_wt_smallRNAseq_76nt_strand_unox.transcripts.gtf.gz

cufflinksOutput=$HOME/scratch/ucsc/Phil.SRA.wt.ox.6wk.testes.raw.xkxh.norm.all0.transcripts.gtf.gz
CAGE=$HOME/scratch/ucsc/mouse_adult_wt_CAGE_PE100nt_strand.peak.bed
PAS=$HOME/scratch/ucsc/mouse_adult_wt_PAS200_100nt.peak.bed

samplename="mouse_adult_wt_smallRNAseq_merge"

[ -d ~/scratch/$samplename ] || mkdir ~/scratch/$samplename; cd ~/scratch/$samplename

###########################################
################# 1.1. merge overlapped transcript from de novo assmelby (e.g. cufflinks)
###########################################

echo "track name=merge1_$samplename description=merge1_$samplename.cufflinks visibility=pack colorByStrand='200,100,0 0,100,200'" > merge1_cufflinks.bed

zcat $cufflinksOutput1 $cufflinksOutput2 | grep -w transcript | grep "transcript_id \"CUFF" | cut -f1,4-7,9 | sed 's/ /\t/g;s/"//g;s/;//g' | awk '{OFS="\t";print $1,$2,$3,$9,$4,$5}' | sort -k1,1 -k2,2n | bedtools merge -s -nms -i stdin | awk '{OFS="\t"; print $1,$2,$3,$4,1000,$5}' > merge1_cufflinks.bed.cp
#zcat $cufflinksOutput2 | grep -w transcript | grep "transcript_id \"CUFF" | cut -f1,4-7,9 | sed 's/ /\t/g;s/"//g;s/;//g' | awk '{OFS="\t";print $1,$2,$3,$9,$4,$5}' > cufflinksOutput2.bed

cat merge1_cufflinks.bed.cp >> merge1_cufflinks.bed

###########################################
################# 1.2. merge contigs into transcript using NorahDesk
###########################################
echo "track name=merge2_$samplename description=merge2_$samplename.NorahDesk visibility=pack colorByStrand='200,100,0 0,100,200'" > merge2_NorahDesk.bed
sort -k1,1 -k2,2n cufflinksOutput.bed ncRNA_prediction.bed | bedtools merge -s -nms -i stdin | awk '{OFS="\t"; print $1,$2,$3,$4,1000,$5}' > merge2_NorahDesk.bed.cp
cat merge2_NorahDesk.bed.cp >> merge2_NorahDesk.bed

###########################################
################# 2. merge gaps overlaping with LTR/LINE repeat elements
###########################################

rm repeatFile.nostrand.bed
awk '{OFS="\t";if($1 !~ /^#/) {print $2,$3,$4,$6,$1,"+"; print $2,$3,$4,$6,$1,"-";} }' $repeatFile | sort -k1,1 -k2,2n | bedtools merge -s -nms -i stdin | sed 's/;/_/g' > repeatFile.nostrand.bed
awk '{OFS="\t";if($1 !~ /^#/) print $2,$3,$4,$6,$1,$5;}' $repeatFile | sort -k1,1 -k2,2n | bedtools merge -nms -i stdin | awk '{OFS="\t";print $1,$2,$3,$4,1000,"+"; print $1,$2,$3,$4,1000,"-"}' | sed 's/;/:/g' > repeatFile.nostrand.bed
sort -k1,1 -k2,2n merge1_cufflinks.bed.cp repeatFile.nostrand.bed | bedtools merge -s -nms -i stdin | awk '{OFS="\t"; if($4 ~ /CUFF/) print $1,$2,$3,$4,1000,$5}' > merge3_repeak.bed.cp
bedtools merge -s -d 50 -nms -i merge3_repeak.bed.cp | awk '{OFS="\t"; print $1,$2,$3,$4,1000,$5}' > merge3_repeak.bed.cp2
echo "track name=merge3_$samplename description=merge3_$samplename.repeak visibility=pack colorByStrand='200,100,0 0,100,200'" > merge3_repeak.bed
cat merge3_repeak.bed.cp2 >> merge3_repeak.bed

###########################################
################# 3. refine the TSS/TTS position using peak submit of CAGE/H3K4me3/PAS
###########################################

echo "track name=merge_$samplename description=merge_$samplename.CAGE visibility=pack colorByStrand='200,100,0 0,100,200'" > merge_CAGE.bed
# convert gtf to bed (Ref: https://lists.soe.ucsc.edu/pipermail/genome/2011-April/025696.html)
#zcat $cufflinksOutput | grep -v track | gtfToGenePred stdin stdout | genePredToBed > transcripts.bed
zcat $cufflinksOutput | grep -v track | perl gtf2bed.pl - > transcripts.bed

more +2 ../data/piRNA.clusters.coordinates.cpg.bed | intersectBed -wb -a stdin -b transcripts.bed -s | cut -f7- | sort -u > piRNA_cluster.bed
more +2 ../data/piRNA.clusters.coordinates.cpg.bed | grep -v "pi-" | intersectBed -wb -a stdin -b transcripts.bed -s | cut -f7- | sort -u | sort -k1,1 -k2,2n > intergenic_piRNA_cluster.bed
bedtools merge -s -nms -d 200 -i piRNA_cluster.bed > piRNA_cluster_merge.bed

# multiple max, choose the closest one
cat intergenic_piRNA_cluster.bed | awk  -v UP=4000 -v DOWN=2000 '{OFS="\t"; s=($6=="+")?($2-UP):($3-DOWN); e=($6=="+")?($2+DOWN):($3+UP); print $1,s,e,$4,$5,$6,$0;}' | intersectBed -a stdin -b $CAGE -wao -s | awk '{OFS="\t"; d=($6=="+")?($8-$20):($9-$21);if(d<0)d=-d;print $0,d;}' | sort -k4,4 -k23,23nr -k26,26n | awk 'BEGIN{id=""}{if($10!=id) {print;id=$10;} }' | awk '{for(i=19;i<25;i++)printf "%s\t", $i; for(i=1;i<19;i++) printf "%s\t", $i; print $25;}' | ./changeTSS | sort -k1,1 -k2,2n > intergenic_piRNA_cluster.merge_CAGE.bed

cut -f1-12 intergenic_piRNA_cluster.merge_CAGE.bed | awk  -v UP=2000 -v DOWN=4000 '{OFS="\t"; s=($6=="+")?($3-UP):($2-DOWN); e=($6=="+")?($3+DOWN):($2+UP); print $1,s,e,$4,$5,$6,$0;}' | intersectBed -a stdin -b $PAS -wao -s | awk '{OFS="\t"; d=($6=="-")?($8-$20):($9-$21);if(d<0)d=-d;print $0,d;}' | sort -k4,4 -k23,23nr -k26,26n | awk 'BEGIN{id=""}{if($10!=id) {print;id=$10;} }' | awk '{for(i=19;i<25;i++)printf "%s\t", $i; for(i=1;i<19;i++) printf "%s\t", $i; print $25;}' | ./changeTTS | sort -k1,1 -k2,2n > intergenic_piRNA_cluster.merge_CAGE_PAS.bed

paste intergenic_piRNA_cluster.merge_CAGE_PAS.bed intergenic_piRNA_cluster.merge_CAGE.bed intergenic_piRNA_cluster.bed | awk '{for(i=27;i<29;i++) printf "%s\t", $i; for(i=13;i<15;i++)printf "%s\t", $i; for(i=1;i<13;i++)printf "%s\t", $i; for(i=15;i<27;i++)printf "%s\t", $i; for(i=29;i<40;i++)printf "%s\t", $i; print $40;}' > intergenic_piRNA_cluster.CAGE.PAS.txt

# UCSC track
echo "track name=cufflinks description=\"intergenic_piRNA_cluster (CUFFLINKS)\" visibility=pack colorByStrand='200,100,0 0,100,200'" > intergenic_piRNA_cluster.cufflinks+CAGE+PAS.bed
cat intergenic_piRNA_cluster.bed >>  !$
echo "track name=cufflinks.CAGE description=\"intergenic_piRNA_cluster (CUFFLINKS + CAGE)\" visibility=pack colorByStrand='200,100,0 0,100,200'" >> !$
cut -f1-12 intergenic_piRNA_cluster.merge_CAGE.bed >> !$
echo "track name=cufflinks.CAGE.PAS description=\"intergenic_piRNA_cluster (CUFFLINKS + CAGE + PAS)\" visibility=pack colorByStrand='200,100,0 0,100,200'" >> !$
cut -f1-12 intergenic_piRNA_cluster.merge_CAGE_PAS.bed >> !$

scp !$ zlab:~/public_html/tracks/RNAseq

echo "JOB DONE";

#after the job done
# scp merge*.bed zlab:~/public_html/tracks/RNAseq
