
# remove the replace/removed ones
more +2 toremove | sed 's/[;,]/\n/g;s/ //g;' | sed '/^$/d' | sort -u | grep -iv -f - transcripts.bed > correction/transcript.bed

# add the added/substituted ones
more +2 toadd | sed 's/\t/\n/g;s/ //g' | sed '/^$/d' | grep -i CUFF | grep -v "(" | sed 's/[;|,]/\n/g;s/ //g' | sort -u | grep -i -f - ../cufflinks_o100/transcripts.gtf | perl $HOME/projects/piRNA/src/gtf2bed.pl - | awk '{OFS="\t";$4=$4"_o100"; print}' >> correction/transcript.bed
more +2 toadd | sed 's/\t/\n/g;s/ //g' | sed '/^$/d' | grep "(-o200" | sed 's/(-o200)//g;s/[;|,]/\n/g;s/ //g' | grep -i -f - ../cufflinks_o200/transcripts.gtf | perl $HOME/projects/piRNA/src/gtf2bed.pl - | awk '{OFS="\t";$4=$4"_o200"; print}' >>  correction/transcript.bed
more +2 toadd | sed 's/\t/\n/g;s/ //g' | sed '/^$/d' | grep "(-j0.2" | sed 's/(-j0.2)//g;s/[;|,]/\n/g;s/ //g' | grep -i -f - ../cufflinks_j0.2/transcripts.gtf | perl $HOME/projects/piRNA/src/gtf2bed.pl - | awk '{OFS="\t";$4=$4"_j0.2"; print}'  >>  correction/transcript.bed
more +2 toadd | sed 's/\t/\n/g;s/ //g' | sed '/^$/d' | grep "(-o350" | sed 's/(-o350)//g;s/[;|,]/\n/g;s/ //g' | grep -i -f - ../cufflinks_o350/transcripts.gtf | perl $HOME/projects/piRNA/src/gtf2bed.pl - | awk '{OFS="\t";$4=$4"_o350"; print}'  >>  correction/transcript.bed
more +2 toadd | sed 's/\t/\n/g;s/ //g' | sed '/^$/d' | grep manual | sed 's/(manual)//g;s/[;|,]/\n/g;s/ //g' | grep -i -f - /home/dongx/projects/piRNA/data/piRNA.clusters.coordinates.bed12 | awk '{OFS="\t";$4=$4"_manual"; print}' >>  correction/transcript.bed

# corret by CAGE/PAS with default setting
correct53byCAGEPAS correction/transcript.bed 1000 500 20 2000 1000 20 | awk '{if($11!~/\-/) print}' > correction/transcript.corrected.bed

# add the one with just coordinate e.g. chr10:86,028,268-86,054,652
more +2 toadd | sed 's/\t/\n/g;s/ //g' | sed '/^$/d' | grep chr | grep -v "(" | sed 's/,//g;s/;/\n/g;s/ //g' | awk '{OFS="\t";split($1,a,":");split(a[2],b,"-");print a[1],b[1],b[2],$1,1000,a[3],b[1],b[2],0,1,b[2]-b[1],0;}' | correct53byCAGEPAS2 - 1000 500 20 2000 1000 20 | awk '{OFS="\t"; $4=$1":"$2"-"$3":"$6; print}' >  correction/newlyadded.transcript.bed
more +2 toadd | sed 's/\t/\n/g;s/ //g' | sed '/^$/d' | grep chr | grep "(" | sed 's/(.*//g;s/,//g;s/;/\n/g;s/ //g' | awk '{OFS="\t";split($1,a,":");split(a[2],b,"-");print a[1],b[1],b[2],$1,1000,a[3],b[1],b[2],0,1,b[2]-b[1],0;}' | correct53byCAGEPAS3 - 1000 500 20 2000 1000 20 | awk '{OFS="\t"; $4=$1":"$2"-"$3":"$6; print}' >>  correction/newlyadded.transcript.bed


# for those need to further correct by CAGE/PAS
> correction/modified.Tx.bed
# case like "CUFF.18534.1	3000,2000	3000,2000"
grep "," tomodify | sed 's/ //g' | grep -v -P "\(|chr" |awk -F"\t" '{split($1,ids,";"); for(i=1;i<=length(ids);i++){if($2~/,/) {split($2,cage,","); cageUp=cage[1]; cageDOWN=cage[2]; cageD=(cage[3]=="")?20:cage[3];} else {cageUp=1000; cageDOWN=500; cageD=20}; if($3~/,/) {split($3,pas,","); pasUp=pas[1]; pasDOWN=pas[2]; pasD=(pas[3]=="")?20:pas[3];} else {pasUp=2000; pasDOWN=1000; pasD=20};  print "cat transcripts.bed | grep -i "ids[i]" | correct53byCAGEPAS - "cageUp" "cageDOWN" "cageD" "pasUp" "pasDOWN" "pasD" >> correction/modified.Tx.bed";}}'

# case like "CUFF.15390.1(-o100)	3000,2000	3000,2000"
grep "," tomodify | sed 's/ //g' | grep -v chr | grep "(-o100" | sed 's/(-o100)//g;' | awk -F"\t" '{split($1,ids,";"); for(i=1;i<=length(ids);i++){if($2~/,/) {split($2,cage,","); cageUp=cage[1]; cageDOWN=cage[2]; cageD=(cage[3]=="")?20:cage[3];} else {cageUp=1000; cageDOWN=500; cageD=20}; if($3~/,/) {split($3,pas,","); pasUp=pas[1]; pasDOWN=pas[2]; pasD=(pas[3]=="")?20:pas[3];} else {pasUp=2000; pasDOWN=1000; pasD=20};  print "grep -i "ids[i]" ../cufflinks_o100/transcripts.gtf | perl $HOME/projects/piRNA/src/gtf2bed.pl - | correct53byCAGEPAS - "cageUp" "cageDOWN" "cageD" "pasUp" "pasDOWN" "pasD" >> correction/modified.Tx.bed";}}'

# for case like "Cuff.1525.1	chr10:62114296-62114934:+	chr10:62167647-62169773:+"
grep "chr" tomodify | sed 's/ //g;s/,//g' | awk -F"\t" '{OFS="\t";split($1,a,";");for(i=1;i<=length(a);i++) print a[i],$2,$3;}' | awk -F"\t" '{OFS="\t";split($2,a,";"); split($3,b,";"); if(length(a)==0) a[1]=$2; if(length(b)==0) b[1]=$3; for(i=1;i<=length(a);i++) for(j=1;j<=length(b);j++) print $1,a[i],b[j];}' | tr "\t" "|" | snip - > correction/modified.bed
perl $HOME/projects/piRNA/src/gtf2bed.pl CUFF.7441.1.gtf >> correction/modified.bed



more +2 tomodify | cut -f1| sort -u | grep -v "(" |sed 's/;/\n/g;s/ //g;/^$/d' | grep -iv -f - correction/transcript.corrected.bed > correction/transcript.corrected.nomodified.bed

cd correction
echo "track name=transcript_CUFFLINKS.CAGE.PAS_ALL visibility=pack colorByStrand='200,100,0 0,100,200'" > transcript.ALL.bed
cat transcript.corrected.nomodified.bed newlyadded.transcript.bed modified.Tx.bed modified.bed >> transcript.ALL.bed
echo "track name=transcript_CUFFLINKS.CAGE.PAS_no visibility=pack colorByStrand='200,100,0 0,100,200'" > transcript.NO.bed
intersectBed -a transcript.ALL.bed -b no.commenttaken.bed -s -wa >> transcript.NO.bed
echo "track name=transcript_CUFFLINKS.CAGE.PAS_yes visibility=pack colorByStrand='200,100,0 0,100,200'" > transcript.YES.bed
intersectBed -a transcript.ALL.bed -b no.commenttaken.bed -s -v -wa >> transcript.YES.bed
#gzip -f transcript.corrected.modified.bed
scp transcript.ALL.bed transcript.NO.bed transcript.YES.bed zlab:~/public_html/tracks/RNAseq

# make nav page for the clusters
annotation=../result/100312.transcript.with.transcripts.removed.new.names.bed
echo "<script src='http://compbio2.csail.mit.edu/sorttable.js'></script>
<a href=old>previous version</a><br><br>
<table border=2 class='sortable'>
<caption>piRNA precursor annotation (v3)</caption>
<tr><td>ID</td><td>Coordinate</td><td><abbr title=\"cluster name\">name</abbr></td><td>#<abbr title=\"isoform number\">Tx</abbr></td><td><abbr title=\"strand\">S</abbr></td><td><abbr title=\"cluster type\">Type</abbr></td><td><abbr title=\"promoter type\">BIP</abbr></td></tr>" > nav.html
sort -k1,1 -k2,2n $annotation | mergeBed -s -n -nms | awk '{OFS="\t"; split($4,a,";"); split("", d); for(i=1;i<=length(a);i++) {gsub(/\.[0-9a-d]+$/, "", a[i]);d[a[i]]=1;} id="";for(j in d) id=j";"id; gsub(/;$/,"",id);$4=id; print}' | sort -k1,1 -k2,2n | awk '{OFMT = "%.2f"; print "<tr><td>cluster_"NR"</td>", "<td><a href=img/region3_"$1"."$2"."$3".png target=png>"$1":"$2"-"$3"</a></td>","<td>"$4"</td>","<td>"$5"</td>", "<td>"$6"</td>", "<td></td><td></td></tr>"; }' >> nav.html
echo "</table>" >> nav.html
scp nav.html zlab:~/public_html/pirna

# generate PDF
sort -k1,1 -k2,2n $annotation | mergeBed -s -n -nms | Rscript UCSC2PDF.R stdin &

# Tx table
snip()
{
    #myfile=$1
    while IFS='|' read -r -a arr
    do
        #echo -e $line
        #arr=$(echo -e $line | sed 's/,//g' | tr "\t" "\n")
        ID0="${arr[0]}"
        [[ $ID0 = *\(* ]] || Txbed=transcripts.bed
        [[ $ID0 = *\(-o100* ]] && Txbed=transcripts_o100.bed
        ID=${ID0/\(*/}
        cage="${arr[1]}"
        pas="${arr[2]}"
        #echo "$ID:$cage:$pas:$Txbed"
        #continue
        [ "$cage" == "" ] && tss=`grep -m 1 -i $ID $Txbed | awk '{OFS="\t"; if($6=="+") print $1,$2-1000,$2+1000,$4,$5,$6; else print $1,$3-1000,$3+1000,$4,$5,$6;}' | intersectBed -a stdin -b $HOME/scratch/ucsc/mouse_adult_wt_CAGE_PE100nt_strand.peak.Dis20.bed -wao -s | sort -k4,4 -k11,11gr -k13,13n | head -n1 | awk '{if($7!=".") {split($10,a,":|-"); print ($12=="+")?a[2]:a[3];} else print ($2+$3)/2;}'`
        [ "$cage" != "" ] && tss=`echo $cage | awk '{OFS="\t"; split($1,a,":"); split(a[2],b,"-"); print a[1],b[1],b[2],".",1000, a[3];}' | intersectBed -a stdin -b $HOME/scratch/ucsc/mouse_adult_wt_CAGE_PE100nt_strand.peak.Dis20.bed -wao -s | sort -k4,4 -k11,11gr -k13,13n | head -n1 | awk '{if($7!=".") {split($10,a,":|-"); print ($12=="+")?a[2]:a[3];} else print ($2+$3)/2;}'`

        [ "$pas" == "" ] && tts=`grep -m 1 -i $ID $Txbed | awk '{OFS="\t"; if($6=="+") print $1,$3-1000,$3+1000,$4,$5,$6; else print $1,$2-1000,$2+1000,$4,$5,$6;}' | intersectBed -a stdin -b $HOME/scratch/ucsc/mouse_adult_wt_PAS200_100nt.peak.Dis20.bed -wao -s | sort -k4,4 -k11,11gr -k13,13n | head -n1 | awk '{if($7!=".") {split($10,a,":|-"); print ($12=="+")?a[3]:a[2];} else print ($2+$3)/2;}'`
        [ "$pas" != "" ] && tts=`echo $pas | awk '{OFS="\t"; split($1,a,":"); split(a[2],b,"-"); print a[1],b[1],b[2],".",1000, a[3];}' | intersectBed -a stdin -b $HOME/scratch/ucsc/mouse_adult_wt_PAS200_100nt.peak.Dis20.bed -wao -s | sort -k4,4 -k11,11gr -k13,13n | head -n1 | awk '{if($7!=".") {split($10,a,":|-"); print ($12=="+")?a[3]:a[2];} else print ($2+$3)/2;}'`

        grep -m 1 -i $ID $Txbed | changeTSSTTS_bed12 -v TSS=$tss -v TTS=$tts | awk -v id=$ID0 '{OFS="\t"; $4=id; print}' | sed 's/[\(\)]//g;'
    done
}

# all mapper
correct53byCAGEPAS3 ()
{
    INPUTBED=$1
    UP=$2
    DOWN=$3
    D=$4

    OUT="$(mktemp)"

    CAGE=$HOME/scratch/ucsc/mouse_adult_wt_CAGE_PE100nt_strand.allpeak.Dis$D.bed
    [ -e $CAGE ] || {
        echo "cluster CAGE with distance cutoff $D"
        cat ~/scratch/mouse_adult_wt_CAGE_PE100nt_strand/mouse_adult_wt_CAGE_PE100nt_strand.plus.bedGraph  | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; for(i=1;i<=length(b);i++) {if(b[i]>max) {max=b[i]; ii=1; delete imax; imax[ii++]=a[i];} else if(b[i]==max) {imax[ii++]=a[i];}} for(j=1;j<ii;j++) print $1,$2,$3,imax[j],max,"+";}' > $CAGE
        cat ~/scratch/mouse_adult_wt_CAGE_PE100nt_strand/mouse_adult_wt_CAGE_PE100nt_strand.minus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; for(i=1;i<=length(b);i++) {if(b[i]<max) {max=b[i]; ii=1; delete imax; imax[ii++]=a[i];} else if(b[i]==max) {imax[ii++]=a[i];}} for(j=1;j<ii;j++) print $1,$2,$3,imax[j],-max,"-";}' >> $CAGE
    }

    cat $INPUTBED | awk -v UP=$UP -v DOWN=$DOWN '{OFS="\t"; s=($6=="+")?($2-UP):((($3-DOWN)>$2)?($3-DOWN):$2); e=($6=="+")?((($2+DOWN)<$3)?($2+DOWN):$3):($3+UP); print $1,(s>0)?s:0,e,$4,$5,$6,$0;}' | intersectBed -a stdin -b $CAGE -wao -s | awk '{OFS="\t"; if($19==".") {$19=$1;$20=$2;$21=$3;$22=($6=="+")?($7":"$8"-"($8+1)):($7":"($9-1)"-"$9);$23=1;$24=$6;} split($22,a,":|-"); $20=a[2]; $21=a[3]; d=($6=="+")?($8-$20):($21-$9); $25=d;if(d<0)d=-d; print $0, $23;}' | sort -k4,4 -k26,26gr -k25,25n | awk 'BEGIN{id=""}{if($10!=id) {print;id=$10;} }' | awk '{OFS="\t"; for(i=19;i<25;i++)printf "%s\t", $i; for(i=1;i<19;i++) printf "%s\t", $i; print $25, $26;}' | $HOME/projects/piRNA/src/changeTSS | cut -f1-12 > $OUT

#debug
#cat $OUT

    UP=$5
    DOWN=$6
    D=$7

    PAS=$HOME/scratch/ucsc/mouse_adult_wt_PAS200_100nt.allpeak.Dis$D.bed
    [ -e $PAS ] || {
        echo "cluster PAS with distance cutoff $D"
        cat ~/scratch/mouse_adult_wt_PAS200_100nt/mouse_adult_wt_PAS200_100nt.plus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; for(i=1;i<=length(b);i++) {if(b[i]>max) {max=b[i]; ii=1; delete imax; imax[ii++]=a[i];} else if(b[i]==max) {imax[ii++]=a[i];}} for(j=1;j<ii;j++) print $1,$2,$3,imax[j],max,"+";}' > $PAS
        cat ~/scratch/mouse_adult_wt_PAS200_100nt/mouse_adult_wt_PAS200_100nt.minus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; for(i=1;i<=length(b);i++) {if(b[i]<max) {max=b[i]; ii=1; delete imax; imax[ii++]=a[i];} else if(b[i]==max) {imax[ii++]=a[i];}} for(j=1;j<ii;j++) print $1,$2,$3,imax[j],-max,"-";}' >> $PAS

    }

    cat $OUT | awk -v UP=$DOWN -v DOWN=$UP '{OFS="\t"; s=($6=="+")?((($3-UP)>$2)?($3-UP):$2):($2-DOWN); e=($6=="+")?($3+DOWN):((($2+UP)<$3)?($2+UP):$3); print $1,(s>0)?s:0,e,$4,$5,$6,$0;}' | intersectBed -a stdin -b $PAS -wao -s | awk '{OFS="\t"; if($19==".") {$19=$1;$20=$2;$21=$3;$22=($6=="-")?($7":"$8"-"($8+1)):($7":"($9-1)"-"$9);$23=1;$24=$6;} split($22,a,":|-"); $20=a[2]; $21=a[3]; d=($6=="-")?($8-$20):($21-$9); $25=d; if(d<0)d=-d;print $0, $23;}' | sort -k4,4 -k26,26gr -k25,25n | awk 'BEGIN{id=""}{if($10!=id) {print;id=$10;} }' | awk '{for(i=19;i<25;i++)printf "%s\t", $i; for(i=1;i<19;i++) printf "%s\t", $i; print $25, $26;}' | $HOME/projects/piRNA/src/changeTTS | cut -f1-12
}



correct53byCAGEPAS2 ()
{
    INPUTBED=$1
    UP=$2
    DOWN=$3
    D=$4

    OUT="$(mktemp)"

    CAGE=$HOME/scratch/ucsc/mouse_adult_wt_CAGE_PE100nt_strand.peak.Dis$D.bed
    [ -e $CAGE ] || {
        echo "cluster CAGE with distance cutoff $D"
        cat ~/scratch/mouse_adult_wt_CAGE_PE100nt_strand/mouse_adult_wt_CAGE_PE100nt_strand.unique.plus.bedGraph  | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; for(i=1;i<=length(b);i++) {if(b[i]>max) {max=b[i]; ii=1; delete imax; imax[ii++]=a[i];} else if(b[i]==max) {imax[ii++]=a[i];}} for(j=1;j<ii;j++) print $1,$2,$3,imax[j],max,"+";}' > $CAGE
        cat ~/scratch/mouse_adult_wt_CAGE_PE100nt_strand/mouse_adult_wt_CAGE_PE100nt_strand.unique.minus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; for(i=1;i<=length(b);i++) {if(b[i]<max) {max=b[i]; ii=1; delete imax; imax[ii++]=a[i];} else if(b[i]==max) {imax[ii++]=a[i];}} for(j=1;j<ii;j++) print $1,$2,$3,imax[j],-max,"-";}' >> $CAGE
    }

    cat $INPUTBED | awk -v UP=$UP -v DOWN=$DOWN '{OFS="\t"; s=($6=="+")?($2-UP):((($3-DOWN)>$2)?($3-DOWN):$2); e=($6=="+")?((($2+DOWN)<$3)?($2+DOWN):$3):($3+UP); print $1,(s>0)?s:0,e,$4,$5,$6,$0;}' | intersectBed -a stdin -b $CAGE -wao -s | awk '{OFS="\t"; if($19==".") {$19=$1;$20=$2;$21=$3;$22=($6=="+")?($7":"$8"-"($8+1)):($7":"($9-1)"-"$9);$23=1;$24=$6;} split($22,a,":|-"); $20=a[2]; $21=a[3]; d=($6=="+")?($8-$20):($21-$9); $25=d;if(d<0)d=-d; print $0, $23;}' | sort -k4,4 -k26,26gr -k25,25n | awk 'BEGIN{id=""}{if($10!=id) {print;id=$10;} }' | awk '{OFS="\t"; for(i=19;i<25;i++)printf "%s\t", $i; for(i=1;i<19;i++) printf "%s\t", $i; print $25, $26;}' | $HOME/projects/piRNA/src/changeTSS | cut -f1-12 > $OUT

#debug
#cat $OUT

    UP=$5
    DOWN=$6
    D=$7

    PAS=$HOME/scratch/ucsc/mouse_adult_wt_PAS200_100nt.peak.Dis$D.bed
    [ -e $PAS ] || {
        echo "cluster PAS with distance cutoff $D"
        cat ~/scratch/mouse_adult_wt_PAS200_100nt/mouse_adult_wt_PAS200_100nt.unique.plus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; for(i=1;i<=length(b);i++) {if(b[i]>max) {max=b[i]; ii=1; delete imax; imax[ii++]=a[i];} else if(b[i]==max) {imax[ii++]=a[i];}} for(j=1;j<ii;j++) print $1,$2,$3,imax[j],max,"+";}' > $PAS
        cat ~/scratch/mouse_adult_wt_PAS200_100nt/mouse_adult_wt_PAS200_100nt.unique.minus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; for(i=1;i<=length(b);i++) {if(b[i]<max) {max=b[i]; ii=1; delete imax; imax[ii++]=a[i];} else if(b[i]==max) {imax[ii++]=a[i];}} for(j=1;j<ii;j++) print $1,$2,$3,imax[j],-max,"-";}' >> $PAS

    }

    cat $OUT | awk -v UP=$DOWN -v DOWN=$UP '{OFS="\t"; s=($6=="+")?((($3-UP)>$2)?($3-UP):$2):($2-DOWN); e=($6=="+")?($3+DOWN):((($2+UP)<$3)?($2+UP):$3); print $1,(s>0)?s:0,e,$4,$5,$6,$0;}' | intersectBed -a stdin -b $PAS -wao -s | awk '{OFS="\t"; if($19==".") {$19=$1;$20=$2;$21=$3;$22=($6=="-")?($7":"$8"-"($8+1)):($7":"($9-1)"-"$9);$23=1;$24=$6;} split($22,a,":|-"); $20=a[2]; $21=a[3]; d=($6=="-")?($8-$20):($21-$9); $25=d; if(d<0)d=-d;print $0, $23;}' | sort -k4,4 -k26,26gr -k25,25n | awk 'BEGIN{id=""}{if($10!=id) {print;id=$10;} }' | awk '{for(i=19;i<25;i++)printf "%s\t", $i; for(i=1;i<19;i++) printf "%s\t", $i; print $25, $26;}' | $HOME/projects/piRNA/src/changeTTS | cut -f1-12
}


correct53byCAGEPAS ()
{
    INPUTBED=$1
    UP=$2
    DOWN=$3
    D=$4

    OUT="$(mktemp)"

    CAGE=$HOME/scratch/ucsc/mouse_adult_wt_CAGE_PE100nt_strand.peak.Dis$D.bed
    [ -e $CAGE ] || {
        echo "cluster CAGE with distance cutoff $D"
        cat ~/scratch/mouse_adult_wt_CAGE_PE100nt_strand/mouse_adult_wt_CAGE_PE100nt_strand.unique.plus.bedGraph  | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; for(i=1;i<=length(b);i++) {if(b[i]>max) {max=b[i]; ii=1; delete imax; imax[ii++]=a[i];} else if(b[i]==max) {imax[ii++]=a[i];}} for(j=1;j<ii;j++) print $1,$2,$3,imax[j],max,"+";}' > $CAGE
        cat ~/scratch/mouse_adult_wt_CAGE_PE100nt_strand/mouse_adult_wt_CAGE_PE100nt_strand.unique.minus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; for(i=1;i<=length(b);i++) {if(b[i]<max) {max=b[i]; ii=1; delete imax; imax[ii++]=a[i];} else if(b[i]==max) {imax[ii++]=a[i];}} for(j=1;j<ii;j++) print $1,$2,$3,imax[j],-max,"-";}' >> $CAGE
    }

    cat $INPUTBED | awk -v UP=$UP -v DOWN=$DOWN '{OFS="\t"; s=($6=="+")?($2-UP):((($3-DOWN)>$2)?($3-DOWN):$2); e=($6=="+")?((($2+DOWN)<$3)?($2+DOWN):$3):($3+UP); print $1,(s>0)?s:0,e,$4,$5,$6,$0;}' | intersectBed -a stdin -b $CAGE -wao -s | awk '{OFS="\t"; if($19==".") {$19=$1;$20=$2;$21=$3;$22=($6=="+")?($7":"$8"-"($8+1)):($7":"($9-1)"-"$9);$23=1;$24=$6;} split($22,a,":|-"); $20=a[2]; $21=a[3]; d=($6=="+")?($8-$20):($21-$9); $25=d;if(d<0)d=-d; print $0, $23/(d+1)/1000;}' | sort -k4,4 -k26,26gr -k25,25n | awk 'BEGIN{id=""}{if($10!=id) {print;id=$10;} }' | awk '{OFS="\t"; for(i=19;i<25;i++)printf "%s\t", $i; for(i=1;i<19;i++) printf "%s\t", $i; print $25, $26;}' | $HOME/projects/piRNA/src/changeTSS | cut -f1-12 > $OUT

#debug
#cat $OUT

    UP=$5
    DOWN=$6
    D=$7

    PAS=$HOME/scratch/ucsc/mouse_adult_wt_PAS200_100nt.peak.Dis$D.bed
    [ -e $PAS ] || {
        echo "cluster PAS with distance cutoff $D"
        cat ~/scratch/mouse_adult_wt_PAS200_100nt/mouse_adult_wt_PAS200_100nt.unique.plus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; for(i=1;i<=length(b);i++) {if(b[i]>max) {max=b[i]; ii=1; delete imax; imax[ii++]=a[i];} else if(b[i]==max) {imax[ii++]=a[i];}} for(j=1;j<ii;j++) print $1,$2,$3,imax[j],max,"+";}' > $PAS
        cat ~/scratch/mouse_adult_wt_PAS200_100nt/mouse_adult_wt_PAS200_100nt.unique.minus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; for(i=1;i<=length(b);i++) {if(b[i]<max) {max=b[i]; ii=1; delete imax; imax[ii++]=a[i];} else if(b[i]==max) {imax[ii++]=a[i];}} for(j=1;j<ii;j++) print $1,$2,$3,imax[j],-max,"-";}' >> $PAS

    }

    cat $OUT | awk -v UP=$DOWN -v DOWN=$UP '{OFS="\t"; s=($6=="+")?((($3-UP)>$2)?($3-UP):$2):($2-DOWN); e=($6=="+")?($3+DOWN):((($2+UP)<$3)?($2+UP):$3); print $1,(s>0)?s:0,e,$4,$5,$6,$0;}' | intersectBed -a stdin -b $PAS -wao -s | awk '{OFS="\t"; if($19==".") {$19=$1;$20=$2;$21=$3;$22=($6=="-")?($7":"$8"-"($8+1)):($7":"($9-1)"-"$9);$23=1;$24=$6;} split($22,a,":|-"); $20=a[2]; $21=a[3]; d=($6=="-")?($8-$20):($21-$9); $25=d; if(d<0)d=-d;print $0, $23/(d+1)/1000;}' | sort -k4,4 -k26,26gr -k25,25n | awk 'BEGIN{id=""}{if($10!=id) {print;id=$10;} }' | awk '{for(i=19;i<25;i++)printf "%s\t", $i; for(i=1;i<19;i++) printf "%s\t", $i; print $25, $26;}' | $HOME/projects/piRNA/src/changeTTS | cut -f1-12

}

correct5byCAGE ()
{
    INPUTBED=$1
    UP=$2
    DOWN=$3
    D=$4

    CAGE=$HOME/scratch/ucsc/mouse_adult_wt_CAGE_PE100nt_strand.peak.Dis$D.bed
    [ -e $CAGE ] || {
        echo "cluster CAGE with distance cutoff $D"
        cat ~/scratch/mouse_adult_wt_CAGE_PE100nt_strand/mouse_adult_wt_CAGE_PE100nt_strand.unique.plus.bedGraph  | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; for(i=1;i<=length(b);i++) {if(b[i]>max) {max=b[i]; ii=1; delete imax; imax[ii++]=a[i];} else if(b[i]==max) {imax[ii++]=a[i];}} for(j=1;j<ii;j++) print $1,$2,$3,imax[j],max,"+";}' > $CAGE
        cat ~/scratch/mouse_adult_wt_CAGE_PE100nt_strand/mouse_adult_wt_CAGE_PE100nt_strand.unique.minus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; for(i=1;i<=length(b);i++) {if(b[i]<max) {max=b[i]; ii=1; delete imax; imax[ii++]=a[i];} else if(b[i]==max) {imax[ii++]=a[i];}} for(j=1;j<ii;j++) print $1,$2,$3,imax[j],-max,"-";}' >> $CAGE
    }

    cat $INPUTBED | awk -v UP=$UP -v DOWN=$DOWN '{OFS="\t"; s=($6=="+")?($2-UP):((($3-DOWN)>$2)?($3-DOWN):$2); e=($6=="+")?((($2+DOWN)<$3)?($2+DOWN):$3):($3+UP); print $1,(s>0)?s:0,e,$4,$5,$6,$0;}' | intersectBed -a stdin -b $CAGE -wao -s | awk '{OFS="\t"; if($19==".") {$19=$1;$20=$2;$21=$3;$22=($6=="+")?($7":"$8"-"($8+1)):($7":"($9-1)"-"$9);$23=1;$24=$6;} split($22,a,":|-"); $20=a[2]; $21=a[3]; d=($6=="+")?($8-$20):($21-$9); $25=d;if(d<0)d=-d; print $0, $23/(d+1)/1000;}' | sort -k4,4 -k26,26gr -k25,25n | awk 'BEGIN{id=""}{if($10!=id) {print;id=$10;} }' | awk '{OFS="\t"; for(i=19;i<25;i++)printf "%s\t", $i; for(i=1;i<19;i++) printf "%s\t", $i; print $25, $26;}' | $HOME/projects/piRNA/src/changeTSS | cut -f1-12
}

correct3byPAS ()
{
    INPUTBED=$1
    D=$4
    UP=$2
    DOWN=$3

    PAS=$HOME/scratch/ucsc/mouse_adult_wt_PAS200_100nt.peak.Dis$D.bed
    [ -e $PAS ] || {
        echo "cluster PAS with distance cutoff $D"
        cat ~/scratch/mouse_adult_wt_PAS200_100nt/mouse_adult_wt_PAS200_100nt.unique.plus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; for(i=1;i<=length(b);i++) {if(b[i]>max) {max=b[i]; ii=1; delete imax; imax[ii++]=a[i];} else if(b[i]==max) {imax[ii++]=a[i];}} for(j=1;j<ii;j++) print $1,$2,$3,imax[j],max,"+";}' > $PAS
        cat ~/scratch/mouse_adult_wt_PAS200_100nt/mouse_adult_wt_PAS200_100nt.unique.minus.bedGraph | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3, $4 }' | mergeBed -d $D -scores collapse -nms | awk '{OFS="\t"; split($4,a,";");split($5,b,","); max=0; for(i=1;i<=length(b);i++) {if(b[i]<max) {max=b[i]; ii=1; delete imax; imax[ii++]=a[i];} else if(b[i]==max) {imax[ii++]=a[i];}} for(j=1;j<ii;j++) print $1,$2,$3,imax[j],-max,"-";}' >> $PAS

    }

    cat $INPUTBED | awk -v UP=$UP -v DOWN=$DOWN '{OFS="\t"; s=($6=="+")?((($3-UP)>$2)?($3-UP):$2):($2-DOWN); e=($6=="+")?($3+DOWN):((($2+UP)<$3)?($2+UP):$3); print $1,(s>0)?s:0,e,$4,$5,$6,$0;}' | intersectBed -a stdin -b $PAS -wao -s | awk '{OFS="\t"; if($19==".") {$19=$1;$20=$2;$21=$3;$22=($6=="-")?($7":"$8"-"($8+1)):($7":"($9-1)"-"$9);$23=1;$24=$6;} split($22,a,":|-"); $20=a[2]; $21=a[3]; d=($6=="-")?($8-$20):($21-$9); $25=d; if(d<0)d=-d;print $0, $23/(d+1)/1000;}' | sort -k4,4 -k26,26gr -k25,25n | awk 'BEGIN{id=""}{if($10!=id) {print;id=$10;} }' | awk '{for(i=19;i<25;i++)printf "%s\t", $i; for(i=1;i<19;i++) printf "%s\t", $i; print $25, $26;}' | $HOME/projects/piRNA/src/changeTTS | cut -f1-12

}
