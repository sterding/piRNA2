#!/bin/sh

# Instead of
# zcat gencode.v7.annotation.gtf.gz | awk 'BEGIN{OFS="\t"}{if($3=="transcript")print $1,$4,$5,$7,$10,$12,$14,$16,$18,$20,$22,$24,$26}' | sed 's/;//g;s/"//g;' | perl ${HOME}/projects/encode/src/get_normalized.CpG.content.pl > gencode_v7_hg19_transcripts.gtf.cpg.tab2 &
# split genome into each chromosome

for chr in {10..19} X Y
do
    grep -P "^chr$chr\t"  $GENOME/mm9/Annotation/Genes/NCBIM37.biomart67.transcripts.tab | perl get_normalized.CpG.content.pl > ../data/NCBIM37.biomart67.transcripts.tab.chr$chr &
done
