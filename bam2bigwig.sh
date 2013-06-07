#!/bin/bash

# script convert BAM to bigwig and bedgraph
# require bedtools and Jim Kent's utility in $PATH

# current usage: bam2bigwig <in.bam> <in.bw>
# optimal usage: bam2bigwig [-sB] -I <mm9> -i <in.bam|sam> -o <out.bw>
# -B: output bedGraph (instead of bigWig)
# -I: index for annnotated genome (e.g. hg19, mm8 etc.)
# -i: input file
# -o: output file

## Log:
#1. change to use "bamToBed -bed12" and "bedItemOverlapCount -bed12", instead of "bamToBed -split"

ANNOTATION=$GENOME/mm9/Annotation/Genes

inputfile=$1

[ -e "$1" ] || {
	echo "
# Script convert BAM to bigwig and bedgraph
# Usage: bam2bw [-split|-nosplit] <in.bam|sam|bed>
# -split is to split input file in +/- strand
# -nosplit is to not split input file by strand
# Right now, it only works for mouse mm9 genome."; exit 0;}

ext=${inputfile##*.}
bname=${inputfile%.*}

case $ext in
    bam|BAM|Bam)
        echo "Input is a BAM file. Converting bam-->sam-->bed ..."
        # use the XS:A for strand
        samtools view $inputfile | sam2bed -v bed12=T -v sCol=NH > $bname.bed;;
        # use the FLAG for strand
        #bamToBed -i $inputfile -bed12 > $bname.bed;;
    sam|SAM|Sam)
        echo "Input is a SAM file. Converting sam->bed ..."
        sam2bed -v bed12=T -v sCol=NH $inputfile > $bname.bed;;
    bed|BED|Bed|bed6|bed12)
        echo "Input is a BED file.";;
    *)
        echo "Unsupported input format. Exit!"
        exit;;
esac

if [ "$2" == "" ] || [ "$2" == "-split" ]
then
    echo "bed-->bw, by strand..."
    echo "bedToBedGraph..."
    > $bname.bed+; > $bname.bed-;
    awk -v OUTPUT=$bname '{print $0 >> OUTPUT".bed"$6}' $bname.bed
    ncol=`head $bname.bed | awk 'END{print NF}'`
    if [ "$ncol" -eq  12 ]
    then
        echo "in bed12 format ..."
        sort -k1,1 $bname.bed+ | bedItemOverlapCount mm9 -bed12 -chromSize=$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n > $bname.plus.bedGraph
        sort -k1,1 $bname.bed- | bedItemOverlapCount mm9 -bed12 -chromSize=$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n | awk '{OFS="\t"; print $1,$2,$3,"-"$4}' > $bname.minus.bedGraph
    else
        echo "in bed6 format ..."
        sort -k1,1 $bname.bed+ | bedItemOverlapCount mm9 -chromSize=$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n > $bname.plus.bedGraph
        sort -k1,1 $bname.bed- | bedItemOverlapCount mm9 -chromSize=$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n | awk '{OFS="\t"; print $1,$2,$3,"-"$4}' > $bname.minus.bedGraph
    fi
    echo "bedGraph2bw..."
    bedGraphToBigWig $bname.plus.bedGraph $ANNOTATION/ChromInfo.txt $bname.plus.bw
    bedGraphToBigWig $bname.minus.bedGraph $ANNOTATION/ChromInfo.txt $bname.minus.bw
fi

if [ "$2" == "-nosplit" ]
then
    echo "bed-->bw, without spliting by strand..."
    echo "bedGraph2bw..."
    ncol=`head $bname.bed | awk 'END{print NF}'`
    if [ "$ncol" -eq  12 ]
    then
        echo "in bed12 format ..."
        sort -k1,1 $bname.bed | bedItemOverlapCount mm9 -bed12 -chromSize=$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n > $bname.bedGraph
    else
        echo "in bed6 format ..."
        sort -k1,1 $bname.bed | bedItemOverlapCount mm9 -chromSize=$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n > $bname.bedGraph
    fi
    echo "bedGraph2bw..."
    bedGraphToBigWig $bname.bedGraph $ANNOTATION/ChromInfo.txt $bname.bw
fi
# remove temp files
#rm ${bamfile/%bam/bed} ${bamfile/%bam/*}bedGraph
