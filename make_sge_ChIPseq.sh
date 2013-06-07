## interactive-version of pipeline for running ChIP-seq mapping, peak calling etc.
## Author: Xianjun Dong (xianjun.dong@umassmed.edu)
## Date: 2012-March-13
## Version: 0.1

#!/bin/sh

#1.Prepare input/parameters for the script

read -p "Pair-end? [No|yes]" pe
pe=${pe:-no}
if [[ "$pe" =~  "YES|Yes|yes|Y|y" ]]; then
    read -p "Please enter the expected inner distance between mate pairs: [300]" mate_inner_dist
    pe="yes"

    read -p "Please enter the input directory: [/home/dongx/nearline/Xin/others/]" inputdir
    inputdir=${inputdir:-/home/dongx/nearline/Xin/others/}
    echo $inputdir
    if test ! -d $inputdir
        then
        echo "Error: $inputdir does not exist!!"
        exit 1
    fi

    cd $inputdir; ls -ogh

    read -p "Please enter the wildcard for left reads files (e.g. 500_mouse_NoIndex_L002_R1*): " wildcard
    wildcard=${wildcard:-"*_R1*"}
    echo $(ls $wildcard) | sed 's/\s/,/g'
    read -p "Is this correct? [Yes|no]" answer
    if [[ "$answer" =~  "YES|Yes|yes|Y|y" ]]; then
        readsfile=$(echo $(ls $wildcard) | sed 's/\s/,/g')
    else
        read -p "Please enter the left reads files (seperated by comma): " readsfile
    fi

    read -p "Please enter the wildcard for right reads files (e.g. 500_mouse_NoIndex_L002_R2*): " wildcard
    wildcard=${wildcard:-"*_R2*"}
    echo $(ls $wildcard) | sed 's/\s/,/g'
    read -p "Is this correct? [Yes|no]" answer
    if [[ "$answer" =~  "YES|Yes|yes|Y|y" ]]; then
        readsfile2=$(echo $(ls $wildcard) | sed 's/\s/,/g')
    else
        read -p "Please enter the left reads files (seperated by comma): " readsfile2
    fi

    #readsfile="-1 $readsfile -2 $readsfile2"

fi
if [[ "$pe" =~  "NO|No|no|N|n" ]]; then
    pe="no"

    read -p "Please enter the input directory: [/home/dongx/nearline/Xin/others/]" inputdir
    inputdir=${inputdir:-/home/dongx/nearline/Xin/others/}
    echo $inputdir
    if test ! -d $inputdir
        then
        echo "Error: $inputdir does not exist!!"
        exit 1
    fi

    cd $inputdir; ls -ogh

    read -p "Please enter the wildcard for right reads files (e.g. *.fastq.gz): " wildcard
    wildcard=${wildcard:-"*.fastq.gz"}
    echo $(ls $wildcard) | sed 's/\s/,/g'
    read -p "Is this correct? [Yes|no]" answer
    answer=${answer:-Y}
    if [[ "$answer" =~  "YES|Yes|yes|Y|y" ]]; then
        readsfile=$(echo $(ls $wildcard) | sed 's/\s/,/g')
    else
        read -p "Please enter the left reads files (seperated by comma): " readsfile
    fi

    #readsfile="-U $readsfile"

fi

echo "Please enter the sample name:[" ${readsfile%.*} "]"
read samplename
samplename=${samplename:-${readsfile%.*}}

read -p "please enter your output dir: [/home/dongx/scratch]" dir
dir=${dir:-/home/dongx/scratch}

if test ! -d $dir/$samplename/fastqc
then
mkdir -p $dir/$samplename/fastqc
fi

read -p "Map to the genome? [mm9]" index
index=${index:-mm9}

if test ! -d "$GENOME/$index/Sequence/Bowtie2Index/"
then
    echo "Error: $GENOME/$index/Sequence/Bowtie2Index/ does not exist!!"
    echo "Please enter the Index directory (full path)"
    read BOWTIE_INDEXES
    if test ! -d $BOWTIE_INDEXES
    then
        echo "Error: $BOWTIE_INDEXES does not exist!!"
        exit 1
    fi
fi

read -p "Predscore encoding? [1:Illumina 1.8+/Sanger; 2:Illumina 1.5+; 3:Solexa]" predscore
predscore=${predscore:-1}
case "$predscore" in
    1)
        bowtie="--phred33-quals"; bowtie2="--phred33"; tophat=""; far="fastq-sanger";
        ;;
    2)
        bowtie="--phred64-quals"; bowtie2="--phred64"; tophat="--solexa1.3-quals"; far="fastq-illumina15";
        ;;
    3)
        bowtie="--solexa-quals"; bowtie2="--solexa-quals"; tophat="--solexa-quals"; far="fastq-illumina13";
        ;;
    *)
        echo "Invalid choice!"
        exit 1
        ;;
esac

read -p "Adaptor? [No|yes]" adaptor
adaptor=${adaptor:-No}
if [[ "$adaptor" =~  "YES|Yes|yes|Y|y" ]]; then
    adaptor="yes"
    read -p "Adaptor at 3' or 5'? [3|5]" trimend
    trimend=${trimend:-3}
    if [[ "$trimend" == "3" ]]; then
        trimend='right'
    else
        trimend='left'
    fi

    read -p "Input the adaptor sequence:[TCGTATGCCGTCTTCTGCTTG]" adaptorsequence
    adaptorsequence=${adaptorsequence:-TCGTATGCCGTCTTCTGCTTG}
else
    adaptor="no"
fi

read -p "Strand-specific? [No|yes]" strand
strand=${strand:-no}
if [[ "$strand" =~  "YES|Yes|yes|Y|y" ]]; then
    strand="yes"
    strandoption=""
fi
if [[ "$strand" =~  "NO|No|no|N|n" ]]; then
    strand="no"
    strandoption="--library-type fr-unstranded"
fi
echo "strand:" $strand

read -p "Your email address for notice: [sterding.hpcc@gmail.com]" EMAIL
EMAIL=${EMAIL:-"sterding.hpcc@gmail.com"}

echo "
## bash script for running chipseq pipeline on cluster
#!/bin/sh
#$ -V
#$ -pe single 4
#$ -cwd
#$ -o \$HOME/sge_jobs_output/sge_job.\$JOB_ID.out -j y
#$ -S /bin/bash
#$ -M $EMAIL
#$ -S /bin/bash
#$ -m e
#$ -l mem_free=2G

###########################################
############## 1. Configuring
###########################################

export BOWTIE_INDEXES=\$GENOME/$index/Sequence/BowtieIndex/
export BOWTIE2_INDEXES=\$GENOME/$index/Sequence/Bowtie2Index/
export ANNOTATION=\$GENOME/$index/Annotation/Genes
export SEQUENCE=\$GENOME/$index/Sequence/WholeGenomeFasta

cd $inputdir

###########################################
###############  2. quality filter: adaptor removal/clip
###########################################

if [[ \"$adaptor\" == \"yes\" ]]; then
    far -s $readsfile -t $readsfile.clipped -f $far -as $adaptorsequence --cut-off 5 --min-overlap 10  --min-readlength 20 --trim-end $trimend --adaptive-overlap yes --nr-threads 8 --max-uncalled 30
    mv $readsfile $readsfile.orgin
    mv $readsfile.clipped $readsfile
fi
#
#############################################
################ 3. QC
############################################
#
fastqc --outdir $dir/$samplename $(echo $readsfile | sed 's/,/ /g')
rm $dir/$samplename/*fastqc.zip
#
#
#############################################
################# 4. mapping to the genome
#############################################
#
## using bowtie2
export BOWTIE2_INDEXES=\$GENOME/$index/Sequence/Bowtie2Index/
if [[ \"$pe\" == \"yes\" ]]; then
    bowtie2 -x genome $bowtie2 -p 8 -1 $readsfile1 -2 $readsfile2 -k 10 -S $dir/$samplename/accepted_hits.sam 2>$dir/$samplename/$samplename.mapping.summary
else
    bowtie2 -x genome $bowtie2 -p 8 -U $readsfile -k 10 -S $dir/$samplename/accepted_hits.sam 2>$dir/$samplename/$samplename.mapping.summary
fi

###########################################
############### 5. post-processing, format converting
###########################################

# sam - bam - sorted - index

cd $dir/$samplename

samtools view -Sbut \$BOWTIE2_INDEXES/genome.fai accepted_hits.sam | samtools sort - accepted_hits.sorted
mv accepted_hits.sorted.bam accepted_hits.bam
samtools index accepted_hits.bam

###########################################
############## 6. prepare for tracks files to display on UCSC
###########################################

[ -d $dir/ucsc ] || mkdir $dir/ucsc

# make index for the (sorted) BAM
cp accepted_hits.bam $dir/ucsc/$samplename.accepted_hits.bam
cp accepted_hits.bam.bai $dir/ucsc/$samplename.accepted_hits.bam.bai

# QC
mv $dir/$samplename/*_fastqc $dir/ucsc

cd $dir/ucsc

# bigwig
bamToBed -i $samplename.accepted_hits.bam -split > $samplename.accepted_hits.bed
if [[ \"$strand\" == \"no\" ]]; then
    sort -k1,1 $samplename.accepted_hits.bed | bedItemOverlapCount $index -chromSize=\$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n > $samplename.accepted_hits.bedGraph
    bedGraphToBigWig $samplename.accepted_hits.bedGraph \$ANNOTATION/ChromInfo.txt $samplename.accepted_hits.bw
fi
if [[ \"$strand\" == \"yes\" ]]; then
    awk '{if(\$6==\"+\") print}' $samplename.accepted_hits.bed | sort -k1,1 | bedItemOverlapCount $index -chromSize=\$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n > $samplename.accepted_hits.+.bedGraph
    awk '{if(\$6==\"-\") print}' $samplename.accepted_hits.bed | sort -k1,1 | bedItemOverlapCount $index -chromSize=\$ANNOTATION/ChromInfo.txt stdin | sort -k1,1 -k2,2n | awk '{OFS=\"\t\"; print \$1,\$2,\$3,\"-\"\$4}' > $samplename.accepted_hits.-.bedGraph
    bedGraphToBigWig $samplename.accepted_hits.+.bedGraph \$ANNOTATION/ChromInfo.txt $samplename.accepted_hits.+.bw
    bedGraphToBigWig $samplename.accepted_hits.-.bedGraph \$ANNOTATION/ChromInfo.txt $samplename.accepted_hits.-.bw
fi
#rm $samplename.accepted_hits.bed $samplename.accepted_hits*.bedGraph

">$dir/$samplename/chiseq_$samplename.sge

cat $dir/$samplename/chiseq_$samplename.sge
#qsub $dir/$samplename/chiseq_$samplename.sge
