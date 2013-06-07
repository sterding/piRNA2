## interactive-version of pipeline for running RNAseq mapping/assmelby on cluster
## Author: Xianjun Dong (xianjun.dong@umassmed.edu)
## Usage: ./make_sge_RNAseq.sh
## Date: 2012-March-13
## Version: 0.1

#!/bin/sh
function die(){ echo $1; exit 1; }

#1.Prepare input/parameters for the script
read -p "Please enter the input directory: [/home/dongx/nearline/Xin/RNAseq]" inputdir
inputdir=${inputdir:-/home/dongx/nearline/Xin/RNAseq}
echo $inputdir
[ -d $inputdir ] || die "Error: $inputdir does not exist!!"

cd $inputdir; ls -ogh

read -p "Please enter the sample name:" samplename
[ -d $samplename ] || die "Error: $samplename does not exist!!"

cd $samplename; ls -ogh

read -p "Pair-end? [Yes|no]" pe
pe=${pe:-yes}
if [[ "$pe" =~  "YES|Yes|yes|Y|y" ]]; then
    pe="yes"
    read -p "Please enter the expected inner distance between mate pairs: [300]" mate_inner_dist
    mate_inner_dist=${mate_inner_dist:-300}
    PE="--mate-inner-dist $mate_inner_dist"

    read -p "Please enter the wildcard for left reads files (e.g. R1.f*q): " wildcard
    wildcard=${wildcard:-"R1.f*q"}
    echo $(ls $wildcard) | sed 's/\s/,/g'
    read -p "Is this correct? [Yes|no]" answer
    answer=${answer:-Yes}
    if [[ "$answer" =~  "YES|Yes|yes|Y|y" ]]; then
        readsfile1=$(echo $(ls $wildcard) | sed 's/\s/,/g')
    else
        read -p "Please enter the left reads files (seperated by comma): " readsfile1
    fi

    read -p "Please enter the wildcard for right reads files (e.g. R2.f*q): " wildcard
    wildcard=${wildcard:-"R2.f*q"}
    echo $(ls $wildcard) | sed 's/\s/,/g'
    read -p "Is this correct? [Yes|no]" answer
    answer=${answer:-Yes}
    if [[ "$answer" =~  "YES|Yes|yes|Y|y" ]]; then
        readsfile2=$(echo $(ls $wildcard) | sed 's/\s/,/g')
    else
        read -p "Please enter the left reads files (seperated by comma): " readsfile2
    fi
fi
if [[ "$pe" =~  "NO|No|no|N|n" ]]; then
    pe="no"
    PE=""
    read -p "Please enter the wildcard for right reads files (e.g. R1.f*q): " wildcard
    wildcard=${wildcard:-"R1.f*q"}
    echo $(ls $wildcard) | sed 's/\s/,/g'
    read -p "Is this correct? [Yes|no]" answer
    answer=${answer:-Yes}
    if [[ "$answer" =~  "YES|Yes|yes|Y|y" ]]; then
        readsfile1=$(echo $(ls $wildcard) | sed 's/\s/,/g')
    else
        read -p "Please enter the left reads files (seperated by comma): " readsfile1
    fi

fi

read -p "Map to the genome? [mm9]" index
index=${index:-mm9}

if test ! -d "$GENOME/$index/Sequence/Bowtie2Index/"
then
    echo "Error: $GENOME/$index/Sequence/Bowtie2Index/ does not exist!!"
    echo "Please enter the Index directory (full path)"
    read BOWTIE2_INDEXES
    [ -d $BOWTIE2_INDEXES ] || die "Error: $BOWTIE2_INDEXES does not exist!!"
fi

read -p "Mismatch? [2|1|0]" mm
mm=${mm:-2}

read -p "Predscore encoding? [1:Phred+33; 2:Phred+64; 3:Solexa+64] " predscore
predscore=${predscore:-1}
case "$predscore" in
    1)
        bowtie="--phred33-quals"; bowtie2="--phred33"; tophat=""; far="fastq-sanger"; fastqmcf="33"; trimmomatic="-phred33"
        ;;
    2)
        bowtie="--phred64-quals"; bowtie2="--phred64"; tophat="--solexa1.3-quals"; far="fastq"; fastqmcf="64"; trimmomatic="-phred64"
        ;;
    3)
        bowtie="--solexa-quals"; bowtie2="--solexa-quals"; tophat="--solexa-quals"; far="fastq-illumina13"; fastqmcf="59"; trimmomatic="-phred59"  # not sure if this is right for solexa when using trimmomatic
        ;;
    *)
        die "Invalid choice!"
        ;;
esac

read -p "Strand-specific? [Yes|no]" strand
strand=${strand:-Yes}
if [[ "$strand" =~  "YES|Yes|yes|Y|y" ]]; then
    strand="yes"
    read -p "Library Type? [fr-firststrand]" libtype
    libtype=${libtype:-"fr-firststrand"}
    strandoption="--library-type $libtype"
fi
if [[ "$strand" =~  "NO|No|no|N|n" ]]; then
    strand="no"
    strandoption="--library-type fr-unstranded"
fi
echo "strand:" $strand

read -p "CPU #: [8]" cpu
cpu=${cpu:-8}

read -p "Your email address for notice: [sterding.hpcc@gmail.com]" EMAIL
EMAIL=${EMAIL:-"sterding.hpcc@gmail.com"}

read -p "please enter your output dir: [$HOME/scratch]" outputdir
outputdir=${outputdir:-$HOME/scratch}

[ -d $outputdir/$samplename ] || mkdir -p $outputdir/$samplename

piRNA_cluster_region_bed=$HOME/projects/piRNA/data/piRNA.clusters.coordinates.bed12
smallRNAreads=$HOME/nearline/Xin/smallRNA/Jia/Phil.SRA.wt.ox.6wk.testes.raw.uniq.reads.fa.gz

echo "
###########################################
## bash script for running RNAseq pipeline on cluster
###########################################
#!/bin/sh
#$ -V
#$ -pe single $cpu
#$ -cwd
#$ -o \$HOME/sge_jobs_output/sge_job.\$JOB_ID.out -j y
#$ -M $EMAIL
#$ -S /bin/bash
#$ -m e
#$ -l mem_free=3G

###########################################
############## 1. Configuring
###########################################

export BOWTIE2_INDEXES=\$GENOME/$index/Sequence/Bowtie2Index/
export ANNOTATION=\$GENOME/$index/Annotation/Genes

[ -d $outputdir/$samplename ] || mkdir -p $outputdir/$samplename
stat_file=$outputdir/$samplename/mapping.stat
ln -fs \$HOME/sge_jobs_output/sge_job.\$JOB_ID.out $outputdir/$samplename/sge.log

cd $inputdir/$samplename

# cat individual files together (DEPRECATED: now we prepare the R1/R2.fastq when we download the data)
#zcat `echo $readsfile1 | sed 's/,/ /g'` > R1.fastq
#[[ \"$pe\" == \"yes\" ]] && zcat `echo $readsfile2 | sed 's/,/ /g'` > R2.fastq

### STAT
echo \"Total number of reads (x2 if PE): \"\`wc -l R1.fastq | awk '{printf(\"%'\\''d\",\$1);}'\` > \$stat_file

############################################
################  2. quality filter: adaptor removal/clip
############################################

## adaptor
#adaptorfile=adaptor.fa
#[ -d \$adaptorfile ] || curl -s http://zlab.umassmed.edu/~dongx/tracks/data/adaptors_list.fa > $inputdir/adaptor.fa
#
###### adaptor removal
#[ -d filtered ] || mkdir filtered
#if [[ \"$pe\" == \"yes\" ]]; then
#    fastq-mcf -o filtered/R1.fastq -o filtered/R2.fastq -l 16 -q 15 -w 4 -x 10 -u -P $fastqmcf $adaptorfile R1.fastq R2.fastq
#else
#    fastq-mcf -o filtered/R1.fastq -l 16 -q 15 -w 4 -x 10 -u -P $fastqmcf $adaptorfile R1.fastq
#fi
#cd filtered
#
### STAT
#echo \"After adaptor removal (x2 if PE): \"\`wc -l R1.fastq | awk '{printf(\"%'\\''d\",\$1);}'\` >> \$stat_file

##############################################
################# 3. QC
#############################################
[ -d $outputdir/$samplename/fastqc ] || mkdir -p $outputdir/$samplename/fastqc
fastqc --outdir $outputdir/$samplename --extract --threads 2 -f fastq \`ls R*.fastq\`
#\`ls R*.fastq\`
rm $outputdir/$samplename/*fastqc.zip

#############################################
################# 4. mapping to the genome
#############################################

# tophat (output accepted_hit.sam, allow multiple hits)
tophat -o $outputdir/$samplename --no-convert-bam -p $cpu --read-mismatches $mm $tophat $PE $strandoption --min-anchor-length 8 --min-intron-length 30 --max-intron-length 50000 --splice-mismatches 1 --max-multihits 100 --no-coverage-search genome \`ls R*.fastq\`

cd $outputdir/$samplename

### STAT
cut -f1 accepted_hits.sam | sort | uniq -c > accepted_hits.sam.IDsorted.temp
echo \"Mapped reads: \"\`wc -l accepted_hits.sam.IDsorted.temp | awk '{printf(\"%'\\''d\",\$1);}'\` >> \$stat_file
echo \"  among which, unique mapper and multi-mappers:\" >> \$stat_file
textHistogram -maxBinCount=2 -noStar accepted_hits.sam.IDsorted.temp >> \$stat_file

###########################################
############### 5. post-processing, format converting
###########################################


cd $outputdir/$samplename

## for strand-specific RNAseq data, we need to manually add the XS:A:+/- tag that cufflinks required for SAM alignment
if [[ \"$strand\" == \"yes\" ]]; then
    ## 5.1 check & fix the XS:A tag (Note: _fixXStag.awk is under ~/bin)
    _fixXStag -v libtype=$libtype -v save_discrepancy_to_file=discrepant_reads.sam accepted_hits.sam > _accepted_hits.sam

    ##4.3: SAM->BAM
    samtools view -Sbut \$BOWTIE2_INDEXES/genome.fai _accepted_hits.sam | samtools sort - accepted_hits.sorted
    mv accepted_hits.sorted.bam accepted_hits.bam
    samtools index accepted_hits.bam
else
    ## sam -> bam -> sorted -> index
    samtools view -Sbut \$BOWTIE2_INDEXES/genome.fai accepted_hits.sam | samtools sort - accepted_hits.sorted
    mv accepted_hits.sorted.bam accepted_hits.bam
    samtools index accepted_hits.bam
fi

###########################################
############## 6. prepare for tracks files to display on UCSC
###########################################

[ -d $outputdir/ucsc ] || mkdir $outputdir/ucsc

# make index for the (sorted) BAM
cp accepted_hits.bam $outputdir/ucsc/$samplename.accepted_hits.bam
cp accepted_hits.bam.bai $outputdir/ucsc/$samplename.accepted_hits.bam.bai

# QC
cp $outputdir/$samplename/*_fastqc $outputdir/ucsc

cd $outputdir/ucsc

# bam -> bigwig
samtools view $samplename.accepted_hits.bam | sam2bed -v bed12=F > $outputdir/ucsc/$samplename.accepted_hits.bed

[[ \"$strand\" == \"no\" ]] && {bam2bw $samplename.accepted_hits.bam -nosplit}
[[ \"$strand\" == \"yes\" ]] && {bam2bw $samplename.accepted_hits.bam}
#rm $outputdir/ucsc/$samplename.accepted_hits.bed $outputdir/ucsc/$samplename.accepted_hits*.bedGraph

###########################################
################# 7. assembly
###########################################
[ -d $outputdir/$samplename/cufflinks ] || mkdir -p $outputdir/$samplename/cufflinks
cufflinks -v --no-update-check $strandoption -o $outputdir/$samplename/cufflinks -p $cpu \$ANNOTATION/genes.gtf -M \$ANNOTATION/chrM.rRNA.tRNA.gtf accepted_hits.bam -A 0.05 --min-intron-length 30 -u -j 0.2 --min-frags-per-transfrag 40 --overlap-radius 250 2&> $outputdir/$samplename/cufflinks/run.log

# gtf of assembly
echo \"track name=$samplename description=$samplename visibility=pack colorByStrand='200,100,0 0,100,200'\" > $outputdir/ucsc/$samplename.transcripts.gtf
cat $outputdir/$samplename/cufflinks/transcripts.gtf >> $outputdir/ucsc/$samplename.transcripts.gtf
gzip -f $outputdir/ucsc/$samplename.transcripts.gtf


echo \"JOBDONE!\"

">$outputdir/$samplename/$samplename.sge

#cat $outputdir/$samplename/$samplename.sge

jobid=`qsub $outputdir/$samplename/$samplename.sge | cut -f3 -d' '`
ln -fs $HOME/sge_jobs_output/sge_job.$jobid.out $outputdir/$samplename/sge.log

echo "
Your job is submitted (jobID: $jobid) !
The SGE script is saved as $outputdir/$samplename/$samplename.sge
The running log is saved as $outputdir/$samplename/sge.log

After the job is done, you can copy the UCSC browser files to webserver for browsing, using command like:

rsync -azv $outputdir/ucsc/$samplename* dongx@zlab.umassmed.edu:~dongx/public_html/tracks/RNAseq

"
