## my Tophat mappers (-k 100)
#piRNA_plus=Phil.SRA.wt.ox.6wk.testes.raw.xkxh.norm.all0.accepted_hits.+.bw
#piRNA_minus=Phil.SRA.wt.ox.6wk.testes.raw.xkxh.norm.all0.accepted_hits.-.bw
## all mapper(Jia)
#piRNA_plus=/home/dongx/nearline/Xin/smallRNA/Jia/Phil.SRA.wt.ox.6wk.testes.raw.xkxh.norm.bed.gz.plus.bw
#piRNA_minus=/home/dongx/nearline/Xin/smallRNA/Jia/Phil.SRA.wt.ox.6wk.testes.raw.xkxh.norm.bed.gz.minus.bw
#piRNA10dpp_plus=/home/dongx/nearline/Xin/smallRNA/Jia/Phil.SRA.10dpp.ox.testes.rep2.raw.xkxh.norm.bed.gz.plus.bw
#piRNA10dpp_minus=/home/dongx/nearline/Xin/smallRNA/Jia/Phil.SRA.10dpp.ox.testes.rep2.raw.xkxh.norm.bed.gz.minus.bw
## uniq mapper (Jia)
piRNA_plus=/home/dongx/nearline/Xin/smallRNA/Jia/Phil.SRA.wt.ox.6wk.testes.raw.uniqmap.xkxh.norm.bed.gz.plus.bw
piRNA_minus=/home/dongx/nearline/Xin/smallRNA/Jia/Phil.SRA.wt.ox.6wk.testes.raw.uniqmap.xkxh.norm.bed.gz.minus.bw
piRNA10dpp_plus=/home/dongx/nearline/Xin/smallRNA/Jia/Phil.SRA.10dpp.ox.testes.rep2.raw.uniqmap.xkxh.norm.bed.gz.plus.bw
piRNA10dpp_minus=/home/dongx/nearline/Xin/smallRNA/Jia/Phil.SRA.10dpp.ox.testes.rep2.raw.uniqmap.xkxh.norm.bed.gz.minus.bw

RNA_plus=/home/dongx/scratch/ucscmouse_adult_wt_RNAseq_PE50nt_strand_R1.accepted_hits.+.bw
RNA_minus=/home/dongx/scratch/ucscmouse_adult_wt_RNAseq_PE50nt_strand_R1.accepted_hits.-.bw
RNA10dpp_plus=/home/dongx/scratch/mouse_10.5dpp_wt_RNAseq_PE50nt_strand/wild.type.10.rep2.rnaseq.accepted_hits_uniqmap.plus.bw
RNA10dpp_minus=/home/dongx/scratch/mouse_10.5dpp_wt_RNAseq_PE50nt_strand/wild.type.10.rep2.rnaseq.accepted_hits_uniqmap.minus.bw
CAGE_plus=/home/dongx/scratch/ucscmouse_adult_wt_CAGE_PE100nt_strand.unique.plus.bw
CAGE_minus=/home/dongx/scratch/ucscmouse_adult_wt_CAGE_PE100nt_strand.unique.minus.bw
CAGE10dpp_plus=/home/dongx/scratch/ucscmouse_10dpp_wt_CAGE_100nt_strand.unique.plus.bw
CAGE10dpp_minus=/home/dongx/scratch/ucscmouse_10dpp_wt_CAGE_100nt_strand.unique.minus.bw
PAS_plus=/home/dongx/scratch/ucscmouse_adult_wt_PAS200_100nt.unique.plus.bw
PAS_minus=/home/dongx/scratch/ucscmouse_adult_wt_PAS200_100nt.unique.minus.bw
#Amyb=/home/dongx/scratch/mouse_adult_wt_InputIP_50nt/input.mouse.by.chris/A-Myb-mouse_R2.gz.k1.sam.sorted.bam.bigwig
#Input=/home/dongx/scratch/mouse_adult_wt_InputIP_50nt/input.mouse.by.chris/Input-mouse_R1.gz.k1.sam.sorted.bam.bigwig
#Pol2=/home/dongx/scratch/mouse_adult_wt_InputIP_50nt/input.mouse.by.chris/pol2.chip.M1.fastq.bw
#H3K4me3=/home/dongx/scratch/mouse_adult_wt_InputIP_50nt/input.mouse.by.chris/H3K4-mouse_R1.gz.k1.sam.sorted.bam.bigwig
Amyb=/home/dongx/scratch/mouse_adult_wt_InputIP_50nt/input.mouse.by.chris/A-Myb-mouse_R2.gz.m1.sam.sorted.bw
Input=/home/dongx/scratch/mouse_adult_wt_InputIP_50nt/input.mouse.by.chris/Input-mouse_R1.gz.m1.sam.sorted.bam.bigwig
Pol2=/home/dongx/scratch/mouse_adult_wt_InputIP_50nt/input.mouse.by.chris/pol2.chip.v3.fastq.sam.sorted.bam.bigWig
H3K4me3=/home/dongx/scratch/mouse_adult_wt_InputIP_50nt/input.mouse.by.chris/H3K4-mouse_R1.gz.m1.sam.sorted.bam.bigwig

rm /tmp/chr*
for i in chr10:66,142,648-66,171,990 chr1:93,294,330-93,307,787 chr4:135,172,729-135,191,732 chr17:49,247,739-49,256,216 chr2:92,361,211-92,457,307 chr7:77,017,882-77,057,749 chr7:77,051,484-77,113,583 chr4:93,942,079-94,005,089 chr11:69,209,105-69,210,391 chr13:23,855,439-23,856,311 chr19:10,056,976-10,059,641 chr13:55,205,999-55,237,283;
## NOT USED FINALLY #for i in chr10:66,140,970-66,172,450 chr1:93,294,330-93,307,787 chr4:135,172,729-135,191,732 chr17:49,246,365-49,257,623 chr2:92,361,211-92,457,307 chr7:77,017,882-77,057,749 chr7:77,051,484-77,113,583 chr4:93,942,079-94,005,089 chr11:69,209,105-69,210,391 chr13:23,855,439-23,856,311 chr19:10,056,976-10,059,641 chr13:55,205,999-55,237,283;
do
    read -a coor <<< $(echo $i | sed 's/,//g;s/[-:]/ /g;')
    chr=${coor[0]}
    start=${coor[1]}
    end=${coor[2]}
    let l=end-start
    # normalization factor
    #echo -e "normalization_factor\t16\t16\t813\t813\t10\t10\t3\t3\t45\t45\t32\t32\t1\t1\t52\t17\t8\t18" > /tmp/$chr-$start-$end.xls
    # header
    echo -e "nt_position\tpiRNA_plus\tpiRNA_minus\tpiRNA10dpp_plus\tpiRNA10dpp_minus\tRNA_plus\tRNA_minus\tRNA10dpp_plus\tRNA10dpp_minus\tCAGE_plus\tCAGE_minus\tCAGE10dpp_plus\tCAGE10dpp_minus\tPAS_plus\tPAS_minus\tH3K4me3\tPol2\tAmyb\tInput" > /tmp/$chr-$start-$end.xls
    seq -f "%.0f" 1 $l > /tmp/$chr-$start-$end.txt
    for bw in $piRNA_plus $piRNA_minus $piRNA10dpp_plus $piRNA10dpp_minus $RNA_plus $RNA_minus $RNA10dpp_plus $RNA10dpp_minus $CAGE_plus $CAGE_minus $CAGE10dpp_plus $CAGE10dpp_minus $PAS_plus $PAS_minus $H3K4me3 $Pol2 $Amyb $Input;
    do
        case "$bw" in
                *SRA.wt.ox.6wk*)
                    let normalizationFactor=16
                    ;;
                *SRA.10dpp.ox.testes*)
                    let normalizationFactor=813
                    ;;
                *mouse_adult_wt_RNAseq_PE50nt_*)
                    let normalizationFactor=10
                    ;;
                *mouse_10.5dpp_wt_RNAseq*)
                    let normalizationFactor=3
                    ;;
                *mouse_adult_wt_CAGE*)
                    let normalizationFactor=45
                    ;;
                *mouse_10dpp_wt_CAGE*)
                    let normalizationFactor=32
                    ;;
                *mouse_adult_wt_PAS*)
                    let normalizationFactor=1
                    ;;
                *A-Myb-mouse*)
                    let normalizationFactor=8
                    ;;
                *Input-mouse*)
                    let normalizationFactor=18
                    ;;
                *pol2.chip*)
                    let normalizationFactor=17
                    ;;
                *H3K4*)
                    let normalizationFactor=52
                    ;;
                *)
                    echo "Invalid value of $bw"
                    exit 1
        esac
        #echo "normalizationFactor" $normalizationFactor;
        bigWigSummary $bw $chr $start $end $l | sed 's/n\/a/0/g;s/\t/\n/g;s/-//g' | awk -v nf=$normalizationFactor '{print $1/nf}' | paste /tmp/$chr-$start-$end.txt - | sed 's/\t\t/\t0\t/g' > /tmp/test.txt
        #bigWigSummary $bw $chr $start $end $l | sed 's/n\/a/0/g;s/\t/\n/g;s/-//g' | paste /tmp/$chr-$start-$end.txt - | sed 's/\t\t/\t0\t/g' > /tmp/test.txt
        cp /tmp/test.txt /tmp/$chr-$start-$end.txt
        echo $bw
    done
    cat /tmp/$chr-$start-$end.txt >> /tmp/$chr-$start-$end.xls
done
