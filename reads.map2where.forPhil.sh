################################################################
# The wrapper script to see where the mapping reads mapped to
# Usage:
# for i in ~/scratch/mouse_*dpp_wt_smallRNA_ox_rep1/allmap.bed; do qsub $HOME/projects/piRNA/src/reads.map2where.forPhil.sh $i; done
# version: 1.0
# date: 2013-05-28
# Author: Xianjun Dong
## Requirement: the input bed files are already sorted (e.g. by sort-bed -)
################################################################

#!/bin/sh

#$ -V
#$ -pe openmpi 8
#$ -cwd
#$ -o $HOME/sge_jobs_output/sge_job.$JOB_ID.out -j y
#$ -S /bin/bash
#$ -l mem_free=5G

export BOWTIE_INDEXES=$GENOME/mm9/Sequence/BowtieIndex/
export ANNOTATION=$GENOME/mm9/Annotation
#ncRNA=$ANNOTATION/ncRNA/ncRNA.wt_miRNA_piRNA.bed12  # ncRNA from fRNAdb (http://www.ncrna.org/frnadb/download); it seems fRNAdb has mixed tRNA with miRNA...
ncRNA=$ANNOTATION/ncRNA/ncRNA.bed  # mergen tRNA from tRNAscan-SE and rRNA/snRNA/snoRNA from biomart
rRNA=$ANNOTATION/ncRNA/rRNA.bed
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
#> $BEDfile_fullpath.mapping.stat.forPhil.txt

#echo "# 1. mapped to rRNA"
#intersectBed -a $BEDfile_fullpath -b $rRNA -s -wo -sorted -f 0.5 > $BEDfile_fullpath"_to_rRNA"
#awk -v type="01.rRNA" -f $HOME/projects/piRNA/src/get_species_reads $BEDfile_fullpath"_to_rRNA" >> $BEDfile_fullpath.mapping.stat.forPhil.txt
#
#echo "# 2. genome mapping reads (-rRNA; +miRNA_hairpin)"
#fgrep -f <(cat $BEDfile_fullpath"_to_rRNA" | awk '{a[$4]=1;}END{for(i in a) print i;}') -wv $BEDfile_fullpath > $BEDfile_fullpath"_to_norRNA"
#awk -v type="02.no_rRNA" -f $HOME/projects/piRNA/src/get_species_reads $BEDfile_fullpath"_to_norRNA" >> $BEDfile_fullpath.mapping.stat.forPhil.txt
#
#echo "# 3. miRNA hairpin reads"
#intersectBed -a $BEDfile_fullpath -b $miRNA -s -wo -sorted -f 0.5 > $BEDfile_fullpath"_to_miRNA"
#awk -v type="03.miRNA" -f $HOME/projects/piRNA/src/get_species_reads $BEDfile_fullpath"_to_miRNA" >> $BEDfile_fullpath.mapping.stat.forPhil.txt
#
#echo "# 4. genome mapping reads (-rRNA; -miRNA_hairpin)"
#fgrep -f <(cat $BEDfile_fullpath"_to_ncRNA" $BEDfile_fullpath"_to_miRNA" | awk '{a[$4]=1;}END{for(i in a) print i;}') -vw $BEDfile_fullpath > $BEDfile_fullpath"_to_norRNA_nomiRNA"
#awk -v type="04.norRNA_nomiRNA" -f $HOME/projects/piRNA/src/get_species_reads $BEDfile_fullpath"_to_norRNA_nomiRNA" >> $BEDfile_fullpath.mapping.stat.forPhil.txt
#
#echo "# 5. genome mapping reads (-rRNA; -miRNA_hairpin; -exon_exon_junction_mapper)"
#awk '{if($10==1) print}' $BEDfile_fullpath"_to_norRNA_nomiRNA" > $BEDfile_fullpath"_to_norRNA_nomiRNA_nojunction"
#awk -v type="05.norRNA_nomiRNA_nojunction" -f $HOME/projects/piRNA/src/get_species_reads $BEDfile_fullpath"_to_norRNA_nomiRNA_nojunction" >> $BEDfile_fullpath.mapping.stat.forPhil.txt
#
#echo "# 6. genome unique mapping reads (-rRNA; -miRNA_hairpi; -exon_exon_junction_mapper)"
#awk -v uniq=$BEDfile_fullpath"_to_norRNA_nomiRNA_nojunction_u" -v multi=$BEDfile_fullpath"_to_norRNA_nomiRNA_nojunction_m" '{if($5==1) print >> uniq; if($5>1) print >> multi;}' $BEDfile_fullpath"_to_norRNA_nomiRNA"
#awk -v type="06.norRNA_nomiRNA_nojunction_uniq" -f $HOME/projects/piRNA/src/get_species_reads $BEDfile_fullpath"_to_norRNA_nomiRNA_nojunction_u" >> $BEDfile_fullpath.mapping.stat.forPhil.txt
#
#echo "# 7. genome multiple mapping reads (-rRNA; -miRNA_hairpi; -exon_exon_junction_mapper)"
#awk -v type="07.norRNA_nomiRNA_nojunction_multi" -f $HOME/projects/piRNA/src/get_species_reads $BEDfile_fullpath"_to_norRNA_nomiRNA_nojunction_m" >> $BEDfile_fullpath.mapping.stat.forPhil.txt
#
#echo "# 8. splice junction reads"
#awk '{if($10==2) print;}' $BEDfile_fullpath > $BEDfile_fullpath"_to_junction"
#awk -v type="08.junction" -f $HOME/projects/piRNA/src/get_species_reads $BEDfile_fullpath"_to_junction" >> $BEDfile_fullpath.mapping.stat.forPhil.txt
#
#echo "# 9. L1 direct mapper reads (independent of genome mapping)"
#intersectBed -a $BEDfile_fullpath -b <(fgrep "L1__LINE" $repeats) -s -wo -sorted -f 0.5 > $BEDfile_fullpath"_to_L1"
#awk -v type="09.L1" -f $HOME/projects/piRNA/src/get_species_reads $BEDfile_fullpath"_to_L1" >> $BEDfile_fullpath.mapping.stat.forPhil.txt
#
#echo "# 10. IAP direct mapper reads (independent of genome mapping)"
#intersectBed -a $BEDfile_fullpath -b <(grep -P "IAP.*__ERVK__LTR" $repeats) -s -wo -sorted -f 0.5 > $BEDfile_fullpath"_to_IAP"
#awk -v type="10.IAP" -f $HOME/projects/piRNA/src/get_species_reads $BEDfile_fullpath"_to_IAP" >> $BEDfile_fullpath.mapping.stat.forPhil.txt

echo "
#awk -v type=\"04.norRNA_nomiRNA\" -f $HOME/projects/piRNA/src/get_species_reads ${BEDfile_fullpath}_to_norRNA_nomiRNA >> $BEDfile_fullpath.mapping.stat.forPhil.txt
#awk -v type=\"05.norRNA_nomiRNA_nojunction\" -f $HOME/projects/piRNA/src/get_species_reads ${BEDfile_fullpath}_to_norRNA_nomiRNA_nojunction >> $BEDfile_fullpath.mapping.stat.forPhil.txt
#awk -v type=\"06.norRNA_nomiRNA_nojunction_uniq\" -f $HOME/projects/piRNA/src/get_species_reads ${BEDfile_fullpath}_to_norRNA_nomiRNA_nojunction_u >> $BEDfile_fullpath.mapping.stat.forPhil.txt
awk -v type=\"07.norRNA_nomiRNA_nojunction_multi\" -f $HOME/projects/piRNA/src/get_species_reads ${BEDfile_fullpath}_to_norRNA_nomiRNA_nojunction_m >> $BEDfile_fullpath.mapping.stat.forPhil.txt
awk '{if($10==2) print;}' $BEDfile_fullpath > ${BEDfile_fullpath}_to_junction
#awk -v type=\"08.junction\" -f $HOME/projects/piRNA/src/get_species_reads ${BEDfile_fullpath}_to_junction >> $BEDfile_fullpath.mapping.stat.forPhil.txt
awk -v type=\"09.L1\" -f $HOME/projects/piRNA/src/get_species_reads ${BEDfile_fullpath}_to_L1 >> $BEDfile_fullpath.mapping.stat.forPhil.txt
awk -v type=\"10.IAP\" -f $HOME/projects/piRNA/src/get_species_reads ${BEDfile_fullpath}_to_IAP >> $BEDfile_fullpath.mapping.stat.forPhil.txt
" > $BEDfile_fullpath.parafly

/home/hanb/nearline/small_RNA_Pipeline/bin/ParaFly -c $BEDfile_fullpath.parafly -CPU 8
