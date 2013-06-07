#!/bin/sh

# script to test Cufflinks parameters

# -j/--pre-mrna-fraction	-a/--junc-alpha	--min-frags-per-transfrag	--trim-3-avgcov-thresh	--trim-3-dropoff-frac	--overlap-radius

pre_mrna_fraction=(0.2 0.3 0.4)
#junc_alpha=(0.01 0.0001 0.00001) --> no change
#min_frags_per_transfrag=(20 40 60) --> 40
trim_3_avgcov_thresh=(15 20 40 60 80)
trim_3_dropoff_frac=(0.2 0.4 0.6 0.8)
#overlap_radius=(200 250 350)  -> 250 !


for j in ${pre_mrna_fraction[@]}
do
    cuffout=/home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks2_j${j}
    [ -d $cuffout ] || mkdir $cuffout
    echo "
    #!/bin/sh
    #$ -V
    #$ -pe single 24
    #$ -cwd
    #$ -o \$HOME/sge_jobs_output/sge_job.\$JOB_ID.out -j y
    #$ -S /bin/bash
    #$ -l mem_free=3G

    export ANNOTATION=\$GENOME/mm9/Annotation/Genes
    cd /home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/merge_reads

    cufflinks -v --no-update-check --frag-len-mean 400 --frag-len-std-dev 100 --library-type fr-firststrand -o $cuffout -p 24 -g \$ANNOTATION/genes.gtf -M \$ANNOTATION/chrM.rRNA.tRNA.gtf accepted_hits.bam -u -j $j --min-frags-per-transfrag 40 --overlap-radius 250 2&> $cuffout/run.log

    echo \"track name=RNAseq_PE50nt_cufflinks2_j${j} description=mouse_adult_wt_RNAseq_PE50nt_strand_cufflinks2_j${j} visibility=pack colorByStrand='200,100,0 0,100,200'\" > /home/dongx/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand.cufflinks2_j${j}.transcripts.gtf
    cat $cuffout/transcripts.gtf >> /home/dongx/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand.cufflinks2_j${j}.transcripts.gtf
    gzip -f /home/dongx/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand.cufflinks2_j${j}.transcripts.gtf

    echo \"JOB DONE!\";
    "  | sed 's/^[ \t]*//g' > $HOME/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks2_j${j}/submit.sge
    qsub $HOME/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks2_j${j}/submit.sge
    #cat $HOME/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks2_j${j}/submit.sge
done


for ta in ${trim_3_avgcov_thresh[@]}
do
    cuffout=/home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks2_ta${ta}
    [ -d $cuffout ] || mkdir $cuffout
    echo "
    #!/bin/sh
    #$ -V
    #$ -pe single 24
    #$ -cwd
    #$ -o \$HOME/sge_jobs_output/sge_job.\$JOB_ID.out -j y
    #$ -S /bin/bash
    #$ -l mem_free=3G

    export ANNOTATION=\$GENOME/mm9/Annotation/Genes
    cd /home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/merge_reads

    cufflinks -v --no-update-check --frag-len-mean 400 --frag-len-std-dev 100 --library-type fr-firststrand -o $cuffout -p 24 -g \$ANNOTATION/genes.gtf -M \$ANNOTATION/chrM.rRNA.tRNA.gtf accepted_hits.bam -u --trim-3-avgcov-thresh $ta --min-frags-per-transfrag 40 --overlap-radius 250 2&> $cuffout/run.log

    echo \"track name=RNAseq_PE50nt_cufflinks2_ta${ta} description=mouse_adult_wt_RNAseq_PE50nt_strand_cufflinks2_ta${ta} visibility=pack colorByStrand='200,100,0 0,100,200'\" > /home/dongx/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand.cufflinks2_ta${ta}.transcripts.gtf
    cat $cuffout/transcripts.gtf >> /home/dongx/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand.cufflinks2_ta${ta}.transcripts.gtf
    gzip -f /home/dongx/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand.cufflinks2_ta${ta}.transcripts.gtf

    echo \"JOB DONE!\";
    "  | sed 's/^[ \t]*//g' > $HOME/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks2_ta${ta}/submit.sge
    qsub $HOME/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks2_ta${ta}/submit.sge
    #cat $HOME/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks2_ta${ta}/submit.sge
done

for td in ${trim_3_dropoff_frac[@]}
do
    cuffout=/home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks2_td${td}
    [ -d $cuffout ] || mkdir $cuffout
    echo "
    #!/bin/sh
    #$ -V
    #$ -pe single 24
    #$ -cwd
    #$ -o \$HOME/sge_jobs_output/sge_job.\$JOB_ID.out -j y
    #$ -S /bin/bash
    #$ -l mem_free=3G

    export ANNOTATION=\$GENOME/mm9/Annotation/Genes
    cd /home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/merge_reads

    cufflinks -v --no-update-check --frag-len-mean 400 --frag-len-std-dev 100 --library-type fr-firststrand -o $cuffout -p 24 -g \$ANNOTATION/genes.gtf -M \$ANNOTATION/chrM.rRNA.tRNA.gtf accepted_hits.bam -u --trim-3-dropoff-frac $td --min-frags-per-transfrag 40 --overlap-radius 250 2&> $cuffout/run.log

    echo \"track name=RNAseq_PE50nt_cufflinks2_td${td} description=mouse_adult_wt_RNAseq_PE50nt_strand_cufflinks2_td${td} visibility=pack colorByStrand='200,100,0 0,100,200'\" > /home/dongx/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand.cufflinks2_td${td}.transcripts.gtf
    cat $cuffout/transcripts.gtf >> /home/dongx/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand.cufflinks2_td${td}.transcripts.gtf
    gzip -f /home/dongx/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand.cufflinks2_td${td}.transcripts.gtf

    echo \"JOB DONE!\";
    "  | sed 's/^[ \t]*//g' > $HOME/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks2_td${td}/submit.sge
    qsub $HOME/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks2_td${td}/submit.sge
    #cat $HOME/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks2_td${td}/submit.sge
done


exit

for a in ${junc_alpha[@]}
do
    cuffout=/home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks2_a${a}
    [ -d $cuffout ] || mkdir $cuffout
    echo "
    #!/bin/sh
    #$ -V
    #$ -pe single 24
    #$ -cwd
    #$ -o \$HOME/sge_jobs_output/sge_job.\$JOB_ID.out -j y
    #$ -S /bin/bash
    #$ -l mem_free=3G

    export ANNOTATION=\$GENOME/mm9/Annotation/Genes
    cd /home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/merge_reads

    cufflinks -v --no-update-check --frag-len-mean 400 --frag-len-std-dev 100 --library-type fr-firststrand -o $cuffout -p 24 -g \$ANNOTATION/genes.gtf -M \$ANNOTATION/chrM.rRNA.tRNA.gtf accepted_hits.bam -u -a $a --overlap-radius 250 2&> $cuffout/run.log

    echo \"track name=RNAseq_PE50nt_cufflinks2_a${a} description=mouse_adult_wt_RNAseq_PE50nt_strand_cufflinks2_a${a} visibility=pack colorByStrand='200,100,0 0,100,200'\" > /home/dongx/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand.cufflinks2_a${a}.transcripts.gtf
    cat $cuffout/transcripts.gtf >> /home/dongx/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand.cufflinks2_a${a}.transcripts.gtf
    gzip -f /home/dongx/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand.cufflinks2_a${a}.transcripts.gtf

    echo \"JOB DONE!\";
    "  | sed 's/^[ \t]*//g' > $HOME/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks2_a${a}/submit.sge
    qsub $HOME/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks2_a${a}/submit.sge
    #cat $HOME/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks2_a${a}/submit.sge
done

for m in ${min_frags_per_transfrag[@]}
do
    cuffout=/home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks2_m${m}
    [ -d $cuffout ] || mkdir $cuffout
    echo "
    #!/bin/sh
    #$ -V
    #$ -pe single 24
    #$ -cwd
    #$ -o \$HOME/sge_jobs_output/sge_job.\$JOB_ID.out -j y
    #$ -S /bin/bash
    #$ -l mem_free=3G

    export ANNOTATION=\$GENOME/mm9/Annotation/Genes
    cd /home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/merge_reads

    cufflinks -v --no-update-check --frag-len-mean 400 --frag-len-std-dev 100 --library-type fr-firststrand -o $cuffout -p 24 -g \$ANNOTATION/genes.gtf -M \$ANNOTATION/chrM.rRNA.tRNA.gtf accepted_hits.bam -u --min-frags-per-transfrag $m 2&> $cuffout/run.log

    echo \"track name=RNAseq_PE50nt_cufflinks2_m${m} description=mouse_adult_wt_RNAseq_PE50nt_strand_cufflinks2_m${m} visibility=pack colorByStrand='200,100,0 0,100,200'\" > /home/dongx/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand.cufflinks2_m${m}.transcripts.gtf
    cat $cuffout/transcripts.gtf >> /home/dongx/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand.cufflinks2_m${m}.transcripts.gtf
    gzip -f /home/dongx/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand.cufflinks2_m${m}.transcripts.gtf

    echo \"JOB DONE!\";
    "  | sed 's/^[ \t]*//g' > $HOME/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks2_m${m}/submit.sge
    qsub $HOME/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks2_m${m}/submit.sge
    #cat $HOME/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks2_m${m}/submit.sge
done

for o in ${overlap_radius[@]}
do
    cuffout=/home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks2_o${o}
    [ -d $cuffout ] || mkdir $cuffout
    echo "
    #!/bin/sh
    #$ -V
    #$ -pe single 24
    #$ -cwd
    #$ -o \$HOME/sge_jobs_output/sge_job.\$JOB_ID.out -j y
    #$ -S /bin/bash
    #$ -l mem_free=3G

    export ANNOTATION=\$GENOME/mm9/Annotation/Genes
    cd /home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/merge_reads

    cufflinks -v --no-update-check --frag-len-mean 400 --frag-len-std-dev 100 --library-type fr-firststrand -o $cuffout -p 24 -g \$ANNOTATION/genes.gtf -M \$ANNOTATION/chrM.rRNA.tRNA.gtf accepted_hits.bam -u --overlap-radius $o 2&> $cuffout/run.log

    echo \"track name=RNAseq_PE50nt_cufflinks2_o${o} description=mouse_adult_wt_RNAseq_PE50nt_strand_cufflinks2_o${o} visibility=pack colorByStrand='200,100,0 0,100,200'\" > /home/dongx/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand.cufflinks2_o${o}.transcripts.gtf
    cat $cuffout/transcripts.gtf >> /home/dongx/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand.cufflinks2_o${o}.transcripts.gtf
    gzip -f /home/dongx/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand.cufflinks2_o${o}.transcripts.gtf

    echo \"JOB DONE!\";
    "  | sed 's/^[ \t]*//g' > $HOME/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks2_o${o}/submit.sge
    qsub $HOME/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks2_o${o}/submit.sge
    #cat $HOME/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks2_o${o}/submit.sge
done


exit;
#### all combinations

for j in ${pre_mrna_fraction[@]}
do
    for a in ${junc_alpha[@]}
    do
        for m in ${min_frags_per_transfrag[@]}
        do
            for ta in ${trim_3_avgcov_thresh[@]}
            do
                for td in ${trim_3_dropoff_frac[@]}
                do
                    for o in ${overlap_radius[@]}
                    do
                        cuffout=/home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks_j${j}_a${a}_m${m}_ta${ta}_td${td}_o${o}
                        [ -d $cuffout ] || mkdir $cuffout
                        echo "
                        #!/bin/sh
                        #$ -V
                        #$ -pe single 24
                        #$ -cwd
                        #$ -o \$HOME/sge_jobs_output/sge_job.\$JOB_ID.out -j y
                        #$ -S /bin/bash
                        #$ -l mem_free=3G

                        export ANNOTATION=\$GENOME/mm9/Annotation/Genes
                        cd /home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/merge_reads

                        cufflinks -v --no-update-check --frag-len-mean 400 --frag-len-std-dev 100 --library-type fr-firststrand -o $cuffout -p 24 -g \$ANNOTATION/genes.gtf -M \$ANNOTATION/chrM.rRNA.tRNA.gtf accepted_hits.bam -u -j $j -a $a --min-frags-per-transfrag $m --trim-3-avgcov-thresh $ta --trim-3-dropoff-frac $td --overlap-radius $o 2&> $cuffout/run.log

                        echo \"track name=RNAseq_PE50nt_cufflinks_j${j}_a${a}_m${m}_ta${ta}_td${td}_o${o} description=mouse_adult_wt_RNAseq_PE50nt_strand_cufflinks_j${j}_a${a}_m${m}_ta${ta}_td${td}_o${o} visibility=pack colorByStrand='200,100,0 0,100,200'\" > /home/dongx/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand.cufflinks_j${j}_a${a}_m${m}_ta${ta}_td${td}_o${o}.transcripts.gtf
                        cat $cuffout/transcripts.gtf >> /home/dongx/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand.cufflinks_j${j}_a${a}_m${m}_ta${ta}_td${td}_o${o}.transcripts.gtf
                        gzip -f /home/dongx/scratch/ucsc/mouse_adult_wt_RNAseq_PE50nt_strand.cufflinks_j${j}_a${a}_m${m}_ta${ta}_td${td}_o${o}.transcripts.gtf

                        echo \"JOB DONE!\";
                        "  | sed 's/^[ \t]*//g' > $HOME/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks_j${j}_a${a}_m${m}_ta${ta}_td${td}_o${o}/submit.sge
                        qsub $HOME/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks_j${j}_a${a}_m${m}_ta${ta}_td${td}_o${o}/submit.sge
                        #cat $HOME/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/cufflinks_j${j}_a${a}_m${m}_ta${ta}_td${td}_o${o}/submit.sge
                    done
                done
            done
        done
    done
done
