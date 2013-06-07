# all junction for adult
junctions=read.table("/home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/junction_analysis_result_forall/junctions2intron.final.tab.withtype")
colnames(junctions)=c('chr', 'start', 'end', 'ID', 'score', 'strand', 'TopTri', 'TophatTrinity_ID', 'Intron.signal', 'Intron.signal_2nt', 'rnaseq.junc.plus.reads', 'rnaseq.junc.minus.reads', 'rnaseq.junc.plus.species', 'rnaseq.junc.minus.species', 'pirna.junc.plus.reads', 'pirna.junc.minus.reads', 'pirna.junc.plus.species', 'pirna.junc.minus.species', 'rnaseq.splicing.plus.reads', 'rnaseq.splicing.minus.reads', 'rnaseq.splicing.plus.species', 'rnaseq.splicing.minus.species', 'pirna.splicing.plus.reads', 'pirna.splicing.minus.reads', 'pirna.splicing.plus.species', 'pirna.splicing.minus.species')

# splicing ratio
junctions2=cbind(junctions[,1:10],
                 rnaseq_junc_reads     = ifelse(junctions$strand=="+", junctions$rnaseq.junc.plus.reads, junctions$rnaseq.junc.minus.reads),
                 rnaseq_junc_species   = ifelse(junctions$strand=="+", junctions$rnaseq.junc.plus.species, junctions$rnaseq.junc.minus.species),
                 pirna_junc_reads      = ifelse(junctions$strand=="+", junctions$pirna.junc.plus.reads, junctions$pirna.junc.minus.reads),
                 pirna_junc_species    = ifelse(junctions$strand=="+", junctions$pirna.junc.plus.species, junctions$pirna.junc.minus.species),
                 rnaseq_splsite_reads     = ifelse(junctions$strand=="+", junctions$rnaseq.splicing.plus.reads, junctions$rnaseq.splicing.minus.reads),
                 rnaseq_splsite_species   = ifelse(junctions$strand=="+", junctions$rnaseq.splicing.plus.species, junctions$rnaseq.splicing.minus.species),
                 pirna_splsite_reads      = ifelse(junctions$strand=="+", junctions$pirna.splicing.plus.reads, junctions$pirna.splicing.minus.reads),
                 pirna_splsite_species    = ifelse(junctions$strand=="+", junctions$pirna.splicing.plus.species, junctions$pirna.splicing.minus.species),
                 rnaseq_splicing_ratio = ifelse(junctions$strand=="+",
                                         junctions$rnaseq.junc.plus.reads / (junctions$rnaseq.splicing.plus.reads + junctions$rnaseq.junc.plus.reads),
                                         junctions$rnaseq.junc.minus.reads / (junctions$rnaseq.splicing.minus.reads + junctions$rnaseq.junc.minus.reads)),
                 pirna_splicing_ratio  = ifelse(junctions$strand=="+",
                                         junctions$pirna.junc.plus.reads / (junctions$pirna.splicing.plus.reads + junctions$pirna.junc.plus.reads),
                                         junctions$pirna.junc.minus.reads / (junctions$pirna.splicing.minus.reads + junctions$pirna.junc.minus.reads))
                 );

write.table(junctions2, "/home/dongx/projects/piRNA/data/mouse_adult_wt_RNAseq_PE50nt_strand.junctions2intron.all.final.merged.tab", sep="\t", row.names=F,col.names=F, quote=F)


# read junctions (10dpp)
junctions=read.table("/home/dongx/scratch/mouse_10.5dpp_wt_RNAseq_PE50nt_strand/junction_analysis_result/junctions2intron.pi.final.tab")
colnames(junctions)=c('chr', 'start', 'end', 'ID', 'score', 'strand', 'Intron.signal', 'Intron.signal_2nt', 'rnaseq.junc.plus.reads', 'rnaseq.junc.minus.reads', 'rnaseq.junc.plus.species', 'rnaseq.junc.minus.species', 'pirna.junc.plus.reads', 'pirna.junc.minus.reads', 'pirna.junc.plus.species', 'pirna.junc.minus.species', 'rnaseq.splicing.plus.reads', 'rnaseq.splicing.minus.reads', 'rnaseq.splicing.plus.species', 'rnaseq.splicing.minus.species', 'pirna.splicing.plus.reads', 'pirna.splicing.minus.reads', 'pirna.splicing.plus.species', 'pirna.splicing.minus.species')

# splicing ratio
junctions2=cbind(junctions[,1:8],
                 rnaseq_junc_reads     = ifelse(junctions$strand=="+", junctions$rnaseq.junc.plus.reads, junctions$rnaseq.junc.minus.reads),
                 rnaseq_junc_species   = ifelse(junctions$strand=="+", junctions$rnaseq.junc.plus.species, junctions$rnaseq.junc.minus.species),
                 pirna_junc_reads      = ifelse(junctions$strand=="+", junctions$pirna.junc.plus.reads, junctions$pirna.junc.minus.reads),
                 pirna_junc_species    = ifelse(junctions$strand=="+", junctions$pirna.junc.plus.species, junctions$pirna.junc.minus.species),
                 rnaseq_splsite_reads     = ifelse(junctions$strand=="+", junctions$rnaseq.splicing.plus.reads, junctions$rnaseq.splicing.minus.reads),
                 rnaseq_splsite_species   = ifelse(junctions$strand=="+", junctions$rnaseq.splicing.plus.species, junctions$rnaseq.splicing.minus.species),
                 pirna_splsite_reads      = ifelse(junctions$strand=="+", junctions$pirna.splicing.plus.reads, junctions$pirna.splicing.minus.reads),
                 pirna_splsite_species    = ifelse(junctions$strand=="+", junctions$pirna.splicing.plus.species, junctions$pirna.splicing.minus.species),
                 rnaseq_splicing_ratio = ifelse(junctions$strand=="+",
                                         junctions$rnaseq.junc.plus.reads / (junctions$rnaseq.splicing.plus.reads + junctions$rnaseq.junc.plus.reads),
                                         junctions$rnaseq.junc.minus.reads / (junctions$rnaseq.splicing.minus.reads + junctions$rnaseq.junc.minus.reads)),
                 pirna_splicing_ratio  = ifelse(junctions$strand=="+",
                                         junctions$pirna.junc.plus.reads / (junctions$pirna.splicing.plus.reads + junctions$pirna.junc.plus.reads),
                                         junctions$pirna.junc.minus.reads / (junctions$pirna.splicing.minus.reads + junctions$pirna.junc.minus.reads)),
                 junctions[,9:ncol(junctions)]);

write.table(junctions2, "/home/dongx/projects/piRNA/data/mouse_10.5dpp_wt_RNAseq_PE50nt_strand.junctions2intron.pi.final.merged.tab", sep="\t", row.names=F, quote=F)


junctions=read.table("/home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/junction_analysis_result/junctions2intron.pi.final.tab")
colnames(junctions)=c('chr', 'start', 'end', 'ID', 'score', 'strand', 'Intron.signal', 'Intron.signal_2nt', 'rnaseq.junc.plus.reads', 'rnaseq.junc.minus.reads', 'rnaseq.junc.plus.species', 'rnaseq.junc.minus.species', 'pirna.junc.plus.reads', 'pirna.junc.minus.reads', 'pirna.junc.plus.species', 'pirna.junc.minus.species', 'rnaseq.splicing.plus.reads', 'rnaseq.splicing.minus.reads', 'rnaseq.splicing.plus.species', 'rnaseq.splicing.minus.species', 'pirna.splicing.plus.reads', 'pirna.splicing.minus.reads', 'pirna.splicing.plus.species', 'pirna.splicing.minus.species')

# splicing ratio
junctions2=cbind(junctions[,1:8],
                 rnaseq_junc_reads     = ifelse(junctions$strand=="+", junctions$rnaseq.junc.plus.reads, junctions$rnaseq.junc.minus.reads),
                 rnaseq_junc_species   = ifelse(junctions$strand=="+", junctions$rnaseq.junc.plus.species, junctions$rnaseq.junc.minus.species),
                 pirna_junc_reads      = ifelse(junctions$strand=="+", junctions$pirna.junc.plus.reads, junctions$pirna.junc.minus.reads),
                 pirna_junc_species    = ifelse(junctions$strand=="+", junctions$pirna.junc.plus.species, junctions$pirna.junc.minus.species),
                 rnaseq_splsite_reads     = ifelse(junctions$strand=="+", junctions$rnaseq.splicing.plus.reads, junctions$rnaseq.splicing.minus.reads),
                 rnaseq_splsite_species   = ifelse(junctions$strand=="+", junctions$rnaseq.splicing.plus.species, junctions$rnaseq.splicing.minus.species),
                 pirna_splsite_reads      = ifelse(junctions$strand=="+", junctions$pirna.splicing.plus.reads, junctions$pirna.splicing.minus.reads),
                 pirna_splsite_species    = ifelse(junctions$strand=="+", junctions$pirna.splicing.plus.species, junctions$pirna.splicing.minus.species),
                 rnaseq_splicing_ratio = ifelse(junctions$strand=="+",
                                         junctions$rnaseq.junc.plus.reads / (junctions$rnaseq.splicing.plus.reads + junctions$rnaseq.junc.plus.reads),
                                         junctions$rnaseq.junc.minus.reads / (junctions$rnaseq.splicing.minus.reads + junctions$rnaseq.junc.minus.reads)),
                 pirna_splicing_ratio  = ifelse(junctions$strand=="+",
                                         junctions$pirna.junc.plus.reads / (junctions$pirna.splicing.plus.reads + junctions$pirna.junc.plus.reads),
                                         junctions$pirna.junc.minus.reads / (junctions$pirna.splicing.minus.reads + junctions$pirna.junc.minus.reads)),
                 junctions[,9:ncol(junctions)]);

write.table(junctions2, "/home/dongx/projects/piRNA/data/mouse_adult_wt_RNAseq_PE50nt_strand.junctions2intron.pi.final.merged.tab", sep="\t", row.names=F, quote=F)



# read junctions mouse_adult_wt_RNAseq_PE50nt_strand
junctions=read.table("/home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/junctions2intron.pi.final.tab")
colnames(junctions)=c('chr', 'start', 'end', 'ID', 'score', 'strand', 'piRNA_cluster', 'TopTri', 'TophatTrinity_ID', 'Intron.signal', 'Intron.signal_2nt', 'rnaseq.junc.plus.reads', 'rnaseq.junc.minus.reads', 'rnaseq.junc.plus.species', 'rnaseq.junc.minus.species', 'pirna.junc.plus.reads', 'pirna.junc.minus.reads', 'pirna.junc.plus.species', 'pirna.junc.minus.species', 'rnaseq.splicing.plus.reads', 'rnaseq.splicing.minus.reads', 'rnaseq.splicing.plus.species', 'rnaseq.splicing.minus.species', 'pirna.splicing.plus.reads', 'pirna.splicing.minus.reads', 'pirna.splicing.plus.species', 'pirna.splicing.minus.species', 'left_exon.rnaseq.reads', 'left_exon.rnaseq.cov', 'intron.rnaseq.reads', 'intron.rnaseq.cov', 'right_exon.rnaseq.reads', 'right_exon.rnaseq.cov', 'left_exon.pirna.reads', 'left_exon.pirna.cov', 'intron.pirna.reads', 'intron.pirna.cov', 'right_exon.pirna.reads', 'right_exon.pirna.cov')

# splicing ratio
junctions2=cbind(junctions[,1:11],
                 rnaseq_junc_reads     = ifelse(junctions$strand=="+", junctions$rnaseq.junc.plus.reads, junctions$rnaseq.junc.minus.reads),
                 rnaseq_junc_species   = ifelse(junctions$strand=="+", junctions$rnaseq.junc.plus.species, junctions$rnaseq.junc.minus.species),
                 pirna_junc_reads      = ifelse(junctions$strand=="+", junctions$pirna.junc.plus.reads, junctions$pirna.junc.minus.reads),
                 pirna_junc_species    = ifelse(junctions$strand=="+", junctions$pirna.junc.plus.species, junctions$pirna.junc.minus.species),
                 rnaseq_splsite_reads     = ifelse(junctions$strand=="+", junctions$rnaseq.splicing.plus.reads, junctions$rnaseq.splicing.minus.reads),
                 rnaseq_splsite_species   = ifelse(junctions$strand=="+", junctions$rnaseq.splicing.plus.species, junctions$rnaseq.splicing.minus.species),
                 pirna_splsite_reads      = ifelse(junctions$strand=="+", junctions$pirna.splicing.plus.reads, junctions$pirna.splicing.minus.reads),
                 pirna_splsite_species    = ifelse(junctions$strand=="+", junctions$pirna.splicing.plus.species, junctions$pirna.splicing.minus.species),
                 rnaseq_splicing_ratio = ifelse(junctions$strand=="+",
                                         junctions$rnaseq.junc.plus.reads / (junctions$rnaseq.splicing.plus.reads + junctions$rnaseq.junc.plus.reads),
                                         junctions$rnaseq.junc.minus.reads / (junctions$rnaseq.splicing.minus.reads + junctions$rnaseq.junc.minus.reads)),
                 pirna_splicing_ratio  = ifelse(junctions$strand=="+",
                                         junctions$pirna.junc.plus.reads / (junctions$pirna.splicing.plus.reads + junctions$pirna.junc.plus.reads),
                                         junctions$pirna.junc.minus.reads / (junctions$pirna.splicing.minus.reads + junctions$pirna.junc.minus.reads)),
                 rnaseq_coverage_diff  = (junctions$left_exon.rnaseq.cov + junctions$right_exon.rnaseq.cov)/2 - junctions$intron.rnaseq.cov,
                 rnaseq_density_diff   = (junctions$left_exon.rnaseq.reads + junctions$right_exon.rnaseq.reads)/400 - junctions$intron.rnaseq.reads / (junctions$end-junctions$start),
                 pirna_coverage_diff   = (junctions$left_exon.pirna.cov + junctions$right_exon.pirna.cov)/2 - junctions$intron.pirna.cov,
                 pirna_density_diff    = (junctions$left_exon.pirna.reads + junctions$right_exon.pirna.reads)/400 - junctions$intron.pirna.reads / (junctions$end-junctions$start),
                 junctions[,12:ncol(junctions)]);

write.table(junctions2, "/home/dongx/projects/piRNA/data/junctions2intron.pi.final.merged.tab", sep="\t", row.names=F, quote=F)

table(junctions2$rnaseq_junc_reads>0, junctions2$pirna_junc_reads>0)
summary(junctions2$TopTri)
library(ggplot2)
ggplot(data=junctions2, aes(x=factor(Intron.signal), y=log(y), fill=factor(TopTri))) + geom_bar() + opts(axis.text.x=theme_text(angle=-90))

junctions2=junctions2[junctions2$rnaseq_junc_reads>0,]

plot(junctions2$rnaseq_splicing_ratio, junctions2$pirna_splicing_ratio, col='#0000ff33')

df=cbind(pirna_splicing_ratio=round(junctions2$pirna_splicing_ratio,1),rnaseq_splicing_ratio=junctions2$rnaseq_splicing_ratio)
boxplot(rnaseq_splicing_ratio~pirna_splicing_ratio, df)

par(mfcol=c(2,2))
hist(junctions2$rnaseq_splicing_ratio, breaks=100, xlim=c(0,1), main="Histogram of RNAseq_splicing_ratio")
hist(junctions2$rnaseq_splicing_ratio[junctions2$pirna_junc_reads>0], xlim=c(0,1), breaks=100, main="")

hist(junctions2$pirna_splicing_ratio, breaks=100, xlim=c(0,1), main="Histogram of piRNA_splicing_ratio")
hist(junctions2$pirna_splicing_ratio[junctions2$pirna_junc_reads>0], xlim=c(0,1), breaks=100, main="")

par(mfcol=c(2,2))
hist(junctions2$rnaseq_splicing_ratio[junctions2$pirna_splicing_ratio>0.9], xlim=c(0,1), breaks=100)
hist(junctions2$pirna_splicing_ratio[junctions2$pirna_splicing_ratio>0.9], xlim=c(0,1), breaks=100)

hist(junctions2$rnaseq_splicing_ratio[!is.na(junctions2$rnaseq_splicing_ratio) & junctions2$rnaseq_splicing_ratio>.9], xlim=c(0,1), breaks=100)
hist(junctions2$pirna_splicing_ratio[!is.na(junctions2$rnaseq_splicing_ratio) & junctions2$rnaseq_splicing_ratio>0.9], xlim=c(0,1), breaks=100)

# junc format for Tophat
junc=junctions2[junctions2$rnaseq_junc_reads>0,c(1,2,3,6)]
junc$start=junc$start-1;
write.table(junc, "/home/dongx/scratch/mouse_adult_wt_RNAseq_PE50nt_strand/junctions2intron.pi.final.merged.junc", sep="\t", row.names=F, col.names=F, quote=F)
