# script to draw aggregation plot of CAGE/Degredome/PAS signal for meta-gene
# input files: genes annotation file
# Usage: for i in all protein-coding pseudogene lincRNA piRNA; do Rscript draw_aggregation_plot.R $i & done
# Rscript draw_aggregation_plot.BIP.R piRNA
# requirement: run bigWigAverageOverBed_81bins.sh before (see detail there)

source('http://zlab.umassmed.edu/~dongx/mylib/mylib.R')

args <- commandArgs(TRUE)
genetype=args[1]
#source('http://zlab.umassmed.edu/~dongx/mylib/mylib.R'); genetype='protein-coding'; type='CAGE'; species='mouse'

#marks=c('H3K4me3','Input','Amyb'); species='rat'
#marks=c('H3K04me1', 'H3K04me3','H3K27ac', 'H3K27me3', 'H3K36me3', 'Pol2', 'Ctcf', 'Input'); species='Encode'
marks=c('CAGE', 'PAS','Degradome.polyA_100nt_strand', 'Degradome.polyA-_100nt_strand', 'H3K4me3', 'Pol2', 'Amyb', 'RNAseq_PE50nt_strand_R1', 'smallRNAseq_76nt_strand_ox'); species='mouse';


# gene file
if(genetype=='piRNA'){
             genes=read.table("../data/piRNA.BIP.cpg.bed", header=T)
             rownames(genes)=genes[,4] # TransID
             pattern="output_10092012.*.[bw|bigWig]"
}
if(genetype=='protein-coding')
{
             genes=read.table("../data/PC.BIP.cpg.bed", header=T)
             rownames(genes)=genes[,4] # TransID
             pattern="output_wgEncodeNhgriBip.*.[bw|bigWig]"
}

# check outliner!! (which can be interesting biologically)
#`%ni%` <- Negate(`%in%`)
#genes=subset(genes, trans_id %ni% c('ENSMUST00000175219','ENSMUST00000175544','ENSMUST00000175227'))

#hcp=rownames(genes)[genes$normalizedCpG>0.4]
#lcp=rownames(genes)[genes$normalizedCpG<=0.4]

QSUBOUPUT="../data"
filenames = paste(QSUBOUPUT, dir(QSUBOUPUT, pattern=pattern), sep="/")  # for de novo piRNA genes

#png(paste("../result/aggregation",species, paste(marks,collapse="_"), genetype,"png",sep="."), width=1200,height=400*length(marks))
#pdf(paste("../result/aggregation","BIP",species, paste(marks,collapse="_"), genetype,"pdf",sep="."), width=12,height=12)
pdf(paste("../result/aggregation","BIP",species, genetype,"pdf",sep="."), width=12,height=12)
par(mfrow=c(3,3))
for(type in marks)
{
             fnames=filenames[grep(paste(species,type,sep=".*"), filenames, ignore.case = T)]
             if(length(fnames)==2)
             {
                          ##  ================= all

                          plus=fnames[grep('plus|\\.\\+\\.', fnames)][1]
                          minus=fnames[grep('minus|\\.-\\.', fnames)][1]

                          PLOTAGGREGATION2(plus, minus, genes, 'All : ', 0)
                          if(genetype=='piRNA'){
                                       PLOTAGGREGATION2(plus, minus, subset(genes, intergenic_genic=='intergenic'), 'Intergenic : ', 0)
                                       PLOTAGGREGATION2(plus, minus, subset(genes, intergenic_genic=='genic'), 'Genic : ', 0)
                          }
             }
             if(length(fnames)==1)
             {
                          PLOTAGGREGATION(fnames, genes, 'All : ', 0)
                          if(genetype=='piRNA'){
                                       #PLOTAGGREGATION(fnames, subset(genes, intergenic_genic=='intergenic'), 'Intergenic : ', 0)
                                       #PLOTAGGREGATION(fnames, subset(genes, intergenic_genic=='genic'), 'Genic : ', 0)
                                       PLOTAGGREGATION3(fnames, genes, 0)
                          }
             }
             print(type)
}
dev.off()
