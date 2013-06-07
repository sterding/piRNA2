# script to piRNA analysis

source('http://zlab.umassmed.edu/~dongx/mylib.R')

# prerequisition
#awk '{OFS="\t";if($10==2 && $5>=3){split($11,a,",");split($12,b,",");print $1,$2+a[1],$2+b[2],$1"_"$2"_"$3"_"$6, $5,$6;}}' ~/scratch/mouse_adult_wt_RNAseq_merge/mouse_adult_wt_RNAseq_merge.junction.bed | sort -k1,1 -k2,2 -k3,3 -k6,6 -k4,4 > ../data/mouse_adult_wt_RNAseq_merge.junctions.size.tab
#awk 'BEGIN{id=""}{OFS="\t";if(id!=$1"."$2"."$3"."$6){if(id!="")print chr,s,e,name,rN,str,sN; delete sArray; rN=1; sN=1; sArray[$4]; chr=$1; s=$2; e=$3 ;str=$6; name=$4; id=$1"."$2"."$3"."$6;} else{rN++; if(!($4 in sArray)){sArray[$4]; sN++;name=name";"$4;} }}END{print chr,s,e,name,rN,str,sN;}'  ../data/mouse_adult_wt_RNAseq_merge.junctions.size.tab > ../data/mouse_adult_wt_RNAseq_merge.junctions.merged.tab &

junctions=read.table("../data/mouse_adult_wt_RNAseq_merge.junctions.merged.tab")
colnames(junctions)=c('chr','start','end','name','readsCount','strand','speciesCount')
junctions=data.frame(junctions,length=junctions$end-junctions$start)

# define positive junction: readCount > 10 & speciesCount>1
junctions=subset(junctions, readsCount>10 & speciesCount>1)


library(ggplot2)
ggplot(junctions[junctions$length<50000,], aes(factor(readsCount),length))+ geom_boxplot()

piRNAs=read.table("../data/piRNA.clusters.coordinates.cpg.bed", header=T)
rownames(piRNAs)=piRNAs[,4] # TransID
piRNAs=piRNAs[,c('chr','start','end','name','normalizedCpG','strand')]

genes=read.table("/home/dongx/nearline/genomes/mm9/Annotation/Genes/NCBIM37.biomart67.transcripts.cpg.tab", header=T)
rownames(genes)=genes[,6] # TransID
genes=genes[,c('chr','start','end','trans_id','normalizedCpG','strand','gene_id', 'gene_type', 'gene_status', 'gene_name', 'trans_type', 'trans_status')]

#piRNA
# piRNA cluster that are not located within a gene
piRNA=bedTools.2in("intersectBed",piRNAs,genes[genes$gene_type=='protein_coding',], "-wa -s -v")
res=bedTools.2in("intersectBed",junctions,piRNAs[grep("^pi",piRNAs$name, invert=T),],"-wa -s -f 0.9")
colnames(res)=c('chr','start','end','name','readsCount','strand','speciesCount', 'length')
res=data.frame(res, type='piRNA')
for(genetype in c('protein_coding','lincRNA')){
             gene=genes[genes$trans_type==genetype,]
             if(genetype!='protein_coding') gene=bedTools.2in("intersectBed",gene,genes[genes$gene_type=='protein_coding',], "-wa -s -v")
             gene=data.frame(bedTools.2in("intersectBed",junctions,gene,"-wa -s -f 0.9"), genetype)
             colnames(gene)=colnames(res)=c('chr','start','end','name','readsCount','strand','speciesCount', 'length','type')
             res=rbind(res, gene)
}

res$type=factor(res$type, levels=c("protein_coding","lincRNA","piRNA"))

write.table(subset(res, readsCount>10 & speciesCount>2), "../result/intronSize_rCount10_sCount2.xls", quote = F, sep ="\t", col.names=NA)


png("../result/intronSize2.png", width=1200, height=900)
par(mfrow=c(3,1), mar=c(4,6,4,4))
for(ty in c("protein_coding","lincRNA","piRNA")){
             res1=subset(res,type==ty & readsCount>10 & speciesCount>8)
             symbols((res1$length), res1$speciesCount, circles=sqrt(res1$readsCount/ pi ), inches=0.35, fg="white", bg="red", cex.lab=2, cex.main=2, cex.axis=2, xlab="intron size (bp)", ylab="species count", main=ty, xlim=range((res$length)))
}
dev.off()

png("../result/intronSize_log.png", width=1200, height=900)
par(mfrow=c(3,1), mar=c(4,6,4,4))
for(ty in c("protein_coding","lincRNA","piRNA")){
             res1=subset(res,type==ty & readsCount>10 & speciesCount>8)
             symbols(log(res1$length), res1$speciesCount, circles=sqrt(res1$readsCount/ pi ), inches=0.35, fg="white", bg="red", cex.lab=2, cex.main=2, cex.axis=2, xlab="intron size (log(bp))", ylab="species count", main=ty, xlim=range(log(res$length)))
}
dev.off()


res2=res[res$readsCount>10,];
png("../result/intronSize.png", width=600, height=600)
boxplot(length ~ type, res2, ylab='Intron size (bp)', outline=F)
text(x=1, y=max(boxplot.stats(res2$length[res2$type=='protein_coding'])$stat), paste("n",sum(res2$type=='protein_coding'),sep="="), pos=3, offset=.2, cex=.8)
text(x=2, y=max(boxplot.stats(res2$length[res2$type=='lincRNA'])$stat), paste("n",sum(res2$type=='lincRNA'),sep="="), pos=3, offset=.2, cex=.8)
text(x=3, y=max(boxplot.stats(res2$length[res2$type=='piRNA'])$stat), paste("n",sum(res2$type=='piRNA'),sep="="), pos=3, offset=.2, cex=.8)
#wilcox.test(length ~ type, data=res[res$readsCount>10,], subset = type %in% c('protein_coding', 'piRNA'))
#ggplot(res2, aes(factor(type), length)) + geom_boxplot() + geom_jitter()
dev.off()
