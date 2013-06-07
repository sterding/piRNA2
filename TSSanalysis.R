
piRNA=read.table("../data/TSS_table_piRNA_adult_mm9.tab", header=T, sep = "\t",fill=TRUE)
rownames(piRNA)=piRNA[,4] # TransID


# intergenic vs. geneic
# CpG score
pdf("piRNA.intergenic.vs.genic.cpg.pdf", width=3, height=6)
boxplot(normalizedCpG ~ intergenic_genic, piRNA, xaxt='n',ylab="normalized CpG score", col=c('gray','white'))
axis(1, at=c(1:2), sapply(levels(piRNA$intergenic_genic), function(x) paste(x, paste("(n=",sum(piRNA$intergenic_genic==x),")",sep=""),sep="\n")), padj=0.5)
dev.off()
# TATA-less percentage
pdf("piRNA.intergenic.vs.genic.TATA.pdf", width=3, height=6)
p1=1-sum(is.na(piRNA$TATA[piRNA$intergenic_genic=="genic"]))/sum(piRNA$intergenic_genic=="genic")
p2=1-sum(is.na(piRNA$TATA[piRNA$intergenic_genic=="intergenic"]))/sum(piRNA$intergenic_genic=="intergenic")
barplot(c(p1,p2)*100, names.arg=c("genic", "intergenic"), col=c('gray','white'), ylab="Percentage of TATA-containing promoter (%)")
dev.off()
# motif

# PAS motif

## BIP: pc vs. piRNA
BIPrest=read.table('../data/x100910.rnaseq.transcripts.ALL.bed.pcpc.BIP.cpg.bed')
BIPpi=read.table("../data/TSS_table_piRNA_adult_mm9.tab.BIP.cpg.bed")
pdf('BIP.piRNA_vs_BIP.pc.pdf', width=5, height=6);
df=rbind(BIPpi[,c(8,9)], BIPrest[,c(8,9)]);
colnames(df)=c('Category', 'normalizedCpG')
df$Category=factor(df$Category, levels=c("pi:pi","pi:nc","pi:pc", "pc:pc"))
boxplot(normalizedCpG~Category, df,xaxt='n', ylab='normalized CpG score');
axis(1, at=c(1:4), sapply(levels(df$Category), function(x) paste(x, paste("(n=",sum(df$Category==x),")",sep=""),sep="\n")), padj=0.5)
dev.off()
## BIP length
pdf('BIP.piRNA_vs_BIP.pc.length.pdf', width=10, height=6);
par(mfrow=c(1,2))
df=rbind(cbind(BIPpi$V3-BIPpi$V2, BIPpi$V8), cbind(BIPrest$V3-BIPrest$V2, BIPrest$V8));
colnames(df)=c('BIP_length','Category')
df$Category=factor(df$Category, levels=c("pi:pi","pi:nc","pi:pc", "pc:pc"))
boxplot(BIP_length~Category, df,xaxt='n', ylab='BIP length (bp)');
axis(1, at=c(1:4), sapply(levels(df$Category), function(x) paste(x, paste("(n=",sum(df$Category==x),")",sep=""),sep="\n")), padj=0.5)
hist(BIPpi[,3], freq=F, breaks=50, xlim=c(0,800),col='red', main='', xlab='BIP length(bp)');
lines(density(BIPpi[,3], bw=25), col='green')
hist(BIPrest[,3], breaks=50,freq=F, add=T, col=rgb(.5, .5, .5, 0.5))
lines(density(BIPrest[,3], bw=25), col='black');
legend('topright', c('piRNA BIP', 'protein-coding BIP'), col=c('red','gray'), pch=15);
dev.off()

# nucleotide frequency
# motif

## T2T


