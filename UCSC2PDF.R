# Script: converting UCSC region (in bed format) to PDF files
# Usage: Rscript UCSC2PDF.R bedfile
source('http://zlab.umassmed.edu/~dongx/mylib/mylib.R')

args <- commandArgs(TRUE)
bedfile=args[1]

# example of controling individual track
#theURL="http://genome.ucsc.edu/cgi-bin/hgTracks?db=mm9&wgRna=hide&cpgIslandExt=pack&ensGene=hide&mrna=hide&intronEst=hide&mgcGenes=hide&hgt.psOutput=on&cons44way=hide&snp130=hide&snpArray=hide&wgEncodeReg=hide&pix=1000&refGene=pack&knownGene=hide&rmsk=hide"
# example of controling via session
theURL="http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=sterding&hgS_otherUserSessionName=mm9_piRNA&hgt.psOutput=on&pix=1000"
library(multicore)

# read regions
toPlot=read.table(bedfile, header=F)

for(i in 1:nrow(toPlot)){
    # display x1.25 region
    screenshotUCSC(theURL, "", as.character(toPlot[i,1]), toPlot[i,2]-round(abs(toPlot[i,3]-toPlot[i,2])/8), toPlot[i,3]+round(abs(toPlot[i,3]-toPlot[i,2])/8), paste("~/scratch/region4_", as.character(toPlot[i,1]),".",toPlot[i,2],".",toPlot[i,3],".pdf", sep=""))

    ## convert to PNG
    #command=paste("convert -density 150 ~/scratch/region3_", as.character(toPlot[i,1]),".",toPlot[i,2],".",toPlot[i,3],".pdf ~/scratch/region3_", as.character(toPlot[i,1]),".",toPlot[i,2],".",toPlot[i,3], ".png", sep="")
    #cat(command,"\n")
    #try(system(command))
    #
    ## scp to zlab
    #command=paste("scp ~/scratch/region3_", as.character(toPlot[i,1]),".",toPlot[i,2],".",toPlot[i,3],".p* zlab:~/public_html/pirna/img", sep="")
    #cat(command,"\n")
    #try(system(command))

    # anti-robot version
    Sys.sleep(5) # Program-driven use of this software is limited to a maximum of one hit every 15 seconds and no more than 5,000 hits per day.
    print(paste(toPlot[i,], sep=" "));
}
