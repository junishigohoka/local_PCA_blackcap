library(RColorBrewer)
library("optparse")
 
option_list = list(
       make_option(c("--input"), type="character", default=NULL, 
              help="Path to local PCA MDS txt file for whole-genome", metavar="character"),
       make_option(c("--thre"), type="character", default=NULL, 
              help="Path to local PCA MDS threshold txt file for whole-genome", metavar="character"),
	make_option(c("--output"), type="character", default="out.txt", 
              help="Output file name", metavar="character"),
	make_option(c("--chrlist"), type="character", default="out.txt", 
              help="Chromosome list ", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

input<-opt$input
output<-opt$output
thre<-read.table(opt$thre,header=T)
chrlist<-read.table(opt$chrlist)
colnames(chrlist)<-c("chr")
chrlist$index<-1:nrow(chrlist)

cols<-c(brewer.pal(7,"Set1")[c(1,2,4)],"gray30")
# transparent
tcols<-adjustcolor(cols,alpha.f = 0.3)


# Merge with chrlist
thre<-merge(chrlist,thre)
thre<-thre[order(thre$index),]

# Sort by length
# thre<-thre[order(thre$len,decreasing = T),]


# Read local PCA results
d<-read.table(input,header=T)
cexout<-1.5


# --------------------------------------------------------------------
# Plot MDS1 vs MDS2
# --------------------------------------------------------------------
n<-3
png(paste0(output,"_MDS1-vs-MDS2_outliers.png"),height=480*n,width=480*n*3,res=350)

par(mfrow=c(1,3),
    oma=c(3,3,3,3),
    mar=c(2,2.5,3,2.5))

for (i in 1:nrow(thre)){
       chr<-thre$chr[i]

       high1<-thre$MDS1.high[i]
       low1<-thre$MDS1.low[i]
       high2<-thre$MDS2.high[i]
       low2<-thre$MDS2.low[i]
       mode1<-thre$MDS1.mode[i]
       mode2<-thre$MDS2.mode[i]

       sub<-d[d$chr==chr,]
       # len<-thre$len[i]

       # outliers
       ol1<-sub[sub$MDS1>high1|sub$MDS1<low1,]
       ol2<-sub[sub$MDS2>high2|sub$MDS2<low2,]
       ol12<-sub[(sub$MDS1>high1|sub$MDS1<low1)&(sub$MDS2>high2|sub$MDS2<low2),]

       # non-outliers
       nol<-sub[sub$MDS1<=high1&sub$MDS1>=low1&sub$MDS2<=high2&sub$MDS2>=low2,]

       plot(nol$MDS1,nol$MDS2,
              xlab="",
              ylab="",
              pch=16,
              col=tcols[4],
              las=1,
              xlim=c(-0.75,0.75),
              ylim=c(-0.75,0.75),
              #xlim=c(min(sub$MDS1),max(sub$MDS1)),
              #ylim=c(min(sub$MDS2),max(sub$MDS2)),
              main=chr)
       points(ol1$MDS1,ol1$MDS2,
              pch=16,
              col=cols[1])
       points(ol2$MDS1,ol2$MDS2,
              pch=16,
              col=cols[2])
       points(ol12$MDS1,ol12$MDS2,
              pch=16,
              col=cols[3])
       # Indicators of modes of MDS1 and MDS2
       abline(v=mode1,
              col=cols[1],
              lty=3)
       abline(h=mode2,
              col=cols[2],
              lty=3)
  }
  
mtext(side=1,
      outer=T,
      text="MDS1",
      line=1,
      cex=cexout)
mtext(side=2,
      outer=T,
      text="MDS2",
      line=1,
      cex=cexout)

dev.off()



