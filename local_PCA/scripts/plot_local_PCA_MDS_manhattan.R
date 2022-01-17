library(RColorBrewer)
library("optparse")
options(show.error.locations = TRUE) 
option_list = list(
       make_option(c("--input"), type="character", default=NULL, 
              help="Path to local PCA MDS txt file for whole-genome", metavar="character"),
       make_option(c("--thre"), type="character", default=NULL, 
              help="Path to local PCA MDS threshold txt file for whole-genome", metavar="character"),
       make_option(c("--bed"), type="character", default=NULL, 
              help="Path to bed file specifying outliers", metavar="character"),
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
sv<-read.table(opt$bed,col.names = c("chr","from","to"))
sv$chr<-factor(sv$chr,levels=levels(thre$chr))
sv$from<-sv$from/10^6
sv$to<-sv$to/10^6

chrlist<-read.table(opt$chrlist)
colnames(chrlist)<-c("chr","len")
chrlist$index<-1:nrow(chrlist)

cols<-brewer.pal(7,"Set1")[c(1,2,4,3)]




# Merge with chrlist
thre<-merge(chrlist,thre)
thre<-thre[order(thre$index),]


# Read local PCA results
d<-read.table(input,header=T)
d$chr<-as.character(d$chr)
d$pos<-d$pos/10^6

cexout<-1.5



# --------------------------------------------------------------------
# Plot MDS1
# --------------------------------------------------------------------
n<-4.5
png(paste0(output,"_MDS1_outliers.png"),height=480*n/3,width=480*n,res=300)

par(mfrow=c(1,3),
    mar=c(3,3,3,1),
    oma=c(2,2,0,0))

for (i in 1:nrow(thre)){
  chr<-thre$chr[i]
  high1<-thre$MDS1.high[i]
  low1<-thre$MDS1.low[i]
  high2<-thre$MDS2.high[i]
  low2<-thre$MDS2.low[i]
  mode1<-thre$MDS1.mode[i]
  mode2<-thre$MDS2.mode[i]

  sub<-d[d$chr==chr,]
  len<-thre$len[i]
  
  svsub<-sv[as.character(sv$chr)==as.character(chr),]
  
  
  
  # outliers
  ol1<-sub[sub$MDS1>high1|sub$MDS1<low1,]
  ol2<-sub[sub$MDS2>high2|sub$MDS2<low2,]
  ol12<-sub[(sub$MDS1>high1|sub$MDS1<low1)&(sub$MDS2>high2|sub$MDS2<low2),]
  
  # non-outliers
  nol<-sub[sub$MDS1<=high1&sub$MDS1>=low1&sub$MDS2<=high2&sub$MDS2>=low2,]
  
  plot(nol$pos,nol$MDS1,
       xlab="",
       ylab="",
       pch=19,
       col="gray",
       xlim=c(0,len/10^6),
       ylim=c(-1,1),
       main=chr)
  if(nrow(svsub)>0){
    for(j in 1:nrow(svsub)){
      f<-svsub$from[j]
      t<-svsub$to[j]
      polygon(c(f,t,t,f),c(-1,-1,1,1),
              col=adjustcolor(cols[4],alpha.f = 0.5),
              border=NA
              )
    }
  }
  points(ol2$pos,ol2$MDS1,
         pch=19,
         col=cols[2])
  points(ol1$pos,ol1$MDS1,
         pch=19,
         col=cols[1])
  points(ol12$pos,ol12$MDS1,
         pch=19,
         col=cols[3])
}
mtext(side=1,
      outer=T,
      text="Position [Mb]"
      )
mtext(side=2,
      outer=T,
      text="MDS1"
      )
      
dev.off()





# --------------------------------------------------------------------
# Plot MDS2
# --------------------------------------------------------------------

png(paste0(output,"_MDS2_outliers.png"),height=480*n/3,width=480*n,res=300)

par(mfrow=c(1,3),
    mar=c(3,3,3,1),
    oma=c(2,2,0,0))

for (i in 1:nrow(thre)){
  chr<-thre$chr[i]
  high1<-thre$MDS1.high[i]
  low1<-thre$MDS1.low[i]
  high2<-thre$MDS2.high[i]
  low2<-thre$MDS2.low[i]
  mode1<-thre$MDS1.mode[i]
  mode2<-thre$MDS2.mode[i]
  
  sub<-d[d$chr==chr,]
  len<-thre$len[i]
  
  svsub<-sv[as.character(sv$chr)==as.character(chr),]
  
  
  # outliers
  ol1<-sub[sub$MDS1>high1|sub$MDS1<low1,]
  ol2<-sub[sub$MDS2>high2|sub$MDS2<low2,]
  ol12<-sub[(sub$MDS1>high1|sub$MDS1<low1)&(sub$MDS2>high2|sub$MDS2<low2),]
  
  # non-outliers
  nol<-sub[sub$MDS1<=high1&sub$MDS1>=low1&sub$MDS2<=high2&sub$MDS2>=low2,]
  
  plot(nol$pos,nol$MDS2,
       xlab="",
       ylab="",
       pch=19,
       col="gray",
       xlim=c(0,len/10^6),
       ylim=c(-1,1),
       main=chr)
  if(nrow(svsub)>0){
    for(j in 1:nrow(svsub)){
      f<-svsub$from[j]
      t<-svsub$to[j]
      polygon(c(f,t,t,f),c(-1,-1,1,1),
              col=adjustcolor(cols[4],alpha.f = 0.5),
              border=NA
      )
    }
  }
  points(ol1$pos,ol1$MDS2,
         pch=19,
         col=cols[1])
  points(ol2$pos,ol2$MDS2,
         pch=19,
         col=cols[2])
  points(ol12$pos,ol12$MDS2,
         pch=19,
         col=cols[3])
}
mtext(side=1,
      outer=T,
      text="Position [Mb]"
      )
mtext(side=2,
      outer=T,
      text="MDS2"
      )

dev.off()


