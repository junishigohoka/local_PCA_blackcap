
library("optparse")
 
option_list = list(
  make_option(c("--input"), type="character", default=NULL, 
              help="Path to local PCA MDS txt file for whole-genome", metavar="character"),
	make_option(c("--output"), type="character", default="out.txt", 
              help="Output file name", metavar="character"),
	make_option(c("--chrlist"), type="character", default="out.txt", 
              help="Chromosome list ", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

input<-opt$input
output<-opt$output
chrlist<-read.table(opt$chrlist)
colnames(chrlist)<-c("chr")


# colnames(chrlist)<-c("chr","len")
# chrlist<-chrlist[order(chrlist$len,decreasing = T),]


# Read whole genome results
d<-read.table(input,header=T)
n<-6



# Plot

# --------------------------------------------------------------------------
# MDS1
# --------------------------------------------------------------------------
png(paste0(output,"_MDS1_distribution.png"),height=480*n/3,width=480*n,res=400)
par(mfrow=c(1,3),mar=c(2,2,3.5,0.5),
    oma=c(3,3,3,3))

# Make data frame to store the outlier thresholds for each scaffold
thre.mds1<-data.frame(matrix(ncol=4,nrow=nrow(chrlist)))
colnames(thre.mds1)<-c("chr","MDS1.mode","MDS1.high","MDS1.low")

for(i in 1:nrow(chrlist)){
  
  chr<-chrlist$chr[i]
  if(chr!="chr_W"&&chr!="chr_w_unkown"){

       sub<-d[d$chr==chr,]
       
       his<-hist(sub$MDS1,breaks=seq(-1,1,length.out = 1001),
              col="gray",
              border="gray",
              main=chr,
              xlab="",
              ylab="")
       # Define mode as the mean of the bin mid positions with the highest frequency
       mode<-mean(his$mids[his$counts==max(his$counts)])
       thre.mds1[i,]<-c(as.character(chr),mode,mode+0.3,mode-0.3)
       
       
       abline(v=mode,
              col="red",
              lty=3)
       abline(v=mode+0.3,
              col="blue",
              lty=3)
       abline(v=mode-0.3,
              col="blue",
              lty=3)
       }
}

# Define mode as the mean of the bin mid positions with the highest frequency
mode<-mean(his$mids[his$counts==max(his$counts)])
abline(v=mode,
       col="red",
       lty=3)
abline(v=mode+0.3,
       col="blue",
       lty=3)
abline(v=mode-0.3,
       col="blue",
       lty=3)

mtext(side=1,
      outer=T,
      text="MDS1"
      )
mtext(side=2,
      outer=T,
      text="Frequency"
)
dev.off()


# --------------------------------------------------------------------------
# MDS2
# --------------------------------------------------------------------------
png(paste0(output,"_MDS2_distribution.png"),height=480*n/3,width=480*n,res=400)
par(mfrow=c(1,3),mar=c(2,2,3.5,0.5),
    oma=c(3,3,3,3))

# Make data frame to store the outlier thresholds for each scaffold
thre.mds2<-data.frame(matrix(ncol=4,nrow=nrow(chrlist)))
colnames(thre.mds2)<-c("chr","MDS2.mode","MDS2.high","MDS2.low")
for(i in 1:nrow(chrlist)){
  chr<-chrlist$chr[i]
  
  if(chr!="chr_W"&&chr!="chr_w_unkown"){

  sub<-d[d$chr==chr,]
  
  his<-hist(sub$MDS2,breaks=seq(-1,1,length.out = 1001),
            col="gray",
            border="gray",
            main=chr,
            xlab="",
            ylab="")
  # Define mode as the mean of the bin mid positions with the highest frequency
  mode<-mean(his$mids[his$counts==max(his$counts)])
  thre.mds2[i,]<-c(as.character(chr),mode,mode+0.3,mode-0.3)
  abline(v=mode,
         col="red",
         lty=3)
  abline(v=mode+0.3,
         col="blue",
         lty=3)
  abline(v=mode-0.3,
         col="blue",
         lty=3)
  }
}

# Define mode as the mean of the bin mid positions with the highest frequency
mode<-mean(his$mids[his$counts==max(his$counts)])
abline(v=mode,
       col="red",
       lty=3)
abline(v=mode+0.3,
       col="blue",
       lty=3)
abline(v=mode-0.3,
       col="blue",
       lty=3)

mtext(side=1,
      outer=T,
      text="MDS2"
)
mtext(side=2,
      outer=T,
      text="Frequency"
)
dev.off()


# -------------------------------------------------------------------------
# Threshold
# -------------------------------------------------------------------------
# Merge data frames
mge<-merge(thre.mds1,thre.mds2)

write.table(mge,paste0(output,"threshold.txt"),
            quote=F,
            row.names = F)
