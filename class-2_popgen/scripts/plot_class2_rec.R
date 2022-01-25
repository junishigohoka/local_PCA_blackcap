library("optparse")
options(show.error.locations = TRUE) 
option_list = list(
       make_option(c("--dirout"), type="character", default=NULL, 
              help="Path to dirout", metavar="character"),
       make_option(c("--dirlist"), type="character", default=NULL, 
              help="Path to dirlist", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dirout<-opt$dirout
dirlist<-opt$dirlist


library("RColorBrewer")



# Define colours
cols<-brewer.pal(8,"Set1")
popcols<-c("black",rainbow(9))

setwd(dirlist)
# Read chromosome length
chrlen<-read.table("chromosomes_length.list",col.names=c("chr","len"))
chr<-"chr_21"
len<-chrlen[chrlen$chr==chr,]$len

# Read local PCA outlier interval
bed<-read.table("local_PCA_MDS_outlier.bed",col.names=c("chr","from","to"))
from<-bed[bed$chr=="chr_21",]$from
to<-bed[bed$chr=="chr_21",]$to




setwd(dirout)
png("class2_rec.png",width=480*4,height=480*6,res=350)

par(mfrow=c(2,1),
    las=1
)

# Comparison of recombination rate
d<-read.table("chr_21_rec_rate_mean_10kb_comp.tsv",header=T)
colnames(d)[4]<-"cont_resident"
colnames(d)[5]<-"cont_medlong"
d<-d[order(d$pos1),]


plot(1,
     type='n',
     xlim=c(0,len)/1e6,
     ylim=c(0,30),
     xlab="Position [Mb]",
     ylab="Recombination rate [cM/Mb]",
     main="chr_21"
)


lines(d$pos1/1e6,d$cont_resident,
      #pch=21,
      #bg=adjustcolor(rainbow(9)[3],alpha.f=0.3),
      type='s',
      col=rainbow(9)[3]
)
lines(d$pos1/1e6,d$cont_medlong,
      #pch=21,
      #bg=adjustcolor(rainbow(9)[3],alpha.f=0.3),
      type='s',
      col=rainbow(9)[1]
)
lines(d$pos1/1e6,d$CapeVerde,
      #pch=21,
      #bg=adjustcolor(rainbow(9)[7],alpha.f=0.3),
      type='s',
      col=rainbow(9)[7]
)
lines(d$pos1/1e6,d$Azores,
      #pch=21,
      #bg=adjustcolor(rainbow(9)[6],alpha.f=0.3),
      type='s',
      col=rainbow(9)[6]
)

abline(v=c(from,to)/1e6,
       lty=2,
       col=cols[8]
)
legend("topright",
       legend=c("Azores","CapeVerde","cont_medlong","cont_resident"),
       col=rainbow(9)[c(6,7,1,3)],
       lty=1,
       #pt.bg=adjustcolor(rainbow(9)[c(6,7,3)],alpha.f=0.3),
       #pch=21,
       cex=0.8
)


plot(1,
     type='n',
     xlim=c(0,len)/1e6,
     ylim=c(-40,15),
     xlab="Position [Mb]",
     ylab="Diff rec rate [cM/Mb]",
     main="chr_21"
)
abline(h=0,
       lty=2,
       col="gray"
)


lines(d$pos1/1e6,d$CapeVerde-d$cont_resident,
      lty=3,
      #bg=adjustcolor(rainbow(9)[7],alpha.f=0.3),
      type='s',
      col=rainbow(9)[7]
)
lines(d$pos1/1e6,d$Azores-d$cont_resident,
      lty=3,
      #bg=adjustcolor(rainbow(9)[6],alpha.f=0.3),
      type='s',
      col=rainbow(9)[6]
)

lines(d$pos1/1e6,d$CapeVerde-d$cont_medlong,
      lty=1,
      #bg=adjustcolor(rainbow(9)[7],alpha.f=0.3),
      type='s',
      col=rainbow(9)[7]
)
lines(d$pos1/1e6,d$Azores-d$cont_medlong,
      lty=1,
      #bg=adjustcolor(rainbow(9)[6],alpha.f=0.3),
      type='s',
      col=rainbow(9)[6]
)

abline(v=c(from,to)/1e6,
       lty=2,
       col=cols[8]
)

legend("bottomright",
       legend=c("Azores - cont_medlong","CapeVerde - cont_medlong","Azores - cont_resident","CapeVerde - cont_resident"),
       col=rainbow(9)[c(6,7,6,7)],
       lty=c(1,1,3,3),
       #pt.bg=adjustcolor(rainbow(9)[c(6,7)],alpha.f=0.3),
       #pch=21,
       cex=0.8
)


dev.off()


