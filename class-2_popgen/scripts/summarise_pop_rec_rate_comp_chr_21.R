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
cols<-brewer.pal(8,"Set1")
setwd(dirlist)

chrlist<-read.table("chromosomes_length.list",col.names=c("chr","len"))
lpca<-read.table("local_PCA_MDS_outlier.bed",col.names=c("chr","pos1","pos2"))

chr<-"chr_21"
pos1<-lpca[lpca$chr==chr,]$pos1
pos2<-lpca[lpca$chr==chr,]$pos2
len<-chrlist[chrlist$chr==chr,]$len

pops<-c("Azores","CapeVerde","ContRes")


setwd(dirout)

d1<-read.table("chr_21_Azores_win.10kb_step.10kb.mean.rec.tab",header=T)
d1<-d1[,c("pos1","cMpMb")]
d2<-read.table("chr_21_CapeVerde_win.10kb_step.10kb.mean.rec.tab",header=T)
d2<-d2[,c("pos1","cMpMb")]
d3<-read.table("chr_21_ContRes_win.10kb_step.10kb.mean.rec.tab",header=T)
d3<-d3[,c("pos1","cMpMb")]
d4<-read.table("chr_21_MeD_SW_Belgium_win.10kb_step.10kb.mean.rec.tab",header=T)
d4<-d4[,c("pos1","cMpMb")]



d<-merge(d1,d2,by="pos1")
d<-merge(d,d3,by="pos1")
d<-merge(d,d4,by="pos1")
colnames(d)<-c("pos1","Azores","CapeVerde","ContRes","medlong")
d$island<-ifelse(d$pos1>=pos1&d$pos1<pos2,"in","out")
d$island<-factor(as.factor(d$island),levels=c("out","in"))
d<-d[order(d$island),]






write.table(d,"chr_21_rec_rate_mean_10kb_comp.tsv",
            sep="\t",
            col.names=T,
            row.names=F,
            quote=F
)

