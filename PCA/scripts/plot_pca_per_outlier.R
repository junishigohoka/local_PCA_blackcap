options <- commandArgs(trailingOnly = TRUE)
library("optparse")
 
option_list = list(
       make_option(c("--dirpca"), type="character", default=NULL, 
              help="Path to PCA output per outlier", metavar="character"),
       make_option(c("--poplist"), type="character", default=NULL, 
              help="Population list", metavar="character"),
       make_option(c("--chrlist"), type="character", default=NULL, 
              help="Chromosome list", metavar="character"),
    	make_option(c("--svlist"), type="character", default=NULL, 
              help="SV list in BED format", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

input<-opt$input
poplist<-opt$poplist
svlist<-read.table(opt$svlist,col.names = c("chr","from","to"))
chrlist<-read.table(opt$chrlist,col.names="chr")
chrlist$index<-1:nrow(chrlist)
dirpca<-opt$dirpca

# Merge svlist with chrlist
svlist<-merge(svlist,chrlist,by="chr")
svlist<-svlist[order(svlist$index,svlist$from),]

# Read the information of population and phenotype
info<-read.table(poplist,header=T,sep="\t")
info<-info[info$species=="Blackcap",]
# info$population<-factor(info$population,levels=c("continent_medlong","continent_short","continent_resident","Azores_resident","CapeVerde_resident","Canary_resident","Madeira_resident","Mallorca_resident","Crete_resident"))
for(i in 2:ncol(info)){
  info[,i]<-as.factor(as.character(info[,i]))
}

info$popforpca<-as.character(info$population)

info$popforpca[info$distance%in%c("med","long")]<-"cont_medlong"
info$popforpca[info$distance=="short"]<-"cont_short"
info$popforpca[info$distance=="resident"&info$land=="continent"]<-"cont_resident"
info$popforpca[info$popforpca=="UK"]<-"cont_medlong"


info$popforpca<-as.factor(info$popforpca)
info$popforpca<-factor(info$popforpca,levels=c("cont_medlong","cont_short","cont_resident","Azores","CapeVerde","Canary","Madeira","Mallorca","Crete"))



# Read the PCA result
setwd(dirpca)
library(RColorBrewer)
# cols<-brewer.pal(7,"Set1")
cols<-rainbow(nlevels(info$popforpca))

#pdf("PCA_localPCA_outlier.pdf",height=5,width=15)
png("PCA_localPCA_outlier.png",height=480*2,width=480*2*3,res=350)
par(mfrow=c(1,3),
    mar=c(4,4,2,2),
    oma=c(0,0,3,0))

for (i in 1:nrow(svlist)){
  chr<-svlist$chr[i]
  f<-svlist$from[i]
  t<-svlist$to[i]
  
  prefix=paste(chr,f,t,sep="_")
  
  
  pca<-read.table(paste(prefix,".eigenvec",sep=""))
  pca<-pca[,-2]
  colnames(pca)<-c("ID",paste("PC",1:20,sep=""))
  
  # Merge
  d<-merge(pca,info,by="ID")
  
  
  for(j in 1){
    for(k in 2){
      plot(d[,j+1],d[,k+1],
          pch=c(1,3)[d$land],
          cex=1.5,
          xlab=paste0("PC",j),
          ylab=paste0("PC",k),
          #col=rainbow(nlevels(d$distance))[d$distance],
          col=cols[d$popforpca],
          main=paste(chr,f,"-",t))
      
    }
} 


}

legend("bottomleft",
       legend=levels(d$popforpca),
       col=cols,
       pch=c(1,1,1,3,3,3,3,3,3),
       bty='n',
       cex=0.75,
       ncol=1
)
    
  
dev.off()


