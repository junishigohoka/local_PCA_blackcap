library("circlize")
library("RColorBrewer")

library("optparse")
 
option_list = list(
  make_option("--chrlen", type="character", default=NULL, 
              help="Path to chromosome length list for the species to compare", metavar="character"),
  make_option("--satsuma", type="character", default=NULL, 
              help="Path to chromosome length list for the species to compare", metavar="character"),
  make_option("--sppname", type="character", default=NULL, 
              help="Species name", metavar="character"),
  make_option("--dirout", type="character", default=NULL, 
              help="output direcotry", metavar="character"),
  make_option("--dirlist", type="character", default=NULL, 
              help="output direcotry", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

spname<-opt$sppname
#pastel<-brewer.pal(9,"Pastel1")
cols<-adjustcolor(brewer.pal(7,"Set1"))


`%notin%`<-Negate(`%in%`)

# Read chromosome length and localPCA bed file
dirlist<-opt$dirlist
setwd(dirlist)
## spp
chrlen<-opt$chrlen
sppchrlist<-read.table(chrlen,col.names=c("chr","len"))
sppchrlist$chr<-paste0("chr_",sppchrlist$chr)
### Remove LGs
sppchrlist$chr<-factor(as.factor(sppchrlist$chr),levels=sppchrlist$chr)


## BC
bcchrlist<-read.table("BC_chromosomes_length.list",col.names=c("chr","len"))
bcchrlist$chr<-factor(as.factor(bcchrlist$chr),levels=bcchrlist$chr)


# Local PCA
bed<-read.table("local_PCA_MDS_outlier.bed",col.names=c("chr","from","to"))
## Refactor the chromosome names in bed
bed$chr<-paste0("BC_",bed$chr)
bed$chr<-factor(as.factor(bed$chr),levels=paste0("BC_",levels(bcchrlist$chr)))

# BC SVs
svbed<-read.table("local_PCA_MDS_class-1.bed",col.names=c("chr","from","to"))
svbed$chr<-paste0("BC_",svbed$chr)
svbed$chr<-factor(as.factor(svbed$chr),levels=paste0("BC_",levels(bcchrlist$chr)))



# Prepare chromosome info
chrinfo<-data.frame(chr=c(paste0("BC_",bcchrlist$chr),rev(paste0(paste0(spname,"_",sppchrlist$chr)))), 
                    start=c(rep(0,nrow(bcchrlist)),rep(0,nrow(sppchrlist))),
                    end=c(bcchrlist$len,rev(sppchrlist$len))
                    )
chrinfo$chr<-factor(chrinfo$chr,levels=chrinfo$chr)


# Read satsuma results
satsuma<-opt$satsuma
satsuma<-read.table(satsuma,col.names = c("spp.chr","spp.pos1","spp.pos2","bc.chr","bc.pos1","bc.pos2"))
satsuma[,1]<-gsub("^",paste0(spname,"_chr_"),satsuma[,1])
satsuma[,4]<-gsub("chr","BC_chr",satsuma[,4])
# filter
satsuma<-satsuma[satsuma[,4]%in%paste0("BC_",bcchrlist$chr),]
satsuma$spp.chr<-factor((satsuma$spp.chr),levels = paste0(paste0(spname,"_"),levels(sppchrlist$chr)))
satsuma$bc.chr<-factor(as.factor(satsuma$bc.chr),paste0("BC_",levels = levels(bcchrlist$chr)))
satsuma<-na.omit(satsuma)

# Subset Satsuma links within the blackcap local PCA outlier intervals
satsuma.sub.bc<-data.frame(matrix(ncol=ncol(satsuma),nrow=0))
colnames(satsuma.sub.bc)<-colnames(satsuma)
for(i in 1:nrow(bed)){
  chr<-bed$chr[i]
  pos1<-bed$from[i]
  pos2<-bed$to[i]
  # Subset only focal chromosome
  #tmp<-satsuma[satsuma$bc.chr==chr&satsuma$bc.pos1>pos1&satsuma$bc.pos2<pos2,]
  tmp <- satsuma[satsuma$bc.chr==chr,]
  tmp<-satsuma[(satsuma$bc.chr==chr&satsuma$bc.pos1>pos1&satsuma$bc.pos2<pos2)|(satsuma$bc.chr==chr&satsuma$bc.pos1<=pos1&satsuma$bc.pos2>pos1)|(satsuma$bc.chr==chr&satsuma$bc.pos2>=pos2&satsuma$bc.pos1<pos2),]
  # Merge proximal links
  tmpmg<-tmp[1,]
  if(nrow(tmp)>0){
          for(j in 2:nrow(tmp)){
            if((abs(tmpmg[nrow(tmpmg),3]-tmp[j,2])<100000)&abs(tmpmg[nrow(tmpmg),6]-tmp[j,5])<100000&(tmpmg[nrow(tmpmg),1]==tmp[j,1])&tmpmg[nrow(tmpmg),4]==tmp[j,4]){
              tmpmg[nrow(tmpmg),2]<-min(tmpmg[nrow(tmpmg),2],tmp[j,2])
              tmpmg[nrow(tmpmg),3]<-max(tmpmg[nrow(tmpmg),3],tmp[j,3])
              tmpmg[nrow(tmpmg),5]<-min(tmpmg[nrow(tmpmg),5],tmp[j,5])
              tmpmg[nrow(tmpmg),6]<-max(tmpmg[nrow(tmpmg),6],tmp[j,6])
            }else{
              tmpmg<-rbind(tmpmg,tmp[j,])
            }
          }
  }
  

  satsuma.sub.bc<-rbind(satsuma.sub.bc,tmpmg) 
  rm(chr,pos1,pos2,tmp,tmpmg)
}


## Subset satsuma.sub.bc further, keeping links with spp inversions
satsuma.bc.spp<-data.frame(matrix(ncol=ncol(satsuma),nrow=0))
for(i in 1:nrow(svbed)){
  chr<-svbed$chr[i]
  pos1<-svbed$from[i]
  pos2<-svbed$to[i]
  # Subset only focal chromosome
  satsuma.bc.spp<-rbind(satsuma.bc.spp,
                       satsuma.sub.bc[(satsuma.sub.bc$bc.chr==chr&satsuma.sub.bc$bc.pos1>pos1&satsuma.sub.bc$bc.pos2<pos2)|(satsuma.sub.bc$bc.chr==chr&satsuma.sub.bc$bc.pos1<=pos1&satsuma.sub.bc$bc.pos2>pos1)|(satsuma.sub.bc$bc.chr==chr&satsuma.sub.bc$bc.pos2>=pos2&satsuma.sub.bc$bc.pos1<pos2),])
}
satsuma.bc.spp<-na.omit(satsuma.bc.spp)

# Randomly sample links from raw satsuma
#satsuma.ran<-satsuma[sample(nrow(satsuma),500),]
satsuma.ran<-satsuma[seq(1,nrow(satsuma),length.out=1000),]

## PLOT #####################################################################################
# Prepare label
chrlab<-chrinfo
# chrlab$lab<-NA
# for(i in 1:nrow(chrlab)){
#   if(chrlab$end[i]>13000000){
#     chrlab$lab[i]<-gsub(".*chr_","",chrlab$chr)
#   }else{
#     chrlab$lab[i]<-" "
#   }
# }
# 

dirout<-opt$dirout
setwd(dirout)

png(paste0("BC_",gsub(" ","_", spname),"_satsuma.png"),height=480*3,width=480*6,res=350)
#pdf("BC_spp_satsuma.pdf",width=10,height=5)
par(mfrow=c(1,2),las=1)

# Plot 1 (synteny only)
circos.par(gap.after=c(rep(1,nrow(bcchrlist)-1),10,rep(1,nrow(sppchrlist)-1),10),
           start.degree=85)
circos.genomicInitialize(chrinfo, plotType = NULL)

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  # circos.axis(h="bottom",labels.cex=0.4,niceFacing=T)
  circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2), 
              facing = "clockwise",
              gsub(".*chr_", "", CELL_META$sector.index),
              cex = 0.6, niceFacing = TRUE)
  }, 
  track.height = mm_h(1), 
  cell.padding = c(0, 0, 0, 0),
  bg.border="gray",
  bg.col="gray"
  # bg.border = NA
  )

# circos.genomicLink(bundle[,1:3],bundle[4:6],
#                    col=rep(cols1,length=nrow(bundle)))
circos.genomicLink(satsuma.ran[,1:3],
                   satsuma.ran[,4:6],
                   col=rainbow(nlevels(satsuma.ran$bc.chr))[satsuma.ran$bc.chr]
                   )

legend("topleft",
       legend=paste0(toupper(substr(spname,1,1)),substr(spname,2,nchar(spname))),
       cex=0.75,
       bty='n')
legend("topright",
       legend="Blackcap",
       cex=0.75,
       bty='n')
circos.clear()


# Plot2 (localPCA)
circos.par(gap.after=c(rep(1,nrow(bcchrlist)-1),10,rep(1,nrow(sppchrlist)-1),10),
           start.degree=85)
circos.genomicInitialize(chrinfo, plotType = NULL)

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  # circos.axis(h="bottom",labels.cex=0.4,niceFacing=T)
  circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2), 
              gsub(".*chr_", "", CELL_META$sector.index), 
              cex = 0.6, 
              facing = "clockwise",
              niceFacing = TRUE)
  }, 
  track.height = mm_h(1), 
  cell.padding = c(0, 0, 0, 0),
  bg.border="gray",
  bg.col="gray"
  # bg.border = NA
  )

# circos.genomicLink(bundle[,1:3],bundle[4:6],
#                    col=rep(cols1,length=nrow(bundle)))
circos.genomicLink(satsuma.sub.bc[,1:3],satsuma.sub.bc[,4:6],
                   col="gray"
                   )

satsuma.bc.spp$bc.chr<-factor(satsuma.bc.spp$bc.chr)
circos.genomicLink(satsuma.bc.spp[,1:3],satsuma.bc.spp[4:6],
                   col=cols[satsuma.bc.spp$bc.chr]
                   )

legend("topleft",
       legend=paste0(toupper(substr(spname,1,1)),substr(spname,2,nchar(spname))),
       cex=0.75,
       bty='n')
legend("topright",
       legend="Blackcap",
       cex=0.75,
       bty='n')
circos.clear()
circos.clear()
dev.off()

