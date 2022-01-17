library("optparse")


option_list = list(
  make_option(c("--dirin"), type="character", default=NULL, 
              help="Path to input bcf directory", metavar="character"),
  make_option(c("--dirsites"), type="character", default=NULL, 
              help="Path to input sites directory", metavar="character"),
  make_option(c("--dirpca"), type="character", default=NULL, 
              help="Path to output directory", metavar="character"),
  make_option(c("--dirmds"), type="character", default=NULL, 
              help="Path to output directory", metavar="character"),
  make_option(c("--dirlist"), type="character", default=NULL, 
              help="Path to list directory", metavar="character"),
  make_option(c("--chr"), type="character", default=NULL, 
              help="chromosome", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dirin <- opt$dirin
dirsites <- opt$dirsites
dirpca <- opt$dirpca
dirmds<-opt$dirmds
dirlist<-opt$dirlist
chr<-opt$chr




library(lostruct)

# Read Super-Scaffold.list
# setwd(dirlist)
# chrlist<-read.table("Super-Scaffold_length.list")
# colnames(chrlist)<-c("superscaffold","length")
# chrlist<-chrlist[order(-chrlist$length),]

# Window size
win<-1000

# loca PCA
# pdf(paste(dirout,"/","local_PCA_MDS.pdf",sep=""),width=9,height=3)
# par(mfrow=c(1,3))

setwd(dirin)

# for (i in 1:nrow(chrlist)){
  # chr<-chrlist$superscaffold[i]
  
  # snps<-vcf_windower(paste(chr,".bcf",sep=""),size=win,type='snp')
  # snps<-read_tped(paste(chr,".tped",sep=""))
  # snps<-read_vcf(paste(chr,".vcf.gz",sep=""))
  snps<-as.matrix(read.table(paste0(chr,".table")))
  sites<-read.table(paste(dirsites,"/",chr,".sites.list",sep=""))
  sites<-as.vector(sites$V1)
  
  # Local PCA
  eigenstuff <- eigen_windows(snps, win=win, k=3)
  
 # detect NA
  nas<-which(is.na(eigenstuff))
  if (length(nas)>0){
    # Remove NA from eigenstuff
    eigenstuff <- eigenstuff[-c(nas),]
    # Remove sites of NA from sites list
    sites<-sites[-c(nas)]
  }
  
  
  
  
  # Dissimilarity matrix of windows
  windist <- pc_dist( eigenstuff, npc=3 )
  
  # MDS on local PCA
  mds <- cmdscale( windist, eig=TRUE, k=3 )
  
  # Coordinates of windows
  coor<-sites[seq(1,length(sites),win)]
  coorright<-sites[seq(win,length(sites),win)]



  if(length(coor)!=length(coorright)){
    coorright<-c(coorright,sites[length(sites)])
  }

  if(length(coor)>nrow(eigenstuff)){
    while(length(coor)>nrow(eigenstuff)){
      coor<-coor[-length(coor)]
      coorright<-coorright[-length(coorright)]
    }
  }

  if(length(coor)<nrow(eigenstuff)){
    while(length(coor)<nrow(eigenstuff)){
      eigenstuff<-eigenstuff[-(nrow(eigenstuff)),]
    }
  }


  # if(length(coor)!=nrow(eigenstuff)){
  #    coor<-coor[-length(coor)]
  #    coorright<-coorright[-length(coorright)]
  #    }
  # Plot
  # pdf(paste(dirout,"/",chr,"_local_PCA_MDS.pdf",sep=""),width=9,height=3)
  # par(mfcol=c(1,3))
  # plot(mds$points[,1],mds$points[,2],
  #    col=rainbow(2*length(coor)),
  #    xlab="MDS coordinate 1",
  #    ylab="MDS coordinate 2",
  #    main=chr)
  # plot(coor,mds$points[,1],
  #    xlab="Position [bp]",
  #    ylab="MDS coordinate 1",
  #    main="MDS coordinate 1",
  #    col=rainbow(2*length(coor)))
  # 
  # plot(coor,mds$points[,2],
  #      xlab="Position [bp]",
  #      ylab="MDS coordinate 2",
  #      main="MDS coordinate 2",
  #      col=rainbow(2*length(coor)))
  #dev.off()
  
  
  # Make output table
  localPCA<-as.data.frame(eigenstuff)
  localPCA$pos<-coor
  localPCA$posright<-coorright
  localPCA$chr<-chr
  localPCA<-localPCA[,c((ncol(localPCA)-c(0,2,1)),1:(ncol(localPCA)-3))]
  
  write.table(localPCA,paste(dirpca,"/localPCA_",chr,".txt",sep=""),
              quote=F,
              row.names = F)

  # MDS output
  out<-mds$points
  out<-as.data.frame(out)
  colnames(out)<-c("MDS1","MDS2","MDS3")
  out<-cbind(localPCA[,1:3],out)
  
  
  write.table(out,paste(dirmds,"/local_PCA_MDS_",chr,".txt",sep=""),
              row.names = F,
              quote=F)
  # write.csv(localPCA,paste(dirout,"/locaPCA_",chr,".csv",sep=""),
  #             quote=F,
  #             row.names = F)

  #rm(chr,snps,sites,eigenstuff,windist,mds,coor,localPCA)
  
# }
# dev.off()
