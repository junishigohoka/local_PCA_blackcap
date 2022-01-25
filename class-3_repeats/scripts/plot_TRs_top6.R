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
cols <- brewer.pal(7,"Set1")

setwd(dirlist)
# Read chromosome length
chrlen <- read.table("chromosomes_length.list",col.names=c("chr","len"))
chrlen$chr <- as.character(chrlen$chr)


# Read local PCA outliers
lpca <- read.table("local_PCA_MDS_outlier.bed",col.names=c("chr","pos1","pos2"))

# Read data
setwd(dirout)
d <- read.table("bSylAtr1.1.filtered.uniqTRs.top6.unmerged.per.chr.txt", col.names=c("chr","pos1","pos2","lrep","nrep","seq"))
d$chr <- as.character(d$chr)
d$mid <- apply(d[,c("pos1","pos2")],1,mean)


# Read supplementary data (info of unique TRs)
dsup <- read.table("bSylAtr1.1.filtered.uniqTRs.top6.per.chr.txt",col.names=c("chr","totalnrep","lrep","seq"))
dsup$chr <- as.character(dsup$chr)


# Make plot
png("bSylAtr1.1.TR.top6.png",height=480*10,width=480*20,res=450)

mat <- matrix(sort(c(rep(1:70,2,each=T),seq(1,70,2))),nrow=7,byrow=T)
layout(mat)
par(oma=c(0.5,0.5,0.5,0.5))
for(i in 1:nrow(chrlen)){
        chr <- chrlen$chr[i]
        len <- chrlen$len[i]
        dsub <- d[d$chr==chr,]
        dsupsub <- dsup[dsup$chr==chr,]
        lpcasub <- lpca[lpca$chr==chr,]
        ylimmax <- 1.1*max(dsub$nrep)
        par(mar=c(3,4.5,3,0.5))
        plot(1,
             type='n',
             xlim=c(0,len)/1e6,
             ylim=c(0,ylimmax),
             #ylim=c(10,max(d$lunitmax)),
               #xlab="Position [Mb]",
               ylab="",
               main="",
               las=1
        )
        title(main=chr,
              adj=0,
              line=0.5
        )
        if(nrow(lpcasub)>0){
                for(j in 1:nrow(lpcasub)){
                        x1 <- lpcasub$pos1[j]/1e6
                        x2 <- lpcasub$pos2[j]/1e6
                        polygon(c(x1,x2,x2,x1),c(0,0,ylimmax,ylimmax),
                                col="gray",
                                border=NA
                        )
                }
        }
        if(nrow(dsub)>0){
                seqs<-unique(dsub$seq)
                for(j in 1:length(seqs)){
                        dsubseq<-dsub[dsub$seq==seqs[j],]
                        points(dsubseq$mid/1e6,
                               dsubseq$nrep,
                               pch=1,
                               col=cols[j],
                               cex=1.5,
                               lwd=2
                        )
                }
        }
        box()
        mtext(side=2,
              cex = 0.75,
              line=3,
              text="TR copy number"
        )
        mtext(side=1,
              cex=0.75,
              line=2.1,
              text="Position [Mb]"
        )
        # Plot bar plot
        par(mar=c(3,3.25,3,2))
        barplot(dsupsub$lrep,
                col=cols,
                border = cols,
                xlab="",
                #ylab="TR monomer size [bp]",
                #xaxt='n',
                main = "",
                las=1,
                names.arg = letters[1:length(seqs)]
        )
        mtext(side=2,
              text="Monomer size [bp]",
              cex=0.75,
              line=2.5
        )

        #mtext(side=1,
        #      cex=0.75,
        #      line=2.1,
        #      text="TRs"
        #)
}

#mtext(side=1,outer=T,
#      text = "Position [Mb]"
#)

#mtext(side=2,outer=T,
#      text = "TR copy number"
#)

dev.off()


