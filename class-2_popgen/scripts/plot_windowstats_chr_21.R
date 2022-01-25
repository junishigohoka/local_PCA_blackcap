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

# Define colours
library(RColorBrewer)
cols<-brewer.pal(8,"Paired")
pale<-cols[c(1,3,5,7)]
solid<-cols[c(2,4,6,8)]

setwd(dirlist) 
# Read chromosome length
chrlen<-read.table("../list/chromosomes_length.list",col.names=c("chr","len"))
chr<-"chr_21"
len<-chrlen[chrlen$chr==chr,]$len

# Read local PCA outlier interval
bed<-read.table("local_PCA_MDS_outlier.bed",col.names=c("chr","from","to"))
from<-bed[bed$chr=="chr_21",]$from
to<-bed[bed$chr=="chr_21",]$to


# Read results
setwd(dirout)
fst<-read.table("blackcap.chr_21_FST_PopGenome10kb.txt",header=T)
dxy<-read.table("blackcap.chr_21_dxy_PopGenome10kb.txt",header=T)
pi<-read.table("blackcap.chr_21_pi_PopGenome10kb.txt",header=T)


# Calculate midpoint
fst$mid<-(fst$start+fst$stop-1)/2
dxy$mid<-(dxy$start+dxy$stop-1)/2
pi$mid<-(pi$start+pi$stop-1)/2

# Remove NAs in fst
fst<-fst[rowSums(is.na(fst))==0, ]
# Keep the same rows as FST
dxy<-dxy[as.numeric(rownames(fst)),]

# Subset for within interval
fstin<-fst[(fst$start>=from&fst$start<to)|(fst$stop>=from&fst$stop<to)|(fst$from>=from&fst$to<to),]
dxyin<-dxy[(dxy$start>=from&dxy$start<to)|(dxy$stop>=from&dxy$stop<to)|(dxy$from>=from&dxy$to<to),]
piin<-pi[(pi$start>=from&pi$start<to)|(pi$stop>=from&pi$stop<to)|(pi$from>=from&pi$to<to),]



# fstin<-fst[fst$mid>=from&fst$mid<to,]
# dxyin<-dxy[dxy$mid>=from&dxy$mid<to,]
# piin<-pi[pi$mid>=from&pi$mid<to,]


# Define order of pop.pairs
pop.pair<-data.frame(matrix(nrow=10,ncol=4))
colnames(pop.pair)<-c("pop1","pop2","pair","idx")
pop.pair$idx<-c(6,5,11,7,13,4,9,12,8,10)
pop.pair$pop1<-c("Azores","Azores","CapeVerde","Azores","CapeVerde","Azores","CapeVerde","Canary","Belgium","Belgium")
pop.pair$pop2<-c("CapeVerde","Canary","Canary","Spain_Caz","Spain_Caz","Belgium","Belgium","Spain_Caz","Canary","Spain_Caz")
pop.pair$pair<-colnames(dxy)[pop.pair$idx]


f1<-8
f2<-12
png(paste0(chr,".windowstats.png"),width=480*f1,height=480*f2,res=350)
par(mar=c(3,4,4,1),
    las=1
)
mat<-matrix(nrow=10,ncol=6)
for(j in 1:nrow(mat)){
        mat[j,]<-c(1,1,2,2,3,4)+(4*(j-1))
}
layout(mat)
lx<-2
ly<-2.5
lm<-0.25
cex1=0.75
cex2=0.75
for(i in 1:nrow(pop.pair)){
        pop1<-pop.pair$pop1[i]
        pop2<-pop.pair$pop2[i]
        pair<-pop.pair$pair[i]
        # Plot dxy
        plot(1,
             type='n',
             xlim=c(0,len)/10^6,
             ylim=c(-0.002,0.0125),
             xlab="",
             ylab="",
             yaxt='n'
        )
        mtext(side=3,
                text=paste0(pop1," vs ",pop2),
                adj=0,
                line=lm
        )
        mtext(side=1,
              text="Position [Mb]",
              cex=0.75,
              line=lx
        )
        mtext(side=2,
              text=expression(paste(pi,", dxy")),
              las=0,
              cex=cex2,
              line=ly
        )

        # Shade interval
        polygon(c(from,to,to,from)/10^6,c(0,0,0.0125,0.0125),
                col=pale[3],
                border=NA
        )

        points(dxy$mid/10^6,
             dxy[,pair],
             col="gray",
             pch=19
        )
        points(dxyin$mid/10^6,
             dxyin[,pair],
             col=solid[3],
             pch=19
        )

        axis(side=2,
             at=c(0,0.005,0.01),
             labels=c(0,"",0.01)
        )
        
        # Plot pi
        lines(pi$start/10^6,
              pi[,pop1],
              type='s',
              col=solid[1]
        )
        lines(pi$start/10^6,
              pi[,pop2],
              type='s',
              col=solid[2]
        )
        # Plot delta pi
        lines(pi$start/10^6,
              pi[,pop1]-pi[,pop2],
              type='s',
              col=solid[4]
        )

        # legend
        legend("topright",
               legend=c(as.expression(bquote(paste(pi [.(pop1)]))),as.expression(bquote(paste(pi [.(pop2)]))),"dxy"),
               lty=c(1,1,NA),
               pch=c(NA,NA,19),
               col=c(solid[1:2],"gray"),
               # bty='n',
               cex=0.75
        )




        # Plot FST
        plot(1,
             type='n',
             xlab="",
             ylab="",
             xlim=c(0,len)/10^6,
             ylim=c(-0.1,0.4),
             pch=19
        )
        # Shade interval
        polygon(c(from,to,to,from)/10^6,c(-0.1,-0.1,0.4,0.4),
                col=pale[3],
                border=NA
        )
        points(fst$mid/10^6,
             fst[,pair],
             col="gray",
             pch=19
        )
        # Add interval points
        points(fstin$mid/10^6,
               fstin[,pair],
               pch=19,
               col=solid[3]
        )


        # mtext(side=3,
        #       text=paste0(pop1," vs ",pop2),
        #       line=lm
        # )
        mtext(side=1,
              text="Position [Mb]",
              cex=0.75,
              line=lx
        )
        mtext(side=2,
              text=expression(F[ST]),
              las=0,
              cex=cex2,
              line=ly
        )
        # Plot FST vs dxy
        plot(x=fst[,pair],
             y=dxy[,pair],
             xlim=c(-0.1,0.4),
             ylim=c(0,0.01),
             xlab="",
             ylab="",
             col=adjustcolor("black",alpha.f=0.1),
             pch=19,
             yaxt='n'
        )
        # Add interval windows
        points(fstin[,pair],
               dxyin[,pair],
               pch=19,
               col=solid[3]
        )

        mtext(side=3,
              text=expression(paste(F[ST] ," vs dxy",sep="")),
              line=lm
        )
        mtext(side=1,
              text=expression(F[ST]),
              cex=0.75,
              line=lx
        )
        mtext(side=2,
              text="dxy",
              las=0,
              cex=cex2,
              line=ly
        )
        axis(side=2,
             at=c(0,0.005,0.01),
             labels=c(0,"",0.01)
        )
        # Plot pi(pop1) vs pi(pop2)
        plot(pi[,pop1],
             pi[,pop2],
             xlim=c(0,0.01),
             ylim=c(0,0.01),
             xaxt='n',
             yaxt='n',
             xlab="",
             ylab="",
             col=adjustcolor("black",alpha.f=0.1),
             pch=19
        )
        # Add interval windows
        points(piin[,pop1],
               piin[,pop2],
               pch=19,
               col=solid[3]
        )

        mtext(side=3,
              text=expression(pi),
              line=lm
        )
        mtext(side=1,
              text=pop1,
              cex=0.75,
              line=lx
        )
        mtext(side=2,
              text=pop2,
              las=0,
              cex=cex2,
              line=ly
        )
        axis(side=1,
             at=c(0,0.005,0.01),
             labels=c(0,"",0.01)
        )
        axis(side=2,
             at=c(0,0.005,0.01),
             labels=c(0,"",0.01)
        )
        
}

dev.off()


