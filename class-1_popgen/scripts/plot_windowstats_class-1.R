library("optparse")
 
option_list = list(
       make_option(c("--dir"), type="character", default=NULL, 
              help="Path to PCA output per outlier", metavar="character"),
       make_option(c("--bedfile"), type="character", default=NULL, 
              help="Path to BED file specifying outlier coordinates", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

dir<-opt$dir
bedfile<-opt$bedfile


# Define colours
library(RColorBrewer)
cols<-brewer.pal(8,"Paired")
pale<-cols[c(1,3,5,7)]
solid<-cols[c(2,4,6,8)]

setwd(dir) 
# Read chromosome length
chrlen<-read.table("../list/chromosomes_length.list",col.names=c("chr","len"))
chrlen$chr<-as.character(chrlen$chr)


# Read local PCA outlier interval
bed<-read.table(bedfile,col.names=c("chr","from","to"))
bed$chr<-as.character(bed$chr)


# Make png
f1<-8
f2<-1.75


png(paste0("chr_12.class-1.windowstats.png"),width=480*f1,height=480*f2,res=350)
par(mar=c(3,4,4,1),
    las=1
)
mat<-matrix(nrow=1,ncol=6)
for(j in 1:nrow(mat)){
        mat[j,]<-c(1,1,2,2,3,4)+(4*(j-1))
}
layout(mat)
lx<-2
ly<-2.5
lm<-0.25
cex1=0.75
cex2=0.75
# Loop over SVs 
for(k in 1:nrow(bed)){
        chr <- bed$chr[k]
        len<-chrlen[chrlen$chr==chr,]$len
        from<-bed[bed$chr==chr,]$from
        to<-bed[bed$chr==chr,]$to

# Read results
        fst<-read.table(paste("blackcap",chr,"AA.BB_FST_PopGenome10kb.txt", sep="."),header=T)
        dxy<-read.table(paste("blackcap",chr,"AA.BB_dxy_PopGenome10kb.txt", sep="."),header=T)
        pi<-read.table(paste("blackcap",chr,"AA.BB_pi_PopGenome10kb.txt", sep="."),header=T)
        if(chr%in%c("chr_6","chr_30")){
                pi[,4:5]<-pi[,5:4]
        }

# Calculate midpoint
        fst$mid<-(fst$start+fst$stop-1)/2
        dxy$mid<-(dxy$start+dxy$stop-1)/2
        pi$mid<-(pi$start+pi$stop-1)/2

# Remove NAs in fst
        fst<-fst[rowSums(is.na(fst))==0, ]
# Keep the same rows as FST
        dxy<-dxy[as.numeric(rownames(fst)),]

# Subset for within interval
        fstin<-fst[(fst$start>=from&fst$start<to)|(fst$stop>=from&fst$stop<to)|(fst$start>=from&fst$stop<to),]
        dxyin<-dxy[(dxy$start>=from&dxy$start<to)|(dxy$stop>=from&dxy$stop<to)|(dxy$start>=from&dxy$stop<to),]
        piin<-pi[(pi$start>=from&pi$start<to)|(pi$stop>=from&pi$stop<to)|(pi$start>=from&pi$stop<to),]

# Two populations here are two genotypes AA and BB
        pop1 <- "AA"
        pop2 <- "BB"
        pair <- 4
# Plot dxy
        plot(1,
             type='n',
             xlim=c(0,len)/10^6,
             ylim=c(0,0.0225),
             xlab="",
             ylab="",
             yaxt='n'
        )
        mtext(side=3,
                #text=paste0(pop1," vs ",pop2),
                text=chr,
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
        polygon(c(from,to,to,from)/10^6,c(0,0,0.0225,0.0225),
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
             at=c(0,0.005,0.01,0.015,0.02),
             labels=c(0,"",0.01,"",0.02)
        )
        
# Plot pi
        lines(pi$start/10^6,
              #pi[,pop1],
              pi[,4],
              type='s',
              col=solid[1]
        )
        lines(pi$start/10^6,
              #pi[,pop2],
              pi[,5],
              type='s',
              col=solid[2]
        )
        #lines(pi$start/1e6,
        #      pi[,4]-pi[,5],
        #      type='s',
        #      col=solid[4]
        #)
        box()

# legend
        legend("topright",
               legend=c(as.expression(bquote(paste(pi [.(pop1)]))),as.expression(bquote(paste(pi [.(pop2)]))),"dxy"),
               lty=c(1,1,NA),
               pch=c(NA,NA,19),
               col=c(solid[c(1,2)],"gray"),
               # bty='n',
               cex=0.75
        )




# Plot FST
        plot(1,
             type='n',
             xlab="",
             ylab="",
             xlim=c(0,len)/10^6,
             ylim=c(-0.1,1),
             pch=19
        )
# Shade interval
        polygon(c(from,to,to,from)/10^6,c(-0.1,-0.1,1,1),
                col=pale[3],
                border=NA
        )
        points(fst$mid/10^6,
             fst[,pair],
             col="gray",
             pch=19
        )
        box()
# Add interval points
        points(fstin$mid/10^6,
               fstin[,pair],
               pch=19,
               col=solid[3]
        )

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
             xlim=c(-0.1,1),
             ylim=c(0,0.02),
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
             at=c(0,0.005,0.01,0.015,0.02),
             labels=c(0,"",0.01,"",0.02)
        )
# Plot pi(pop1) vs pi(pop2)
        plot(pi[,4],
             pi[,5],
             xlim=c(0,0.02),
             ylim=c(0,0.02),
             xaxt='n',
             yaxt='n',
             xlab="",
             ylab="",
             col=adjustcolor("black",alpha.f=0.1),
             pch=19
        )
# Add interval windows
        points(piin[,4],
               piin[,5],
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
             at=c(0,0.005,0.01,0.015,0.02),
             labels=c(0,"",0.01,"",0.02)
        )
        axis(side=2,
             at=c(0,0.005,0.01,0.015,0.02),
             labels=c(0,"",0.01,"",0.02)
        )
}
dev.off()

