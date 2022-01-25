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

setwd(dirlist)

# Read chromosome lengths
chrlen <- read.table("chromosomes_length.list", col.names = c("chr","len"))
chrlen$chr <- as.character(chrlen$chr)


# Read local PCA outliers
lpca<-read.table("local_PCA_MDS_outlier.bed")
colnames(lpca)<-c("chr","pos1","pos2")
lpca$chr<-as.character(lpca$chr)



# Read TR summary
setwd(dirout)

d <- read.table("bSylAtr1.1.filtered.TRs.summary.wholegenome.txt",header=T)
thre<-200
d[d$nrep>thre,]$nrep<-thre

# Make colour palette
#redcol <- brewer.pal(3,"Set1")[1]
colfunc <- colorRampPalette(c("blue","purple","red"))
cols <- colfunc(51)[-1]
d$col <- cols[round(d$nrep/thre*49)+1]


# Plot
png("bSylAtr1.1.TR.heatmap.png",height=480*10,width=480*10,res=450)
#m<-2
#pdf("bSylAtr1.1.TR.heatmap.pdf",height=7*m,width=7*m)
par(mfrow=c(7,5),
    mar=c(3,3,3,1),
    oma=c(2,2,0,0))
for(i in 1:nrow(chrlen)){
        chr <- chrlen$chr[i]
        len <- chrlen$len[i]
        dsub <- d[d$chr==chr,]
        lpcasub <- lpca[lpca$chr==chr,]
        plot(1,
             type='n',
             xlim=c(0,len)/1e6,
             ylim=c(0,500),
             #ylim=c(10,max(d$lunitmax)),
               xlab="",
               ylab="",
               main=chr,
               las=1
        )
        if(nrow(lpcasub)>0){
                for(j in 1:nrow(lpcasub)){
                        x1 <- lpcasub$pos1[j]/1e6
                        x2 <- lpcasub$pos2[j]/1e6
                        polygon(c(x1,x2,x2,x1),c(0,0,500,500),
                                col="gray",
                                border=NA
                        )
                }
        }
        if(nrow(dsub)>0){
                for(j in 1:nrow(dsub)){
                        x1 <- dsub$pos1[j]/1e6
                        x2 <- dsub$pos2[j]/1e6
                        y1 <- dsub$lunitmin[j]
                        y2 <- dsub$lunitmax[j]
                        col <- dsub$col[j]
                        polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),
                                col = col,
                                border = col,
                                lwd=0.25
                        )
                }
                box()
        }
}
plot(1,type='n',
     xlim = c(0,10),
     ylim = c(0,thre+1),
     xlab = "",
     ylab = "",
     xaxt = 'n',
     yaxt = 'n',
     bty = 'n'
)

for(j in 1:thre){
        if(j>10){
                polygon(c(0,1,1,0),c(j,j,j+1,j+1),
                        col=cols[round(j/thre*49)+1],
                        border=cols[round(j/thre*49)+1]
                )
        }
}

axis(side=2,
     las=1,
     at=c(10,seq(50,thre,by=50)),
     label=c(10,seq(50,thre,by=50))
)

title("total n.rep",
      adj=0
)

mtext(side=1,
      outer=T,
      text="Position [Mb]"
)

mtext(side=2,
      outer=T,
      text="TR monomer size"
)
dev.off()


