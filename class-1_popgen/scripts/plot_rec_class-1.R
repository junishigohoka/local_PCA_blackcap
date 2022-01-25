library("optparse")
 
option_list = list(
       make_option(c("--dirlist"), type="character", default=NULL, 
              help="Path to dirlist", metavar="character"),
       make_option(c("--dirout"), type="character", default=NULL, 
              help="Path to dirout", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

dirlist=opt$dirlist
dirout=opt$dirout


library("RColorBrewer")
cols<-brewer.pal(7,"Set1")[c(2,1,3)]


setwd(dirlist)
chrlist<-read.table("chromosomes_length.list",col.names=c("chr","len"))
lpca<-read.table("class-1.chr_12.bed",col.names=c("chr","POS1","POS2"))

setwd(dirout)
png("rec_rate_class-1.png",height=480*2,width=480*3,res=350)


for(i in 1:nrow(lpca)){
        chr<-lpca$chr[i]
        pos1<-lpca$POS1[i]
        pos2<-lpca$POS2[i]
        len<-chrlist[chrlist$chr==chr,]$len[1]
        d<-read.table(paste0(chr,".rec.10kb.tab"),header=T)
        plot(1,
             xlim=c(0,len)/1e6,
             ylim=c(0,50),
             type='n',
             xlab="Position [Mb]",
             ylab="Recombination rate [cM/Mb]",
             main=chr,
             las=1
        )
        polygon(c(pos1,pos2,pos2,pos1)/1e6,c(0,0,50,50),
                col="gray",
                border="gray"
        )
        for(j in 1:3){
                lines(d$pos1/1e6,d[,5+j],
                      type='s',
                      col=cols[j]
                )
        }
        if(i==5){
                legend("topright",
                       lty=1,
                       col=cols,
                       legend=c("AA","AB","BB")

                )
        }

}


dev.off()
