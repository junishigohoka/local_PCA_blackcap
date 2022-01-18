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


options(scipen = 999)
library(RColorBrewer)
cols<-brewer.pal(7,"Set1")[c(2,1,3)]
#library("beeswarm")

setwd(dirlist)
chrlen<-read.table("chromosomes_length.list",col.names = c("chr","len"))
chrlen$chr<-as.character(chrlen$chr)

prefixlist<-chrlen
chrft<-read.table("class-1.chr_12.bed",col.names = c("chr","from","to"))
chrft$chr<-as.character(chrft$chr)

prefixlist<-merge(prefixlist,chrft)
genos<-c("AA","AB","BB")

prefixlist<-prefixlist[order(-prefixlist$len),]

# Make empty data frame to put stats results

stats<-data.frame(matrix(ncol=5,nrow=0))
colnames(stats) <- c("prefix","interval","KW","MW_AA_AB","MW_BB_AB")


setwd(dirout)
pdf("chr_12.tmrca.pdf",width=3.5*5,height=3.5*nrow(prefixlist))
m<-2.5
par(mfrow=c(nrow(prefixlist),5),
    las=1
    )
for(i in 1:nrow(prefixlist)){
    chr<-prefixlist$chr[i]
    len<-prefixlist$len[i]
    from<-prefixlist$from[i]
    to<-prefixlist$to[i]

    for(k in 1:3){
            GENO=genos[k]
            if(chr%in%c("chr_6","chr_30")){
                    k<-4-k
            }
        geno<-genos[k]
        files<-list.files(".",paste0(chr,".",geno,".*tmrca.txt"))
        # Make empty plot
        plot(1,type='n',
            xlab="Position [Mb]",
            ylab="TMRCA",
            xlim=c(0,len/1000000),
            ylim=c(0,32),
            main=paste(chr,GENO)
            )
            if(chr%in%c("chr_6","chr_30")){
                    k<-4-k
            }
        # Add lines at the interval boundaries
        polygon(c(from,to,to,from)/1000000,c(-10,-10,40,40),
            col=cols[k],
            border=F
            )
        if(length(files)>0){
          # Read TMRCA data into a single data frame
            for(j in 1:length(files)){
                tmp<-read.table(files[j],col.names = c("pos",paste0("tmrca",j)))
                if(j==1){
                    d<-tmp
                }else{
                    d<-merge(d,tmp)
                }
            }
            # I removed consecutive windows with exactly same TMRCA in all samples because they must be under LD
            # The first window stays
            for(l in 2:nrow(d)){
                if(all(as.vector(d[l,2:ncol(d)])==as.vector(d[l-1,2:ncol(d)]))){
                    d[l,1]<-NA
                } 
            }
            d <- na.omit(d)
            # calculate median TMRCA across 4 samples for each genomic location
            d$medtmrca<-NA
            for(l in 1:nrow(d)){
                d$medtmrca[l]<-round(median(as.numeric(d[l,2:(length(files)+1)])))
            }
            # Plot median TMRCA along the chromosome
            lines(d$pos/1000000,d$medtmrca)
            # Make dsum
            d$geno<-genos[k]
            if(k==1){
              dsum<-d[,c("pos","medtmrca","geno")]
            }else{
              dsum<-rbind(dsum,d[,c("pos","medtmrca","geno")])
            }
        }
        # Write title
        mtext(side=3,outer=T,
            text=chr
            )
    }
    # Perform Mann-Whitney test on TMRCA within the interval (AB-AA, AB-BB), and out of the interval
    dsum$interval<-"out"
    dsum[dsum$pos>=from&dsum$pos<to,"interval"]<-"in"
    write.table(dsum,paste0(chr,".tmrca.txt"),
                row.names = F,
                quote=F)

    for(interval in c("in","out")){
        if(interval=="in"){
            colbox<-cols[1:3]
        }else{
            colbox<-NA
        }
        dsub<-dsum[dsum$interval==interval,]
        dsub$geno<-factor(dsub$geno,levels=genos)


        # Plot 
        boxplot(dsub$medtmrca~dsub$geno,
                ylab="TMRCA",
                ylim=c(0,32),
                main=interval,
                col=colbox)
        #beeswarm(dsub$medtmrca~dsub$geno,
        #    add=T,
        #    method="square",
        #    cex=min(1,300/nrow(dsub))
        #    # corral="omit"
        #    )
        legend("topleft",
                col=cols[1:3],
                legend=paste(genos,c(nrow(dsub[dsub$geno=="AA",]),
                        nrow(dsub[dsub$geno=="AB",]),
                        nrow(dsub[dsub$geno=="BB",])),"windows"),
                        bty='n')
        # Subset for AA-AB and BB-AB
        # If the number of windows is greater than 100, I sampled random 100 windows, to balance the power of test across within-interval and out of interval


        dsubAA<-dsub[dsub$geno!="BB",]
        dsubAA$geno<-factor(dsubAA$geno,levels=c("AA","AB"))
        if(nrow(dsubAA)>100){
            dsubAA<-dsubAA[sample(nrow(dsubAA),100),]
        }
        dsubBB<-dsub[dsub$geno!="AA",]
        dsubBB$geno<-factor(dsubBB$geno,levels=c("BB","AB"))
        if(nrow(dsubBB)>100){
            dsubBB<-dsubBB[sample(nrow(dsubBB),100),]
        }
        # Kruskal-Wallis test

        stats[nrow(stats)+1,]<-NA
        stats$prefix[nrow(stats)]<-chr
        stats$interval[nrow(stats)]<-interval
        stats$KW[nrow(stats)]<-round(kruskal.test(dsub$medtmrca~dsub$geno)$p.value,3)
        # Mann-Whitney test
        if(stats$KW[nrow(stats)]<=0.05&nrow(dsub[dsub$geno=="AB",])>0){
            if(nrow(dsubAA)>0){
                stats$MW_AA_AB[nrow(stats)]<- round(wilcox.test(dsubAA$medtmrca~dsubAA$geno)$p.value,3)
                }
            if(nrow(dsubBB)>0){
                stats$MW_BB_AB[nrow(stats)]<- round(wilcox.test(dsubBB$medtmrca~dsubBB$geno)$p.value,3)
                }
        }
        rm(dsubAA,dsubBB)
    }
    # dev.off()
}
dev.off()

write.table(stats,"tmrca.stats.txt",
    quote=F,
    row.names = F
    )
