library("optparse")
 
option_list = list(
       make_option(c("--dirpca"), type="character", default=NULL, 
              help="Path to PCA output per outlier", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

dirpca<-opt$dirpca




library("RColorBrewer")
cols<-c(brewer.pal(7,"Set2")[1:3])


setwd(dirpca)

d <- read.table("eigenvalues.txt",header=T)

d$pc1pc2<-d$PC1/d$PC2


d$class<-3
d[d$chr%in%paste0("chr_",c(6,12,14,28,30)),]$class <- 1
d[d$chr%in%paste0("chr_",c(21)),]$class <- 2



png("eigenval.png",height=480*3,width=480*3*3,res=350)

par(
    mfrow=c(1,3)
)



for(i in 1:nrow(d)){
        barplot(as.vector(t(d[i,4:23])),
                names.arg=paste0("PC",1:20),
                main=paste0(d$chr[i],"\n",d$pos1[i]," - ",d$pos2[i]),
                col=cols[d$class[i]],
                border=cols[d$class[i]],
                las=2,
                ylab="eigenvalue"
        )
}


dev.off()

