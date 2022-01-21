library("optparse")
options(show.error.locations = TRUE) 
option_list = list(
       make_option(c("--dirout"), type="character", default=NULL, 
              help="Path to dirout", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dirout<-opt$dirout


setwd(dirout)

d<-read.table("ZF_synteny_BC_chr_12.14126710.22227355.eigenvec",header=F)
d<-d[,-2]
colnames(d)<-c("ID",paste0("PC",1:19))

png("ZF_synteny_BC_chr_12.14126710.22227355.PCA.png",height=480*3,width=480*3,res=350)
plot(d$PC1,d$PC2,
     xlab="PC1",
     ylab="PC2",
     main="zebra finch"
)

dev.off()
