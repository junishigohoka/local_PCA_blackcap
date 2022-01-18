library("optparse")
 
option_list = list(
       make_option(c("--dir"), type="character", default=NULL, 
              help="Path to PCA output per outlier", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

dir<-opt$dir
setwd(dir)

library("RColorBrewer")
cols <- brewer.pal(7,"Set1")[c(2,1,3)]

png("chr_12.het.png",height=480*4,width=480*6,res=350)
par(las=1)
d <- read.table("chr_12.het.txt",header=T)

        boxplot(d$het~d$geno,
                ylab="Heterozygosity",
                xlab="Genotype",
                pch = NA,
                col = cols[1:3],
                main="chr_12"
        )
        points(jitter(as.numeric(d$geno)),d$het)
dev.off()


