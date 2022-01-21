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

library("ape")

d1<-read.tree("RAxML_bipartitions.BC.chr_12.background.bootstrap.final.txt")
#d11<-rotate(d1,7);d11<-rotate(d11,8);d11<-rotate(d11,9)
d11<-rotate(d1,9);d11<-rotate(d11,12);d11<-rotate(d11,11);#d11<-rotate(d11,8)

d2<-read.tree("RAxML_bipartitions.BC.chr_12.inversion.bootstrap.final.txt")
d21<-rotate(d2,11)


#pdf("chr_12.raxml.pdf",width=18)
png("chr_12.raxml.png",width=480*6,height=480*3,res=300)
par(mfrow=c(1,2),
    mar=c(1,1,3,1)
)
# Plot background
plot.phylo(d11,
           show.node.label=T
)

add.scale.bar()
mtext(side=3,
      text="Chromosome 12 (background)",
      adj=0,
      cex=1.25
)


# Plot inversion
plot.phylo(d21,
           show.node.label=T
)
add.scale.bar()
mtext(side=3,
      text="Chromosome 12 (class-1 genomic island)",
      adj=0,
      cex=1.25
)
dev.off()

