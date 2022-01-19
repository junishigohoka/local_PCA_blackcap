library("optparse")
 
option_list = list(
       make_option(c("--dirout"), type="character", default=NULL, 
              help="Path to dirout", metavar="character"),
       make_option(c("--dirfig"), type="character", default=NULL, 
              help="Path to dirfig", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

dirout<-opt$dirout
dirfig<-opt$dirfig


library("RColorBrewer")
cols<-brewer.pal(7,"Set1")[c(2,1,3:7)]

setwd(dirout)
d<-read.table("model_id_gen_tmrca.txt",header=T)

d$model<-gsub("FD","frequency-dependent selection",gsub("OD","overdominance",gsub("_",", ",d$model)))
d$model<-factor(as.factor(d$model),levels=c("neutral, neutral","neutral, overdominance","neutral, frequency-dependent selection","mixed, neutral","mixed, overdominance","mixed, frequency-dependent selection","deleterious, neutral","deleterious, overdominance","deleterious, frequency-dependent selection"))

d$gen<-d$gen-4000

setwd(dirfig)
png("tmrca_in.png",width=480*8,height=480*5,res=350)
par(mfrow=c(3,3),
    mar=c(4,4,3,1)
)
deltax=100
cexpt=0.75
for(i in 1:nlevels(d$model)){
        model<-levels(d$model)[i]
        sub<-d[d$model==model,]
        plot(jitter(sub$gen,factor=0.4)-deltax,sub$NNin,
             col=cols[1],
             pch=1,
             cex=cexpt,
             xlab="Time (generations)",
             ylab="TMRCA",
             xaxt='n',
             main=model,
             xlim=c(0,4000),
             ylim=c(0,32),
             las=1
        )
        points(jitter(sub$gen,factor=0.4),sub$NIin,
               col=cols[2],
               cex=cexpt,
               pch=1
        )
        points(jitter(sub$gen,factor=0.4)+deltax,sub$IIin,
               col=cols[3],
               cex=cexpt,
               pch=1
        )
        axis(side=1,
             at=c(100,500,1000,2000,4000)
        )
}
legend("topright",
       col=cols[1:3],
       pch=1,
       pt.cex=cexpt,
       legend=c("NN","NI","II")
)
dev.off()


png("tmrca_out.png",width=480*8,height=480*5,res=350)
par(mfrow=c(3,3),
    mar=c(4,4,3,1)
)
deltax=100
for(i in 1:nlevels(d$model)){
        model<-levels(d$model)[i]
        sub<-d[d$model==model,]
        plot(jitter(sub$gen,factor=0.4)-deltax,sub$NNout,
             col=cols[1],
             pch=1,
             cex=cexpt,
             xlab="Time (generations)",
             ylab="TMRCA",
             xaxt='n',
             main=model,
             xlim=c(0,4000),
             ylim=c(0,32),
             las=1
        )
        points(jitter(sub$gen,factor=0.4),sub$NIout,
               col=cols[2],
               cex=cexpt,
               pch=1
        )
        points(jitter(sub$gen,factor=0.4)+deltax,sub$IIout,
               col=cols[3],
               cex=cexpt,
               pch=1
        )
        axis(side=1,
             at=c(100,500,1000,2000,4000)
        )
}
legend("topright",
       col=cols[1:3],
       pch=1,
       pt.cex=cexpt,
       legend=c("NN","NI","II")
)
dev.off()

