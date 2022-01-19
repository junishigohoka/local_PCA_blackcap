library("optparse")
 
option_list = list(
       make_option(c("--dirout"), type="character", default=NULL, 
              help="Path to output", metavar="character"),
       make_option(c("--dirfig"), type="character", default=NULL, 
              help="Path to figures", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

dirout<-opt$dirout
dirfig<-opt$dirfig

setwd(dirout)

d<-read.table("model_gen_id_gf_af.txt",col.names=c("model","gen","id","NN","NI","II","N","I"))
# Remove ODDF
d<-d[grepl("ODFD",as.character(d$model),fixed=T)==F,]
d$model<-factor(as.factor(d$model),levels=c("neutral_neutral","neutral_OD","neutral_FD","mixed_neutral","mixed_OD","mixed_FD","deleterious_neutral","deleterious_OD","deleterious_FD"))
d$gen<-d$gen-4000

setwd(dirfig)
png("gt_af.png",width=480*8,height=480*5,res=350)
par(mfrow=c(3,3),
    mar=c(4,4,3,1)
)
for(i in 1:nlevels(d$model)){
        model<-levels(d$model)[i]
        sub<-d[d$model==model,]
        plot(jitter(sub$gen),sub$I/2000,
             xlab="Time (generations)",
             ylab="Inversion frequency",
             xaxt='n',
             main=model,
             xlim=c(0,4000),
             ylim=c(0,1),
             las=1
        )
        axis(side=1,
             at=c(100,500,1000,2000,4000)
        )
        for(gen in c(100,500,1000,2000,4000)){
                subsub<-sub[sub$gen==gen,]
                m<-median(subsub$I/2000)
                #s<-sd(subsub$I/2000)
                s<-quantile(subsub$I/2000)
                points(rep(gen,3),s[2:4],
                       col="red",
                       pch="-",
                       cex=2
                )
                segments(gen,s[2],gen,s[4],
                         col="red",
                         lty=2,
                         cex=2
                )
                text(gen,1,
                     labels=nrow(subsub)
                )
        }
}
dev.off()

