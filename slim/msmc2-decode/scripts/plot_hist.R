library("optparse")
 
option_list = list(
       make_option(c("--dirbase"), type="character", default=NULL, 
              help="Path to base", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

dirbase<-opt$dirbase




setwd(paste0(dirbase))

d<-read.table("output/log/model_id_lastgen.txt",col.names=c("model","id","gen"))
# Remove ODDF
d<-d[grepl("ODFD",as.character(d$model),fixed=T)==F,]
d$model<-factor(as.factor(d$model),levels=c("neutral_neutral","neutral_OD","neutral_FD","mixed_neutral","mixed_OD","mixed_FD","deleterious_neutral","deleterious_OD","deleterious_FD"))

png("figures/gen.png",width=480*5,height=480*5,res=350)
par(mfrow=c(3,3))
for(i in 1:nlevels(d$model)){
        model<-levels(d$model)[i]
        sub<-d[d$model==model,]
        hist(sub$gen,main=model,
             breaks = seq(0,4000,100),
             col="gray",
             border=NA,
             xlim=c(0,4000),
             ylim=c(0,100),
             xlab="Time (generations)",
             las=1
        )
}
dev.off()


