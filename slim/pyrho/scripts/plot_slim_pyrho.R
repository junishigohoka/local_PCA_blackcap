library("optparse")
 
option_list = list(
       make_option(c("--dirout"), type="character", default=NULL, 
              help="Path to dirout", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

dirout<-opt$dirout


library("RColorBrewer")
cols<-brewer.pal(7,"Set1")[c(2,1,3)]

setwd(dirout)

png("slim.pyrho.png",height=480*7.5,width=480*5,res=350)
par(mfrow=c(6,2),
    mar=c(3,3,3,1),
    oma=c(2,2,0,0)
)

for(i in 1:6){
        model<-paste0("model",i)
        for(j in 1:2){
                chr<-paste0("chr",j)
                dNN<-read.table(paste(model,"NN",chr,"rmap.cMMb.txt",sep="."),col.names=c("pos","cM.Mb","cM"))
                dNI<-read.table(paste(model,"NI",chr,"rmap.cMMb.txt",sep="."),col.names=c("pos","cM.Mb","cM"))
                dII<-read.table(paste(model,"II",chr,"rmap.cMMb.txt",sep="."),col.names=c("pos","cM.Mb","cM"))
                # Plot
                plot(1,type='n',
                     xlim=c(0,5),
                     ylim=c(0,70),
                     xlab="",
                     ylab="",
                     main=paste0(model,", ",chr),
                     las=1
                )
                lines(dNN$pos/1e6,
                      dNN$cM.Mb,
                      type='s',
                      col=cols[1]
                )
                lines(dNI$pos/1e6,
                      dNI$cM.Mb,
                      type='s',
                      col=cols[2]
                )
                lines(dII$pos/1e6,
                      dII$cM.Mb,
                      type='s',
                      col=cols[3]
                )
                mtext(side=3,
                      line=1,
                      adj=-0.1,
                      text=LETTERS[c(2*i-1,2*i)][j],
                      font=2
                )
                if(j==2&i==1){
                        legend("bottomright",
                               lty=1,
                               col=cols,
                               legend=c("NN","NI","II")
                        )
                }
        }
}

mtext(side=1,
      outer=T,
      text="Position [Mb]"
)

mtext(side=2,
      outer=T,
      text="Recombination rate [cM/Mb]"
)

dev.off()



