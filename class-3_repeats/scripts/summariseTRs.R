library("optparse")
options(show.error.locations = TRUE) 
option_list = list(
       make_option(c("--dirout"), type="character", default=NULL, 
              help="Path to dirout", metavar="character"),
       make_option(c("--dirlist"), type="character", default=NULL, 
              help="Path to dirlist", metavar="character"),
       make_option(c("--chr"), type="character", default=NULL, 
              help="Chromosome", metavar="character")

); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
dirout<-opt$dirout
dirlist<-opt$dirlist
chr<-opt$chr

# Read chr list
setwd(dirlist)
chrlen <- read.table("chromosomes_length.list",col.names=c("chr","len"))
chrlen$chr <- as.character(chrlen$chr)


options(scipen = 999)
library("RColorBrewer")

setwd(dirout)


# Read TSV
d <- read.table("bSylAtr1.1.filtered.TRs.tsv", col.names = c("chr","pos1","pos2","lunit","nrep","seq"),header)
d$chr <- as.character(d$chr)
# Calculate mid position
d$mid <-apply(d[,c("pos1","pos2")],1,mean)

# Define bin size
win <- 100000

# Make empty data frame for output
smr <- data.frame(matrix(nrow=0,ncol=6))
colnames(smr) <- c("chr","pos1","pos2","lunitmin","lunitmax","nrep")

# Loop over chromosomes
len <- chrlen[chrlen==chr,]$len
dsub <- d[d$chr==chr,]
pos1s <- seq(1,len,by=win)
pos2s <- pos1s + win -1
# Loop over unit length
for(j in seq(10,max(d$lunit),10)){
        minl <- j
        maxl <- j+10
        dsublen <- dsub[dsub$lunit>=minl & dsub$lunit<maxl,]    
        # Loop over positions
        for(k in 1:length(pos1s)){
                print(paste(chr,j,k/length(pos1s)*100))
                pos1 <- pos1s[k]
                pos2 <- pos2s[k]
                dsublenwin <- dsublen[dsublen$mid>=pos1&dsublen$mid<pos2,]
                tmp <- data.frame(matrix(nrow=0,ncol=6))
                colnames(tmp) <- c("chr","pos1","pos2","lunitmin","lunitmax","nrep")
                if(nrow(dsublenwin)>0){
                        tmp[1,]<-c(chr,pos1,pos2,minl,maxl,sum(dsublenwin$nrep))
                        smr <- rbind(smr,tmp)
                }
        }
}


write.table(smr,paste0("bSylAtr1.1.filtered.TRs.summary.",chr,".txt"),
            row.names = F,
            quote = F
)

