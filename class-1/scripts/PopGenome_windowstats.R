library("optparse")

dir<-getwd()

option_list = list(
  make_option(c("--vcf"), type="character", default=NULL, 
              help="bgzipped vcf file. Also needs to be tabixed", metavar="character"),
  make_option(c("--chr_list"), type="character", default=NULL, 
              help="list of chromosomes and their lengths", metavar="character"),
  make_option(c("--chr"), type="character", default=NULL, 
              help="chromosome", metavar="character"),
  make_option(c("--pop_list"), type="character", default=NULL, 
              help="list of individuals and populations. \n
              1st column:ID, 2nd column:population \n
              header is required, with \"population\" for the 2nd" , metavar="character"),
  make_option(c("--window_size"), type="numeric", default=10000, 
              help="window size [bp] [default= %default]", metavar="numeric"),
  make_option(c("--window_jump"), type="numeric", default=10000, 
              help="window step size [bp] [default= %default]", metavar="numeric"),
  make_option(c("--out"), type="character", default=NULL,
              help="output prefix [default= %default] \n
              window size is automatically written after the prefix", metavar="character")
  
); 



opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);




if (is.null(opt$vcf)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).", call=FALSE)
}


vcf<-opt$vcf;print(paste("vcf",vcf))
chr.list<-opt$chr_list;print(paste("chr_list",chr.list))
chr<-opt$chr;print(paste("chr",chr))
pop.list<-opt$pop_list;print(paste("pop_list",pop.list))
window_size<-as.numeric(opt$window_size);print(paste("window_size",window_size))
window_jump<-as.numeric(opt$window_jump);print(paste("window_jump",window_jump))


out<-opt$out;print(paste("out",out))





# Load libaries
library(PopGenome)
library(tidyverse)


# Read chromosome length list
chrlen<-read.table(chr.list)
colnames(chrlen)<-c("chr","len")
head(chrlen)

# Read list of ids and populations
poplist<-read.table(pop.list,header=T)
# Omit NA population
poplist<-poplist[is.na(poplist$population)==F,]
poplist$population<-as.factor(poplist$population)


pops<-list()
for(i in 1:nlevels(poplist$population)){
  pops[[i]]<-poplist[poplist$population==levels(poplist$population)[i],]$sample
}
names(pops)<-levels(poplist$population)
# PopGenome automatycally makes list of individuals ind1,ind2,..., indn, ind1.2, ind2.2,..., indn.n for diploid VCF
# I will make the same list with populations assigned
pops.tmp<-lapply(pops,paste0,".2")
pops<-mapply(c,pops,pops.tmp)
rm(pops.tmp)





# Take the length of the focal chromosome 
len<-chrlen[chrlen$chr==chr,]$len




# Load VCF
d<-readVCF(filename = vcf,
           numcols=100000, # just to speed up loading
           tid=chr,
           frompos = 1, 
           topos=len)
# Define populations
d<-set.populations(d,pops)





# set window size and window jump
# window_size <- 10000
# window_jump <- 10000

# use seq to find the start points of each window
window_start <- seq(from = 1, to = len, by = window_jump)
# add the size of the window to each start point 
window_stop <- window_start + window_size -1

# Remove the last incomplete window (if it exists)
window_start <- window_start[window_stop < len]
window_stop <- window_stop[window_stop < len]


# save sliding windows as a data.frame
windows <- data.frame(chr=chr,start = window_start, stop = window_stop)

# Put the sliding windows into PopGenome data d
d.win <- sliding.window.transform(d, width = window_size , jump = window_jump, type = 2)

# Calculate pi
# d.win<-diversity.stats(d.win,pi=T)

# Calculate FST
d.win<-F_ST.stats(d.win,mode="nucleotide")
fst<-t(d.win@nuc.F_ST.pairwise)

# Calculate pi
pi<-d.win@nuc.diversity.within/window_size
colnames(pi)<-levels(poplist$population)
pi<-cbind(windows,pi)


# Calculate dxy
dxy<-t(d.win@nuc.diversity.between/window_size)





# Edit the colnames of fst and dxt, to show which pair of populations are focused for each column
# Make functino to do it
make.pop.pairs<-function(data){
  poppairs<-c()
  count<-0
  for(i in 1:(length(data)-1)){
    for(j in (i+1):length(data)){
      # print(paste0("pop",i,"/pop",j))
      count<-count+1
      #print(count)
      #print(paste0(levels(poplist$population)[i],"/",levels(poplist$population)[j]))
      poppairs<-c(poppairs,paste0(data[i],"/",data[j]))
    }
  }
  return(poppairs)
}

colnames(fst)<-make.pop.pairs(data=levels(poplist$population))
fst<-cbind(windows,fst)

colnames(dxy)<-make.pop.pairs(data=levels(poplist$population))
dxy<-cbind(windows,dxy)


# Output tables
windkb<-window_size/1000

write.table(fst,paste0(out,"_FST_PopGenome",windkb,"kb.txt"),
            row.names = F,
            quote=F)
write.table(dxy,paste0(out,"_dxy_PopGenome",windkb,"kb.txt"),
            row.names = F,
            quote=F)
write.table(pi,paste0(out,"_pi_PopGenome",windkb,"kb.txt"),
            row.names = F,
            quote=F)
