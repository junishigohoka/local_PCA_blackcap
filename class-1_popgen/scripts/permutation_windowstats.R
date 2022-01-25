library("optparse")
#n<-10000
#chr<-"chr_28"
#pos1<-917037
#pos2<-1154843
#pos1<-3000000
#pos2<-3500000
#
#pifile<-"/home/ishigohoka/projects/Miriam/PhD/PhD/PopGenome/output/blackcap.chr_28.917037.1154843.AA.BB_pi_PopGenome10kb.txt"
#dxyfile<-"/home/ishigohoka/projects/Miriam/PhD/PhD/PopGenome/output/blackcap.chr_28.917037.1154843.AA.BB_dxy_PopGenome10kb.txt"
#fstfile<-"/home/ishigohoka/projects/Miriam/PhD/PhD/PopGenome/output/blackcap.chr_28.917037.1154843.AA.BB_FST_PopGenome10kb.txt"
#
#
#
#

option_list = list(
  make_option(c("--chr"), type="character", default=NULL, 
              help="Chromosome name", metavar="character"),
  make_option(c("--pos1"), type="numeric", default=NULL, 
              help="Starting position of the interval", metavar="numeric"),
  make_option(c("--pos2"), type="numeric", default=NULL, 
              help="End position of the interval", metavar="numeric"),
  make_option(c("--n"), type="numeric", default=NULL, 
              help="Number of sampling sessions",metavar="numeric"),
  make_option(c("--pifile"), type="character", default=10000, 
              help="Path to pi window stats file", metavar="character"),
  make_option(c("--dxyfile"), type="character", default=NULL, 
              help="Path to dxy window stats file", metavar="character"),
  make_option(c("--fstfile"), type="character", default=NULL,
              help="Path to FST window stats file", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


chr<-opt$chr
pos1<-opt$pos1
pos2<-opt$pos2
n<-opt$n
pifile<-opt$pifile
dxyfile<-opt$dxyfile
fstfile<-opt$fstfile


# Read pi
pi<-read.table(pifile,header=T)
colnames(pi)[4:5]<-c("pi1","pi2")


# Read dxy
dxy<-read.table(dxyfile,header=T)
colnames(dxy)[4]<-"dxy"


# Read fst
fst<-read.table(fstfile,header=T)
colnames(fst)[4]<-"fst"


# Bind
d<-cbind(pi,dxy$dxy);d<-cbind(d,fst$fst)
colnames(d)[6:7]<-c("dxy","fst")


# Add column for pi-pi
d$pi1pi2<-d$pi1-d$pi2
d<-na.omit(d)

# Extract windows in the interval
dsub<-d[(d$start<=pos1&pos1<d$stop)|(pos1<=d$start&d$stop<pos2)|(d$start<=pos2&pos2<d$stop),]
obs.stats<-sapply(dsub[,4:8],mean)

# Permutation
for(i in 1:(n-1)){
        tmp<-d
        tmp[,2:3]<-tmp[sample(1:nrow(tmp)),2:3]
        # Shuffle each column
        for(j in 4:8){
                tmp[,j]<-tmp[sample(1:nrow(tmp)),j]
        }
        tmpsub<-tmp[(tmp$start<=pos1&pos1<tmp$stop)|(pos1<=tmp$start&tmp$stop<pos2)|(tmp$start<=pos2&pos2<tmp$stop),]
        if(i==1){
                nulldist<-as.data.frame(t(sapply(tmpsub[,4:8],mean)))
        }else{
                nulldist<-rbind(nulldist,as.data.frame(t(sapply(tmpsub[,4:8],mean))))
        }
}


# Calculate p.value for each stat
pval<-c()

for(i in 1:length(obs.stats)){
        if(i<=2){ # Left side test for pi
                pval<-c(pval,(sum(nulldist[,i]<=obs.stats[i])+1)/n)
        }else if(i<=4){ # Right side test for dxy and fst
                pval<-c(pval,(sum(nulldist[,i]>obs.stats[i])+1)/n)
        }else{ # Two side test for pi1-pi2
                pval<-c(pval,2*min(sum((nulldist[,i])>=(obs.stats[i]))+1,sum((nulldist[,i])<(obs.stats[i]))+1)/n)
        }
        #if(mean(nulldist[,i])<=obs.stats[i]){
        #        pval<-c(pval,(sum(nulldist[,i]>=obs.stats[i])+1)/n)
        #}else{
        #        pval<-c(pval,(sum(nulldist[,i]<obs.stats[i])+1)/n)
        #}
}
cat(chr,pos1,pos2,pval,"\n")

