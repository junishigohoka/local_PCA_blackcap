
library("optparse")
 
option_list = list(
       make_option(c("--input"), type="character", default=NULL, 
              help="Path to local PCA MDS txt file for whole-genome", metavar="character"),
       make_option(c("--thre"), type="character", default=NULL, 
              help="Path to local PCA MDS threshold txt file for whole-genome", metavar="character"),
	make_option(c("--output"), type="character", default="out.txt", 
              help="Output file name", metavar="character"),
	make_option(c("--chrlist"), type="character", default="out.txt", 
              help="Chromosome list ", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

input<-opt$input
output<-opt$output
thre<-read.table(opt$thre,header=T)
chrlist<-read.table(opt$chrlist)
colnames(chrlist)<-c("chr")
chrlist$index<-1:nrow(chrlist)


# Read local PCA results for the whole genome
d<-read.table(input,header=T)




# Merge with chrlist
thre<-merge(chrlist,thre)
thre<-thre[order(thre$index),]

# Prepare result data frame
res<-data.frame(matrix(ncol=3,nrow=0))
colnames(res)<-c("chr","from","to")


for (i in 1:nrow(thre)){
  chr<-thre$chr[i]
  high1<-thre$MDS1.high[i]
  low1<-thre$MDS1.low[i]
  high2<-thre$MDS2.high[i]
  low2<-thre$MDS2.low[i]
  mode1<-thre$MDS1.mode[i]
  mode2<-thre$MDS2.mode[i]
  
  sub<-d[d$chr==chr,]
  len<-thre$len[i]
  
  # outliers
  ol1<-sub[sub$MDS1>high1|sub$MDS1<low1&(sub$MDS2<=high2&sub$MDS2>=low2),]
  ol2<-sub[(sub$MDS1<=high1&sub$MDS1>=low1)&(sub$MDS2>high2|sub$MDS2<low2),]
  ol12<-sub[(sub$MDS1>high1|sub$MDS1<low1)&(sub$MDS2>high2|sub$MDS2<low2),]
  
  # non-outliers
  nol<-sub[sub$MDS1<=high1&sub$MDS1>=low1&sub$MDS2<=high2&sub$MDS2>=low2,]
  
  
  # Put outliers in one list
  ols<-list(ol1=ol1,ol2=ol2,ol12=ol12)
  
  # Step 1. Remove outlier group with less than 10 windows
  ols<-ols[lapply(ols,nrow)>=10]
  
  if(length(ols)>0){
    # Step 2. Calculate range of positions
    ## Make function
    iqr<-function(x){
      #tmp<-quantile(x,probs=c(0.05,0.95))
      tmp<-range(x)
      tmp<-(x[x>=tmp[1]&x<=tmp[2]])
      return(tmp[c(1,length(tmp))])
    }
    ## Get inter 5 percentile range in bed-like format
    olsbed<-as.data.frame(t(sapply(ols,function(x){iqr(c(x[,"pos"],x[,"posright"]))})))
    colnames(olsbed)<-c("from","to")
    ## sort
    olsbed<-olsbed[order(olsbed$from,olsbed$to),]
    
    # Step 3. Join overlaps
    if(nrow(olsbed)>1){
      for (j in 1:(nrow(olsbed)-1)){
        if(olsbed$to[j]>olsbed$from[j+1]){
          olsbed$from[j+1]<-olsbed$from[j]
          olsbed$to[j+1]<-max(olsbed$to[j:(j+1)])
          olsbed[j,]<-c(NA)
        }
      }
    }
    olsbed<-na.omit(olsbed)
    
    # add chr column
    olsbed<-data.frame(chr=chr,olsbed)
    
    # bind to the res data.frame
    res<-rbind(res,olsbed)
    
    # remove
    
  }
}

# Sort res
res<-res[order(res$chr),]

write.table(res,paste0(output,"_outlier.bed"),
            quote=F,
            col.names = F,
            row.names = F
            )

