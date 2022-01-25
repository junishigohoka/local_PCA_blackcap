# Install optparse if not installed yet
if(!require("optparse",quietly=T)){
        install.packages("optparse",repo="https://cloud.r-project.org")
}
library("optparse")
options(scipen=999)

#mapname<-"/home/ishigohoka/projects/Miriam/PhD/PhD/pyrho/data/pyrhoout/chr_31_ContRes_n38_Pen20_W50.rmap"
#winsize<-10000
#winstep<-5000
#chrlen<-559822
#chr<-"chr_31"
#outprefix<-"/home/ishigohoka/projects/Miriam/PhD/PhD/pyrho/output/test"

option_list = list(
  make_option("--map", action="store", default=NA, type='character',
              help="Path to pyrho output"),
  make_option( "--winsize", action="store", default=10000, type='numeric',
              help="Window size. [default %default]"),
  make_option( "--winstep", action="store", default=10000, type='numeric',
              help="Window step size [default %default]"),
  make_option("--chrlen", action="store", default=NA, type='numeric',
              help="Chromosome length"),
  make_option("--chr", action="store", default=NA,type='character',
              help="Chromosome name"),
  make_option("--output", action="store", default=NA, type='character',
              help="Prefix of output file")
)
opt = parse_args(OptionParser(option_list=option_list))


mapname<-opt$map
winsize<-opt$winsize
winstep<-opt$winstep
chrlen<-opt$chrlen
chr<-opt$chr
outprefix<-opt$output

# Read pyrho output
pyrho<-read.table(mapname, col.names=c("pos1","pos2","r"))



# Add physical length of each segment [bp]
pyrho$phylen<-pyrho$pos2-pyrho$pos1

# Add columns for map length (cM for the focal segment) and map positions (cM)
pyrho$mappos1<-NA
pyrho$mappos2<-NA
pyrho$maplen<-NA


mappre<-0 # Initial map position is 0 cM

for(i in 1:nrow(pyrho)){
        # Map position of pos1
        pyrho$mappos1[i]<-mappre
        # Map length of i-th segment
        pyrho$maplen[i]<-pyrho$r[i]*pyrho$phylen[i]*100
        # Map position of pos2
        pyrho$mappos2[i]<-pyrho$mappos1[i]+pyrho$maplen[i]
        # Update previous map position
        mappre<-pyrho$mappos2[i]
}





# Prepare a dataframe for window average

d<-data.frame(chr=chr,
              pos1=seq(0,chrlen,winstep),
              pos2=seq(0,chrlen,winstep)+winsize,
              mappos1=NA,
              mappos2=NA,
              cMpMb=NA
)

# Remove the first and last windows that are only partially covered by pyrho output
d<-d[d$pos1>=pyrho$pos1[1]&d$pos2<pyrho$pos2[nrow(pyrho)],]


# Calculate mean recombination rate per window
for(i in 1:nrow(d)){
        # Get positions of the i-th window borders
        pos1<-d$pos1[i]
        pos2<-d$pos2[i]
        mappospre<-pyrho$mappos1[pyrho$pos1<=pos1&pos1<pyrho$pos2] # map position of left border of pyrho row overlapping pos1
        phypospre<-pyrho$pos1[pyrho$pos1<=pos1&pos1<pyrho$pos2] # physical position of left border of pyrho row overlapping pos1
        rpre<-pyrho$r[pyrho$pos1<=pos1&pos1<pyrho$pos2] # recombination rate in this pyrho row
        d$mappos1[i]<-mappospre+rpre*(pos1-phypospre)*100
        # Calculate map position [cM] for pos2
        mappospre<-pyrho$mappos1[pyrho$pos1<=pos2&pos2<pyrho$pos2]
        phypospre<-pyrho$pos1[pyrho$pos1<=pos2&pos2<pyrho$pos2]
        rpre<-pyrho$r[pyrho$pos1<=pos2&pos2<pyrho$pos2]
        d$mappos2[i]<-mappospre+rpre*(pos2-phypospre)*100

        # Calculate mean recombination rate of the i-th window
        d$cMpMb[i]<-(d$mappos2[i]-d$mappos1[i])/winsize*1e6
}


# Output
stepkb<-winstep/1000
sizekb<-winsize/1000
write.table(d,paste0(outprefix,"_win.",sizekb,"kb_step.",stepkb,"kb.mean.rec.tab"),
            row.names= F,
            quote=F,
            sep="\t"
)

