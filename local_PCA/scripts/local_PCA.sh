#!/bin/bash
#
#  example submit script for a serial job
#  submit by  sbatch serial-job.sh
#
#  specify the job name
#SBATCH --job-name=localPCA
#  how many cpus are requested
#SBATCH --ntasks=1
#  run on one node, importand if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime, here 10min
#SBATCH --time=48:00:00
#  maximum requested memory
#SBATCH --mem=100G
#  write std out and std error to these files
#SBATCH --error=/home/ishigohoka/stdout/locPCA%J.err
#SBATCH --output=/home/ishigohoka/stdout/locPCA%J.out
#  send a mail for job start, end, fail, etc.
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ishigohoka@evolbio.mpg.de
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=standard

#  add your code here:
dirtped=$1
dirsites=$2
dirlocalpca=$3
dirlocalpcamds=$4
dirlist=$5
chr=$6
dirscripts=$7

module load R/3.5.3

Rscript $dirscripts/local_PCA.R --dirin $dirtped --dirsites $dirsites --dirpca $dirlocalpca --dirmds $dirlocalpcamds --dirlist $dirlist --chr $chr 

