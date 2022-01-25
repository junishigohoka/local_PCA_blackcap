#!/bin/bash
#
#  example submit script for a serial job
#  submit by  sbatch serial-job.sh
#
#  specify the job name
#SBATCH --job-name=PopGenome
#  how many cpus are requested
#SBATCH --ntasks=1
#  run on one node, importand if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime, here 10min
#SBATCH --time=03:00:00
#  maximum requested memory
#SBATCH --mem=12G
#  write std out and std error to these files
#SBATCH --error=/home/ishigohoka/stdout/PopGenome%J.err
#SBATCH --output=/home/ishigohoka/stdout/PopGenome%J.out
#  send a mail for job start, end, fail, etc.
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ishigohoka@evolbio.mpg.de
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=standard

#  add your code here:



vcf=$1
chrlist=$2
chr=$3
poplist=$4
win=$5
out=$6
dirscripts=$7

module load R/3.5.3
Rscript $dirscripts/PopGenome_windowstats.R \
--vcf $vcf \
--chr_list $chrlist \
--chr $chr \
--pop_list $poplist \
--window_size $win \
--window_jump $win \
--out $out

