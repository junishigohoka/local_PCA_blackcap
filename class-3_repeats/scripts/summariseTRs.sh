#!/bin/bash
#
#  example submit script for a serial job
#  submit by  sbatch serial-job.sh
#
#  specify the job name
#SBATCH --job-name=sumTRs
#  how many cpus are requested
#SBATCH --ntasks=1
#  run on one node, importand if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime, here 10min
#SBATCH --time=00:05:00
#  maximum requested memory
#SBATCH --mem=12G
#  write std out and std error to these files
#SBATCH --error=/home/ishigohoka/stdout/summariseTRs%J.err
#SBATCH --output=/home/ishigohoka/stdout/summariseTRs%J.out
#  send a mail for job start, end, fail, etc.
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ishigohoka@evolbio.mpg.de
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=testing

#  add your code here:
module load R/3.5.3

chr=$1
dirout=$2
dirlist=$3

Rscript summariseTRs.R --chr $chr --dirout $dirout --dirlist $dirlist

