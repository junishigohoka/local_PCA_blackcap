#!/bin/bash
#
#  example submit script for a serial job
#  submit by  sbatch serial-job.sh
#
#  specify the job name
#SBATCH --job-name=MSMC2
#  how many cpus are requested
#SBATCH --ntasks=1
#  run on one node, importand if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime, here 10min
#SBATCH --time=72:00:00
#  maximum requested memory
#SBATCH --mem=100G
#  write std out and std error to these files
#SBATCH --error=/home/ishigohoka/stdout/MSMC2_per_pop_%J.err
#SBATCH --output=/home/ishigohoka/stdout/MSMC2_per_pop_%J.out
#  send a mail for job start, end, fail, etc.
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ishigohoka@evolbio.mpg.de
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=standard

#  add your code here:

#  directories
dir_multihetsep=
dir_out=

cd $dir_out


msmc2 -t 4 -I -o  $dir_multihetsep/INPUT
