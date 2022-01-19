#!/bin/bash
#
#  example submit script for a serial job
#  submit by  sbatch serial-job.sh
#
#  specify the job name
#SBATCH --job-name=slim
#  how many cpus are requested
#SBATCH --ntasks=1
#  run on one node, importand if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime, here 10min
#SBATCH --time=00:30:00
#  maximum requested memory
#SBATCH --mem=25G
#  write std out and std error to these files
#SBATCH --error=/home/ishigohoka/stdout/slim%J.err
#SBATCH --output=/home/ishigohoka/stdout/slim%J.out
#  send a mail for job start, end, fail, etc.
#SBATCH --mail-type=fail
#SBATCH --mail-user=ishigohoka@evolbio.mpg.de
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=standard

#  add your code here:

dirbase=$1
dirscripts=$dirbase/scripts
dirlist=$dirbase/list
dirout=$dirbase/output
        dirlog=$dirout/log
        dirvcf=$dirout/vcf
        dirmultihetsep=$dirout/multihetsep
        dirdecode=$dirout/msmc-decode

id=$2

while read f0 f1 s0 h0 s1 h1 s2 h2 FD model
do
        slim -d "OUTPUT='${model}_$id'" $dirscripts/$model.slim  > $dirlog/${model}_$id.log 
done<$dirlist/parameters.list


