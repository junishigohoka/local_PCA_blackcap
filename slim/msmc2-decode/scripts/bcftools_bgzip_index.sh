#!/bin/bash
#
#  example submit script for a serial job
#  submit by  sbatch serial-job.sh
#
#  specify the job name
#SBATCH --job-name=bcftools
#  how many cpus are requested
#SBATCH --ntasks=1
#  run on one node, importand if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime, here 10min
#SBATCH --time=03:00:00
#  maximum requested memory
#SBATCH --mem=25G
#  write std out and std error to these files
#SBATCH --error=/home/ishigohoka/stdout/bcftools%J.err
#SBATCH --output=/home/ishigohoka/stdout/bcftools%J.out
#  send a mail for job start, end, fail, etc.
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ishigohoka@evolbio.mpg.de
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=standard

#  add your code here:

model=$1
gen=$2
dirvcf=$3
dirlist=$4

while read id
do
        bgzip -f $dirvcf/${model}_${id}_gen.$gen.vcf 
        bcftools index $dirvcf/${model}_${id}_gen.$gen.vcf.gz
done<$dirlist/vcf_id/${model}_gen.$gen.list


