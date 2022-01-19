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

dirbase=$1
dirscripts=$dirbase/scripts
dirlist=$dirbase/list
dirout=$dirbase/output
        dirlog=$dirout/log
        dirsum=$dirout/summary
        dirvcf=$dirout/vcf
        dirvcfind=$dirout/vcfind
        dirmultihetsep=$dirout/multihetsep
        dirdecode=$dirout/msmc-decode
dirfigures=$dirbase/figures

model=$2
gen=$3

while read id
do
        for geno in NN NI II
        do
                while read ID
                do
                        bcftools view -s $ID $dirvcf/${model}_${id}_gen.$gen.vcf.gz | bgzip > $dirvcfind/${model}_${id}_gen.$gen.$ID.vcf.gz
                        bcftools index $dirvcfind/${model}_${id}_gen.$gen.$ID.vcf.gz
                done<$dirlist/genotype/${model}_${id}_gen.$gen.$geno.4samples.list
        done
done<$dirlist/vcf_id/${model}_gen.$gen.list





