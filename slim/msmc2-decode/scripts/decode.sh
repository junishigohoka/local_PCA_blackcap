#!/bin/bash
#
#  example submit script for a serial job
#  submit by  sbatch serial-job.sh
#
#  specify the job name
#SBATCH --job-name=decode
#  how many cpus are requested
#SBATCH --ntasks=1
#  run on one node, importand if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime, here 10min
#SBATCH --time=00:05:00
#  maximum requested memory
#SBATCH --mem=12G
#  write std out and std error to these files
#SBATCH --error=/home/ishigohoka/stdout/decode%J.err
#SBATCH --output=/home/ishigohoka/stdout/decode%J.out
#  send a mail for job start, end, fail, etc.
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ishigohoka@evolbio.mpg.de
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=testing

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


mu=4e-4
rec=4e-3

model=$2
gen=$3
id=$4


for geno in NN NI II
do
        while read ID idx1 idx2
        do
                decode -m $mu -r $rec -I $idx1,$idx2 -t 32 -s 10000 $dirmultihetsep/${model}_${id}_gen.$gen.$geno.multihetsep.txt > $dirdecode/${model}_${id}_gen.$gen.$geno.$ID.posterior.txt
        done < $dirlist/genotype/${model}_${id}_gen.$gen.$geno.4samples.idx.list
done


