#!/bin/bash
#
#  example submit script for a serial job
#  submit by  sbatch serial-job.sh
#
#  specify the job name
#SBATCH --job-name=multihetsep
#  how many cpus are requested
#SBATCH --ntasks=1
#  run on one node, importand if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime, here 10min
#SBATCH --time=03:00:00
#  maximum requested memory
#SBATCH --mem=25G
#  write std out and std error to these files
#SBATCH --error=/home/ishigohoka/stdout/multihetsep_%J.err
#SBATCH --output=/home/ishigohoka/stdout/multihetsep_%J.out
#SBATCH --mail-user=ishigohoka@evolbio.mpg.de
#SBATCH --mail-type=FAIL
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



module load python/3.5.3

while read id
do
        for geno in NN NI II
        do
                vcfs=`cat $dirlist/genotype/${model}_${id}_gen.$gen.$geno.4samples.list | awk -v ORS=" " -v prefix=$dirvcfind/${model}_${id}_gen.$gen '{print prefix"."$1".vcf.gz"}'`
                generate_multihetsep.py --chr 1 --mask $dirlist/mask.bed.gz $vcfs  > $dirmultihetsep/${model}_${id}_gen.$gen.$geno.multihetsep.txt
        done
done<$dirlist/vcf_id/${model}_gen.$gen.list


