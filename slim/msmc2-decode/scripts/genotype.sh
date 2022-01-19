#!/bin/bash
#
#  example submit script for a serial job
#  submit by  sbatch serial-job.sh
#
#  specify the job name
#SBATCH --job-name=genotype
#  how many cpus are requested
#SBATCH --ntasks=1
#  run on one node, importand if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime, here 10min
#SBATCH --time=03:00:00
#  maximum requested memory
#SBATCH --mem=25G
#  write std out and std error to these files
#SBATCH --error=/home/ishigohoka/stdout/genotype%J.err
#SBATCH --output=/home/ishigohoka/stdout/genotype%J.out
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
        dirmultihetsep=$dirout/multihetsep
        dirdecode=$dirout/msmc-decode
dirfigures=$dirbase/figures

model=$2
gen=$3

while read id
do
        bcftools query -r 1:100001  -f '[%GT ]' $dirvcf/${model}_${id}_gen.$gen.vcf.gz | awk '{for(i=1;i<=NF;i++){print "i"i-1,$(i)}}' > $dirlist/genotype/${model}_${id}_gen.$gen.ID.gt.list
        awk '$2=="0|0"{print $1}'  $dirlist/genotype/${model}_${id}_gen.$gen.ID.gt.list > $dirlist/genotype/${model}_${id}_gen.$gen.NN.list
        awk '$2=="0|1"||$2=="1|0"{print $1}' $dirlist/genotype/${model}_${id}_gen.$gen.ID.gt.list > $dirlist/genotype/${model}_${id}_gen.$gen.NI.list
        awk '$2=="1|1"{print $1}' $dirlist/genotype/${model}_${id}_gen.$gen.ID.gt.list > $dirlist/genotype/${model}_${id}_gen.$gen.II.list
done<$dirlist/vcf_id/${model}_gen.$gen.list

