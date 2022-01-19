#!/bin/bash
#
#  example submit script for a serial job
#  submit by  sbatch serial-job.sh
#
#  specify the job name
#SBATCH --job-name=sumdecode
#  how many cpus are requested
#SBATCH --ntasks=1
#  run on one node, importand if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime, here 10min
#SBATCH --time=00:20:00
#  maximum requested memory
#SBATCH --mem=12G
#  write std out and std error to these files
#SBATCH --error=/home/ishigohoka/stdout/sumdecode%J.err
#SBATCH --output=/home/ishigohoka/stdout/sumdecode%J.out
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
        dirdecodesum=$dirout/msmc-decode_summary
dirfigures=$dirbase/figures

model=$2
gen=$3
id=$4

function mjoin() {
    tmp1=.tmp.$RANDOM
    tmp2=.tmp.$RANDOM
    cat $1 | awk '{print $1}' > $tmp1
    for f in $*
    do 
            join $tmp1 $f > $tmp2
            mv $tmp2 $tmp1
    done
    cat $tmp1
    rm $tmp1
}


for geno in NN NI II
do
        #if [[ ! -f $dirdecodesum/${model}_${id}_gen.$gen.$geno.decode-sum.txt ]]; then
                while read ID
                do
                        awk '{p=0;k=0;for(i=2;i<=NF;i++){if($i>p){p=$i;k=i}};print $1,k}' $dirdecode/${model}_${id}_gen.$gen.$geno.$ID.posterior.txt > $dirdecodesum/${model}_${id}_gen.$gen.$geno.$ID.decode-sum.txt.tmp
                done<$dirlist/genotype/${model}_${id}_gen.$gen.$geno.4samples.list
                files=`cat $dirlist/genotype/${model}_${id}_gen.$gen.$geno.4samples.list | awk -v prefix=$dirdecodesum/${model}_${id}_gen.$gen.$geno -v ORS=" " '{print prefix"."$1".decode-sum.txt.tmp"}'`
                mjoin $files > $dirdecodesum/${model}_${id}_gen.$gen.$geno.decode-sum_ind.txt
                rm $files
                awk '{s=0;for(i=2;i<=NF;i++){s+=$i};print $1,s/(NF-1)}' $dirdecodesum/${model}_${id}_gen.$gen.$geno.decode-sum_ind.txt > $dirdecodesum/${model}_${id}_gen.$gen.$geno.decode-sum.txt
        #fi
done

