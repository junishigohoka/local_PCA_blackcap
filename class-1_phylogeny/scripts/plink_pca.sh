#!/bin/bash
#
#  example submit script for a serial job
#  submit by  sbatch serial-job.sh
#
#  specify the job name
#SBATCH --job-name=pca
#  how many cpus are requested
#SBATCH --ntasks=1
#  run on one node, importand if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime, here 10min
#SBATCH --time=03:00:00
#  maximum requested memory
#SBATCH --mem=15G
#  write std out and std error to these files
#SBATCH --error=/home/ishigohoka/stdout/pca%J.err
#SBATCH --output=/home/ishigohoka/stdout/pca%J.out
#  send a mail for job start, end, fail, etc.
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ishigohoka@evolbio.mpg.de
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=standard

#  add your code here:

vcfin=$1
prefixout=$2



plink --vcf $vcfin --mind 0.2 --not-chr chr_Z --double-id --geno 0.9 --maf 0.001 --pca var-wts  --allow-no-sex --allow-extra-chr --out $prefixout  # 2> $dir_log/${prefix}_${chr}_plink_pca.log

