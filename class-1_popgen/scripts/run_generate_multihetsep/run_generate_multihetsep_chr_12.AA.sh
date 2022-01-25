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
#SBATCH --time=00:20:00
#  maximum requested memory
#SBATCH --mem=5G
#  write std out and std error to these files
#SBATCH --error=/home/ishigohoka/stdout/multihetsep_%J.err
#SBATCH --output=/home/ishigohoka/stdout/multihetsep_%J.out
#SBATCH --mail-user=ishigohoka@evolbio.mpg.de
#SBATCH --mail-type=FAIL
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=standard

#  add your code here:


dir_mask=$1
dir_vcf=$2
dir_multihetsep=$3

module load python/3.5.3

generate_multihetsep.py --chr chr_12 --mask $dir_mask/O065107_chr_12.mask.bed.gz --mask $dir_mask/O065118_chr_12.mask.bed.gz --mask $dir_mask/KZ97910_chr_12.mask.bed.gz --mask $dir_mask/G05_chr_12.mask.bed.gz   $dir_vcf/O065107.chr_12.vcf.gz $dir_vcf/O065118.chr_12.vcf.gz $dir_vcf/KZ97910.chr_12.vcf.gz $dir_vcf/G05.chr_12.vcf.gz  > $dir_multihetsep/chr_12.AA.multihetsep.txt

 