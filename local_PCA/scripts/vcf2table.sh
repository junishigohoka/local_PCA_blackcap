#!/bin/bash
#
#  example submit script for a serial job
#  submit by  sbatch serial-job.sh
#
#  specify the job name
#SBATCH --job-name=vcf2tab
#  how many cpus are requested
#SBATCH --ntasks=1
#  run on one node, importand if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime, here 10min
#SBATCH --time=12:00:00
#  maximum requested memory
#SBATCH --mem=25G
#  write std out and std error to these files
#SBATCH --error=/home/ishigohoka/stdout/vcf2tab%J.err
#SBATCH --output=/home/ishigohoka/stdout/vcf2tab%J.out
#  send a mail for job start, end, fail, etc.
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ishigohoka@evolbio.mpg.de
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=standard

#  add your code here:

vcf=$1
tab=$2

bcftools view -m2 -M2 $vcf |bcftools query -f '[%GT ]\n' - | sed 's@0/0@0@g; s@0|0@0@g; s@1/1@2@g; s@1|1@2@g; s@0/1@1@g; s@0|1@1@g; s@1/0@1@g; s@1|0@1@g; s@./.@NA@g; s@.|.@NA@g'> $tab
