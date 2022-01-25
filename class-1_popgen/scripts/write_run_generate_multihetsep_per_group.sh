#!/bin/bash

# This script generate a bash file run_generate_multihetsep_<Super_scaffold_XX>.sh


dirbase=$1
dirlist=$dirbase/list
dirin=$dirbase/input
dirout=$dirbase/output
dirscripts=$dirbase/scripts
# dirmask=$dirin/mask

dirmask=$dirin/mask

dirvcfind=$dirin/vcf

#mkdir $dirbashout

cd $dirmask

# Make one "run_generate_multihetsep_<Super_scaffold_XX>.sh" per Super-Scaffold



if [ ! -d $dirlist/vcflist ];then
	mkdir $dirlist/vcflist
fi

if [ ! -d $dirlist/masklist ];then
	mkdir $dirlist/masklist
fi

if [ ! -d $dirscripts/run_generate_multihetsep ];then
	mkdir $dirscripts/run_generate_multihetsep
fi

for chr in chr_12
do
	for geno in AA AB BB	
	do
		# VCF list
		cat $dirlist/$chr.$geno.4samples.list  | awk -v chr=$chr '{print $1"."chr".vcf.gz"}' | awk '{ print "\\$dir_vcf\\ /"$1 }' | sed 's/ //2' | sed 's/\\ //' | awk '{ ORS=" " }{ print }' > $dirlist/vcflist/$chr.$geno.vcf.list
		vcflist=`cat $dirlist/vcflist/$chr.$geno.vcf.list`

		# Mask list
		cat $dirlist/$chr.$geno.4samples.list | awk -v chr=$chr '{print $1"_"chr".mask.bed.gz"}' | awk '{ print "--mask \\$dir_mask\\ /"$1 }'| sed 's/ //2' | sed 's/\\ //' | awk '{ ORS=" " }{ print }' > $dirlist/masklist/$chr.$geno.mask.list
		masklist=`cat $dirlist/masklist/$chr.$geno.mask.list`
			
		sed -e "s/chromosome/$chr/" $dirscripts/run_generate_multihetsep_template.sh | sed -e "s@--mask@$masklist@" | sed "s/outmultihetsep/$chr.$geno.multihetsep.txt/" | sed "s@inputvcf@$vcflist@" | sed '' > $dirscripts/run_generate_multihetsep/run_generate_multihetsep_$chr.$geno.sh
		chmod +x $dirscripts/run_generate_multihetsep/run_generate_multihetsep_$chr.$geno.sh
	done
done 

