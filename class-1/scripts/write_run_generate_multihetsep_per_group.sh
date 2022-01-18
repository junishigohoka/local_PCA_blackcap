#!/bin/bash

# This script generate a bash file run_generate_multihetsep_<Super_scaffold_XX>.sh


dirbase=$1
dirlist=$dirbase/list
dirin=$dirbase/input
dirout=$dirbase/output
dirscripts=$dirbase/scripts
# dirmask=$dirin/mask

dirmask=/home/ishigohoka/projects/Miriam/PhD/PhD/callability_mask/output/mask

dirvcfind=$dirin/vcf

#mkdir $dirbashout

cd $dirmask

# Make one "run_generate_multihetsep_<Super_scaffold_XX>.sh" per Super-Scaffold



# # Create a list of mask files for one Sper-Scaffold, which will directly used as an arguement in sed to edit the run_generate_multihetsep.sh
# cat $dirlist/$pop.list  | awk '{print $1".mask.bed.gz"}' | awk '{ print "--mask \\$dir_mask\\ /" $1 }' | sed 's/ //2' | sed 's/ //2' | awk '{ ORS=" " }{ print }' >$dirlist/$pop.mask.list
# masklist=`cat $dirlist/$pop.mask.list`

# # Create a list of vcf files for one chromosome, which will be directly used as an arguement in sed to edit the run_generate_multihetsep.sh

if [ ! -d $dirlist/vcflist ];then
	mkdir $dirlist/vcflist
fi

if [ ! -d $dirlist/masklist ];then
	mkdir $dirlist/masklist
fi

if [ ! -d $dirscripts/run_generate_multihetsep ];then
	mkdir $dirscripts/run_generate_multihetsep
fi

# while read ID
# do
while read chr prefix pc left right
do
	for geno in AA AB BB	
	do
		# VCF list
		cat $dirlist/$prefix.$geno.callable.list  | awk -v chr=$chr '{print $1"."chr".vcf.gz"}' | awk '{ print "\\$dir_vcf\\ /"$1 }' | sed 's/ //2' | sed 's/\\ //' | awk '{ ORS=" " }{ print }' > $dirlist/vcflist/$prefix.$geno.vcf.list
		vcflist=`cat $dirlist/vcflist/$prefix.$geno.vcf.list`

		# Mask list
		cat $dirlist/$prefix.$geno.callable.list | awk -v chr=$chr '{print $1"_"chr".mask.bed.gz"}' | awk '{ print "--mask \\$dir_mask\\ /"$1 }'| sed 's/ //2' | sed 's/\\ //' | awk '{ ORS=" " }{ print }' > $dirlist/masklist/$prefix.$geno.mask.list
		masklist=`cat $dirlist/masklist/$prefix.$geno.mask.list`
			
		# done
		# echo ""
		# # Edit the run_gemerate_multihetsep_<Super-Scaffold_XX>.sh
		#echo $masklist
		#echo $vcflist
		# for chr in chr1 chr2
		# do
		sed -e "s/chromosome/$chr/" $dirscripts/run_generate_multihetsep_template.sh | sed -e "s@--mask@$masklist@" | sed "s/outmultihetsep/$prefix.$geno.multihetsep.txt/" | sed "s@inputvcf@$vcflist@" | sed '' > $dirscripts/run_generate_multihetsep/run_generate_multihetsep_$prefix.$geno.sh
		chmod +x $dirscripts/run_generate_multihetsep/run_generate_multihetsep_$prefix.$geno.sh
	done
done < $dirlist/chr_interval_PC_boundary.list

