#!/bin/bash
dirbase=/home/ishigohoka/projects/Miriam/PhD/PhD/PSMCprime_blackcap
dir_bash=$dirbase/scripts
dir_multihetsep=$dirbase/input/multihetsep
dir_list=$dirbase/list
dir_msmc2=$dirbase/output
dir_out=$dirbase/scripts/PSMCprime


if [ ! -d ${dir_out} ];then
    mkdir $dir_out
fi

echo $dirbase	
# Edit the MSMC2_per_pop_<population>.sh
while read chr prefix pc left right
do
	for geno in AA AB BB
	do
		while read id idx1 idx2
		do
            sed -e "s/-I/-I $idx1,$idx2/g" $dir_bash/MSMC2_per_group_template.sh | sed -e "s/-o /-o $prefix.$geno.$id/g" |sed 's/name=MSMC2/name=PSMCp/g' | sed "s@dir_multihetsep=@dir_multihetsep=$dir_multihetsep@g" | sed "s@dir_out=@dir_out=$dir_msmc2@g" | sed "s/INPUT/$prefix.$geno.multihetsep.txt/g"> $dir_out/PSMCprime_$prefix.$geno.$id.sh
		done < $dir_list/$prefix.$geno.callable.idx.list
	done
done < $dir_list/chr_interval_PC_boundary.list

# while read spp
# do
#     while read ID idx1 idx2
#     do
#         sed -e "s/-I/-I $idx1,$idx2/g" $dir_bash/MSMC2_per_group_template.sh | sed -e "s/-o /-o $ID/g" |sed 's/name=MSMC2/name=PSMCp/g' | sed "s@dir_multihetsep=@dir_multihetsep=$dir_multihetsep@g" | sed "s@dir_out=@dir_out=$dir_msmc2@g" | sed "s/INPUT/\*.multihetsep.txt/g"> $dir_out/PSMCprime_${spp}_$ID.sh
#     done < ${dir_list}/$spp.idx.list
# done < ${dir_list}/out.species.list


