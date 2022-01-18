# Pipeline for "Recombination suppression and selection affect local ancestries in genomes of a migratory songbird"

## Local PCA

Here I show how local PCA was run with three example chromosomes, 12, 20, and 21, which are listed in `local_PCA/list/chromosomes/list`

```bash
dirbase=$PWD/local_PCA
dirin=$dirbase/input
dirvcf=$dirin/vcf
dirtab=$dirin/table
dirlist=$dirbase/list
dirout=$dirbase/output
dirlocalpca=$dirout/local_PCA
dirlocalpcamds=$dirout/local_PCA_MDS
dirscripts=$dirbase/scripts

```

Download filtered VCF files for chromosomes 12, 20, 21.

```bash


```


Make list of sites from VCF.
```bash
while read chr
do
        bcftools query -f '%POS\n' $dirvcf/$chr.vcf.gz > $dirtab/$chr.sites.list
done <$dirlist/chromosomes.list

```

Make genotype table from VCF.
```bash

while read chr
do
        sbatch $dirscripts/vcf2table.sh $dirvcf/$chr.vcf.gz $dirtab/$chr.table
done <$dirlist/chromosomes.list


```


Run local PCA.
`local_PCA/scripts/local_PCA.sh` submits `local_PCA/scripts/local_PCA.R` via slurm.
Check `local_PCA/scripts/local_PCA.sh` and  `local_PCA/scripts/local_PCA.R` for detail.

```bash
module load R/3.5.3

while read chr 
do
        sbatch $dirscripts/local_PCA.sh $dirtab $dirtab $dirlocalpca $dirlocalpcamds $dirlist $chr $dirscripts
done <$dirlist/chromosomes.list


```


Concatenate the output for the three chromosomes.
```bash

while read chr
do
        cat $dirlocalpcamds/local_PCA_MDS_$chr.txt
done <$dirlist/chromosomes.list | awk 'NR==1{print $0}NR>1{if($1!="chr")print $0}' > $dirlocalpcamds/local_PCA_MDS_3_chromosomes.txt


```


Get threshold of MDS values and plot MDS distribution.
```bash

Rscript $dirscripts/plot_local_PCA_MDS_distribution.R --input $dirlocalpcamds/local_PCA_MDS_3_chromosomes.txt --output $dirlocalpcamds/local_PCA_MDS_3_chromosomes --chrlist $dirlist/chromosomes.list

```
![](local_PCA/output/local_PCA_MDS/local_PCA_MDS_3_chromosomes_MDS1_distribution.png)

![](local_PCA/output/local_PCA_MDS/local_PCA_MDS_3_chromosomes_MDS2_distribution.png)


Plot MDS1 vs MDS2.

```bash

Rscript $dirscripts/plot_local_PCA_MDS_MDS1-vs-MDS2.R --input $dirlocalpcamds/local_PCA_MDS_3_chromosomes.txt --thre $dirlocalpcamds/local_PCA_MDS_3_chromosomesthreshold.txt --output $dirlocalpcamds/local_PCA_MDS --chrlist $dirlist/chromosomes.list

```
![](local_PCA/output/local_PCA_MDS/local_PCA_MDS_MDS1-vs-MDS2_outliers.png)



Get coordinates of outlier windows.
```bash
Rscript $dirscripts/getOutliers.R --input $dirlocalpcamds/local_PCA_MDS_3_chromosomes.txt --thre $dirlocalpcamds/local_PCA_MDS_3_chromosomesthreshold.txt --output $dirlocalpcamds/local_PCA_MDS --chrlist $dirlist/chromosomes.list

```

Make Manhattan plots.
```bash
Rscript $dirscripts/plot_local_PCA_MDS_manhattan.R --input $dirlocalpcamds/local_PCA_MDS_3_chromosomes.txt --thre $dirlocalpcamds/local_PCA_MDS_3_chromosomesthreshold.txt --bed $dirlocalpcamds/local_PCA_MDS_outlier.bed --output $dirlocalpcamds/local_PCA_MDS --chrlist $dirlist/chromosomes_length.list

```

![](local_PCA/output/local_PCA_MDS/local_PCA_MDS_MDS1_outliers.png)

![](local_PCA/output/local_PCA_MDS/local_PCA_MDS_MDS2_outliers.png)



## PCA in local PCA outliers


Here I show how PCA was run for three local PCA outlier regions of chromosomes 12, 20, and 21, which represent class-1, 2, and 3 outliers.

```bash
dirbase=$PWD/PCA
dirin=$dirbase/input
dirvcf=$dirin/vcf
dirlist=$dirbase/list
dirout=$dirbase/output
dirscripts=$dirbase/scripts

```

Link VCF files from `local_PCA/input/vcf` and index them.

```bash
ln $dirbase/../local_PCA/input/vcf/*vcf.gz $dirvcf

while read chr
do
        echo $chr
        bcftools index $dirvcf/$chr.vcf.gz
done<$dirlist/chromosomes.list


```

Copy BED file for coordinates of local PCA outliers of chromosomes 12, 20, and 21.

```bash
ln $dirbase/../local_PCA/output/local_PCA_MDS/local_PCA_MDS_outlier.bed $dirlist

```

Extract SNPs within the local PCA outlier regions.

```bash

while read chr from to
do
        bcftools view -O z -r $chr:$from-$to -S $dirlist/blackcap_id.list $dirvcf/$chr.vcf.gz | vcftools --gzvcf - --max-missing 0.9 --recode --recode-INFO-all -c |bgzip > $dirvcf/${chr}_${from}_$to.vcf.gz
done < $dirlist/local_PCA_MDS_outlier.bed

```

Run `PLINK` for PCA.
```bash

while read chr from to
do
        sbatch $dirscripts/plink_pca.sh  $dirvcf/${chr}_${from}_$to $dirout/${chr}_${from}_$to
done < $dirlist/local_PCA_MDS_outlier.bed

```


Plot PCA results.
```bash

Rscript $dirscripts/plot_pca_per_outlier.R --dirpca $dirout --poplist $dirlist/id_spp_pop_site_pheno.tsv --outlierlist $dirlist/local_PCA_MDS_outlier.bed --chrlist $dirlist/chromosomes.list

```

![](PCA/output/PCA_localPCA_outlier.png)



Get Eigenvalues.

```bash

while read chr pos1 pos2
do
        awk -v chr=$chr -v pos1=$pos1 -v pos2=$pos2 'BEGIN{printf "%s %s %s ", chr,pos1,pos2}{printf "%s ",$1}END{print ""}' $dirout/${chr}_${pos1}_${pos2}.eigenval
done < $dirlist/local_PCA_MDS_outlier.bed | awk 'NR==1{printf "%s %s %s ","chr","pos1","pos2";for(i=4;i<=NF;i++){printf "%s ","PC"i-3};print ""}{print $0}'> $dirout/eigenvalues.txt


```

Plot Eigenvalues.

```bash

Rscript $dirscripts/plot_eigenval.R --dirpca $dirout 

```

![](PCA/output/eigenval.png)



## Class-1 genomic islands

```bash

dirbase=$PWD/class-1
dirin=$dirbase/input
dirvcf=$dirin/vcf
dirmask=$dirin/mask
dirlist=$dirbase/list
dirout=$dirbase/output
dirscripts=$dirbase/scripts

```

### Heterozygosity


Make link to VCF of chromosome 12, as an example chromosome harbouring a class-1 genomic island.
```bash
ln $dirbase/../PCA/input/vcf/chr_12.vcf.gz  $dirvcf
ln $dirbase/../PCA/input/vcf/chr_12.vcf.gz.csi  $dirvcf

```

Get coordinates of class-1 genomic island of chromosome 12.

```bash

awk -v OFS="\t" '$1=="chr_12"{$1=$1;print $0}' $dirbase/../local_PCA/output/local_PCA_MDS/local_PCA_MDS_outlier.bed > $dirlist/class-1.chr_12.bed

```

Make list of individuals with AA, AB, and BB. 

```bash

awk '$3<0{print $1}' $dirbase/../PCA/output/chr_12_14126710_22227355.eigenvec > $dirlist/chr_12.AA.list
awk '$3>0.15{print $1}' $dirbase/../PCA/output/chr_12_14126710_22227355.eigenvec > $dirlist/chr_12.BB.list
awk '$3>0&&$3<0.15{print $1}' $dirbase/../PCA/output/chr_12_14126710_22227355.eigenvec > $dirlist/chr_12.AB.list

```


Compute heterozygosity.
```bash

while read chr pos1 pos2
do
        for geno in AA AB BB
        do
                while read id
                do
                        # echo $id $chr $geno
                        bcftools query -f '[%GT ]\n' -r $chr:$pos1-$pos2 -s $id $dirvcf/$chr.vcf.gz |sed 's@0/0@0@g;s@0|0@0@g;s@1/1@0@g;s@1|1@0@g;s@0/1@1@g;s@0|1@1@g;s@1|0@1@g;s@\./\.@@g;s@\.|\.@@g' | awk -v chr=$chr -v geno=$geno -v id=$id '{i++;s+=$1}END{print chr,geno,id,s,i,s/i}'
                done<$dirlist/${chr}.$geno.list
        done 
done < $dirlist/class-1.chr_12.bed | awk 'BEGIN{print "chr","geno","id","n.het","n.sites","het"}{print $0}' > $dirout/chr_12.het.txt


```

Plot heterozygosity for AA, AB, and BB.

```bash

Rscript $dirscripts/plot_het.R --dir $dirout

```

![](class-1/output/chr_12.het.png)

### F<sub>ST</sub>, d<sub>XY</sub> and π

Copy `chromosomes_length.list`.

```bash
cp $dirbase/../local_PCA/list/chromosomes_length.list $dirlist

```


Make `chr_12_IDgeno.list`.

```bash

while read chr pos1 pos2
do
        for geno in AA BB
        do
                awk -v geno=${chr}.${geno} '{print $1,geno}' $dirlist/${chr}.$geno.list 
        done |awk 'BEGIN{print "sample","population"}{print $0}' > $dirlist/${chr}_IDgeno.list
done<$dirlist/class-1.chr_12.bed

```

Tabix input VCF.
```bash
tabix $dirvcf/chr_12.vcf.gz

```


Run PopGenome.
`class-1/scripts/PopGenome_windowstats.sh` submits `class-1/scripts/PopGenome_windowstats.R` via slurm.
Check `class-1/scripts/PopGenome_windowstats.sh` and `class-1/scripts/PopGenome_windowstats.R` for detail.

```bash

win=10000
while read chr pos1 pos2
do
sbatch $dirscripts/PopGenome_windowstats.sh \
$dirvcf/$chr.vcf.gz \
$dirlist/chromosomes_length.list \
$chr \
$dirlist/${chr}_IDgeno.list \
$win \
$dirout/blackcap.$chr.AA.BB
done<$dirlist/class-1.chr_12.bed


```


Plot the results.
```bash
module load R/3.5.3
Rscript $dirscripts/plot_windowstats_class-1.R --dir $dirout --bedfile $dirlist/class-1.chr_12.bed

```

![](class-1/output/chr_12.class-1.windowstats.png)


Run permutation tests.

```bash

while read chr pos1 pos2
do
        if [ $chr != "chr_6" ];then
                Rscript $dirscripts/permutation_windowstats.R --chr $chr --pos1 $pos1 --pos2 $pos2 --n 10000 --pifile $dirout/blackcap.$chr.AA.BB_pi_PopGenome10kb.txt --dxyfile $dirout/blackcap.$chr.AA.BB_dxy_PopGenome10kb.txt --fstfile $dirout/blackcap.$chr.AA.BB_FST_PopGenome10kb.txt
        fi
done<$dirlist/class-1.chr_12.bed | awk 'BEGIN{print "chr","pos.from","pos.to","p.val_piAA","p.val_pi.BB","p.val_dxy","p.val_FST","p.val_pi1-pi2"}{print $0}' | sed 's/ /,/g' > $dirout/permutation_windowstats_class1.csv

```

The result is written in [`class-1/output/permutation_windowstats_class1.csv`](class-1/output/permutation_windowstats_class1.csv).




### Coalescent time (empirical)

Samples used for coalescent time analysis are listed in [`class-1/list/chr_12.AA.4samples.list`](class-1/list/chr_12.AA.4samples.list), [`class-1/list/chr_12.AB.4samples.list`](class-1/list/chr_12.AB.4samples.list) and [`class-1/list/chr_12.BB.4samples.list`](class-1/list/chr_12.BB.4samples.list).


Mask files were created using `generate_multihetsep.py` following <https://github.com/stschiff/msmc-tools> from alignment files (BAM) and reference file (FASTA).
Precomputed mask files are found in [`class-1/input/mask/`](class-1/input/mask/).


Split VCF file into individuals.

```bash

for geno in AA AB BB
do
        while read id
        do
                echo $geno $id
                bcftools view -s $id -M2 -m2 $dirvcf/chr_12.vcf.gz | vcftools --vcf - --mac 1 --recode -c | bgzip > $dirvcf/$id.chr_12.vcf.gz
        done < $dirlist/chr_12.$geno.4samples.list 
done

```


Make `run_generate_multihetsep.sh` files.
```bash
$dirscripts/write_run_generate_multihetsep_per_group.sh $dirbase

```






