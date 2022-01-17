# Methods for "Recombination suppression and selection affect local ancestries in genomes of a migratory songbird"

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



## 



