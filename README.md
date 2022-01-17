# Methods for "Recombination suppression and selection affect local ancestries in genomes of a migratory songbird"

# Local PCA

Here I show how local PCA was run with three example chromosomes, 12, 20, and 21, which are listed in `local_PCA/list/chromosomes/list`

```bash
dirbase=local_PCA
dirin=$dirbase/input
dirvcf=$dirin/vcf
dirtab=$dirin/table
dirlist=$dirbase/list
dirout=$dirbase/output

```


```bash
while read chr
do
        echo $chr
        bcftools query -f '%POS\n' $dirvcf/$chr.vcf.gz > $dirtab/$chr.sites.list
done <$dirlist/chromosomes.list

```



# 



