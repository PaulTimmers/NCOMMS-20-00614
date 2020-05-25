#!/bin/bash

#----
# Setup environment
#----

MAIN=${ROOT}/projects/prj_108_lifegen_phase3/p06_triple_manova
cd $MAIN

#----
# START
#----

mkdir -p ${MAIN}/q12_iron_mr
mkdir -p ${MAIN}/q12_iron_mr/z012_mr_results
cd ${MAIN}/q12_iron_mr

#----
# Format multivariate
#----

Rscript ../scripts/a012_01_format_manova.R


#----
# Format iron
#----

# https://www.ncbi.nlm.nih.gov/pubmed/25352340

wget "https://genepi.qimr.edu.au/contents/p/studies/2013/GIS%20Consortium%20-%20Iron%20data%20-%20Benyamin%20et%20al%202014.zip" -O st012_02_iron.zip
unzip st012_02_iron.zip && rm st012_02_iron.zip && rm *Replication.txt && rm *.doc && mv Iron_Benyamin_NatCom_2014.txt st012_02_iron.tsv && mv Log10Ferritin_Benyamin_NatCom_2014.txt st012_02_ferritin.tsv && mv Saturation_Benyamin_NatCom_2014.txt st012_02_saturation.tsv && mv Transferrin_Benyamin_NatCom_2014.txt st012_02_transferrin.tsv

for iron_trait in iron ferritin saturation transferrin
do
awk -v OFS="\t" \
'NR == 1 {print "rsid","snpid","chr","pos","a1","a0","n","freq1","beta1","se","p"; next}
{print $1,$2"_"$3,$2,$3,$4,$5,48972,$6,$7,$8,$9}' st012_02_${iron_trait}.tsv | gzip -9 > st012_02_${iron_trait}.tsv.gz \
&& rm st012_02_${iron_trait}.tsv
done



#----
# Get lead SNPs
#----

iron_traits="iron ferritin saturation transferrin"


for haem_trait in $iron_traits
do
(
# Look for overlap between manova and haem trait
awk 'NR == FNR {rsid[$1]++; snpid[$2]++; next} FNR == 1 || ($11 < 5e-8 && (rsid[$1] > 0 || snpid[$2] > 0))' \
<(zcat st012_01_manova.tsv.gz) <(zcat st012_02_${haem_trait}.tsv.gz) > st012_03_${haem_trait}_hits.tmp

# Find lead SNPs by distance
gwas_extract_hits_and_lz.sh -m 0 -q -p 5e-08 st012_03_${haem_trait}_hits.tmp -e st012_03_${haem_trait}_hits.tsv \
&& rm st012_03_${haem_trait}_hits.tmp
) &
done

wait


# Merge lead SNPs

cut -f1-11 st012_03_*_hits.tsv | sort -ur > st012_04_all_hits.tsv # combine all genome-wide significant hits
zgrep -wFf <(cut -f1 st012_04_all_hits.tsv | sort -ur) st012_01_manova.tsv.gz > st012_04_manova_all_univariate_hits.tsv


# Run univariate MR

iron_traits="iron ferritin saturation transferrin"
Rscript ../scripts/a012_04_two_sample_mr_univariate.R $iron_traits &> log_a012_04_univariate.log



#----
# Merge lead SNPs
#----


# By distance

mv_iron_traits="iron ferritin saturation transferrin"
truncate -s 0 st012_05_all_hits.tsv
for mv_iron_trait in ${mv_iron_traits}
do
    cut -f1-11 st012_03_${mv_iron_trait}_hits.tsv >> st012_05_all_hits.tsv
    sort -ur st012_05_all_hits.tsv > st012_05_all_hits.tmp && mv st012_05_all_hits.tmp st012_05_all_hits.tsv
done

gwas_extract_hits_and_lz.sh -m 0 -q st012_05_all_hits.tsv -e st012_05_all_hits_clumped.tsv # subset to independent loci




#----
# Extract SNPs
#----

# Manova 
zgrep -wFf <(cut -f1 st012_05_all_hits_clumped.tsv | sort -ur) st012_01_manova.tsv.gz  > st012_06_manova_all_hits.tsv


for haem_trait in $mv_iron_traits 
do
zgrep -wFf <(cut -f1 st012_06_manova_all_hits.tsv | sort -ur) st012_02_${haem_trait}.tsv.gz > st012_06_${haem_trait}_all_hits.tsv &
done
wait



#----
# Perform MVMR analysis
#----

mv_iron_traits="iron ferritin saturation transferrin"
Rscript ../scripts/a012_07_two_sample_mr_multivariate.R $mv_iron_traits &> log_a012_07_multivariate.log

