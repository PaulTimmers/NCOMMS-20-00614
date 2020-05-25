#!/bin/bash

#$ -N run_rg
#$ -l h_vmem=8G
#$ -l h_rt=00:30:00
#$ -j y
#$ -o ${ROOT}/projects/prj_108_lifegen_phase3/p06_triple_manova/q05_ldsc/ldhub/t006_rg.$TASK_ID.log
#$ -cwd
#$ -V

truncate -s 0 $SGE_STDOUT_PATH

#----
# Setup environment
#----

MAIN=${ROOT}/projects/prj_108_lifegen_phase3/p06_triple_manova
cd ${MAIN}/q05_ldsc/ldhub

source ../venv/bin/activate

weights=${ROOT}/apps_by_us_full_stuff/ldsc/data/eur_w_ld_chr/
reference=${ROOT}/apps_by_us_full_stuff/ldsc/data/eur_w_ld_chr/
snp_list=${ROOT}/apps_by_us_full_stuff/ldsc/data/w_hm3.snplist

n_traits=`wc -l < t006_traits.tsv`
trait=`awk -v id=$SGE_TASK_ID 'NR == id {print $1}' t006_traits.tsv`
pmid=`awk -v id=$SGE_TASK_ID 'NR == id {print $2}' t006_traits.tsv`
category=`awk -v id=$SGE_TASK_ID 'NR == id {print $3}' t006_traits.tsv`


#----
# Perform rg regression
#----

rg_list=`awk -v trait=$trait '$1 != trait {print "t005_"$1".sumstats.gz"}' t006_traits.tsv | xargs printf ",%s"`

ldsc.py \
--rg t005_${trait}.sumstats.gz${rg_list} \
--ref-ld-chr $reference \
--w-ld-chr $weights \
--out t006_${trait}


#----
# Format results
#----

h2=`grep -A 2 "Heritability of phenotype 1" t006_${trait}.log | awk 'NR == 3 {print $(NF-1)}'`
h2_se=`grep -A 2 "Heritability of phenotype 1" t006_${trait}.log | awk -v FS="[ |(|)]" 'NR == 3 {print $(NF-1)}'`

cat t006_${trait}.log | grep "p1" -A `echo $n_traits-1 | bc` | sed 's/^[ ]*//g; s/t005_//g; s/.sumstats.gz//g' | tr -s ' ' '\t' \
| cut -f1-6 | awk -v category=$category -v pmid=$pmid -v h2=$h2 -v h2_se=$h2_se -v OFS="\t" 'NR == 1 {print "category","pmid","trait1","trait2","rg","se","z","p","h2","h2_se"; next} {print category,pmid,$0,h2,h2_se,n_snp}' > t007_${trait}.tsv
awk -v OFS="\t" 'FNR == NR {n_snp[FNR]=$0; next} {print $0,n_snp[FNR]}' <(grep "SNPs with valid alleles" t006_${trait}.log | awk 'BEGIN{print "n_snp"} {print $1}') t007_${trait}.tsv | tr -s ' ' '\t' > ${trait}.tmp && mv ${trait}.tmp t007_${trait}.tsv

#----
# Log
#----

mv $SGE_STDOUT_PATH t006_${trait}.log