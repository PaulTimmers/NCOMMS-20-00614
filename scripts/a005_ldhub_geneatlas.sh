#!/bin/bash

#$ -N GeneAtlas
#$ -l h_vmem=8G
#$ -l h_rt=00:30:00
#$ -j y
#$ -o ${ROOT}/projects/prj_108_lifegen_phase3/p06_triple_manova/q05_ldsc/ldhub/t006_geneatlas.$TASK_ID.log
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

download=`awk -v id=$SGE_TASK_ID 'NR == id+1 {print $1}' t000_geneatlas_traits.txt`
trait=`awk -v id=$SGE_TASK_ID 'NR == id+1 {print $2}' t000_geneatlas_traits.txt`


#----
# Download data
#----

wget -nv "ftp://ftp.igmm.ed.ac.uk/pub/GeneATLAS/${download}.v2.tar" -O ${download}.v2.tar
mkdir -p ${download}
tar --overwrite -xvf ${download}.v2.tar -C ${download} && rm ${download}.v2.tar

# INFO > 0.9

echo -e "SNP\tA1\tBETA\tSE\tP" > t001_${trait}.tsv
for CHR in {1..22}
do
zcat ${download}/results/${download}/imputed.allWhites.${download}.*chr${CHR}.csv.gz | \
awk -v OFS="\t" 'NR == 1 {next} $3 > 0.9 {print $1,$2,$4,$5,$6}' >> t001_${trait}.tsv
done

rm -r ${download}


#----
# Subset SNPs 
#----


# HapMap3

awk -v OFS="\t" \
'NR == 1 {print "SNP","A1","A2","N","Z","P"; next} 
 NR == FNR {a1[$1]=$2; a2[$1]=$3; next} 
 a1[$1] == $2 || a2[$1] == $2 {printf "%s\t%s\t%s\t%i\t%.6g\t%.6g\n",$1,$2,(a1[$1]==$2 ? a2[$1] : a1[$1]),452264,$3/$4,$5}' \
 <(zcat t005_healthspan.sumstats.zip) t001_${trait}.tsv > t002_${trait}.tsv

# Remove MHC

grep -wvFf mhc_snps.tsv t002_${trait}.tsv | gzip -9 -f > t005_${trait}.sumstats.gz && rm t002_${trait}.tsv


#----
# Perform rg regression
#----

#ldsc.py \
#--rg t005_${trait}.sumstats.gz,t005_healthspan.sumstats.gz,t005_lifespan.sumstats.gz,t005_longevity.sumstats.gz \
#--ref-ld-chr $reference \
#--w-ld-chr $weights \
#--out t006_${trait}


#----
# Format results
#----

#cat t006_${trait}.log | grep "p1" -A 3 | sed 's/^[ ]*//g' | tr -s ' ' '\t' | \
#awk -v OFS="," -v id=$SGE_TASK_ID \
#'BEGIN{traits[1]="header"; traits[2]="healthspan"; traits[3]="lifespan"; traits[4]="longevity"}
# NR == 1 {next}
# {gsub("t005_","",$1); gsub(".sumstats.gz","",$1)}
# {print traits[NR],$1,30349118,"GeneATLAS","European","SNPs from the MHC (chr6 26M~34M) region was removed for this traits",$3,$4,$5,$6,$7,$8,$9,$10,$11,$12 >> "t006_"traits[NR]"_geneatlas."id".csv"}'


