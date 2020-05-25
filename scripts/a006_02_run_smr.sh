#!/bin/bash

#$ -cwd
#$ -l h_vmem=32G
#$ -l h_rt=01:00:00

truncate -s 0 $SGE_STDOUT_PATH

#----
# Setup environment
#----

. /etc/profile.d/modules.sh
module load igmm/apps/R/3.6.0
module load igmm/apps/plink/1.90b4

trait=$1
eqtl_data=$2
eqtl_data_type=$3
rsid=`awk -v i=$SGE_TASK_ID 'NR == i+1 {print $1}' ../q03_analyse/st03_03_triple_manova_res.tsv`
chr=` awk -v i=$SGE_TASK_ID 'NR == i+1 {print $2}' ../q03_analyse/st03_03_triple_manova_res.tsv`

eqtl_name=`echo $eqtl_data | sed 's/_.*$//'`

if [[ "${eqtl_data_type}" == "trans" ]]; then
    trans="--trans"
else
    trans=""
fi


#----
# START
#----

sumstats=${trait}/st006_01_${rsid}_${trait}.tsv
	
echo -e "\n\n#-- $trait $rsid ----------------------------------------------------------------\n"
    
echo $sumstats

smr_script=./smr_Linux

$smr_script \
--bfile ${ROOT}/data/reference_data/1000genomes/phase3/plink/1000genomes_chr${chr} \
--keep EUR_id_list.txt \
--gwas-summary $sumstats \
--beqtl-summary ${eqtl_data} \
$trans \
--out st006_02_${rsid}_${trait}_${eqtl_name}${trans/--/-} \
--thread-num 1


echo -e "\n#--------------------------------------------------------------------\n\n"





