#!/bin/bash

#----
# Setup environment
#----

MAIN=${ROOT}/projects/prj_108_lifegen_phase3/p06_triple_manova
cd $MAIN

mkdir -p q04_diseases
cd q04_diseases


#----
# Get list of close proxies
#----

mkdir -p proxy_info

tail -n+2 ../q03_analyse/st03_03_triple_manova_res.tsv | cut -f1 | while read RSID
do
echo -e "\n\n===== ${RSID} ====="
curl --silent -k -X GET "https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?var=${RSID}&pop=EUR&r2_d=r2&token=2183e2d5c529" | \
awk -v rsid=$RSID -v OFS="\t" -v FS="\t|,|=" 'NR == 1 {print "rsid","a1","a0","proxy","p.a1","p.a0","r2"; next} $1 ~ /rs/ && $7 >= 0.6 {print rsid,$8,$10,$1,$9,$11,$7}' | gzip -c | tee proxy_info/st04_01_${RSID}_info.tsv.gz | zcat | column -t
done

zcat proxy_info/st04_01_*_info.tsv.gz | cut -f4 | grep -v "proxy" | awk '{i=int((NR-1)/100)+1} {print $0 > "proxy_info/st04_00_rsids."i".txt"}' 


#----
# Get phenoscanner information
#----


Rscript ../scripts/a004_get_phenoscanner_data.R &> log_a004_get_phenoscanner_data.log


#----
# Create table
#----

Rscript ../scripts/a004_parse_gwas_catalog.R &> log_a004_format_gwas_catalog.log &


#----
# Create heatmap
#----

Rscript ../scripts/a004_05_disease_heatmap.R