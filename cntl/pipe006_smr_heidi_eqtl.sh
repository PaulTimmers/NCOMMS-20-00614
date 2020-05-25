#!/bin/bash

#----
# Setup environment
#----

MAIN=${ROOT}/projects/prj_108_lifegen_phase3/p06_triple_manova
cd $MAIN

mkdir -p q06_smr_heidi
cd q06_smr_heidi


#----
# Format GWAS data
#----

TRIPLEMANOVA=../q02_manova/st02_01_triple_manova.tsv.gz


mkdir -p manova

awk -v flank=500000 -v OFS="\t" \
'BEGIN{OUT[0]="init"; OUT[1]="health"; OUT[2]="life"; OUT[3]="longevity"; OUT[4]="manova"}
 BEGINFILE{print OUT[ARGIND-1],(ARGIND-1); chr=0}
 
 NR == 1 {next}
 NR == FNR {locus[$2][$3] = $1; next}
 
 (ARGIND-1) == 1 && $3 in locus {for (i in locus[$3]) { if($4 >= i-flank && $4 <= i+flank) {health_a1[$2]=$5; health_beta1[$2]=$9; health_se[$2]=$10; health_n[$2]=$7}}}
 (ARGIND-1) == 1 {next}

 (ARGIND-1) == 2 && $3 in locus {for (i in locus[$3]) { if($4 >= i-flank && $4 <= i+flank) {life_a1[$2]=$5; life_beta1[$2]=$9; life_se[$2]=$10; life_n[$2]=$7}}}
 (ARGIND-1) == 2 {next}

 (ARGIND-1) == 3 && $3 in locus {for (i in locus[$3]) { if($4 >= i-flank && $4 <= i+flank) {long_a1[$2]=$5; long_beta1[$2]=$9; long_se[$2]=$10; long_n[$2]=$7}}}
 (ARGIND-1) == 3 {next}


 FNR == 1 {for (i in locus) {for (j in locus[i]) {print "SNP","CHR","POS","A1","A2","freq","b","se","p","health_a1","health_beta1","health_se","health_n","life_a1","life_beta1","life_se","life_n","long_a1","long_beta1","long_se","long_n" > OUT[ARGIND-1]"/st006_01_"locus[i][j]"_"OUT[ARGIND-1]".tsv"}}}
 FNR == 1 {next}

 $3 in locus {for (i in locus[$3]) { if($4 >= i-flank && $4 <= i+flank) {print $1,$3,$4,$5,$6,$8,$9,$10,$11,health_a1[$2],health_beta1[$2],health_se[$2],health_n[$2],life_a1[$2],life_beta1[$2],life_se[$2],life_n[$2],long_a1[$2],long_beta1[$2],long_se[$2],long_n[$2] > OUT[ARGIND-1]"/st006_01_"locus[$3][i]"_"OUT[ARGIND-1]".tsv"}}}
 FNR % 100000 == 0 {printf "%02i\t%i\r",$3,FNR}
 ENDFILE{if((ARGIND-1) > 3){print $3,FNR}}' \
../q03_analyse/st03_03_triple_manova_res.tsv  <(zcat $TRIPLEMANOVA)

Rscript ../scripts/a006_01_format_smr.R





#----
# Gather eQTL datasets
#----

if [[ ! -f smr_Linux ]]; then
    wget http://cnsgenomics.com/software/smr/download/smr_Linux.zip
    unzip smr_Linux.zip && rm smr_Linux.zip
fi


if [[ ! -f vosa_eqtl_full.besd ]]; then
    wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/cis-eQTLs_full_20180905.txt.gz
    wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz

    n_lines=`zcat cis-eQTLs_full_20180905.txt.gz | wc -l`

    awk -v n_lines=$n_lines -v OFS="\t" \
    'NR == FNR {alleleA[$1]=$4; alleleB[$1]=$5; maf[$1]=$9; next} 
     FNR == 1 {print "SNP", "Chr", "BP", "A1", "A2", "Freq", "Probe", "Probe_Chr", "Probe_bp", "Gene", "Orientation", "b", "se", "p" > "vosa_eqtl_full.txt"; next}
     maf[$2] > 0 && maf[$2] != "NA" && (alleleA[$2] == $6 || alleleB[$2] == $6) && (alleleA[$2] == $7 || alleleB[$2] == $7) {{if(alleleB[$2] != $6) {maf[$2]=1-maf[$2]; alleleB[$2]=$6}}; print $2, $3, $4, $6, $7, maf[$2], $8, $10, $11, $9, "-",$5 / sqrt( 2*maf[$2]*(1-maf[$2])*($13+$5^2)), 1 / sqrt(2*maf[$2]*(2-maf[$2])*($13 + $5^2)), $1 > "vosa_eqtl_full.txt"}
     FNR % 100000 == 0 {printf "%.1f%\r",FNR/n_lines*100}
     END{print "100.0%"}' \
     <(zcat 2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz) <(zcat cis-eQTLs_full_20180905.txt.gz)

    (./smr_Linux --qfile vosa_eqtl_full.txt --make-besd --out vosa_eqtl_full &> log_make_besd_vosa_eqtl_full.log && rm vosa_eqtl_full.txt) &

    rm cis-eQTLs_full_20180905.txt.gz 2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz
fi


if [[ ! -f vosa_trans_full.besd ]]; then
    wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/trans-eqtl/trans-eQTLs_full_20180905.txt.gz
    wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz

    n_lines=`zcat trans-eQTLs_full_20180905.txt.gz | wc -l`

    awk -v n_lines=$n_lines -v OFS="\t" \
    'NR == FNR {alleleA[$1]=$4; alleleB[$1]=$5; maf[$1]=$9; next} 
     FNR == 1 {print "SNP", "Chr", "BP", "A1", "A2", "Freq", "Probe", "Probe_Chr", "Probe_bp", "Gene", "Orientation", "b", "se", "p" > "vosa_trans_full.txt"; next}
     maf[$2] > 0 && maf[$2] != "NA" && (alleleA[$2] == $6 || alleleB[$2] == $6) && (alleleA[$2] == $7 || alleleB[$2] == $7) {{if(alleleB[$2] != $6) {maf[$2]=1-maf[$2]; alleleB[$2]=$6}}; print $2, $3, $4, $6, $7, maf[$2], $8, $10, $11, $9, "-",$5 / sqrt( 2*maf[$2]*(1-maf[$2])*($13+$5^2)), 1 / sqrt(2*maf[$2]*(2-maf[$2])*($13 + $5^2)), $1 > "vosa_trans_full.txt"}
     FNR % 100000 == 0 {printf "%.1f%\r",FNR/n_lines*100}
     END{print "100.0%"}' \
    <(zcat 2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz) <(zcat trans-eQTLs_full_20180905.txt.gz)

    (./smr_Linux --qfile vosa_trans_full.txt --make-besd --out vosa_trans_full &> log_make_besd_vosa_trans_full.log && rm vosa_trans_full.txt) &

    rm trans-eQTLs_full_20180905.txt.gz 2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz
fi


if [[ ! -f westra_eqtl_hg19.besd ]]; then
    wget http://cnsgenomics.com/data/SMR/westra_eqtl_hg19.zip -O st006_cis_eqtl_westra.zip
    unzip st006_cis_eqtl_westra.zip && rm st006_cis_eqtl_westra.zip
    awk -v OFS="\t" '$5 ~/C*ORF/ {sub("ORF","orf",$5)} 1' westra_eqtl_hg19.epi > tmp && mv tmp westra_eqtl_hg19.epi # Coordinate open reading frame transcript names with other datasets
fi


if [[ ! -f cage_eqtl.summary ]]; then
    wget http://cnsgenomics.com/data/SMR/cage_eqtl_data_hg19.tgz -O st006_cis_eqtl_cage.tar.gz
    utar st006_cis_eqtl_cage.tar.gz
    for ext in besd epi esi summary
    do
         mv cage_eqtl_data/CAGE.sparse.${ext} cage_eqtl_hg19.${ext}
    done
    rm -r cage_eqtl_data
    awk -v OFS="\t" '$5 == "NA" {$5=$2} 1' cage_eqtl_hg19.epi > cage_eqtl_hg19.tmp.epi && mv cage_eqtl_hg19.tmp.epi cage_eqtl_hg19.epi
fi


if [[ ! -d gtex_eqtl_5pct ]]; then
    wget "https://cnsgenomics.com/data/SMR/GTEx_V7_cis_eqtl_summary_lite.tar.gz" -O getx_eqtl_5pct.tar.gz
    utar getx_eqtl_5pct.tar.gz
    mv GTEx_V7_cis_eqtl_summary_lite gtex_eqtl_5pct

    for file in `echo gtex_eqtl_5pct/*.summary`
    do
    file=`echo $file | sed 's/.summary//'`
    tissue=`echo $file | sed 's=^.*/==; s/_1e-05//; s/_/-/g' | tr '[:upper:]' '[:lower:]'`

    for suffix in besd epi esi summary
    do
    ln -sf $file.${suffix} gtex-${tissue}.${suffix}
    done

    done


fi



#----
# Sensitivity
#----

# Perform the following analysis once
# Move results (st006_0*.smr) to effective/
# Rename manova/*_manova.tsv.simple to manova/*_manova.tsv
# Redo the analysis below with the new files
# Move results to simple/


#----
# Run SMR
#----

rm -rf logs
mkdir -p logs

awk '{print $1,$1}' ${ROOT}/data/reference_data/1000genomes/meta_data/1000genomes_europeans.tsv > EUR_id_list.txt
n_rsid=`awk 'END{print NR-1}' ../q03_analyse/st03_03_triple_manova_res.tsv`

for eqtl_data in  cage_eqtl_hg19 `echo gtex-*.summary | sed 's/.summary//g'` vosa_eqtl_full vosa_trans_full westra_eqtl_hg19 
do
    echo "$eqtl_data" | grep -q "vosa_trans" && eqtl_data_type="trans" || eqtl_data_type="cis"
    for trait in manova 
    do
        qsub -N ${eqtl_data:0:4}.${eqtl_data_type:0:1}.smr.${trait} -j y -o logs/log_a006_02_run_smr_${trait}.${eqtl_data:0:4}.${eqtl_data_type}.\$TASK_ID.log -t 1-${n_rsid} ../scripts/a006_02_run_smr.sh ${trait} ${eqtl_data} ${eqtl_data_type}
    done
done


# Review and remove frequency mismatch files

for file in st006_02_*_*_*.snp_failed_freq_ck.list
do
    n=`wc -l < $file`
    if [[ $n -gt 1 ]]; then
        awk '{print $0,(NR == 1 ? "file" : FILENAME)}' $file | column -t && echo ""
    fi
    rm $file
done




#----
# Format results
#----

awk -v OFS="\t" \
'BEGINFILE{N=split(FILENAME,INFO,"[_|.]")}
 INFO[(N-1)] ~ /trans/ && (NR == 1 || ($22 < 1 && $23 > 0)) {print $3,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,(NR == 1 ? "rsid" : INFO[3]),(NR == 1 ? "trait" : INFO[4]),(NR == 1 ? "eqtls" : INFO[5]); next} 
 INFO[(N-1)] !~ /trans/ && (NR == 1 || ($19 < 1 && $20 > 0)) {print $3,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,(NR == 1 ? "rsid" : INFO[3]),(NR == 1 ? "trait" : INFO[4]),(NR == 1 ? "eqtls" : INFO[5])}' st006_02_rs*.smr | sort -k10r,10 -k7g,7 > st006_02_smr_results.txt

Rscript ../scripts/a006_03_smr_format_results.R &> log_st006_03_smr_format_results.log

