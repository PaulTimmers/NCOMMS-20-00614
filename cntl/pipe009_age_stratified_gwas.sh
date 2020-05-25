#!/bin/bash

#----
# Setup environment
#----

MAIN=${ROOT}/projects/prj_108_lifegen_phase3/p06_triple_manova


#----
# Run GWAS
#----

mkdir -p ${MAIN}/q09_age_stratified_gwas
cd ${MAIN}/q09_age_stratified_gwas

# Setup GWAS pipeline
# Edit bespoke.R and analysis_plan.txt
# Run all steps in pipeline

cd ${MAIN}/q09_age_stratified_gwas/p02_qc_phenotypes

Rscript ../../scripts/a009_01_sample_descriptives.R


cd ${MAIN}/q09_age_stratified_gwas/p03_residuals


R --quiet --no-save <<"Rcode"
library(data.table)
resids <- fread("../p03_residuals/st03_01_resids.txt")
cor1 <- cor(resids[,-1], use="pairwise.complete.obs")
cor1 <- cor1[grepl("fath",colnames(cor1)), grepl("moth",rownames(cor1))]
fwrite(data.table(x=rownames(cor1),cor1), "st03_02_resids_cor.txt", sep="\t", quote=F, na="NA")
Rcode


#----
# Meta-analyse 
#----

MAIN=${ROOT}/projects/prj_108_lifegen_phase3/p06_triple_manova

mkdir -p ${MAIN}/q09_age_stratified_gwas/q01_combined/p05_comb_chr
cd ${MAIN}/q09_age_stratified_gwas/q01_combined/p05_comb_chr



# Meta-analyse

age=(4060 6080 80120 e)
for i in {1..4}
do
    age_range=${age[i-1]}
    fath_stats="../../p05_comb_chr/fath_lif${age_range}.tsv.gz"
    moth_stats="../../p05_comb_chr/moth_lif${age_range}.tsv.gz"
    (metal2wilson.sh $fath_stats $moth_stats -o "bothpl_lif${age_range}_ungcd.tsv" &> log_bothpl_lif${age_range}.log) &
done

wait


# Correct for correlation

age=(4060 6080 80120 e)
for i in {1..4}
do
    age_range=${age[i-1]}
    GC=`awk -v n=$i 'NR == n+1 {print 1+$(n+1)}' ../../p03_residuals/st03_02_resids_cor.txt` # Get correlation between residuals
    (metal2wilson.sh bothpl_lif${age_range}_ungcd.tsv -m "GENOMICCONTROL ${GC}\nMINMAXFREQ OFF" -o "bothpl_lif${age_range}.tsv" &>> log_bothpl_lif${age_range}.log && rm bothpl_lif${age_range}_ungcd.tsv) &
done

wait


# Remove duplicates

age=(4060 6080 80120 e)

for i in {1..4}
do
    age_range=${age[i-1]}
    bothpl_stats="../p05_comb_chr/bothpl_lif${age_range}.tsv"
    awk '{rsid[$1]++} rsid[$1] > 1 {print $1}' ${bothpl_stats} > ${bothpl_stats}.duplicated
    (awk -v OFS="\t" 'NR == FNR {rsid[$1]++; next} rsid[$1] == 0 {print}' ${bothpl_stats}.duplicated ${bothpl_stats} | gzip -9 > ${bothpl_stats}.gz && rm ${bothpl_stats}) &
done


# Visualise results

mkdir -p ${MAIN}/q09_age_stratified_gwas/q01_combined/p06_analyse
cd ${MAIN}/q09_age_stratified_gwas/q01_combined/p06_analyse


age=(4060 6080 80120 e)

for i in {1..4}
do
    age_range=${age[i-1]}
    bothpl_stats="../p05_comb_chr/bothpl_lif${age_range}.tsv.gz"
echo -e "gwas_results=${bothpl_stats}
output_file_name=mh_both_lif${age_range}.png
do_thinning=TRUE
p_value_threshold=0.01
maf_threshold=0.0005
sample_rate=100
chromosomes=1-22" > parameter_file_bothpl_lif${age_range}.txt
    plot_manhattan_and_QQ_to_png.R parameter_file_bothpl_lif${age_range}.txt && rm parameter_file_bothpl_lif${age_range}.txt
done

age=(4060 6080 80120 e)
for i in {1..4}
do
    age_range=${age[i-1]}
    bothpl_stats="../p05_comb_chr/bothpl_lif${age_range}.tsv.gz"
    maf_threshold=`zcat $bothpl_stats | head -n 500000 | awk 'NR > 1 && $7 > N {N=$7} END{print 200/N}'`
    gwas_extract_hits_and_lz.sh -p 5e-9 -m $maf_threshold -q $bothpl_stats -e hits_bothpl_lif${age_range}.tsv
    paste <(cut -f1-11,13 hits_bothpl_lif${age_range}.tsv) <(snp2gene.sh -s hits_bothpl_lif${age_range}.tsv | awk '{print $NF}') > tmp && mv tmp hits_bothpl_lif${age_range}.tsv

done



#----
# Phenotypic correlations
#----

MAIN=${ROOT}/projects/prj_108_lifegen_phase3/p06_triple_manova

mkdir -p ${MAIN}/q09_age_stratified_gwas/q02_correlations
cd ${MAIN}/q09_age_stratified_gwas/q02_correlations

Rscript ../../scripts/a009_02_pheno_correlations.R



#----
# Genetic correlations
#----


source ../../q05_ldsc_intercept/venv/bin/activate
# If you need to reinstall the virtual environment do:
#virtualenv-2.7 venv
#source venv/bin/activate
#pip install --upgrade pip
#pip install pandas
#pip install scipy
#pip install bitarray

echo "../../q01_sumstats/st01_01_healthspan.tsv.gz
../q01_combined/p05_comb_chr/bothpl_lif4060.tsv.gz
../q01_combined/p05_comb_chr/bothpl_lif6080.tsv.gz
../q01_combined/p05_comb_chr/bothpl_lif80120.tsv.gz
../../q01_sumstats/st01_01_longevity_90_pct.tsv.gz" > files.txt


weights=${ROOT}/apps_by_us_full_stuff/ldsc/data/eur_w_ld_chr/
reference=${ROOT}/apps_by_us_full_stuff/ldsc/data/eur_w_ld_chr/
snp_list=${ROOT}/apps_by_us_full_stuff/ldsc/data/w_hm3.snplist

for SUMSTATS in `cat files.txt`
do
(outname=`basename ${SUMSTATS} | sed 's/.tsv.*//'`
zcat $SUMSTATS | awk -v snp_list=$snp_list 'NR ==1 {print; next} {while (getline < snp_list) {snp[$1]++; a1[$1]=$2; a0[$1]=$3; next}} snp[$1] > 0 && ((a1[$1] == $5 && a0[$1] == $6) || (a1[$1] == $6 && a2[$1] == $5))' | cut -f1-11,14 > ${outname}.tmp


munge_sumstats.py \
--sumstats ${outname}.tmp \
--out $outname \
--signed-sumstats beta1,0 \
--snp rsid \
--a1 a1 \
--a2 a0 \
--ignore snpid &> ${outname}.log && rm ${outname}.tmp )  &
done

wait

n=`wc -l < files.txt`
for i in `eval echo {1..$n}`
do
ldsc.py \
--rg `awk -v ORS="," -v i=$i '{gsub("^.*/","",$0); gsub("tsv","sumstats",$0)} NR == FNR && FNR == i {print; next} NR != FNR && FNR != i {print}' files.txt files.txt | sed 's/,$//'` \
--ref-ld-chr $reference \
--w-ld-chr $weights \
--out st02_genetic_correlations.${i}
done

echo -e "p1\tp2\trg\tse\tz\tp\th2_obs\th2_obs_se\th2_int\th2_int_se\tgcov_int\tgcov_int_se" > st02_genetic_correlations.tsv
grep -h p1 -A `wc -l < files.txt` --no-group-separator st02_genetic_correlations.*.log | sed 's/^[ ]*//; /^$/d; /p1/d' | tr -s ' ' '\t' >> st02_genetic_correlations.tsv


Rscript ../../scripts/a009_03_plot_correlations.R