#!/bin/bash

#----
# Setup environment
#----

MAIN=${ROOT}/projects/prj_108_lifegen_phase3/p06_triple_manova
mkdir -p ${MAIN}/q10_age_stratified_vars
cd ${MAIN}/q10_age_stratified_vars

#----
# START
#----

# Extract loci that are nominally significant in all datasets but genome-wide significant in none

awk -v FS="," 'NR == FNR && $8 > 5e-8 && $11 > 5e-8 && $14 > 5e-8 {gene[$1]++; next}
               FNR ==1 || gene[$1] > 0 && $10 < 0.05 && $13 < 0.05 && $16 < 0.05 {print}' \
../q03_analyse/st03_06_triple_manova_region_table.csv  ../q03_analyse/st03_06_triple_manova_table.csv | sed 's/,/\t/g' > st10_01_novel_loci.tsv


# Extract age band stats for SNPs of interest 

zgrep -wFf <(cut -f2 st10_01_novel_loci.tsv) ../q09_age_stratified_gwas/q01_combined/p05_comb_chr/bothpl_lif4060.tsv.gz | cut -f1-11,13 > st10_02_bothpl_lif4060_stats.tsv &
zgrep -wFf <(cut -f2 st10_01_novel_loci.tsv) ../q09_age_stratified_gwas/q01_combined/p05_comb_chr/bothpl_lif6080.tsv.gz | cut -f1-11,13 > st10_02_bothpl_lif6080_stats.tsv &
zgrep -wFf <(cut -f2 st10_01_novel_loci.tsv) ../q09_age_stratified_gwas/q01_combined/p05_comb_chr/bothpl_lif80120.tsv.gz | cut -f1-11,13 > st10_02_bothpl_lif80120_stats.tsv &

for par in fath moth
do
    zgrep -wFf <(cut -f2 st10_01_novel_loci.tsv) ../q09_age_stratified_gwas/p05_comb_chr/${par}_lif4060.tsv.gz | cut -f1-11,13 > st10_02_${par}_lif4060_stats.tsv &
    zgrep -wFf <(cut -f2 st10_01_novel_loci.tsv) ../q09_age_stratified_gwas/p05_comb_chr/${par}_lif6080.tsv.gz | cut -f1-11,13 > st10_02_${par}_lif6080_stats.tsv &
    zgrep -wFf <(cut -f2 st10_01_novel_loci.tsv) ../q09_age_stratified_gwas/p05_comb_chr/${par}_lif80120.tsv.gz | cut -f1-11,13 > st10_02_${par}_lif80120_stats.tsv &

done


wait

# Extract all loci

cat ../q03_analyse/st03_06_triple_manova_table.csv | sed 's/,/\t/g' | awk 'NR ==1 || $10 < 0.05 && $13 < 0.05 && $16 < 0.05' > st10_04_all_loci.tsv

(zcat ../q09_age_stratified_gwas/q01_combined/p05_comb_chr/bothpl_lif4060.tsv.gz | head -1 
zgrep -wHFf <(cut -f2 st10_04_all_loci.tsv | tail -n+2) ../q09_age_stratified_gwas/q01_combined/p05_comb_chr/bothpl_life.tsv.gz &&
zgrep -wHFf <(cut -f2 st10_04_all_loci.tsv | tail -n+2) ../q09_age_stratified_gwas/q01_combined/p05_comb_chr/bothpl_lif4060.tsv.gz &&
zgrep -wHFf <(cut -f2 st10_04_all_loci.tsv | tail -n+2) ../q09_age_stratified_gwas/q01_combined/p05_comb_chr/bothpl_lif6080.tsv.gz &&
zgrep -wHFf <(cut -f2 st10_04_all_loci.tsv | tail -n+2) ../q09_age_stratified_gwas/q01_combined/p05_comb_chr/bothpl_lif80120.tsv.gz) | cut -f1-11,13 > st10_05_bothpl_all.tsv

awk -v OFS="\t" 'NR == 1 {print $0,"par","band"; next}
     {par="bothpl"; band=$1; sub(".tsv.gz.*","",band); sub(".*lif","",band);
     sub(".*:","",$1); print $0, par, substr(band,0,2)"_"substr(band,3,9)}' st10_05_bothpl_all.tsv > bothpl.tmp && mv bothpl.tmp st10_05_bothpl_all.tsv

for par in fath moth
do
(zcat ../q09_age_stratified_gwas/p05_comb_chr/${par}_lif4060.tsv.gz | head -1 && 
    zgrep -wHFf <(cut -f2 st10_04_all_loci.tsv | tail -n+2) ../q09_age_stratified_gwas/p05_comb_chr/${par}_life.tsv.gz &&
    zgrep -wHFf <(cut -f2 st10_04_all_loci.tsv | tail -n+2) ../q09_age_stratified_gwas/p05_comb_chr/${par}_lif4060.tsv.gz &&
    zgrep -wHFf <(cut -f2 st10_04_all_loci.tsv | tail -n+2) ../q09_age_stratified_gwas/p05_comb_chr/${par}_lif6080.tsv.gz && 
    zgrep -wHFf <(cut -f2 st10_04_all_loci.tsv | tail -n+2) ../q09_age_stratified_gwas/p05_comb_chr/${par}_lif80120.tsv.gz) | cut -f1-11,14 > st10_05_${par}_all.tsv
awk -v OFS="\t" -v par=$par 'NR == 1 {print $0,"par","band"; next} {band=$1; sub(".tsv.gz.*","",band); sub(".*lif","",band); sub(".*:","",$1); print $0, par, substr(band,0,2)"_"substr(band,3,9)}' st10_05_${par}_all.tsv > ${par}.tmp && mv ${par}.tmp st10_05_${par}_all.tsv
done

Rscript ../scripts/a010_06_all_age_bands_rel.R |& tee log_a010_06_all_age_bands_rel.log
