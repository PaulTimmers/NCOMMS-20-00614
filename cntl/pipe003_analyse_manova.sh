#!/bin/bash

#----
# Setup environment
#----

MAIN=${ROOT}/projects/prj_108_lifegen_phase3/p06_triple_manova
cd $MAIN


#----
# START
#----

mkdir -p q03_analyse
cd q03_analyse

# Clump by LD

HEALTHSPAN="../q01_sumstats/st01_01_healthspan.tsv.gz"
TRIPLEMANOVA=../q02_manova/st02_01_triple_manova.tsv.gz
LONGEVITY="../q01_sumstats/st01_01_longevity_90_pct.tsv.gz"
LIFESPAN="../q01_sumstats/st01_01_lifespan.tsv.gz"

p_thresh=5e-8
flank=500kb
r2=0.02

gwas_extract_hits_and_lz.sh -q $flank $HEALTHSPAN -e st03_01_healthspan_hits.tsv
gwas_extract_hits_and_lz.sh -q $flank $LIFESPAN -e st03_01_lifespan_hits.tsv
gwas_extract_hits_and_lz.sh -q $flank $LONGEVITY -e st03_01_longevity_hits.tsv
gwas_extract_hits_and_lz.sh -q $flank $TRIPLEMANOVA -e st03_01_triple_manova_hits.tsv

for trait in healthspan lifespan longevity triple_manova
do
bash ../scripts/a003_01_ld_clump.sh st03_01_${trait}_hits.tsv $p_thresh $r2 st03_01_${trait}_hits_r2 &> log_a003_01_ld_clump_${trait}.log
done


# Find gene names
snp2gene.sh -s st03_01_healthspan_hits_r2.tsv > st03_01_healthspan_genes.tsv
snp2gene.sh -s st03_01_lifespan_hits_r2.tsv > st03_01_lifespan_genes.tsv
snp2gene.sh -s st03_01_longevity_hits_r2.tsv > st03_01_longevity_genes.tsv
snp2gene.sh -s st03_01_triple_manova_hits_r2.tsv > st03_02_triple_manova_genes.tsv

awk -v OFS="\t" 'FNR == NR {gene[$1]=$NF; next} {print $1,$3,$4,$11,gene[$1]}' st03_01_healthspan_genes.tsv st03_01_healthspan_hits_r2.tsv > st03_03_healthspan_res.tsv
awk -v OFS="\t" 'FNR == NR {gene[$1]=$NF; next} {print $1,$3,$4,$11,gene[$1]}' st03_01_lifespan_genes.tsv st03_01_lifespan_hits_r2.tsv > st03_03_lifespan_res.tsv
awk -v OFS="\t" 'FNR == NR {gene[$1]=$NF; next} {print $1,$3,$4,$11,gene[$1]}' st03_01_longevity_genes.tsv st03_01_longevity_hits_r2.tsv > st03_03_longevity_res.tsv
awk -v OFS="\t" 'FNR == NR {gene[$1]=$NF; next} {print $1,$3,$4,$11,gene[$1]}' st03_02_triple_manova_genes.tsv st03_01_triple_manova_hits_r2.tsv > st03_03_triple_manova_res.tsv


# Draw manhattan

Rscript ../scripts/a003_04_manhattan.R &> log_manhattan.log &


# Make LocusZoom plots

mkdir -p lz
cd lz

n_snps=`tail -n+2 ../st03_03_triple_manova_res.tsv | wc -l`
qsub -t 1-$n_snps ../../scripts/a003_04_locuszoom.sh
cat `eval echo log_lz_st03_04_manova.{1..$n_snps}.log` > log_lz_st03_04_manova.log && rm log_lz_st03_04_manova.*.log
cd ..



# Extract statistics for table

gwas_extract_hits_and_lz.sh -s st03_01_triple_manova_hits_r2.tsv $HEALTHSPAN -e st03_05_manova_healthspan_stats.tsv
gwas_extract_hits_and_lz.sh -s st03_01_triple_manova_hits_r2.tsv $LIFESPAN -e st03_05_manova_lifespan_stats.tsv
gwas_extract_hits_and_lz.sh -s st03_01_triple_manova_hits_r2.tsv $LONGEVITY -e st03_05_manova_longevity_stats.tsv


wait 

cut -f1-12 st03_05_manova_healthspan_stats.tsv | gzip -9 -f > st03_05_manova_healthspan_stats.tsv.gz && rm st03_05_manova_healthspan_stats.tsv
cut -f1-12 st03_05_manova_lifespan_stats.tsv | gzip -9 -f > st03_05_manova_lifespan_stats.tsv.gz && rm st03_05_manova_lifespan_stats.tsv
cut -f1-12 st03_05_manova_longevity_stats.tsv | gzip -9 -f > st03_05_manova_longevity_stats.tsv.gz && rm st03_05_manova_longevity_stats.tsv


# Extract sex-specific effects

FATH_LIFESPAN="../../../prj_065b_lifegen_phase2/p005_make_alldr/st005_03_cleaned_fath_alldr.tsv"
MOTH_LIFESPAN="../../../prj_065b_lifegen_phase2/p005_make_alldr/st005_03_cleaned_moth_alldr.tsv"

gwas_extract_hits_and_lz.sh -s st03_01_triple_manova_hits_r2.tsv $FATH_LIFESPAN -e st03_05_manova_fath_lifespan_stats.tsv &
gwas_extract_hits_and_lz.sh -s st03_01_triple_manova_hits_r2.tsv $MOTH_LIFESPAN -e st03_05_manova_moth_lifespan_stats.tsv &
wait
cut -f1-11,13 st03_05_manova_fath_lifespan_stats.tsv | gzip -9 -f > st03_05_manova_fath_lifespan_stats.tsv.gz && rm st03_05_manova_fath_lifespan_stats.tsv
cut -f1-11,13 st03_05_manova_moth_lifespan_stats.tsv | gzip -9 -f > st03_05_manova_moth_lifespan_stats.tsv.gz && rm st03_05_manova_moth_lifespan_stats.tsv


# Extract regional statistics for supplementary table

extract_region() {
    flank=250000
    
    gzfread() {
        file -b `readlink -f $1` | grep gzip -q && zcat $1 || cat $1
    }

    awk -v OFS="\t" \
    -v gene=$1 \
    -v rsid=$2 \
    -v chr=$3  \
    -v pos=$4 \
    -v flank=$flank \
    'NR==1 {p=1; next} $3 == chr && $4 >= pos-flank && $4 <= pos+flank && $11 < p \
    {p=$11; answ=sprintf("%s\t%s\t%i\t%i\t%s\t%s\t%i\t%.4g\t%.6g\t%.6g\t%.6g",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11)} END{print gene,rsid,answ}' <(gzfread $5)
}
export -f extract_region

HEALTHSPAN="../q01_sumstats/st01_01_healthspan.tsv.gz"
LONGEVITY="../q01_sumstats/st01_01_longevity_90_pct.tsv.gz"
LIFESPAN="../q01_sumstats/st01_01_lifespan.tsv.gz"
TRIPLEMANOVA=../q02_manova/st02_01_triple_manova.tsv.gz

echo -e "gene\trsid_lead\trsid\tsnpid\tchr\tpos\ta1\ta0\tn\tfreq1\tbeta1\tse\tp" | gzip -9 > st03_01_healthspan_region_res.tsv.gz
tail -n+2 st03_03_triple_manova_res.tsv | parallel --no-notice -j 22 --col-sep="\t" extract_region {5} {1} {2} {3} $HEALTHSPAN | sort -k5g,5 -k6g,6 | gzip -9 >> st03_01_healthspan_region_res.tsv.gz

echo -e "gene\trsid_lead\trsid\tsnpid\tchr\tpos\ta1\ta0\tn\tfreq1\tbeta1\tse\tp" | gzip -9 > st03_01_lifespan_region_res.tsv.gz
tail -n+2 st03_03_triple_manova_res.tsv | parallel --no-notice -j 22 --col-sep="\t" extract_region {5} {1} {2} {3} $LIFESPAN | sort -k5g,5 -k6g,6 | gzip -9 >> st03_01_lifespan_region_res.tsv.gz

echo -e "gene\trsid_lead\trsid\tsnpid\tchr\tpos\ta1\ta0\tn\tfreq1\tbeta1\tse\tp" | gzip -9 > st03_01_longevity_region_res.tsv.gz
tail -n+2 st03_03_triple_manova_res.tsv | parallel --no-notice -j 22 --col-sep="\t" extract_region {5} {1} {2} {3} $LONGEVITY | sort -k5g,5 -k6g,6 | gzip -9 >> st03_01_longevity_region_res.tsv.gz

echo -e "gene\trsid_lead\trsid\tsnpid\tchr\tpos\ta1\ta0\tn\tfreq1\tbeta1\tse\tp" | gzip -9 > st03_01_triple_manova_region_res.tsv.gz
tail -n+2 st03_03_triple_manova_res.tsv | parallel --no-notice -j 22 --col-sep="\t" extract_region {5} {1} {2} {3} $TRIPLEMANOVA | sort -k5g,5 -k6g,6 | gzip -9 >> st03_01_triple_manova_region_res.tsv.gz



#----
# Create table
#----

Rscript ../scripts/a003_06_create_table.R |& tee log_a003_06_create_table.log
Rscript ../scripts/a003_06_create_region_table.R |& tee log_a003_06_create_table.log



# Extract loci that are nominally significant in all datasets but genome-wide significant in none

awk -v FS="," 'NR == FNR && $8 > 5e-8 && $11 > 5e-8 && $14 > 5e-8 {gene[$1]++; next}
               FNR ==1 || gene[$1] > 0 && $10 < 0.05 && $13 < 0.05 && $16 < 0.05 {print}' \
st03_06_triple_manova_region_table.csv  st03_06_triple_manova_table.csv | sed 's/,/\t/g' > st03_08_novel_loci.tsv



# Make upset plot

Rscript ../scripts/a003_08_upset_plot.R &> log_upset_plot.log &