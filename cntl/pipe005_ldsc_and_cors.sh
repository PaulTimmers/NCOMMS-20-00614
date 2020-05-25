#!/bin/bash

#----
# Setup environment
#----

MAIN=${ROOT}/projects/prj_108_lifegen_phase3/p06_triple_manova
cd $MAIN

mkdir -p ${MAIN}/q05_ldsc
cd ${MAIN}/q05_ldsc


# Install virtual environment

if [[ ! -d venv ]]; then
	virtualenv-2.7 venv
	source venv/bin/activate
	pip install --upgrade pip
	pip install pandas
	pip install scipy
	pip install bitarray
fi


weights=${ROOT}/apps_by_us_full_stuff/ldsc/data/eur_w_ld_chr/
reference=${ROOT}/apps_by_us_full_stuff/ldsc/data/eur_w_ld_chr/
snp_list=${ROOT}/apps_by_us_full_stuff/ldsc/data/w_hm3.snplist

TRIPLEMANOVA=../q02_manova/st02_01_triple_manova.tsv.gz


#~~~
#~~~~~~
# Step 1: Lifespan trait correlations
#~~~~~~
#~~~


#----
# START
#----


source venv/bin/activate


munge_sumstats.py \
--sumstats $TRIPLEMANOVA \
--out st005_00_triple_manova \
--signed-sumstats beta1,0 \
--snp rsid \
--a1 a1 \
--a2 a0 \
--ignore snpid


ldsc.py \
--h2 st005_00_triple_manova.sumstats.gz \
--ref-ld-chr ${reference} \
--w-ld-chr ${weights} \
--out st005_01_triple_manova

rm st005_00_triple_manova.sumstats.gz


#----
# Calculate genetic correlations
#----

# Basic three-way correlation

echo "../q01_sumstats/st01_01_healthspan.tsv.gz
../q01_sumstats/st01_01_lifespan.tsv.gz
../q01_sumstats/st01_01_longevity_90_pct.tsv.gz" > st005_03_sumstat_files.txt


for SUMSTATS in `cat st005_03_sumstat_files.txt`
do
(outname=`basename ${SUMSTATS} | sed 's/.tsv.*//; s/st01_01_/st005_04_/'`

# Merge alleles with SNP list

zcat $SUMSTATS | head -1 | cut -f1-11,14 > ${outname}.tmp
awk -v OFS="\t" \
'NR == FNR {a1[$1]=$2; a2[$1]=$3; next} 
a1[$1] > 0 && ((a1[$1] == $5 && a2[$1] == $6) || (a1[$1] == $6 && a2[$1] == $5)) {print $1,$2,$3,$4,a1[$1],a2[$1],$7,(a1[$1]==$5 ? $8 : 1-$8),(a1[$1]==$5 ? $9 : -$9),$10,$11,$14}' $snp_list <(zcat $SUMSTATS) >> ${outname}.tmp


# Remove MHC
mhc=`wc -l < ${outname}.tmp`
awk 'NR == 1 || !($3 == 6 && $4 >= 28477797 && $4 <= 33448354)' ${outname}.tmp > ${outname}.tmp2
echo ${outname}.tmp: $(echo `wc -l < ${outname}.tmp`-`wc -l < ${outname}.tmp2` | bc -l) MHC SNPs removed
mv ${outname}.tmp2 ${outname}.tmp

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

n=`wc -l < st005_03_sumstat_files.txt`
for i in `eval echo {1..$n}`
do
ldsc.py \
--rg `awk -v ORS="," -v i=$i '{gsub("^.*/","",$0); gsub("tsv","sumstats",$0); gsub("st01_01_","st005_04_",$0)} NR == FNR && FNR == i {print; next} NR != FNR && FNR != i {print}' st005_03_sumstat_files.txt st005_03_sumstat_files.txt | sed 's/,$//'` \
--ref-ld-chr $reference \
--w-ld-chr $weights \
--out st005_05_genetic_correlations.${i}
done

echo -e "p1\tp2\trg\tse\tz\tp\th2_obs\th2_obs_se\th2_int\th2_int_se\tgcov_int\tgcov_int_se" > st005_05_genetic_correlations.tsv
grep -h p1 -A `wc -l < st005_03_sumstat_files.txt` --no-group-separator st005_05_genetic_correlations.*.log | sed 's/^[ ]*//; /^$/d; /p1/d' | tr -s ' ' '\t' >> st005_05_genetic_correlations.tsv
cat st005_05_genetic_correlations.*.log > st005_05_genetic_correlations.log && rm st005_05_genetic_correlations.*.log && rm st005_04_*.sumstats.gz




# Age-stratified 5-way correlation

echo "../q01_sumstats/st01_01_healthspan.tsv.gz
../q09_age_stratified_gwas/q01_combined/p05_comb_chr/bothpl_lif4060.tsv.gz
../q09_age_stratified_gwas/q01_combined/p05_comb_chr/bothpl_lif6080.tsv.gz
../q09_age_stratified_gwas/q01_combined/p05_comb_chr/bothpl_lif80120.tsv.gz
../q01_sumstats/st01_01_longevity_90_pct.tsv.gz" > st005_06_age_sumstat_files.txt


for SUMSTATS in `cat st005_06_age_sumstat_files.txt`
do
(outname=st005_07_`basename ${SUMSTATS} | sed 's/.tsv.*//; s/st01_01_//'`

# Merge alleles with SNP list
zcat $SUMSTATS | head -1 | cut -f1-11,14 > ${outname}.tmp
awk -v OFS="\t" \
'NR == FNR {a1[$1]=$2; a2[$1]=$3; next} 
a1[$1] > 0 && ((a1[$1] == $5 && a2[$1] == $6) || (a1[$1] == $6 && a2[$1] == $5)) {print $1,$2,$3,$4,a1[$1],a2[$1],$7,(a1[$1]==$5 ? $8 : 1-$8),(a1[$1]==$5 ? $9 : -$9),$10,$11,$14}' $snp_list <(zcat $SUMSTATS) >> ${outname}.tmp

# Remove MHC
mhc=`wc -l < ${outname}.tmp`
awk 'NR == 1 || !($3 == 6 && $4 >= 28477797 && $4 <= 33448354)' ${outname}.tmp > ${outname}.tmp2
echo ${outname}.tmp: $(echo `wc -l < ${outname}.tmp`-`wc -l < ${outname}.tmp2` | bc -l) MHC SNPs removed
mv ${outname}.tmp2 ${outname}.tmp

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

n=`wc -l < st005_06_age_sumstat_files.txt`
for i in `eval echo {1..$n}`
do
ldsc.py \
--rg `awk -v ORS="," -v i=$i '{gsub("^.*/","",$0); gsub("tsv","sumstats",$0); gsub("st01_01_","")} NR == FNR && FNR == i {print "st005_07_"$0; next} NR != FNR && FNR != i {print "st005_07_"$0}' st005_06_age_sumstat_files.txt st005_06_age_sumstat_files.txt | sed 's/,$//'` \
--ref-ld-chr $reference \
--w-ld-chr $weights \
--out st005_08_age_genetic_correlations.${i}
done

echo -e "p1\tp2\trg\tse\tz\tp\th2_obs\th2_obs_se\th2_int\th2_int_se\tgcov_int\tgcov_int_se" > st005_08_age_genetic_correlations.tsv
grep -h p1 -A `wc -l < st005_06_age_sumstat_files.txt` --no-group-separator st005_08_age_genetic_correlations.*.log | sed 's/^[ ]*//; /^$/d; /p1/d' | tr -s ' ' '\t' >> st005_08_age_genetic_correlations.tsv
cat st005_08_age_genetic_correlations.*.log > st005_08_age_genetic_correlations.log && rm st005_08_age_genetic_correlations.*.log && rm st005_07_*.sumstats.gz


#----
# Create plots
#----


Rscript ../scripts/a005_09_plot_correlations.R




#~~~
#~~~~~~
# Step 2: Setup environment for disease correlations
#~~~~~~
#~~~

MAIN=${ROOT}/projects/prj_108_lifegen_phase3/p06_triple_manova
cd $MAIN

mkdir -p ${MAIN}/q05_ldsc
cd ${MAIN}/q05_ldsc


#----
# Format for LD Hub
#----

source venv/bin/activate

mkdir -p ${MAIN}/q05_ldsc/ldhub
cd ${MAIN}/q05_ldsc/ldhub

ln -s ../../q01_sumstats/st01_01_healthspan.tsv.gz t001_healthspan.tsv.gz
ln -s ../../q01_sumstats/st01_01_lifegen.tsv.gz t001_lifespan.tsv.gz
ln -s ../../q01_sumstats/st01_01_longevity_90_pct.tsv.gz t001_longevity.tsv.gz
ln -s ../../q02_manova/st02_01_triple_manova.tsv.gz t001_manova.tsv.gz

awk 'NR == FNR {snpid[$2]++; next} snpid[$2] > 0' <(zcat t001_manova.tsv.gz) <(zcat t001_healthspan.tsv.gz) > t002_healthspan.tsv &
awk 'NR == FNR {snpid[$2]++; next} snpid[$2] > 0' <(zcat t001_manova.tsv.gz) <(zcat t001_lifespan.tsv.gz) > t002_lifespan.tsv &
awk 'NR == FNR {snpid[$2]++; next} snpid[$2] > 0' <(zcat t001_manova.tsv.gz) <(zcat t001_longevity.tsv.gz) > t002_longevity.tsv &


R --quiet --no-save <<"Rcode"
options(width=200)
library(data.table)

health <- fread("t002_healthspan.tsv", select=c("snpid","a1","a0","beta1","se"), key="snpid")
health[,c("z","beta1","se"):=.(beta1/se, NULL, NULL)]


life <- fread("t002_lifespan.tsv", select=c("snpid","a1","a0","beta1","se"), key="snpid")
life[,c("z","beta1","se"):=.(beta1/se, NULL, NULL)]

long <- fread("t002_longevity.tsv", select=c("snpid","a1","a0","beta1","se"), key="snpid")
long[,c("z","beta1","se"):=.(beta1/se, NULL, NULL)]

manova <- fread("t001_manova.tsv.gz", select=c("rsid","snpid","chr","pos","a1","a0","n","freq1","beta1","se","p"), key="snpid")
manova[,c("z","beta1","se"):=.(abs(qnorm(p/2)),NULL,NULL)]


m1 <- merge(manova,health,on="snpid", suffix=c("","_health"))
m1[a1_health == a0, c("a1_health","a0_health","z_health"):=.(a0_health, a1_health, -z_health)]
m1[,c("a1_health","a0_health"):=NULL]

m2 <- merge(m1,life,on="snpid", suffix=c("","_life"))
m2[a1_life == a0, c("a1_life","a0_life","z_life"):=.(a0_life, a0_life, -z_life)]
m2[,c("a1_life","a0_life"):=NULL]

m3 <- merge(m2,long,on="snpid", suffix=c("","_long"))
m3[a1_long == a0, c("a1_long","a0_long","z_long"):=.(a0_long, a0_long, -z_long)]
m3[,c("a1_long","a0_long"):=NULL]

export <- m3[,.(rsid,snpid,chr,pos,a1,a0,n,freq1,z=z * sign(z_health + z_life + z_long),p)]

fwrite(export, "t002_manova.tsv", sep="\t", quote=FALSE, na="NA")
Rcode

rm t002_healthspan.tsv t002_lifespan.tsv t002_longevity.tsv

for trait in healthspan lifespan longevity
do
munge_sumstats.py \
--sumstats t001_${trait}.tsv.gz \
--out t003_${trait} \
--signed-sumstats beta1,0 \
--snp rsid \
--a1 a1 \
--a2 a0 \
--ignore snpid &> t003_${trait}.log &
done


munge_sumstats.py \
--sumstats t002_manova.tsv \
--out t003_manova \
--signed-sumstats z,0 \
--snp rsid \
--a1 a1 \
--a2 a0 \
--ignore snpid &> t003_manova.log &




for trait in healthspan lifespan longevity manova
do
(awk -v OFS="\t" 'ARGIND == 1 {a1[$1]=$2; a2[$1]=$3; next} a1[$1] > 0 && (a1[$1] == $2 || a2[$1] == $2) {print $1,a1[$1],a2[$1],$4,(a1[$1]==$2 ? $5 : -$5)}' $snp_list <(zcat t003_${trait}.sumstats.gz) > t004_${trait}.sumstats && rm t003_${trait}.sumstats.gz) &
done

for trait in healthspan lifespan longevity manova
do
(awk -v OFS="\t" \
'NR == 1 {print "SNP","A1","A2","N","Z","P"; next} 
ARGIND == 1 {p[$1]=$11; next} p[$1] > 0 {print $1,$2,$3,$4,$5,p[$1]}' <(zcat t001_${trait}.tsv.gz) t004_${trait}.sumstats > t005_${trait}.sumstats
zip -m t005_${trait}.sumstats.zip t005_${trait}.sumstats && rm t004_${trait}.sumstats) &
done


#----
# "LD Hub"
#----

#----
# Format healthspan/lifespan/longevity
#----

zcat t005_healthspan.sumstats.zip | gzip -9 -f > t005_healthspan.sumstats.gz
zcat t005_lifespan.sumstats.zip | gzip -9 -f > t005_lifespan.sumstats.gz
zcat t005_longevity.sumstats.zip | gzip -9 -f > t005_longevity.sumstats.gz


# MHC snps

wget "http://static.geneatlas.roslin.ed.ac.uk/gwas/allWhites/snps/extended/snps.imputed.chr6.csv.gz" -O chr6_snps.csv.gz
awk '$2 >= 28477797 && $2 <= 33448354 {print $1}' <(zcat chr6_snps.csv.gz) > mhc_snps.tsv && rm chr6_snps.csv.gz

zcat t001_healthspan.tsv.gz | awk '$3 == 6 && $4 >= 28477797 && $4 <= 33448354 {print $1}' >> mhc_snps.tsv
sort -u mhc_snps.tsv > tmp && mv tmp mhc_snps.tsv


#~~~
#~~~~~~
# Step 3: Format GeneATLAS stats
#~~~~~~
#~~~


echo -e "GENEATLAS TRAIT" > t000_geneatlas_traits.txt
echo -e "cancer_c_Block_C15-C26 cancer_digestive" >> t000_geneatlas_traits.txt
echo -e "cancer_c_Block_C30-C39 cancer_respiratory" >> t000_geneatlas_traits.txt
echo -e "cancer_c_Block_C43-C44 cancer_melanoma" >> t000_geneatlas_traits.txt
echo -e "cancer_c_Block_C50-C50 cancer_breast" >> t000_geneatlas_traits.txt
echo -e "cancer_c_Block_C60-C63 cancer_genital_male" >> t000_geneatlas_traits.txt
echo -e "cancer_c_Block_C64-C68 cancer_urinary" >> t000_geneatlas_traits.txt
echo -e "cancer_c_Block_C81-C96 cancer_lymphoid" >> t000_geneatlas_traits.txt
echo -e "clinical_c_Block_K40-K46 hernia" >> t000_geneatlas_traits.txt
echo -e "clinical_c_Block_J40-J47 copd" >> t000_geneatlas_traits.txt


#----
# Run GeneATLAS formatting
#----

n_traits=`tail -n+2 t000_geneatlas_traits.txt | wc -l`
qsub -t 1-$n_traits ../../scripts/a005_ldhub_geneatlas.sh
qsub ../../scripts/a005_ldhub_geneatlas_merge.sh




#~~~
#~~~~~~
# Step 4: Format other sumstats
#~~~~~~
#~~~


source ../venv/bin/activate

weights=${ROOT}/apps_by_us_full_stuff/ldsc/data/eur_w_ld_chr/
reference=${ROOT}/apps_by_us_full_stuff/ldsc/data/eur_w_ld_chr/
snp_list=${ROOT}/apps_by_us_full_stuff/ldsc/data/w_hm3.snplist
 


#----
# Colorectal cancer
#----

wget "ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/ZhouW_30104761_GCST008372/PheCode_153_SAIGE_MACge20.txt.vcf.gz" -O t001_cancer_colorectal.tsv.gz

awk -v OFS="\t" 'ARGIND == 1 {a1[$1]=$2; a2[$1]=$3; next} 
FNR==1 {print "SNP","CHR","POS","A1","A2","N","Z","P"; next} 
a1[$3] > 0 && (a1[$3] == $5 || a2[$3] == $5) {print $3,$1,$2,$5,(a1[$3]==$5 ? a2[$3] : a1[$3]),2/(1/$8+1/$9),$10/$11,$13}' <(zcat t005_healthspan.sumstats.zip) <(zcat t001_cancer_colorectal.tsv.gz) > t002_cancer_colorectal.tsv

munge_sumstats.py \
--sumstats t002_cancer_colorectal.tsv \
--out t003_cancer_colorectal

awk -v OFS="\t" 'NR == FNR {(NR == 1 ? p[$1]="P": p[$1]=$NF); next} p[$1] > 0 {print $0,p[$1]}' t002_cancer_colorectal.tsv <(zcat t003_cancer_colorectal.sumstats.gz) > t004_cancer_colorectal.sumstats
grep -wvFf mhc_snps.tsv t004_cancer_colorectal.sumstats | gzip -9 -f > t005_cancer_colorectal.sumstats.gz && rm t004_cancer_colorectal.sumstats t003_cancer_colorectal.sumstats.gz t002_cancer_colorectal.tsv


#----
# Male Pattern Baldness
#----


mv UKBB_imputed_discovery_FULL.txt.gz t001_baldness.tsv.gz

awk -v OFS="\t" 'ARGIND == 1 {a1[$1]=$2; a2[$1]=$3; next} 
FNR==1 {print "SNP","CHR","POS","A1","A2","Z","P"; next} 
a1[$1] > 0 && (a1[$1] == $5 || a2[$1] == $5) {print $1,$3,$4,$5,(a1[$1]==$5 ? a2[$1] : a1[$1]),$7/$8,$9}' <(zcat t005_healthspan.sumstats.zip) <(zcat t001_baldness.tsv.gz | tr -s ' ' '\t') > t002_baldness.tsv

munge_sumstats.py \
--sumstats t002_baldness.tsv \
--out t003_baldness \
--N `echo "2 / (1/25662 + 1/17928)" | bc -l`

awk -v OFS="\t" 'NR == FNR {(NR == 1 ? p[$1]="P": p[$1]=$7); next} p[$1] > 0 {print $0,p[$1]}' t002_baldness.tsv <(zcat t003_baldness.sumstats.gz) > t004_baldness.sumstats
grep -wvFf mhc_snps.tsv t004_baldness.sumstats | gzip -9 -f > t005_baldness.sumstats.gz && rm t004_baldness.sumstats t003_baldness.sumstats.gz t002_baldness.tsv


#----
# Depression
#----

wget "https://datashare.is.ed.ac.uk/bitstream/handle/10283/3203/PGC_UKB_depression_genome-wide.txt?sequence=3&isAllowed=y" -O t001_depression.tsv

awk -v OFS="\t" 'ARGIND == 1 {a1[$1]=$2; a2[$1]=$3; next} 
FNR==1 {print "SNP","A1","A2","Z","P"; next} 
a1[$1] > 0 && (a1[$1] == toupper($2) || a2[$1] == toupper($2)) {print $1,toupper($2),(a1[$1]==toupper($2) ? a2[$1] : a1[$1]),$5/$6,$7}' <(zcat t005_healthspan.sumstats.zip) <(cat t001_depression.tsv | tr -s ' ' '\t') > t002_depression.tsv

munge_sumstats.py \
--sumstats t002_depression.tsv \
--out t003_depression \
--N `echo "2 / (1/(127552 + 43204) + 1/(233763 + 95680))" | bc -l`

awk -v OFS="\t" 'NR == FNR {(NR == 1 ? p[$1]="P": p[$1]=$NF); next} p[$1] > 0 {print $0,p[$1]}' t002_depression.tsv <(zcat t003_depression.sumstats.gz) > t004_depression.sumstats
grep -wvFf mhc_snps.tsv t004_depression.sumstats | gzip -9 -f > t005_depression.sumstats.gz && rm t004_depression.sumstats t003_depression.sumstats.gz t002_depression.tsv


#----
# Alzheimers
#----

# ncbi.nlm.nih.gov/pmc/articles/PMC5959890/
wget "https://datashare.is.ed.ac.uk/bitstream/handle/10283/3364/4_UKB_IGAP_AD_meta_summary_output_June2019.txt?sequence=5&isAllowed=y" -O t001_alzheimer.tsv

awk -v OFS="\t" 'ARGIND == 1 {a1[$1]=$2; a2[$1]=$3; next} 
FNR==1 {print "SNP","CHR","POS","A1","A2","Z","P"; next} 
a1[$3] > 0 && (a1[$3] == toupper($4) || a2[$3] == toupper($4)) {print $3,$1,$2,toupper($4),(a1[$3]==toupper($4) ? a2[$3] : a1[$3]),$6/$7,$8}' <(zcat t005_healthspan.sumstats.zip) <(cat t001_alzheimer.tsv | tr -s ' ' '\t') > t002_alzheimer.tsv

munge_sumstats.py \
--sumstats t002_alzheimer.tsv \
--out t003_alzheimer \
--N `echo "2 / (1/(25580+(27696+14338)/4) + 1/(48466 + (260980+245941)/4))" | bc -l`

awk -v OFS="\t" 'NR == FNR {(NR == 1 ? p[$1]="P": p[$1]=$NF); next} p[$1] > 0 {print $0,p[$1]}' t002_alzheimer.tsv <(zcat t003_alzheimer.sumstats.gz) > t004_alzheimer.sumstats
grep -wvFf mhc_snps.tsv t004_alzheimer.sumstats | gzip -9 -f > t005_alzheimer.sumstats.gz && rm t004_alzheimer.sumstats t003_alzheimer.sumstats.gz t002_alzheimer.tsv


#----
# Risk Taking
#----

wget "https://rp.mrc-epid.cam.ac.uk/files/Risk_SumStats_Clifton_2018.csv.gz" -O t001_risk_taking.csv.gz


awk -v OFS="\t" 'ARGIND == 1 {a1[$1]=$2; a2[$1]=$3; next} 
FNR==1 {print "SNP","A1","A2","Z","P"; next} 
a1[$1] > 0 && (a1[$1] == toupper($2) || a2[$1] == toupper($2)) {print $1,toupper($2),(a1[$1]==toupper($2) ? a2[$1] : a1[$1]),$5/$6,$7}' <(zcat t005_healthspan.sumstats.zip) <(zcat t001_risk_taking.csv.gz | tr -s ',' '\t') > t002_risk_taking.tsv

munge_sumstats.py \
--sumstats t002_risk_taking.tsv \
--out t003_risk_taking \
--N `echo "2 / (1/129877 + 1/352296)" | bc -l`


awk -v OFS="\t" 'NR == FNR {(NR == 1 ? p[$1]="P": p[$1]=$NF); next} p[$1] > 0 {print $0,p[$1]}' t002_risk_taking.tsv <(zcat t003_risk_taking.sumstats.gz) > t004_risk_taking.sumstats
grep -wvFf mhc_snps.tsv t004_risk_taking.sumstats | gzip -9 -f > t005_risk_taking.sumstats.gz && rm t004_risk_taking.sumstats t003_risk_taking.sumstats.gz t002_risk_taking.tsv


#----
# Pubertal Growth
#----

# http://egg-consortium.org/pubertal-growth.html
wget "http://egg-consortium.org/downloads/Pubertal_growth_10F_12M_combined.txt.gz" -O t001_pubertal_growth.tsv.gz

awk -v OFS="\t" 'ARGIND == 1 {a1[$1]=$2; a2[$1]=$3; next} 
FNR==1 {print "SNP","CHR","POS","A1","A2","N","Z","P"; next} 
a1[$1] > 0 && (a1[$1] == toupper($4) || a2[$1] == toupper($4)) {print $1,$2,$3,toupper($4),(a1[$1]==toupper($4) ? a2[$1] : a1[$1]),$9,$6/$7,$8}' <(zcat t005_healthspan.sumstats.zip) <(zcat t001_pubertal_growth.tsv.gz) > t002_pubertal_growth.tsv

munge_sumstats.py \
--sumstats t002_pubertal_growth.tsv \
--out t003_pubertal_growth 

awk -v OFS="\t" 'NR == FNR {(NR == 1 ? p[$1]="P": p[$1]=$NF); next} p[$1] > 0 {print $0,p[$1]}' t002_pubertal_growth.tsv <(zcat t003_pubertal_growth.sumstats.gz) > t004_pubertal_growth.sumstats
grep -wvFf mhc_snps.tsv t004_pubertal_growth.sumstats | gzip -9 -f > t005_pubertal_growth.sumstats.gz && rm t004_pubertal_growth.sumstats t003_pubertal_growth.sumstats.gz t002_pubertal_growth.tsv


#----
# Smoking
#----

# UKBB_Ben_Neale Ever Smoked
wget "https://www.dropbox.com/s/t0svj30x1hugqh4/20160.assoc.tsv.gz?dl=0%20-O%2020160.assoc.tsv.gz" -O t001_smoking.tsv.gz

awk -v OFS="\t" 'ARGIND == 1 {a1[$1]=$2; a2[$1]=$3; next} 
FNR==1 {print "SNP","A1","A2","N","Z","P"; next} 
{gsub(".*:","",$1)}
a1[$2] > 0 && (a1[$2] == toupper($1) || a2[$2] == toupper($1)) {print $2,toupper($1),(a1[$2]==toupper($1) ? a2[$2] : a1[$2]),$3,$8,$9}' <(zcat t005_healthspan.sumstats.zip) <(zcat t001_smoking.tsv.gz) > t002_smoking.tsv

munge_sumstats.py \
--sumstats t002_smoking.tsv \
--out t003_smoking

awk -v OFS="\t" 'NR == FNR {(NR == 1 ? p[$1]="P": p[$1]=$NF); next} p[$1] > 0 {print $0,p[$1]}' t002_smoking.tsv <(zcat t003_smoking.sumstats.gz) > t004_smoking.sumstats
grep -wvFf mhc_snps.tsv t004_smoking.sumstats | gzip -9 -f > t005_smoking.sumstats.gz && rm t004_smoking.sumstats t003_smoking.sumstats.gz t002_smoking.tsv


#----
# Drinking
#----

# https://www.ncbi.nlm.nih.gov/pubmed/31358974
wget "ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/EvangelouE_31358974_GCST008757/EvangelouE_31358974.txt" -O t001_drinking.tsv

# Add rsids
sed -i 's/:/_/' t001_drinking.tsv
awk 'NR == FNR {rsid[$2]=$1; next} FNR == 1 || rsid[$1] > 0 {print (FNR == 1 ? "rsid" : rsid[$1]),$0}' <(zcat t001_healthspan.tsv.gz) t001_drinking.tsv > tmp && mv tmp t001_drinking.tsv

awk -v OFS="\t" 'ARGIND == 1 {a1[$1]=$2; a2[$1]=$3; next} 
FNR==1 {print "SNP","A1","A2","N","Z","P"; next} 
a1[$1] > 0 && (a1[$1] == toupper($3) || a2[$1] == toupper($3)) {print $1,toupper($3),(a1[$1]==toupper($3) ? a2[$1] : a1[$1]),$11,$6/$7,$8}' <(zcat t005_healthspan.sumstats.zip) <(cat t001_drinking.tsv | tr -s ' ' '\t') > t002_drinking.tsv

munge_sumstats.py \
--sumstats t002_drinking.tsv \
--out t003_drinking

awk -v OFS="\t" 'NR == FNR {(NR == 1 ? p[$1]="P": p[$1]=$NF); next} p[$1] > 0 {print $0,p[$1]}' t002_drinking.tsv <(zcat t003_drinking.sumstats.gz) > t004_drinking.sumstats
grep -wvFf mhc_snps.tsv t004_drinking.sumstats | gzip -9 -f > t005_drinking.sumstats.gz && rm t004_drinking.sumstats t003_drinking.sumstats.gz t002_drinking.tsv


#----
# Coronary Artery Disease
#----

wget "http://www.cardiogramplusc4d.org/media/cardiogramplusc4d-consortium/data-downloads/cardiogram_gwas_results.zip" -O t001_coronary_artery.tsv.gz

awk -v OFS="\t" 'ARGIND == 1 {a1[$1]=$2; a2[$1]=$3; next} 
FNR==1 {print "SNP","A1","A2","N","Z","P"; next} 
a1[$1] > 0 && (a1[$1] == toupper($3) || a2[$1] == toupper($3)) {print $1,toupper($3),(a1[$1]==toupper($3) ? a2[$1] : a1[$1]),2/(1/$10 + 1/$11),$8/$9,$6}' <(zcat t005_healthspan.sumstats.zip) <(zcat t001_coronary_artery.tsv.gz) > t002_coronary_artery.tsv

munge_sumstats.py \
--sumstats t002_coronary_artery.tsv \
--out t003_coronary_artery

awk -v OFS="\t" 'NR == FNR {(NR == 1 ? p[$1]="P": p[$1]=$NF); next} p[$1] > 0 {print $0,p[$1]}' t002_coronary_artery.tsv <(zcat t003_coronary_artery.sumstats.gz) > t004_coronary_artery.sumstats
grep -wvFf mhc_snps.tsv t004_coronary_artery.sumstats | gzip -9 -f > t005_coronary_artery.sumstats.gz && rm t004_coronary_artery.sumstats t003_coronary_artery.sumstats.gz t002_coronary_artery.tsv


#----
# Educational Attainment
#----

wget "https://www.dropbox.com/s/ho58e9jmytmpaf8/GWAS_EA_excl23andMe.txt?dl=0" -O t001_years_of_schooling.1.tsv
wget "http://ssgac.org/documents/GWAS_EA.to10K.txt" -O t001_years_of_schooling.2.tsv

awk -v OFS="\t" 'NR == 1 {print $0,"N"; next} NR == FNR {snp[$1]=$0; next} snp[$1] > 0 {print snp[$1],1131881; next} FNR > 1 {print $0,766345} ' t001_years_of_schooling.2.tsv t001_years_of_schooling.1.tsv \
> t001_years_of_schooling.tsv && rm t001_years_of_schooling.{1..2}.tsv

awk -v OFS="\t" 'ARGIND == 1 {a1[$1]=$2; a2[$1]=$3; next} 
FNR==1 {print "SNP","CHR","POS","A1","A2","N","Z","P"; next} 
a1[$1] > 0 && (a1[$1] == toupper($4) || a2[$1] == toupper($4)) {print $1,$2,$3,toupper($4),(a1[$1]==toupper($4) ? a2[$1] : a1[$1]),$10,$7/$8,$9}' <(zcat t005_healthspan.sumstats.zip) <(cat t001_years_of_schooling.tsv) > t002_years_of_schooling.tsv

munge_sumstats.py \
--sumstats t002_years_of_schooling.tsv \
--out t003_years_of_schooling

awk -v OFS="\t" 'NR == FNR {(NR == 1 ? p[$1]="P": p[$1]=$NF); next} p[$1] > 0 {print $0,p[$1]}' t002_years_of_schooling.tsv <(zcat t003_years_of_schooling.sumstats.gz) > t004_years_of_schooling.sumstats
grep -wvFf mhc_snps.tsv t004_years_of_schooling.sumstats | gzip -9 -f > t005_years_of_schooling.sumstats.gz && rm t004_years_of_schooling.sumstats t003_years_of_schooling.sumstats.gz t002_years_of_schooling.tsv


#----
# Rheumatoid Arthritis
#----

wget "http://plaza.umin.ac.jp/~yokada/datasource/files/GWASMetaResults/RA_GWASmeta_European_v2.txt.gz" -O t001_rheumatoid_arthritis.tsv.gz

# Add rsids
awk -v OFS="\t" 'NR == FNR {rsid[$2]=$1; next} FNR == 1 || rsid[$2"_"$3] > 0 {print (FNR == 1 ? "rsid" : rsid[$2"_"$3]),$0}' <(zcat t001_healthspan.tsv.gz) <(zcat t001_rheumatoid_arthritis.tsv.gz) > t001_rheumatoid_arthritis.tsv

awk -v OFS="\t" 'ARGIND == 1 {a1[$1]=$2; a2[$1]=$3; next} 
FNR==1 {print "SNP","CHR","POS","A1","A2","Z","P"; next} 
a1[$1] > 0 && (a1[$1] == toupper($5) || a2[$1] == toupper($5)) {print $1,$3,$4,toupper($5),(a1[$1]==toupper($5) ? a2[$1] : a1[$1]), log($7) / ((log($7)-log($8)) / 1.959964),$10}' <(zcat t005_healthspan.sumstats.zip) <(cat t001_rheumatoid_arthritis.tsv) > t002_rheumatoid_arthritis.tsv

munge_sumstats.py \
--sumstats t002_rheumatoid_arthritis.tsv \
--out t003_rheumatoid_arthritis \
--N `echo "2 / (1/4361 + 1/43923)" | bc -l`

awk -v OFS="\t" 'NR == FNR {(NR == 1 ? p[$1]="P": p[$1]=$NF); next} p[$1] > 0 {print $0,p[$1]}' t002_rheumatoid_arthritis.tsv <(zcat t003_rheumatoid_arthritis.sumstats.gz) > t004_rheumatoid_arthritis.sumstats
grep -wvFf mhc_snps.tsv t004_rheumatoid_arthritis.sumstats | gzip -9 -f > t005_rheumatoid_arthritis.sumstats.gz && rm t004_rheumatoid_arthritis.sumstats t003_rheumatoid_arthritis.sumstats.gz t002_rheumatoid_arthritis.tsv


#----
# Inflammatory Bowel Disease
#----

wget "ftp://ftp.sanger.ac.uk/pub/consortia/ibdgenetics/iibdgc-trans-ancestry-filtered-summary-stats.tgz"
utar iibdgc-trans-ancestry-filtered-summary-stats.tgz && \
rm *trans_ethnic_association*
mv  EUR.IBD.gwas_info03_filtered.assoc t001_inflammatory_bowel.tsv
rm EUR.* README

awk -v OFS="\t" 'ARGIND == 1 {a1[$1]=$2; a2[$1]=$3; next} 
FNR==1 {print "SNP","CHR","POS","A1","A2","Z","INFO","P"; next} 
a1[$2] > 0 && (a1[$2] == toupper($4) || a2[$2] == toupper($4)) {print $2,$1,$3,toupper($4),(a1[$2]==toupper($4) ? a2[$2] : a1[$2]), log($9) / $10,$8,$11}' <(zcat t005_healthspan.sumstats.zip) <(cat t001_inflammatory_bowel.tsv) > t002_inflammatory_bowel.tsv
 
munge_sumstats.py \
--sumstats t002_inflammatory_bowel.tsv \
--out t003_inflammatory_bowel \
--N `echo "2 / (1/12882 + 1/21770)" | bc -l`

awk -v OFS="\t" 'NR == FNR {(NR == 1 ? p[$1]="P": p[$1]=$NF); next} p[$1] > 0 {print $0,p[$1]}' t002_inflammatory_bowel.tsv <(zcat t003_inflammatory_bowel.sumstats.gz) > t004_inflammatory_bowel.sumstats
grep -wvFf mhc_snps.tsv t004_inflammatory_bowel.sumstats | gzip -9 -f > t005_inflammatory_bowel.sumstats.gz && rm t004_inflammatory_bowel.sumstats t003_inflammatory_bowel.sumstats.gz t002_inflammatory_bowel.tsv


#----
# Allergic disease
#----

wget "ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/FerreiraMA_29083406_GCST005038/SHARE-without23andMe.LDSCORE-GC.SE-META.v0.gz" -O t001_allergic_disease.tsv.gz

awk -v OFS="\t" 'ARGIND == 1 {a1[$1]=$2; a2[$1]=$3; next} 
FNR==1 {print "SNP","CHR","POS","A1","A2","N","Z","P"; next} 
a1[$10] > 0 && (a1[$10] == toupper($4) || a2[$10] == toupper($4)) {print $10,$2,$3,toupper($4),(a1[$10]==toupper($4) ? a2[$10] : a1[$10]), $13,$6/$7,$8}' <(zcat t005_healthspan.sumstats.zip) <(zcat t001_allergic_disease.tsv.gz | tr -s ' ' '\t') > t002_allergic_disease.tsv

munge_sumstats.py \
--sumstats t002_allergic_disease.tsv \
--out t003_allergic_disease

awk -v OFS="\t" 'NR == FNR {(NR == 1 ? p[$1]="P": p[$1]=$NF); next} p[$1] > 0 {print $0,p[$1]}' t002_allergic_disease.tsv <(zcat t003_allergic_disease.sumstats.gz) > t004_allergic_disease.sumstats
grep -wvFf mhc_snps.tsv t004_allergic_disease.sumstats | gzip -9 -f > t005_allergic_disease.sumstats.gz && rm t004_allergic_disease.sumstats t003_allergic_disease.sumstats.gz t002_allergic_disease.tsv


#----
# Type 2 diabetes
#----

wget "http://cnsgenomics.com/data/t2d/Xue_et_al_T2D_META_Nat_Commun_2018.gz" -O t001_type_2_diabetes.tsv.gz

awk -v OFS="\t" 'ARGIND == 1 {a1[$1]=$2; a2[$1]=$3; next} 
FNR==1 {print "SNP","CHR","POS","A1","A2","N","Z","P"; next} 
a1[$3] > 0 && (a1[$3] == toupper($4) || a2[$3] == toupper($4)) {print $3,$1,$2,toupper($4),(a1[$3]==toupper($4) ? a2[$3] : a1[$3]), $10, $7/$8,$9}' <(zcat t005_healthspan.sumstats.zip) <(zcat t001_type_2_diabetes.tsv.gz | tr -s ' ' '\t') > t002_type_2_diabetes.tsv

munge_sumstats.py \
--sumstats t002_type_2_diabetes.tsv \
--out t003_type_2_diabetes

awk -v OFS="\t" 'NR == FNR {(NR == 1 ? p[$1]="P": p[$1]=$NF); next} p[$1] > 0 {print $0,p[$1]}' t002_type_2_diabetes.tsv <(zcat t003_type_2_diabetes.sumstats.gz) > t004_type_2_diabetes.sumstats
grep -wvFf mhc_snps.tsv t004_type_2_diabetes.sumstats | gzip -9 -f > t005_type_2_diabetes.sumstats.gz && rm t004_type_2_diabetes.sumstats t003_type_2_diabetes.sumstats.gz t002_type_2_diabetes.tsv


#----
# Age at menarche
#----

wget "https://www.reprogen.org/Menarche_1KG_NatGen2017_WebsiteUpload.zip" -O t001_age_at_menarche.zip
unzip t001_age_at_menarche.zip && rm README.txt t001_age_at_menarche.zip
mv Menarche_1KG_NatGen2017_WebsiteUpload.txt t001_age_at_menarche.tsv

# Add rsids
sed -i 's/chr//; s/:/_/' t001_age_at_menarche.tsv
awk -v OFS="\t" 'NR == FNR {rsid[$2]=$1; next} FNR == 1 || rsid[$1] > 0 {print (FNR == 1 ? "rsid" : rsid[$1]),$0}' <(zcat t001_healthspan.tsv.gz) t001_age_at_menarche.tsv > tmp && mv tmp t001_age_at_menarche.tsv

awk -v OFS="\t" 'ARGIND == 1 {a1[$1]=$2; a2[$1]=$3; next} 
FNR==1 {print "SNP","A1","A2","BETA","P"; next} 
a1[$1] > 0 && (a1[$1] == toupper($3) || a2[$1] == toupper($3)) {print $1,toupper($3),(a1[$1]==toupper($3) ? a2[$1] : a1[$1]),$5,$6}' <(zcat t005_healthspan.sumstats.zip) <(cat t001_age_at_menarche.tsv) > t002_age_at_menarche.tsv

# Calculate Z score
R --quiet --no-save <<"Rcode"
library(data.table); data <- fread("t002_age_at_menarche.tsv"); fwrite(data[,.(SNP,A1,A2,Z=sign(BETA) * -qnorm(P/2), P)], "t002_age_at_menarche.tsv", quote=F, sep="\t", na="NA")
Rcode

munge_sumstats.py \
--sumstats t002_age_at_menarche.tsv \
--out t003_age_at_menarche \
--N 329345

awk -v OFS="\t" 'NR == FNR {(NR == 1 ? p[$1]="P": p[$1]=$NF); next} p[$1] > 0 {print $0,p[$1]}' t002_age_at_menarche.tsv <(zcat t003_age_at_menarche.sumstats.gz) > t004_age_at_menarche.sumstats
grep -wvFf mhc_snps.tsv t004_age_at_menarche.sumstats | gzip -9 -f > t005_age_at_menarche.sumstats.gz && rm t004_age_at_menarche.sumstats t003_age_at_menarche.sumstats.gz t002_age_at_menarche.tsv


#----
# Age at menopause
#----

wget "https://www.reprogen.org/Menopause_HapMap2_DayNG2015_18112015.txt.gz" -O t001_age_at_menopause.tsv.gz

awk -v OFS="\t" 'ARGIND == 1 {a1[$1]=$2; a2[$1]=$3; next} 
FNR==1 {print "SNP","A1","A2","Z","P"; next} 
a1[$1] > 0 && (a1[$1] == toupper($2) || a2[$1] == toupper($2)) {print $1,toupper($2),(a1[$1]==toupper($2) ? a2[$1] : a1[$1]), $5/$6,$7}' <(zcat t005_healthspan.sumstats.zip) <(zcat t001_age_at_menopause.tsv.gz) > t002_age_at_menopause.tsv

munge_sumstats.py \
--sumstats t002_age_at_menopause.tsv \
--out t003_age_at_menopause \
--N 69360

awk -v OFS="\t" 'NR == FNR {(NR == 1 ? p[$1]="P": p[$1]=$NF); next} p[$1] > 0 {print $0,p[$1]}' t002_age_at_menopause.tsv <(zcat t003_age_at_menopause.sumstats.gz) > t004_age_at_menopause.sumstats
grep -wvFf mhc_snps.tsv t004_age_at_menopause.sumstats | gzip -9 -f > t005_age_at_menopause.sumstats.gz && rm t004_age_at_menopause.sumstats t003_age_at_menopause.sumstats.gz t002_age_at_menopause.tsv


#----
# Stroke
#----

wget "ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/MalikR_29531354_GCST006906/MEGASTROKE.1.AS.EUR.out" -O t001_stroke.tsv

awk -v OFS="\t" 'ARGIND == 1 {a1[$1]=$2; a2[$1]=$3; next} 
FNR==1 {print "SNP","A1","A2","Z","P"; next} 
a1[$1] > 0 && (a1[$1] == toupper($2) || a2[$1] == toupper($2)) {print $1,toupper($2),(a1[$1]==toupper($2) ? a2[$1] : a1[$1]), $5/$6,$7}' <(zcat t005_healthspan.sumstats.zip) <(cat t001_stroke.tsv) > t002_stroke.tsv

munge_sumstats.py \
--sumstats t002_stroke.tsv \
--out t003_stroke \
--N 446696

awk -v OFS="\t" 'NR == FNR {(NR == 1 ? p[$1]="P": p[$1]=$NF); next} p[$1] > 0 {print $0,p[$1]}' t002_stroke.tsv <(zcat t003_stroke.sumstats.gz) > t004_stroke.sumstats
grep -wvFf mhc_snps.tsv t004_stroke.sumstats | gzip -9 -f > t005_stroke.sumstats.gz && rm t004_stroke.sumstats t003_stroke.sumstats.gz t002_stroke.tsv


#----
# BMI
#----

wget "https://zenodo.org/record/1251813/files/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz?download=1" -O t001_bmi.tsv.gz



awk -v OFS="\t" 'ARGIND == 1 {a1[$1]=$2; a2[$1]=$3; next} 
FNR==1 {print "SNP",CHR,"A1","A2","N","INFO","Z","P"; next}
{gsub(":.*$","",$3)}
a1[$3] > 0 && (a1[$3] == toupper($4) || a2[$3] == toupper($4)) {print $3,toupper($4),(a1[$3]==toupper($4) ? a2[$3] : a1[$3]),$10,$11,$7/$8,$9}' <(zcat t005_healthspan.sumstats.zip) <(zcat t001_bmi.tsv.gz | tr -s ' ' '\t') > t002_bmi.tsv

munge_sumstats.py \
--sumstats t002_bmi.tsv \
--out t003_bmi

awk -v OFS="\t" 'NR == FNR {(NR == 1 ? p[$1]="P": p[$1]=$NF); next} p[$1] > 0 {print $0,p[$1]}' t002_bmi.tsv <(zcat t003_bmi.sumstats.gz) > t004_bmi.sumstats
grep -wvFf mhc_snps.tsv t004_bmi.sumstats | gzip -9 -f > t005_bmi.sumstats.gz && rm t004_bmi.sumstats t003_bmi.sumstats.gz t002_bmi.tsv


#~~~
#~~~~~~
# Step 5: Run genetic correlations
#~~~~~~
#~~~

rm t001_* t002_*.tsv

R --quiet --no-save <<"Rcode"
library(yaml); traits <- data.frame(pmid=unlist(yaml.load_file("t000_traits.yml")))
traits$category <- gsub("\\..*$","",rownames(traits))
rownames(traits) <- gsub("^.*\\.","",rownames(traits))
write.table(traits, "t006_traits.tsv", col.names=F, sep="\t", quote=F)
Rcode


n_traits=`wc -l < t006_traits.tsv`
qsub -t 1-$n_traits ../../scripts/a005_ldhub_run_rg.sh



#----
# Combine results
#----

head -1 t007_healthspan.tsv > t008_rg.tsv
cut -f1 t006_traits.tsv | while read trait
do
	tail -n+2 t007_${trait}.tsv >> t008_rg.tsv
done



#----
# Plot
#----

Rscript ../../scripts/a005_ldhub_plot_correlations.R  &> log_a005_ldhub_plot_correlations.log
