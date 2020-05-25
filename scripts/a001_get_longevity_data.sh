#!/bin/bash

#----
# Setup environment
#----

# Define variables

set -e
LIFEGEN2=${ROOT}/projects/prj_065b_lifegen_phase2/p007_gather_res/bothpl_alldr.tsv
REFERENCE="${ROOT}/data/reference_data/1000genomes/phase3/plink/1000genomes_chr"

#----
# START
#----

echo -e "--- START ---\n"`date`"\n\n"

wget --no-verbose "https://storage.googleapis.com/charge-longevity/Results_90th_percentile.txt.gz" -O raw_dl_longevity_90_pct.tsv.gz


# Print descriptives

echo -e "\nFile 'raw_dl_longevity_90_pct.tsv.gz' downloaded from 'https://storage.googleapis.com/charge-longevity/Results_90th_percentile.txt.gz'" 

awk -v FS="\t" \
'NR == 1 {print "Header:\n"; min_maf=1; max_maf=0} 
 NR == 1 {print; next}
 
 $6 > max_maf {max_maf=$6}
 $6 < min_maf {min_maf=$6}

 $6 < 0.05 || $6 > 0.95 {rare++}
 $6 >= 0.05 && $6 <= 0.95 {common++}

 NR <= 10 {print} 

 END{printf "\n%s: %'"'"'d\n%s: %'"'"'d\n%s: %'"'"'d\n%s: %.6g-%.6g\n\n",
            "Number of lines",NR,"Common SNPs (MAF >= 0.05)",common,"Rare SNPs (MAF < 0.05)",rare,"MAF range",min_maf,max_maf}' \
<(zcat raw_dl_longevity_90_pct.tsv.gz)


# Reformat data

echo -n "Reformatting to wilson format... "

awk -v OFS="\t" 'NR == 1 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11; next} 
                 NR == FNR {rsid[$2]=$1; next} 
                 FNR > 1 {print (rsid[$2"_"$3] > 0) ? rsid[$2"_"$3] : "NA",
                 $2"_"$3,
                 $2, $3,
                 toupper($4), toupper($5),
                 $10, $6,
                 $7, $8, $9}' $LIFEGEN2 <(zcat raw_dl_longevity_90_pct.tsv.gz) \
                 > st01_01_longevity_90_pct.tsv && rm raw_dl_longevity_90_pct.tsv.gz

echo "done."

echo -e "\nHeader:\n"
head st01_01_longevity_90_pct.tsv
echo -e "\n\n"


# Fill out missing rsids using 1000G data


echo -n "Filling out missing rsids... "

CHRS=`awk '$1 == "NA" {chr[$3]++} END{ for (i in chr) {print i} }' st01_01_longevity_90_pct.tsv`


for CHR in $CHRS
do
echo -n "$CHR "
awk -v CHR=$CHR -v OFS="\t" 'NR == FNR {rsid[$1"_"$4]=$2; next}
                 $1 == "NA" && $3==CHR {$1=(rsid[$2] > 0) ? rsid[$2] : $3":"$4}
                 {print}' ${REFERENCE}${CHR}.bim st01_01_longevity_90_pct.tsv \
                  > .tmp.$CHR && mv .tmp.$CHR st01_01_longevity_90_pct.tsv
done

echo "done."

echo -e "\nHeader:\n"
head st01_01_longevity_90_pct.tsv
echo -e "\n\n"



# Remove duplicates

echo -n "Removing all SNPs with non-unique rsids... "

awk 'NR==1 {next} NR == FNR {rsid[$1]++; next} rsid[$1] < 2 {print $0 > ".tmp.longevity.no.duplicates"}
rsid[$1] > 1 {i++} END {printf " (%s: %'"'"'d) ","N", i}' st01_01_longevity_90_pct.tsv st01_01_longevity_90_pct.tsv \
&& mv .tmp.longevity.no.duplicates st01_01_longevity_90_pct.tsv

echo "done."


echo -n "Removing all SNPs with non-unique snpids..."

awk 'NR==1 {next} NR == FNR {snpid[$2]++; next} snpid[$2] < 2 {print $0 > ".tmp.longevity.no.duplicates"}
snpid[$1] > 1 {i++} END {printf " (%s: %'"'"'d) ","N", i}' st01_01_longevity_90_pct.tsv st01_01_longevity_90_pct.tsv \
&& mv .tmp.longevity.no.duplicates st01_01_longevity_90_pct.tsv

echo "done."


# Print more descriptives


awk -v FS="\t" \
'NR == 1 {print "Header:\n"; min_maf=1; max_maf=0} 
 NR == 1 {print; next}
 
 $8 > max_maf {max_maf=$8}
 $8 < min_maf {min_maf=$8}

 $8 < 0.05 || $8 > 0.95 {rare++}
 $8 >= 0.05 && $8 <= 0.95 {common++}

 NR <= 10 {print} 

 END{printf "\n%s: %'"'"'d\n%s: %'"'"'d\n%s: %'"'"'d\n%s: %.6g-%.6g\n\n",
            "Number of lines",NR,"Common SNPs (MAF >= 0.05)",common,"Rare SNPs (MAF < 0.05)",rare,"MAF range",min_maf,max_maf}' \
st01_01_longevity_90_pct.tsv


# Gzip results

echo -n "Gzipping results... "
gzip -9 -f st01_01_longevity_90_pct.tsv
echo -e "done.\n\n"

date
echo -e "--- DONE ---\n\n"
