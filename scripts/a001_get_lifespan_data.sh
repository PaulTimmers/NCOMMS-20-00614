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

wget --no-verbose "https://datashare.is.ed.ac.uk/bitstream/handle/10283/3209/lifegen_phase2_bothpl_alldr_2017_09_18.tsv.gz?sequence=1&isAllowed=y" -O raw_dl_lifespan.tsv.gz

# Print descriptives

echo -e "\nFile 'raw_dl_lifespan.tsv.gz' downloaded from 'https://datashare.is.ed.ac.uk/bitstream/handle/10283/3209/lifegen_phase2_bothpl_alldr_2017_09_18.tsv.gz?sequence=1&isAllowed=y'" 

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
<(zcat raw_dl_lifespan.tsv.gz)



# Reformat data

echo -n "Reformatting to wilson format... "

zcat raw_dl_lifespan.tsv.gz | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$13}' > st01_01_lifespan.tsv && rm raw_dl_lifespan.tsv.gz
echo "done."

echo -e "\nHeader:\n"
head st01_01_lifespan.tsv
echo -e "\n\n"


# Remove duplicates

echo -n "Removing all SNPs with non-unique rsids..."

awk 'NR==1 {next} NR == FNR {rsid[$1]++; next} rsid[$1] < 2 {print $0 > ".tmp.lifespan.no.duplicates"}
rsid[$1] > 1 {i++} END {printf " (%s: %'"'"'d) ","N", i}' st01_01_lifespan.tsv st01_01_lifespan.tsv \
&& mv .tmp.lifespan.no.duplicates st01_01_lifespan.tsv

echo "done."


echo -n "Removing all SNPs with non-unique snpids..."

awk 'NR==1 {next} NR == FNR {snpid[$2]++; next} snpid[$2] < 2 {print $0 > ".tmp.lifespan.no.duplicates"}
snpid[$1] > 1 {i++} END {printf " (%s: %'"'"'d) ","N", i}' st01_01_lifespan.tsv st01_01_lifespan.tsv \
&& mv .tmp.lifespan.no.duplicates st01_01_lifespan.tsv

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
st01_01_lifespan.tsv


# Retrieve correct N

awk -v OFS="\t" 'NR == FNR {N[$1]=$7; next} N[$1] > 0 {$7=N[$1]; print $0 > ".tmp.lifespan.fixed.n"}' $LIFEGEN2 st01_01_lifespan.tsv \
&& mv .tmp.lifespan.fixed.n st01_01_lifespan.tsv


# Gzip results

echo -n "Gzipping results... "
gzip -9 -f st01_01_lifespan.tsv
echo -e "done.\n\n"

date
echo -e "--- DONE ---\n\n"
