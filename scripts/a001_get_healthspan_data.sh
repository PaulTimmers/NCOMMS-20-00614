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

wget --no-verbose "https://zenodo.org/record/1302861/files/healthspan_summary.csv.gz?download=1" -O raw_dl_healthspan.csv.gz

# Print descriptives

echo -e "\nFile 'raw_dl_healthspan.csv.gz' downloaded from 'https://zenodo.org/record/1302861/files/healthspan_summary.csv.gz?download=1'" 

awk -v FS="," \
'NR == 1 {print "Header:\n"; min_maf=1; max_maf=0} 
 NR == 1 {print; next}
 
 $6 > max_maf {max_maf=$6}
 $6 < min_maf {min_maf=$6}

 $6 < 0.05 || $6 > 0.95 {rare++}
 $6 >= 0.05 && $6 <= 0.95 {common++}

 NR <= 10 {print} 

 END{printf "\n%s: %'"'"'d\n%s: %'"'"'d\n%s: %'"'"'d\n%s: %.6g-%.6g\n\n",
            "Number of lines",NR,"Common SNPs (MAF >= 0.05)",common,"Rare SNPs (MAF < 0.05)",rare,"MAF range",min_maf,max_maf}' \
<(zcat raw_dl_healthspan.csv.gz)



# Reformat data

echo -n "Reformatting to wilson format... "

zcat raw_dl_healthspan.csv.gz | awk -v FS="," -v OFS="\t" 'BEGIN {print "rsid","snpid","chr","pos","a1","a0","n","freq1","beta1","se","p"}
NR > 1 {print $1,$2"_"$3,$2,$3,$4,$5,300447,$6,-$7,$8,10^-$10}' > st01_01_healthspan.tsv && rm raw_dl_healthspan.csv.gz # Convert lnHR to -lnHR

echo "done."

echo -e "\nHeader:\n"
head st01_01_healthspan.tsv
echo -e "\n\n"


# Remove duplicates

echo -n "Removing all SNPs with non-unique rsids..."

awk 'NR==1 {next} NR == FNR {rsid[$1]++; next} rsid[$1] < 2 {print $0 > ".tmp.healthspan.no.duplicates"}
rsid[$1] > 1 {i++} END {printf " (%s: %'"'"'d) ","N", i}' st01_01_healthspan.tsv st01_01_healthspan.tsv \
&& mv .tmp.healthspan.no.duplicates st01_01_healthspan.tsv

echo "done."


echo -n "Removing all SNPs with non-unique snpids..."

awk 'NR==1 {next} NR == FNR {snpid[$2]++; next} snpid[$2] < 2 {print $0 > ".tmp.healthspan.no.duplicates"}
snpid[$1] > 1 {i++} END {printf " (%s: %'"'"'d) ","N", i}' st01_01_healthspan.tsv st01_01_healthspan.tsv \
&& mv .tmp.healthspan.no.duplicates st01_01_healthspan.tsv

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
st01_01_healthspan.tsv



# Gzip results

echo -n "Gzipping results... "
gzip -9 -f st01_01_healthspan.tsv
echo -e "done.\n\n"

date
echo -e "--- DONE ---\n\n"
