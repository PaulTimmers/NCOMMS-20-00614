#!/bin/bash

#----
# Setup environment
#----

sumstats=$1
p=$2
r2=$3
outfile=$4

function gzfread(){
    for x in "$@"; do [[ -f $x ]] && break; done
    file `readlink -f $x` | grep -q gzip && zcat "$@" || cat "$@"
}


#================ Define clump function ================================================


function ld_clump(){

SUMSTATS=$1
P=$2
R2=$3
OUTFILE=$4
CHR=$5

echo -e "\n\n========================================================"
echo "Sumstats: $SUMSTATS"
echo "Chromosome: $CHR"
echo "========================================================\n"

    function gzfread(){
        for x in "$@"; do [[ -f $x ]] && break; done
        file `readlink -f $x` | grep -q gzip && zcat "$@" || cat "$@"
    }


# Inner join sumstat SNPs x LD reference SNPs

bfile=${ROOT}/data/reference_data/1000genomes/phase3/plink/1000genomes_chr${CHR}
awk -v OFS="\t" 'NR == 1 {print "SNP","P"; next} 
                 NR == FNR {var_id[$1] = $11; next}
                 var_id[$2] > 0 {print $2, var_id[$2]}' \
<(gzfread $SUMSTATS) $bfile.bim > ${OUTFILE}.snps.${CHR}.tsv

N_SNPS=`wc -l < ${OUTFILE}.snps.${CHR}.tsv`
echo -e "Number of SNPs in common between LD reference and sumstats: $N_SNPS\n"

# Stop if no SNPs

if [[ $N_SNPS -eq 0 ]]; then
    rm ${OUTFILE}.snps.${CHR}.tsv
    echo "No clumping necessary."
    echo -e "\n---- DONE ----\n\n"
    exit 0
fi


# Remove duplicates

awk 'NR == FNR {var_id[$1]++; next} var_id[$1] == 1 {print}' ${OUTFILE}.snps.${CHR}.tsv ${OUTFILE}.snps.${CHR}.tsv \
> ${OUTFILE}.snps.unique.${CHR}.tsv && mv ${OUTFILE}.snps.unique.${CHR}.tsv ${OUTFILE}.snps.${CHR}.tsv

N_SNPS_UNIQ=`wc -l < ${OUTFILE}.snps.${CHR}.tsv`

if [[ $N_SNPS -ne $N_SNPS_UNIQ ]]; then
    echo -e "Removed `echo "(${N_SNPS}-${N_SNPS_UNIQ})/2" | bc` duplicated SNPs.\n"
fi

# Clump

plink \
--bfile $bfile \
--no-sex \
--no-pheno \
--extract ${OUTFILE}.snps.${CHR}.tsv \
--clump ${OUTFILE}.snps.${CHR}.tsv \
--clump-field P \
--clump-p1 $P \
--clump-p2 1 \
--clump-r2 $R2 \
--clump-kb 1000 \
--out ${OUTFILE}.${CHR} && rm ${OUTFILE}.snps.${CHR}.tsv ${OUTFILE}.${CHR}.nosex ${OUTFILE}.${CHR}.log || rm ${OUTFILE}.snps.${CHR}.tsv ${OUTFILE}.${CHR}.nosex

# Format output

if [[ ! -f ${OUTFILE}.${CHR}.clumped ]]; then
    echo -e "\nNumber of clumps in this chromosome: 0"
    exit 0
fi

N_SNPS_CLUMP=`wc -l < ${OUTFILE}.${CHR}.clumped`

if [[ ! -z $N_SNPS_CLUMP  ]]; then
    cat ${OUTFILE}.${CHR}.clumped | awk -v OFS="\t" '$0 ~/rs/ {print $3,$1"_"$4,$5}' > tmp && mv tmp ${OUTFILE}.${CHR}.clumped
else
    N_SNPS_CLUMP=0
fi

echo -e "\nNumber of clumps in this chromosome: $N_SNPS_CLUMP"
echo -e "\n---- DONE ----\n\n"

}

export -f ld_clump


#================ End clump function ================================================



#----
# START
#----

chrs=`gzfread $sumstats | awk 'NR > 1 && chrs[$3] == 0 {chrs[$3]++; next} END {for (i in chrs) {print i}} '`
parallel --no-notice -j 11 -k ld_clump $sumstats $p $r2 $outfile ::: $chrs


#----
# Combine and format
#----

awk 'ARGIND < ARGC-1 {hits[$2]++; next} FNR == 1 || hits[$2] > 0' ${outfile}.*.clumped <(gzfread $sumstats) \
> ${outfile}.tsv && rm ${outfile}.*.clumped


