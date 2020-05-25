#!/bin/bash

#----
# Setup environment
#----

MAIN=${ROOT}/projects/prj_108_lifegen_phase3/p06_triple_manova
cd $MAIN

#----
# START
#----

mkdir -p ${MAIN}/q11_fuma
cd ${MAIN}/q11_fuma


#----
# Get datasets
#----

mkdir -p datasets
# Download http://software.broadinstitute.org/gsea/downloads.jsp#msigdb GO HALLMARK as hallmark.entrez.gmt
# Download http://software.broadinstitute.org/gsea/downloads.jsp#msigdb GO BIOLOGICAL PROCESSES as go_bioproc.entrez.gmt


cat datasets/c5.bp.v7.0.symbols.gmt | tr -s '\t' '\n' | sed '/http/d' | awk '$1 ~ /GO_/ {col++; row=1} {data[row][col]=$1; row++} END{ {for (row in data) { for(col in data[row]) {printf "%s\t",data[row][col]} {print "\n"}}   }  }'
awk 'FNR == 1 {print; exit}' datasets/c5.bp.v7.0.symbols.gmt | tr -s '\t' '\n' | sed '/http/d'


# Magic AWK code which creates a heatmap out of the horrible gmt format

for collection in hallmark go_bioproc
do
cat ${collection}.entrez.gmt | awk '{for (i=3; i<=NF; i++) {all_genes[$i]++; this_set[NR]=$1; this_set_genes[NR][$i]++}} END{ {printf "%s\t","gene_set"} {for (a in all_genes) {printf "%s\t",a}} {printf "\n",1} {for (s in this_set) {{ printf "%s\t",this_set[s]} {for (a in all_genes) {printf "%s\t",(this_set_genes[s][a] > 0 ? 1 : 0)} {printf "\n",1}}}}}' | gzip -9 -f > ${collection}.entrez.heat.gz
done


#----
# Calculate enrichment
#----

Rscript ../scripts/a011_gene_set_enrichment.R



#----
# Plot results
#----

Rscript ../scripts/a011_pretty_plot_fuma.R