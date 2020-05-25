#!/bin/bash

#----
# Setup environment
#----

MAIN=${ROOT}/projects/prj_108_lifegen_phase3/p06_triple_manova
cd $MAIN

#----
# START
#----

mkdir -p q02_manova
cd q02_manova

HEALTHSPAN="../q01_sumstats/st01_01_healthspan.tsv.gz"
LONGEVITY="../q01_sumstats/st01_01_longevity_90_pct.tsv.gz"
LIFESPAN="../q01_sumstats/st01_01_lifespan.tsv.gz"

Rscript ../scripts/a002_run_metabel.R $HEALTHSPAN $LONGEVITY $LIFESPAN st02_01_triple_manova.tsv &> log_a002_run_metabel.log &
gzip -9 -f st02_01_triple_manova.tsv

