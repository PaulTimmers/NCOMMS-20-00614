#!/bin/bash

#----
# Setup environment
#----

MAIN=${ROOT}/projects/prj_108_lifegen_phase3/p06_triple_manova
cd $MAIN


#----
# START
#----

mkdir -p ${MAIN}/q01_sumstats
cd ${MAIN}/q01_sumstats

nohup bash ../scripts/a001_get_healthspan_data.sh &> log_st01_get_healthspan_data.log &
nohup bash ../scripts/a001_get_lifespan_data.sh &> log_st01_get_lifespan_data.log &
nohup bash ../scripts/a001_get_longevity_data.sh &> log_st01_get_longevity_data.log &

