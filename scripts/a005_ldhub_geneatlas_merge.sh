#!/bin/bash

#$ -N merge.GeneAtlas
#$ -hold_jid GeneAtlas
#$ -l h_vmem=1G
#$ -l h_rt=00:10:00
#$ -j y
#$ -o ${ROOT}/projects/prj_108_lifegen_phase3/p06_triple_manova/q05_ldsc/ldhub/t006_geneatlas.log
#$ -cwd
#$ -V

truncate -s 0 $SGE_STDOUT_PATH

#----
# Setup environment
#----

MAIN=${ROOT}/projects/prj_108_lifegen_phase3/p06_triple_manova
cd ${MAIN}/q05_ldsc/ldhub

#----
# Merge
#----

n_traits=`tail -n+2 t000_geneatlas_traits.txt | wc -l`

echo "trait1,trait2,PMID,Category,ethnicity,note,rg,se,z,p,h2_obs,h2_obs_se,h2_int,h2_int_se,gcov_int,gcov_int_se" > t006_healthspan_geneatlas.csv
echo "trait1,trait2,PMID,Category,ethnicity,note,rg,se,z,p,h2_obs,h2_obs_se,h2_int,h2_int_se,gcov_int,gcov_int_se" > t006_lifespan_geneatlas.csv
echo "trait1,trait2,PMID,Category,ethnicity,note,rg,se,z,p,h2_obs,h2_obs_se,h2_int,h2_int_se,gcov_int,gcov_int_se" > t006_longevity_geneatlas.csv

for i in `eval echo {1..${n_traits}}`
do
cat t006_geneatlas.${i}.log && rm t006_geneatlas.${i}.log
cat t006_healthspan_geneatlas.${i}.csv >> t006_healthspan_geneatlas.csv && rm t006_healthspan_geneatlas.${i}.csv
cat t006_lifespan_geneatlas.${i}.csv >> t006_lifespan_geneatlas.csv && rm t006_lifespan_geneatlas.${i}.csv
cat t006_longevity_geneatlas.${i}.csv >> t006_longevity_geneatlas.csv && rm t006_longevity_geneatlas.${i}.csv
done