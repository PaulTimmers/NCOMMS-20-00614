#$ -N lz_manova
#$ -l h_vmem=8G
#$ -l h_rt=00:15:00
#$ -j y
#$ -cwd
#$ -o log_lz_st03_04_manova.$TASK_ID.log
#$ -V

truncate -s 0 $SGE_STDOUT_PATH

#----
# Setup environment
#----

module load igmm/apps/R/3.6.0

snp=`awk -v i=$SGE_TASK_ID 'NR == i+1 {print $1}' ../st03_03_triple_manova_res.tsv`
gene=`awk -v i=$SGE_TASK_ID 'NR == i+1 {gsub("/","_",$NF); print tolower($NF)}' ../st03_03_triple_manova_res.tsv`


#----
# START
#----

locuszoom \
--build hg19 \
--source 1000G_Nov2014 \
--plotonly \
--pvalcol p \
--no-date \
--markercol rsid \
--metal ../../q02_manova/st02_01_triple_manova.tsv.gz \
--refsnp $snp \
--flank 500kb \
--pop EUR \
--prefix lz_st03_03_manova_${gene} \
--rundir ./ \
--cache 'None' \
--delim tab  
