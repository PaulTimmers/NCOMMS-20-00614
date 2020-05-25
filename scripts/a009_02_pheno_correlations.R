#!/usr/bin/env R

#----
# Setup environment
#----

options(width=200)
library(data.table)
load("../../data/d001_indep.snps.RData")

#----
# Load data
#----

healthspan <- fread(cmd="zcat ../../q01_sumstats/st01_01_healthspan.tsv.gz", select=c("rsid","snpid","a1","a0","freq1","beta1","se","p"))
life4060 <- fread(cmd="zcat ../q01_combined/p05_comb_chr/bothpl_lif4060.tsv.gz", select=c("rsid","snpid","a1","a0","freq1","beta1","se","p"))
life6080 <- fread(cmd="zcat ../q01_combined/p05_comb_chr/bothpl_lif6080.tsv.gz", select=c("rsid","snpid","a1","a0","freq1","beta1","se","p"))
lifegen <- fread(cmd="zcat ../../q01_sumstats/st01_01_lifegen.tsv.gz", select=c("rsid","snpid","a1","a0","freq1","beta1","se","p"))
life80120 <- fread(cmd="zcat ../q01_combined/p05_comb_chr/bothpl_lif80120.tsv.gz", select=c("rsid","snpid","a1","a0","freq1","beta1","se","p"))
longevity <- fread(cmd="zcat ../../q01_sumstats/st01_01_longevity_90_pct.tsv.gz", select=c("rsid","snpid","a1","a0","freq1","beta1","se","p"))

bothpl <- fread(cmd="zcat ../q01_combined/p05_comb_chr/bothpl_life.tsv.gz", select=c("rsid","snpid","a1","a0","freq1","beta1","se","p"))

# Subset

healthspan <- healthspan[rsid %in% indep.snps]
life4060 <- life4060[rsid %in% indep.snps]
life6080 <- life6080[rsid %in% indep.snps]
lifegen <- lifegen[rsid %in% indep.snps]
life80120 <- life80120[rsid %in% indep.snps]
longevity <- longevity[rsid %in% indep.snps]

bothpl <- bothpl[rsid %in% indep.snps]


#----
# Age-stratified correlations 
#----

m1 <- merge(healthspan, life4060, by="snpid", suffix=c("_health","_4060"))
m1[a1_health==a0_4060,c("a1_health","a0_health","freq1_health","beta1_health"):=.(a0_health, a1_health,1-freq1_health,-beta1_health)]
m1 <- m1[a1_health==a1_4060,.(rsid=rsid_health,snpid,a1=a1_health,a0=a0_health,
	freq1_health,freq1_4060,
	beta1_health,se_health,p_health,
	beta1_4060,se_4060,p_4060)]

m2 <- merge(m1, life6080, by="snpid", suffix=c("","_6080"))
m2[a1_6080==a0,c("a1_6080","a0_6080","freq1","beta1"):=.(a0_6080, a1_6080,1-freq1,-beta1)]
m2 <- m2[a1_6080==a1,.(rsid,snpid,a1,a0,
	freq1_health,freq1_4060,freq1_6080=freq1,
	beta1_health,se_health,p_health,
	beta1_4060,se_4060,p_4060,
	beta1_6080=beta1,se_6080=se,p_6080=p)]

m3 <- merge(m2, lifegen, by="snpid", suffix=c("","_life"))
m3[a1_life==a0,c("a1_life","a0_life","freq1","beta1"):=.(a0_life, a1_life,1-freq1,-beta1)]
m3 <- m3[a1_life==a1,.(rsid,snpid,a1,a0,
	freq1_health,freq1_4060,freq1_6080,freq1_life=freq1,
	beta1_health,se_health,p_health,
	beta1_4060,se_4060,p_4060,
	beta1_6080,se_6080,p_6080,
	beta1_life=beta1,se_life=se,p_life=p)]

m4 <- merge(m3, life80120, by="snpid", suffix=c("","_80120"))
m4[a1_80120==a0,c("a1_80120","a0_80120","freq1","beta1"):=.(a0_80120, a1_80120,1-freq1,-beta1)]
m4 <- m4[a1_80120==a1,.(rsid,snpid,a1,a0,
	freq1_health,freq1_4060,freq1_6080,freq1_life,freq1_80120=freq1,
	beta1_health,se_health,p_health,
	beta1_4060,se_4060,p_4060,
	beta1_6080,se_6080,p_6080,
	beta1_life,se_life,p_life,
	beta1_80120=beta1,se_80120=se,p_80120=p)]

m5 <- merge(m4, longevity, by="snpid", suffix=c("","_long"))
m5[a1_long==a0,c("a1_long","a0_long","freq1","beta1"):=.(a0_long, a1_long,1-freq1,-beta1)]
m5 <- m5[a1_long==a1,.(rsid,snpid,a1,a0,
	freq1_health,freq1_4060,freq1_6080,freq1_life,freq1_80120,freq1_long=freq1,
	beta1_health,se_health,p_health,
	beta1_4060,se_4060,p_4060,
	beta1_6080,se_6080,p_6080,
	beta1_life,se_life,p_life,
	beta1_80120,se_80120,p_80120,
	beta1_long=beta1,se_long=se,p_long=p)]

cor <- cor(m5[,.(beta1_health, beta1_4060, beta1_6080, beta1_life, beta1_80120, beta1_long)])

cor2 <- cbind(data.frame(x=row.names(cor)),cor)
write.table(cor2, "st02_stratified_pheno_correlations.tsv", quote=FALSE, sep="\t", row.names=FALSE)


#----
# Full correlations 
#----

m1 <- merge(healthspan, lifegen, by="snpid", suffix=c("_health","_lifegen"))
m1[a1_health==a0_lifegen,c("a1_lifegen","a0_lifegen","freq1_lifegen","beta1_lifegen"):=.(a0_lifegen, a1_lifegen, 1-freq1_lifegen,-beta1_lifegen)]
m1 <- m1[a1_health==a1_lifegen,.(rsid=rsid_health,snpid,a1=a1_health,a0=a0_health,
	freq1_health,freq1_lifegen,
	beta1_health,se_health,p_health,
	beta1_lifegen,se_lifegen,p_lifegen)]

m2 <- merge(m1, bothpl, by="snpid", suffix=c("","_bothpl"))
m2[a1_bothpl==a0,c("a1_bothpl","a0_bothpl","freq1","beta1"):=.(a0_bothpl, a1_bothpl,1-freq1,-beta1)]
m2 <- m2[a1_bothpl==a1,.(rsid,snpid,a1,a0,
	freq1_health,freq1_lifegen,freq1_bothpl=freq1,
	beta1_health,se_health,p_health,
	beta1_lifegen,se_lifegen,p_lifegen,
	beta1_bothpl=beta1,se_bothpl=se,p_bothpl=p)]

m3 <- merge(m2, longevity, by="snpid", suffix=c("","_long"))
m3[a1_long==a0,c("a1_long","a0_long","freq1","beta1"):=.(a0_long, a1_long,1-freq1,-beta1)]
m3 <- m3[a1_long==a1,.(rsid,snpid,a1,a0,
	freq1_health,freq1_lifegen,freq1_bothpl,freq1_long=freq1,
	beta1_health,se_health,p_health,
	beta1_lifegen,se_lifegen,p_lifegen,
	beta1_bothpl,se_bothpl,p_bothpl,
	beta1_long=beta1,se_long=se,p_long=p)]

cor <- cor(m3[,.(beta1_health, beta1_lifegen, beta1_bothpl, beta1_long)])

cor2 <- cbind(data.frame(x=row.names(cor)),cor)
write.table(cor2, "st02_full_pheno_correlations.tsv", quote=FALSE, sep="\t", row.names=FALSE)
