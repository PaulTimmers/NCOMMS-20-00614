#!/usr/bin/env Rscript

#----
# Setup environment
#----

options(width=200)
library(data.table)
library(dplyr)

#----
# START
#----

# Load data

manova <- fread("../q02_manova/st02_01_triple_manova.tsv.gz", select=1:12)
healthspan <- fread("../q01_sumstats/st01_01_healthspan.tsv.gz", select=c(2,5,6,9,10))
lifespan <- fread("../q01_sumstats/st01_01_lifespan.tsv.gz", select=c(2,5,6,9,10))
longevity <- fread("../q01_sumstats/st01_01_longevity_90_pct.tsv.gz", select=c(2,5,6,9,10))


# Merge data

df1 <- healthspan[manova,.(rsid, snpid, chr, pos, a1=i.a1, a0=i.a0, n, freq1, beta1=i.beta1, se=i.se, p, info, beta1_healthspan=ifelse(a1==i.a1, beta1, -beta1), se_healthspan=se),on="snpid"]
df2 <- lifespan[df1,.(rsid, snpid, chr, pos, a1=i.a1, a0=i.a0, n, freq1, beta1=i.beta1, se=i.se, p, info, beta1_healthspan, se_healthspan, beta1_lifespan=ifelse(a1==i.a1, beta1, -beta1), se_lifespan=se),on="snpid"]
df3 <- longevity[df2,.(rsid, snpid, chr, pos, a1=i.a1, a0=i.a0, n, freq1, beta1=i.beta1, se=i.se, p, info, beta1_healthspan, se_healthspan, beta1_lifespan, se_lifespan, beta1_longevity=ifelse(a1==i.a1, beta1, -beta1), se_longevity=se),on="snpid"]

# Calculate Z scores

df3[,c("z1_healthspan","z1_lifespan","z1_longevity"):=.(beta1_healthspan/se_healthspan, beta1_lifespan/se_lifespan, beta1_longevity/se_longevity)]
df3[,z_add := z1_healthspan + z1_lifespan + z1_longevity]
df3[,z1 := sign(z_add) * abs(qnorm(p/2))]


# Effective sample size based on events

df3[, n := round(84949 + 609139 + 2/(1/11262 + 1/25483))]

# Calculate betas and SE based on sample size and allele frequency

df3[, se :=  1 / sqrt( 2 * freq1 * (1 - freq1) * (n + z1^2)) ]
df3[, beta1 := z1 * se]

out1 <- df3[,.SD,.SDcols=!grepl("z", names(df3))]


out2 <- out1 %>% mutate_at(grep("beta|se",names(out1),value=T), round, digits=6)

fwrite(out2, "st012_01_manova.tsv.gz", quote=F, sep="\t", na="NA")