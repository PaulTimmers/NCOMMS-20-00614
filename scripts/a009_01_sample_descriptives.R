#!/usr/bin/env Rscript

#----
# Setup environment
#----

options(width=200)
library(data.table)
library(survival)

#----
# Get stats
#----

stats <- fread("st02_02_table_1_data.txt")

counts_n <- stats[stat=="count" & sex=="ALL", .SD, .SDcols=grepl("age",colnames(stats))]
counts_dead <- stats[stat=="sum" & sex=="ALL", .SD, .SDcols=grepl("dead",colnames(stats))]

age_stats <- stats[sex=="ALL" & stat %in% c("mean","sd"), .SD, .SDcols=c("stat",grep("age",colnames(stats), value=TRUE))]

surv <- fread("st02_01_qcd_phenotypes.txt")
surv <- surv[!is.na(g_bri) & !is.na(unrelated) & !is.na(is_not_withdrawn_20181016),]


#----
# Format
#----

# Counts

mcounts_n <- melt(counts_n, value.name="n", variable.name="group")
mcounts_n[,group:=sub("_age","",group)]
mcounts_dead <- melt(counts_dead, value.name="dead", variable.name="group")
mcounts_dead[,group:=sub("_dead","",group)]


# Age distribution

mage_stats <- melt(age_stats, id.var="stat", value.name="age", variable.name="group")
mage_stats[,group:=sub("_age","",group)]
dage_stats <- dcast(mage_stats, group~stat, value.var="age")
names(dage_stats)[-1] <- paste0(names(dage_stats)[-1], "_age")



# Median survival

fath <- surv[,.(iid=paste0(iid,"FATH"), age=fath_age, dead=fath_dead, par="fath")]
moth <- surv[,.(iid=paste0(iid,"MOTH"), age=moth_age, dead=moth_dead, par="moth")]
ph <- rbind(fath, moth)
ph[,band:=ifelse(age > 80, "80_120", ifelse(age > 60, "60_80", ifelse(age > 40, "40_60", NA)))]
m0 <- survfit(Surv(age, dead) ~ par, data=ph); q0 <- quantile(m0, probs=0.5, conf.int=T)
m1 <- survfit(Surv(age, dead) ~ par + band, data=ph); q1 <- quantile(m1, probs=0.5, conf.int=T)

medians <- data.table(strata=c(rownames(q0$quantile),rownames(q1$quantile)), median_survival=c(q0$quantile,q1$quantile))
medians[,parents:=c("Fathers","Mothers",rep(c("Fathers","Mothers"),each=3))]
medians[,age_band:=c("ALL","ALL",rep(c("40-60","60-80","80+"), times=2))]


#----
# Merge
#----

mcounts <- mcounts_n[mcounts_dead,,on="group"]
all_stats <- mcounts[dage_stats,,on="group"]

# Format names
nlist <- c(fath="ALL", moth="ALL", `40_60`="40-60", `60_80`="60-80", `80_120`="80+")
all_stats[,c("parents","age_band","group"):=.(ifelse(grepl("fath",group),"Fathers","Mothers"), nlist[gsub(".*th_","",group)],NULL)]

all_stats1 <- all_stats[medians,,on=c("parents","age_band")]
all_stats2 <- all_stats1[,.(parents, age_band, n, dead, mean_age, sd_age, median_survival)]

fwrite(all_stats2,"st02_03_n_dead_stats.txt", sep="\t", na="NA", quote=FALSE)