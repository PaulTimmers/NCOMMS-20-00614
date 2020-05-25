#!/usr/bin/Rscript

#----
# Setup environment
#----

options(width=200)
library(data.table)
library(ggplot2)
library(ggpubr)

#----
# Load data
#----

cat("\n\n==========\nLoading Data\n==========\n\n")

stats_manova <- fread("zcat st03_01_triple_manova_region_res.tsv.gz", select=c("gene","rsid_lead","rsid","chr","pos","p"))
stats_health <- fread(cmd="zcat st03_01_healthspan_region_res.tsv.gz", select=c("gene","rsid_lead","rsid","chr","pos","p"))
stats_long <- fread(cmd="zcat st03_01_longevity_region_res.tsv.gz", select=c("gene","rsid_lead","rsid","chr","pos","p"))
stats_life <- fread(cmd="zcat st03_01_lifespan_region_res.tsv.gz", select=c("gene","rsid_lead","rsid","chr","pos","p"))

cat("Done!")


#----
# Combine all regional stats
#----

cat("\n\n==========\nCombining Long/Health/Life region stats\n==========\n\n")

stats1 <- merge(stats_manova, stats_long, by=c("gene","rsid_lead"), suffixes=c("_manova","_long"))
stats2 <- merge(stats1, stats_life, by=c("gene","rsid_lead"), suffixes=c("","_life"))
stats3 <- merge(stats2, stats_health, by=c("gene","rsid_lead"), suffixes=c("","_health"))


#----
# Create table
#----

cat("\n\n==========\nCreating tables\n==========\n\n")

table <- stats3[,.(gene, rsid_manova, chr_manova, pos_manova, p_manova, rsid_health, pos_health, p_health,
	                   rsid_life=rsid, pos_life=pos, p_life=p, rsid_long,  pos_long, p_long)]

table[rsid_manova=="rs140570886", gene:="LPAL2/LPA"]
table[rsid_manova=="rs8042849", gene:="CHRNA3/5"]
table[gene=="HLA-DRB1/HLA-DQA1",gene:="HLA-DRB1/DQA1"]
table[gene=="SLC22A2/SLC22A3", gene:="SLC22A2/3"]
table[gene=="NOL4L/NOL4L-DT", gene:="NOL4L"]

table <- table[order(chr_manova,pos_manova)]

print(table)

#----
# Export
#----

fwrite(table, "st03_06_triple_manova_region_table.csv", sep=",",quote=F,na="NA" )