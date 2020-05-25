#!/usr/bin/Rscript

#----
# Setup environment
#----

options(width=200)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggthemes)

#----
# Load data
#----

cat("\n\n==========\nLoading Data\n==========\n\n")

genes <- fread("st03_03_triple_manova_res.tsv")

stats_health <- fread(cmd="zcat st03_05_manova_healthspan_stats.tsv.gz")
stats_long <- fread(cmd="zcat st03_05_manova_longevity_stats.tsv.gz")
stats_life <- fread(cmd="zcat st03_05_manova_lifespan_stats.tsv.gz")

stats_fath <- fread(cmd="zcat st03_05_manova_fath_lifespan_stats.tsv.gz")
stats_moth <- fread(cmd="zcat st03_05_manova_moth_lifespan_stats.tsv.gz")

cat("Done!")


#----
# Combine sex-specific stats
#----

cat("\n\n==========\nCombining Fath/Moth stats\n==========\n\n")

sex_stats <- merge(stats_fath, stats_moth, by=c("rsid","snpid","chr","pos"), suffixes=c("_fath","_moth"))

# Correct and verify allele mismatches

sex_stats[a1_fath==a0_moth, c("a1_fath","a0_fath","freq1_fath","beta1_fath") := .(a0_fath, a1_fath, 1-freq1_fath, -beta1_fath)]
freq_mismatch <- sex_stats[abs(freq1_fath-freq1_moth) > 0.05, ]
cat(freq_mismatch[,.N],"SNPs have allele frequency differences > 0.05\n\n")
if(freq_mismatch[,.N] > 0 ) {print(freq_mismatch); cat("\n\n")}

# Format


sex_stats <- sex_stats[,.(rsid,snpid,chr,pos,a1=a1_fath, a0=a0_fath, n_fath, n_moth, freq1=round((freq1_fath+freq1_moth)/2,4), beta1_fath, se_fath, p_fath, beta1_moth, se_moth, p_moth)]
print(genes[,.(rsid,gene=summary)][sex_stats,,on="rsid"])


#----
# Combine all other stats
#----

cat("\n\n==========\nCombining Long/Health/Life stats\n==========\n\n")

stats_long_health <- merge(stats_long, stats_health, by=c("rsid","snpid","chr","pos"), suffixes=c("_long","_health"))
stats_lhl <- merge(stats_long_health, stats_life, by=c("rsid","snpid","chr","pos"), suffixes=c("","_life"))

# Correct and verify allele mismatches

stats_lhl[a1_long==a0, c("a1_long","a0_long","freq1_long","beta1_long") := .(a0_long, a1_long, 1-freq1_long, -beta1_long)]
stats_lhl[a1_health==a0, c("a1_health","a0_health","freq1_health","beta1_health") := .(a0_health, a1_health, 1-freq1_health, -beta1_health)]
stats_lhl <- stats_lhl[,.(rsid, snpid, chr, pos, a1, a0, n_health, n_life=n, n_long, freq1_health=round(freq1_health,2), freq1_life=round(freq1,2), freq1_long, beta1_health, se_health, p_health, beta1_life=beta1, se_life=se, p_life=p, beta1_long, se_long, p_long)]

freq_mismatch <- stats_lhl[abs(freq1_health-freq1_life) > 0.05 | abs(freq1_long-freq1_life) > 0.05, ]
cat(freq_mismatch[,.N],"SNPs have allele frequency differences > 0.05\n\n")
if(freq_mismatch[,.N] > 0 ) {print(genes[,.(rsid,gene=summary)][freq_mismatch,,on="rsid"]); cat("\n\n")}

freq_caution <- stats_lhl[abs(freq1_health - 0.5) < 0.05 | abs(freq1_life - 0.5) < 0.05 | abs(freq1_long - 0.5) < 0.05,]
cat(freq_caution[,.N],"SNPs have allele frequency differences close to 0.50 -- please check alleles are OK.\n\n")
if(freq_caution[,.N] > 0 ) {print(genes[,.(rsid,gene=summary)][freq_caution,,on="rsid"]); cat("\n\n")}

# Format

stats_lhl <- stats_lhl[,.(rsid,snpid,chr,pos,a1,a0,freq1=round((freq1_health + freq1_life + freq1_long)/3,4),beta1_health,se_health,p_health,beta1_life,se_life,p_life,beta1_long, se_long, p_long)]
print(stats_lhl)



#----
# Merge sex-specific statistics with other stats
#----

cat("\n\n==========\nMerging sex-specific with other stats\n==========\n\n")

stats <- merge(sex_stats, stats_lhl, by=c("rsid","snpid","chr","pos"), suffixes=c("_par",""))
stats[a1_par==a0, c("a1_par","a0_par","freq1_par","beta1_fath","beta1_moth") := .(a0_par, a1_par, 1-freq1_par, -beta1_fath, -beta1_moth)]

freq_mismatch <- stats[abs(freq1_par-freq1) > 0.05, ]
cat(freq_mismatch[,.N],"SNPs have allele frequency differences > 0.05\n\n")
if(freq_mismatch[,.N] > 0 ) {print(freq_mismatch); cat("\n\n")}


print(stats)

#----
# Create table
#----

cat("\n\n==========\nCreating tables\n==========\n\n")

table <- genes[stats,.(gene=summary,rsid,chr,pos,a1,a0,freq1,
	                   beta1_health,se_health,p_health,
	                   beta1_life,se_life,p_life,
	                   beta1_fath_life=beta1_fath,se_fath_life=se_fath, p_fath_life=p_fath,
	                   beta1_moth_life=beta1_moth,se_moth_life=se_moth, p_moth_life=p_moth,
	                   beta1_long, se_long, p_long, 
	                   p_manova=p),on="rsid"]

# Flip alleles to protective
table[(beta1_health/se_health) + (beta1_life/se_life) + (beta1_long/se_long)  < 0, c("a1","a0","freq1","beta1_health","beta1_life","beta1_fath_life","beta1_moth_life","beta1_long"):=.(a0,a1,1-freq1,-beta1_health,-beta1_life,-beta1_fath_life,-beta1_moth_life,-beta1_long)]

table[rsid=="rs8042849", gene:="CHRNA3/5"]
table[gene=="HLA-DRB1/HLA-DQA1",gene:="HLA-DRB1/DQA1"]
table[gene=="SLC22A2/SLC22A3", gene:="SLC22A2/3"]
table[gene=="NOL4L/NOL4L-DT", gene:="NOL4L/-DT"]

table <- table[order(chr,pos)]

print(table)




cat("\n\n==========\nHealthspan sample size comparison\n==========\n\n")

ss <- table[stats_health[,.(rsid,n)], , on="rsid"][, .(gene, p_manova=sprintf("%.2g", p_manova), p_health=sprintf("%.2g", p_health), n_ratio=(qnorm(p_manova/2)/(qnorm(p_health/2)))^2)]
print(ss)


cat("\nLoci of interest sample size equivalent:\n\n")

n_health = 300477
ss1 <- ss[gene %in% c("SLC4A7","LINC02513","FOXO3","LPA","CDKN2B-AS1","ZW10","FGD6","TOX3","LDLR","APOE")]
ss2 <- ss1[, .(gene, n_health, n_equiv=round(n_ratio * n_health), n_increase= round(n_health * (n_ratio - 1)))]

print(ss2)
cat("\n")

# https://www.ucl.ac.uk/child-health/short-courses-events/about-statistical-courses/research-methods-and-statistics/chapter-8-content-8
n_equiv <- sort(ss2$n_equiv)
median <- median(ss2$n_equiv)
lci <- n_equiv[round(ss2[,.N]/2 - qnorm(0.975)*sqrt(ss2[,.N])/2)]
uci <- n_equiv[1 + round(ss2[,.N]/2 + qnorm(0.975)*sqrt(ss2[,.N])/2)]

print(
	data.frame(` `="<MEDIAN>", n_health, n_equiv=median, n_equiv_lci=lci, n_equiv_uci=uci, 
	n_increase=round(median-n_health), n_increase_lci=round(lci-n_health), n_increase_uci=round(uci-n_health), 
	pct_increase=round((median-n_health)/n_health,2), pct_increase_lci=round((lci-n_health)/n_health,2), pct_increase_uci=round((uci-n_health)/n_health,2),
	check.names=F),
row.names=F)



cat("\n\n==========\nLifespan sample size comparison\n==========\n\n")

ss <- table[stats_life[,.(rsid,n)], , on="rsid"][, .(gene, p_manova=sprintf("%.2g", p_manova), p_life=sprintf("%.2g", p_life), n_ratio=(qnorm(p_manova/2)/(qnorm(p_life/2)))^2)]
print(ss)


cat("\nLoci of interest sample size equivalent:\n\n")

n_life = 1012240
ss1 <- ss[gene %in% c("SLC4A7","LINC02513","FOXO3","LPA","CDKN2B-AS1","ZW10","FGD6","TOX3","LDLR","APOE")]
ss2 <- ss1[, .(gene, n_life, n_equiv=round(n_ratio * n_life), n_increase= round(n_life * (n_ratio - 1)))]

print(ss2)
cat("\n")

# https://www.ucl.ac.uk/child-health/short-courses-events/about-statistical-courses/research-methods-and-statistics/chapter-8-content-8
n_equiv <- sort(ss2$n_equiv)
median <- median(ss2$n_equiv)
lci <- n_equiv[round(ss2[,.N]/2 - qnorm(0.975)*sqrt(ss2[,.N])/2)]
uci <- n_equiv[1 + round(ss2[,.N]/2 + qnorm(0.975)*sqrt(ss2[,.N])/2)]

print(
	data.frame(` `="<MEDIAN>", n_life, n_equiv=median, n_equiv_lci=lci, n_equiv_uci=uci, 
	n_increase=round(median-n_life), n_increase_lci=round(lci-n_life), n_increase_uci=round(uci-n_life), 
	pct_increase=round((median-n_life)/n_life,2), pct_increase_lci=round((lci-n_life)/n_life,2), pct_increase_uci=round((uci-n_life)/n_life,2),
	check.names=F),
row.names=F)



cat("\n\n==========\nLongevity sample size comparison\n==========\n\n")

ss <- table[stats_long[,.(rsid,n)], , on="rsid"][, .(gene, p_manova=sprintf("%.2g", p_manova), p_long=sprintf("%.2g", p_long), n_ratio=(qnorm(p_manova/2)/(qnorm(p_long/2)))^2)]
print(ss)


cat("\nLoci of interest sample size equivalent:\n\n")

n_long = 11262
ss1 <- ss[gene %in% c("SLC4A7","LINC02513","FOXO3","LPA","CDKN2B-AS1","ZW10","FGD6","TOX3","LDLR","APOE")]
ss2 <- ss1[, .(gene, n_long, n_equiv=round(n_ratio * n_long), n_increase= round(n_long * (n_ratio - 1)))]

print(ss2)
cat("\n")

# https://www.ucl.ac.uk/child-health/short-courses-events/about-statistical-courses/research-methods-and-statistics/chapter-8-content-8
n_equiv <- sort(ss2$n_equiv)
median <- median(ss2$n_equiv)
lci <- n_equiv[round(ss2[,.N]/2 - qnorm(0.975)*sqrt(ss2[,.N])/2)]
uci <- n_equiv[1 + round(ss2[,.N]/2 + qnorm(0.975)*sqrt(ss2[,.N])/2)]

print(
	data.frame(` `="<MEDIAN>", n_long, n_equiv=median, n_equiv_lci=lci, n_equiv_uci=uci, 
	n_increase=round(median-n_long), n_increase_lci=round(lci-n_long), n_increase_uci=round(uci-n_long), 
	pct_increase=round((median-n_long)/n_long,2), pct_increase_lci=round((lci-n_long)/n_long,2), pct_increase_uci=round((uci-n_long)/n_long,2),
	check.names=F),
row.names=F)




#----
# Export
#----


table1 <- table[,.SD,.SDcols=names(table)[!grepl("fath|moth",names(table))]]
fwrite(table1, "st03_06_triple_manova_table.csv", sep=",",quote=F,na="NA" )

table2 <- table[,.SD,.SDcols=names(table)[!grepl("[^th]_life",names(table))]]
table2[,p_sex:=2*pnorm(-abs((beta1_fath_life-beta1_moth_life)/sqrt(se_fath_life^2 + se_moth_life^2)))]

fwrite(table2, "st03_06_triple_manova_table_sex.csv", sep=",",quote=F,na="NA" )


#----
# Plot sex differences
#----

cat("\n\n==========\nCreating sex-difference plot\n==========\n\n")

t1 <- melt(table2, id.vars=c("gene","rsid","a1","p_sex"), measure.vars=c("beta1_fath_life","beta1_moth_life","se_fath_life","se_moth_life"))
t1[,c("par","variable"):=.(ifelse(grepl("fath",variable),"fath","moth"),gsub("_..th_life","",variable))]
t1[, label:=paste0(gene," (",rsid,"_",a1,")")]
t2 <- dcast(t1, label+par+p_sex~variable, value.var="value")
t2[,c("ymin","ymax"):=.(beta1 - qnorm(0.975) * se, beta1 + qnorm(0.975) * se)]

ord <- t2[order(p_sex),unique(label)]
t2[,label:=factor(label, levels=rev(ord))]

star <- table2[,.(label=paste0(gene," (",rsid,"_",a1,")"), beta1=(beta1_fath_life+beta1_moth_life)/2, 
	star=ifelse(p_sex <= 0.0005,"***",ifelse(p_sex < 0.005, "**", ifelse(p_sex < 0.05, "*", ""))))]


pdf("st03_06_triple_manova_sex.1.pdf", width=8, height=7)

p1 <- ggplot(t2, aes(x=label, y=beta1, colour=par)) +
geom_hline(yintercept=0, linetype=3) + geom_linerange(aes(ymin=ymin, ymax=ymax), size=1.1, alpha=0.5) + geom_point(size=3) + coord_flip() +
geom_text(data=star, aes(x=label, y=beta1, label=star), col="black", size=8 ) +
ylab("Lifespan (-lnHR)") + xlab("") +
scale_colour_manual(values=c("#2166ac","#b2182b"),breaks=c("moth","fath"),labels=c("Mothers","Fathers"), name="") + 
theme_pubr()  + theme(legend.position="right")
print(p1)

dev.off()

ggsave("st03_06_triple_manova_sex.1.svg", p1, width=8, height=7)

# Loci of interest


table2[,group:=ifelse(p_sex < 0.05, "A", ifelse(p_sex < 0.5, "B", "C"))]


loci <- table1[p_health < 0.05 & p_life < 0.05 & p_long < 0.05, gene]

t1 <- melt(table2[gene %in% loci], id.vars=c("gene","rsid","a1","p_sex","group"), measure.vars=c("beta1_fath_life","beta1_moth_life","se_fath_life","se_moth_life"))
t1[,c("par","variable"):=.(ifelse(grepl("fath",variable),"fath","moth"),gsub("_..th_life","",variable))]
t1[, label:=paste0(gene," (",rsid,"_",a1,")")]
t2 <- dcast(t1, label+par+p_sex+group~variable, value.var="value")
t2[,c("ymin","ymax"):=.(beta1 - qnorm(0.975) * se, beta1 + qnorm(0.975) * se)]

ord <- t2[order(p_sex),unique(label)]
t2[,label:=factor(label, levels=rev(ord))]




#lims <- data.table(x=c(0,0.15), y=t2$label, group=t2$group)
#lims[group!="C",x:=NA]

pdf("st03_06_triple_manova_sex.2.pdf", width=8, height=4)

p2 <- ggplot(t2, aes(x=label, y=beta1, colour=par)) +
#geom_hline(yintercept=seq(0, 0.15, 0.05), colour="grey90") +
#geom_linerange(aes(x=label), ymin=-Inf, ymax=Inf, colour="grey90") +
#geom_rangeframe(data=lims, aes(x=y, y=x), colour="black") +
geom_hline(yintercept=0, linetype=3) + 
geom_linerange(aes(ymin=ymin, ymax=ymax), size=1.1, alpha=0.5) + geom_point(size=3) + coord_flip() +
#facet_grid(group~., scales="free_y", space="free") +
ylab("Lifespan (-lnHR)") + xlab("") +
scale_colour_manual(values=c("#2166ac","#b2182b"),breaks=c("moth","fath"),labels=c("Mothers","Fathers"), name="") + 
theme_bw() #+
#theme(legend.position="right", panel.grid.minor=element_blank(), panel.border=element_blank(), panel.spacing = unit(-0.5, "lines"), panel.grid=element_blank(), strip.text=element_blank(), strip.background=element_blank())
print(p2)

dev.off()


ggsave("st03_06_triple_manova_sex.2.svg", p2, width=8, height=4)

cat("--- Done ---\n")