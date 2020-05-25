#!/usr/bin/env Rscript

#----
# Setup environment
#----

options(width=200)
library(data.table)
library(survival)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggthemes)
library(meta)

#----
# Get GWAS data
#----

df <- data.table()

loci <- fread("st10_04_all_loci.tsv")

stats_bothpl <- fread("st10_05_bothpl_all.tsv", select=c("rsid","a1","a0","beta1","se","p","par","band"))
stats_fath <- fread("st10_05_fath_all.tsv", select=c("rsid","a1","a0","beta1","se","p","par","band"))
stats_moth <- fread("st10_05_moth_all.tsv", select=c("rsid","a1","a0","beta1","se","p","par","band"))

stats <- rbind(stats_bothpl,stats_fath,stats_moth)
stats <- stats[order(rsid)]

# Align alleles
stats <- loci[,.(rsid,a1,a0)][stats,,on="rsid"]
stats[i.a1==a0,c("i.a1","i.a0","beta1"):=.(i.a0,i.a1,-beta1)]
stats[,c("i.a1","i.a0"):=NULL]

stats[,ratio:=beta1/.SD[band=="e_",beta1] - 1,by=c("rsid","par")] # Calculate ratios
stats[,se_ratio:=sqrt((se/.SD[band=="e_",beta1])^2 + (beta1^2 * .SD[band=="e_",se]^2)/(.SD[band=="e_",beta1]^2)^2),by=c("rsid","par")]

stats <- stats[!band=="e_",]
stats[,c("beta1","se","ratio","se_ratio"):=.(ratio, se_ratio, NULL, NULL)]



#----
# Get Median Lifespan
#----

pheno <- fread("../q09_age_stratified_gwas/p02_qc_phenotypes/st02_01_qcd_phenotypes.txt")
pheno <- pheno[!is.na(g_bri) & !is.na(unrelated) & !is.na(is_not_withdrawn_20181016),]

fath <- pheno[,.(iid=paste0(iid,"FATH"), age=fath_age, dead=fath_dead, par="fath")]
moth <- pheno[,.(iid=paste0(iid,"MOTH"), age=moth_age, dead=moth_dead, par="moth")]

ph <- rbind(fath, moth)
ph[,band:=ifelse(age > 80, "80_120", ifelse(age > 60, "60_80", ifelse(age > 40, "40_60", NA)))]

m0 <- survfit(Surv(age, dead) ~ par + band, data=ph)
quantile(m0,probs=0.5, conf.int=T)


s0 <- summary(m0)
t0 <- cbind(data.table(s0$table), par=c(rep("fath",3),rep("moth",3)), band=rep(c("40_60","60_80","80_120"),2))

medians <- t0[,.(par,band, median)]



#----
# Plot
#----

# Configure dataset
stats1 <- stats[medians,,on=c("par","band")]

plot_df <- loci[stats1, , on="rsid"]
plot_df[i.a1==a0, c("beta1", "i.a1", "i.a0"):=.(-beta1, i.a0, i.a1)]
plot_df[,label:=paste0(gene, "\n(",rsid,"_",i.a1,")")]

# Calculate regression by locus
newx <- seq(40,100,1)

predict_df <- data.table(	
	plot_df %>% group_by(label) %>% do({
		m1 <- lm(beta1 ~ median, weight=1/se^2, data=.)	
		p1 <- predict(m1,newdata=list(median=newx), se.fit=T)
		beta <- summary(m1)$coefficients[2,1]
		se <- summary(m1)$coefficients[2,2] 
		z <- summary(m1)$coefficients[2,3] 
		p <- summary(m1)$coefficients[2,4]
		data.table(median=newx, beta1=p1$fit, se=p1$se.fit, mbeta=beta, mse=se, mz=z, mp=p, label=unique(.$label))
		})
	)

pvalues <- predict_df[,.(beta=unique(mbeta), se=unique(mse), z=unique(mz), p=unique(mp)),by="label"][order(p),]
pvalues[,label:=factor(label, levels=label)]
predict_df[,c("mbeta","mse","mz","mp"):=NULL]


plot_df[,label:=factor(label, levels=pvalues[,label])]
predict_df[,label:=factor(label, levels=pvalues[,label])]


# Create dataframe for pretty axes
panels <- levels(plot_df[,label])
lims <- data.table(x=rep(c(40,100),length(panels)), y=rep(c(-4,4, rep(NA,(length(panels)-2))),2), label=factor(rep(panels,each=2), levels=panels))



# Draw actual plot
pdf("st10_06_all_age_bands_rel.pdf", width=8, height=4)
p1 <- ggplot(plot_df, aes(x=median, y=beta1, colour=par)) +
geom_hline(yintercept=seq(-4, 4, 2), colour="grey90") +
geom_hline(yintercept=0, colour="grey80") +
geom_vline(xintercept=seq(40,100, 20), colour="grey90") +
geom_ribbon(data=predict_df, aes(x=median, ymin=beta1-qnorm(0.975)*se, ymax=beta1+qnorm(0.975)*se), alpha=0.1, inherit.aes=FALSE) +
geom_linerange(data=plot_df, aes(ymin=beta1 -qnorm(0.975) * se, ymax=beta1 + qnorm(0.975) * se), alpha=0.5) +
geom_point() + 
geom_text(data=pvalues, aes(x=70, y=Inf, label=ifelse(p < 0.001, sprintf("P[age]: %.0e",p), sprintf("P[age]: %.3f",round(p,4)))), hjust=0.5, vjust=1.5, size=3, colour="black", parse=T) +
facet_wrap(~label, nrow=2, scales="free_x") + geom_rangeframe(data=lims, aes(x=x,y=y), colour="black") +
scale_y_continuous(breaks=seq(-10,10,2)) +
scale_x_continuous(breaks=seq(40,100,20),lim=c(40,100), expand=c(0,0)) +
scale_colour_manual(values=c("#2166ac","#b2182b"), name="", breaks=c("moth","fath"), labels=c("Mothers","Fathers")) +
theme_bw() + labs(x="Median survival of parental age band", y=c("Relative change in effect size")) +
coord_cartesian(ylim=c(-5,5)) +
theme(legend.position="right", panel.grid.minor=element_blank(), panel.border=element_blank(), panel.spacing = unit(1.5, "lines"), panel.grid=element_blank())
print(p1)
dev.off()

ggsave("st10_06_all_age_bands_rel.svg", p1, width=8, height=4, dpi='retina')

# Overall trends

m1 <- metagen(beta,se,studlab=gsub("\n.*$","",label), data=pvalues)
print(m1)
i1 <- metainf(m1, pooled="fixed")
print(i1)

t1 <- data.table(mbeta=i1$TE[length(i1$TE)], mse=i1$seTE[length(i1$seTE)])
t1[,c("mz","mp"):=.(mbeta/mse, 2*pnorm(-abs(mbeta/mse)))]

loci[,label:=paste0(gene, "\n(",rsid,"_",a1,")")]

p_table <- loci[pvalues,.(gene,rsid,a1,a0,beta,se,z,p),on="label"]
pt_table <- rbind(p_table,t1[,.(gene="<TOTAL>",rsid=NA,a1=NA,a0=NA,beta=mbeta,se=mse,z=mz,p=mp)])


# APOE sensitivity

t2 <- pvalues[!grepl("APOE",label),.(mbeta=sum(beta/se^2)/sum(1/se^2), mse=sqrt(1/sum(1/se^2)))]
t2[,c("mz","mp"):=.(mbeta/mse, 2*pnorm(-abs(mbeta/mse)))]


# Remaining loci

t3 <- pvalues[!grepl("APOE|FOXO|CDKN2B",label),.(mbeta=sum(beta/se^2)/sum(1/se^2), mse=sqrt(1/sum(1/se^2)))]
t3[,c("mz","mp"):=.(mbeta/mse, 2*pnorm(-abs(mbeta/mse)))]


pt_table <- rbind(pt_table, 
	t2[,.(gene="<TOTAL w/o APOE>",rsid=NA,a1=NA,a0=NA,beta=mbeta,se=mse,z=mz,p=mp)],
	t3[,.(gene="<TOTAL w/o APOE, FOXO3, CDKN2B-AS1>",rsid=NA,a1=NA,a0=NA,beta=mbeta,se=mse,z=mz,p=mp)])


# Export

fwrite(format(pt_table,digits=6), "st10_06_all_age_bands_rel_meta.tsv", sep="\t", na="NA", quote=F)
