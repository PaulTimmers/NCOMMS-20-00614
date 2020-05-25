#!/usr/bin/env Rscript

#$ -l h_vmem=32G

#----
# Setup environment
#----

options(width=200)
library(data.table)
library(TwoSampleMR)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(MendelianRandomization)

haem_traits <- commandArgs(T)

# Pretty names

labels <- c(iron="Serum iron", transferrin="Transferrin", saturation="Transferrin saturation", ferritin="Ferritin")
units <- c(iron="µmol/L", transferrin="g/L", saturation="percentage", ferritin="log µg/L")

#----
# Multivariable MR
#----

all_res <- data.table()

for (outcome in c("manova","healthspan","lifespan","longevity")) {


writeLines(paste0("\n\n\n=========================== ", outcome," ===========================\n"))
# Load outcome stats

sink("/dev/null")
outcome_dat <- read_outcome_data(
	filename=paste0("st012_06_manova_all_hits.tsv"),
	sep="\t",
	snp_col="rsid",
	beta_col=ifelse(outcome=="manova", "beta1", paste0("beta1_",outcome)),
	se_col=ifelse(outcome=="manova", "se", paste0("se_",outcome)),
	effect_allele_col="a1",
	other_allele_col="a0",
	eaf_col="freq1",
	pval_col="p",
	samplesize_col="n")

outcome_dat$outcome <- outcome
sink()


plots <- list()
exposure_dat <- data.table()

for (haem_trait in haem_traits) {

# Load exposure stats

sink("/dev/null")
exposure_dat1 <- read_exposure_data(
	filename=paste0("st012_06_",haem_trait,"_all_hits.tsv"),
	sep="\t",
	snp_col="rsid",
	beta_col="beta1",
	se_col="se",
	effect_allele_col="a1",
	other_allele_col="a0",
	eaf_col="freq1",
	pval_col="p",
	samplesize_col="n")

exposure_dat1$exposure <- haem_trait
sink()

exposure_dat <- rbind(exposure_dat, exposure_dat1)
}



#----
# Harmonise data
#----


# Get LD correlations

x <- harmonise_data(exposure_dat1, outcome_dat)
ld <- ld_matrix(unique(x$SNP))
out <- TwoSampleMR:::harmonise_ld_dat(x, ld)


# Harmonize multivariate data

mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)

dat_to_MVMRInput <- function(mvdat, correlation=matrix()) {
	snps <- rownames(mvdat$exposure_beta)
	MRMVInputObject <- MendelianRandomization::mr_mvinput(
		snps=snps,
		exposure=mvdat$expname[match(colnames(mvdat$exposure_beta), mvdat$expname$id.exposure),exposure],
		outcome=mvdat$outname$outcome,
		bx = mvdat$exposure_beta,
		bxse = mvdat$exposure_se,
		by = mvdat$outcome_beta,
		byse = mvdat$outcome_se,
		correlation=correlation)
	
	return(MRMVInputObject)
}


mvdat2 <- dat_to_MVMRInput(mvdat, correlation=out$ld)


# Perform MR


if(outcome=="manova") { 
	
	loci <- fread(cmd=paste0("snp2gene.sh -s ",paste0(mvdat2@snps,collapse=" ")))
	loci[,summary:=gsub("/.*","",summary)]
	res_loo <- data.table()

	for(snp in c("None",mvdat2@snps)) {
		x_loo <- mv_harmonise_data(exposure_dat=exposure_dat[exposure_dat$SNP != snp,], outcome_dat=outcome_dat[outcome_dat$SNP != snp,])
		x2_loo <- dat_to_MVMRInput(x_loo, correlation=out$ld[rownames(out$ld)!=snp,colnames(out$ld)!=snp])
		m2_loo <- mr_mvivw(x2_loo, model = "default", correl = TRUE, distribution = "normal", alpha = 0.05)
		res1_loo <- data.table(exposure=m2_loo@Exposure, excluding_snp=snp, snp_label=loci[rsid==snp, summary], beta=m2_loo@Estimate, se=m2_loo@StdError, p=m2_loo@Pvalue)
		res_loo <- rbind(res_loo, res1_loo, fill=T)
	}

	fwrite(format(res_loo,digits=4), "z012_mr_results/st012_03_multivariate_mr_tech_stats.tsv", sep="\t", na="NA", quote=F)

}

m1 <- mr_mvivw(mvdat2, model = "default", correl = TRUE, distribution = "normal", alpha = 0.05)
print(m1)




m1_egger <- mr_mvegger(mvdat2, correl = TRUE, distribution = "normal", alpha = 0.05)
print(m1_egger)

res1 <- data.table(exposure=m1@Exposure, outcome=m1@Outcome, nsnp=m1@SNPs, beta=m1@Estimate, se=m1@StdError, p=m1@Pvalue, p_intercept=m1_egger@Pvalue.Int)
res1[,c("fdr","fdr_intercept"):=.(p.adjust(p,"fdr",n=8),p.adjust(p_intercept,"fdr",n=8))]
print(res1)
print(res1)

all_res <- rbind(all_res, res1)

}


t1 <- all_res[,.(exposure,outcome,nsnp,beta,se,p,fdr,p_intercept,fdr_intercept)]

manova_stats <- t1[outcome=="manova",]
trait_betas <- dcast(t1[outcome!="manova",], exposure ~ outcome, value.var=c("beta","se"))
t2 <- manova_stats[trait_betas,,on="exposure"][order(p)]


bse <- function(beta, digits=2, se=NA, se_digits=NA) {
	se_digits <- ifelse(is.na(se_digits), digits, se_digits)
	ifelse(is.na(se), return(format(round(beta,digits), nsmall=digits)), return(paste0(format(round(beta,digits),nsmall=digits), " (",format(round(se,digits), nsmall=se_digits),")")))
}

t3 <- t2[,.(exposure=labels[exposure], nsnp, beta, se, p, fdr, 
	healthspan=bse(beta_healthspan,2,se_healthspan,2),
	lifespan=bse(beta_lifespan,2,se_lifespan,2),
	longevity=bse(beta_longevity,2,se_longevity,2))]

writeLines("\n\n")
print(t3)
fwrite(t3, "z012_mr_results/st012_03_multivariate_mr_stats.tsv", sep="\t", na="NA", quote=F)
