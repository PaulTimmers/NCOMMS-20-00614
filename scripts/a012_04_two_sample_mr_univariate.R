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
library(scales)
library(ggrepel)

haem_traits <- commandArgs(T)


# Pretty names

labels <- c(iron="Serum iron", transferrin="Transferrin", saturation="Transferrin saturation", ferritin="Ferritin")


#----
# Multivariable MR
#----


# Load manova stats

outcome_dat <- read_outcome_data(
	filename=paste0("st012_04_manova_all_univariate_hits.tsv"),
	sep="\t",
	snp_col="rsid",
	beta_col="beta1",
	se_col="se",
	effect_allele_col="a1",
	other_allele_col="a0",
	eaf_col="freq1",
	pval_col="p",
	samplesize_col="n")

outcome_dat$outcome <- "healthspan/lifespan/longevity"


# Run 

plots <- list()
out <- data.table()
technical_res <- data.table()

for (haem_trait in haem_traits) {

cat("\n\n\n---------------------",haem_trait,"---------------------\n\n")

# Load exposure stats

exposure_dat <- read_exposure_data(
	filename=paste0("st012_03_",haem_trait,"_hits.tsv"),
	sep="\t",
	snp_col="rsid",
	beta_col="beta1",
	se_col="se",
	effect_allele_col="a1",
	other_allele_col="a0",
	eaf_col="freq1",
	pval_col="p",
	samplesize_col="n")

exposure_dat$exposure <- haem_trait


# Harmonise data

dat <- harmonise_data(
    exposure_dat = exposure_dat, 
    outcome_dat = outcome_dat
)


# Perform MR


res <- mr(dat)
res_loo <- mr_leaveoneout(dat)


temp <- mr_egger_regression(dat$beta.exposure, 
                  dat$beta.outcome, dat$se.exposure, dat$se.outcome, 
                  default_parameters())

temp$dat$SNP <- dat$SNP
loci <- fread(cmd=paste0("snp2gene.sh -s ",paste0(dat$SNP,collapse=" "))) # Get locus names
temp$dat$label <- gsub("/.*","",loci$summary)
technical_res <- rbind(technical_res, data.table(res_loo)[,.(exposure, excluding_snp=SNP, snp_label=loci[match(res_loo$SNP, rsid), summary], beta=b, se, p)])



plot_res <- cbind(res, data.table(a=ifelse(res$method == "MR Egger", temp$b_i, 0)))
plot_res <- plot_res[plot_res$method %in% c("MR Egger", "Inverse variance weighted"),]


lims <- data.table(temp$dat)[,.(b_exp=pretty(c(b_exp-se_exp, 0, b_exp+se_exp)), b_out=pretty(c(b_out-se_out, 0, b_out+se_out)))]
lines <- data.table(plot_res)[,.(b_exp=lims$b_exp, b_out=a + b * lims$b_exp), by="method"]


plot <- ggplot(data = temp$dat, aes(x = b_exp, y = b_out)) +
geom_rangeframe(data=lims) +
geom_smooth(method="lm", formula= y ~ x, aes(weight=1/se_out^2), colour=NA, fill="#a6cee3", fullrange=TRUE) +
geom_smooth(method="lm", formula= y ~ 0 + x, aes(weight=1/se_out^2), colour=NA, fill="#1f78b4", fullrange=TRUE) +
geom_line(data=lines, aes(colour=method), show.legend=TRUE) +
geom_errorbar(aes(ymin = b_out - se_out, ymax = b_out + se_out), colour = "grey20", width = 0) +
geom_errorbarh(aes(xmin = b_exp - se_exp, xmax = b_exp + se_exp), colour = "grey20", height = 0) + 
geom_point() +
geom_text_repel(aes(label=label), size=2, fontface="italic", point.padding=0.15, segment.color="black") +
scale_y_continuous(breaks=lims$b_out) +
scale_x_continuous(breaks=lims$b_exp) +
coord_cartesian(xlim=c(min(0,min(lims$b_exp)),max(0,max(lims$b_exp))), ylim=c(min(0,min(lims$b_out)),max(0,max(lims$b_out)))) +
labs(colour = "MR Test", x = labels[dat$exposure[1]], y = "Beta") + 
scale_colour_manual(values = c("#1f78b4","#a6cee3")) + 
theme_minimal() +
theme(legend.position = "top", legend.direction = "vertical", panel.grid=element_blank(), axis.ticks=element_line()) + 
guides(colour = guide_legend(ncol = 2))



plots[[haem_trait]] <- plot

fwrite(res, paste0("z012_mr_results/st012_01_",haem_trait,"_univariate_mr.tsv"), sep="\t", na="NA", quote=F)



res <- data.table(res)
out <- rbind(out, data.table(res[method=="Inverse variance weighted",.(exposure, nsnp, beta=b, se, p=pval)], p_egger=res[method=="MR Egger",pval]))
}


technical_res[excluding_snp=="All",excluding_snp:="None"]
fwrite(format(technical_res,digits=4), "z012_mr_results/st012_02_univariate_mr_tech_stats.tsv", sep="\t", na="NA", quote=F)

mout <- melt(out, id.var=c("exposure","nsnp","beta","se"), value.name="p")
mout[,fdr:=p.adjust(p,"fdr")]
dmout <- dcast(mout, exposure + nsnp + beta + se ~ variable, value.var=c("p","fdr"))
out <- dmout[,.(exposure,nsnp,beta,se,p=p_p, fdr=fdr_p, p_egger=p_p_egger, fdr_egger=fdr_p_egger)]

out <- out[order(fdr),]
fwrite(format(out,digits=4), "z012_mr_results/st012_02_univariate_mr_stats.tsv", sep="\t", na="NA", quote=F)
print(out)

g1 <- ggarrange(plotlist=plots, labels="", common.legend=T)
ggsave("z012_mr_results/st012_02_univariate_mr_plot.svg", g1, width=8, height=6)
ggsave("z012_mr_results/st012_02_univariate_mr_plot.pdf", g1, width=8, height=6)
