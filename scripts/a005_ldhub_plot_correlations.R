#!/usr/bin/env Rscript

#----
# Setup environment
#----

options(width=200)
library(data.table)
library(ggplot2)
library(meta)
library(dplyr)
library(yaml)

cmyk = function(c, m, y, k=0, alpha, maxColorValue=1){
  if(maxColorValue != 1) { c <- c/maxColorValue; m <- m/maxColorValue; y <- y/maxColorValue }
  c <- (c * (1-k) + k)
  m <- (m * (1-k) + k)  
  y <- (y * (1-k) + k)

  c <- 1-c; m <- 1-m; y <- 1-y
  hex <- function(v) substring(rgb(v,0,0),2,3)
  if(!missing(alpha)) alpha <- hex(alpha) else alpha <- ''
  paste0('#',hex(c), hex(m), hex(y), alpha)
}


#----
# Load data
#----

rg <- fread("t008_rg.tsv")
rg_list <- fread("t006_traits.tsv", col.names=c("trait","pmid","category"))


rg1 <- rg[trait2 %in% c("healthspan", "lifespan", "longevity") & !trait1 %in% c("healthspan", "lifespan", "longevity"),]

#----
# Descriptives
#----

heterogeneity <- data.table(
rg1 %>% 
group_by(trait1) %>% 
do({
	m1 <- metagen(TE=rg, seTE=se, data=., studlab=trait2)
	i1 <- metainf(m1)
  i2 <- i1$I2[1:3]; names(i2) <- paste0("iomit_",m1$studlab)
  data.table(p_het=m1$pval.Q, t(i2))
	})
)

rg1 <- heterogeneity[rg1,,on="trait1"]
rg1[,q:=p.adjust(p,"fdr")]

# Signifant in all, no heterogeneity
print(rg1[p < 0.05/.N & p_het > 0.05,][,.(N=length(unique(trait2)),rg_min=min(rg),se_min=.SD[which.min(rg),se],rg_max=max(rg),se_max=.SD[which.max(rg),se],trait_min=.SD[which.min(rg),trait2],trait_max=.SD[which.max(rg),trait2]),by="trait1"][N == 3][order(pmin(-abs(rg_max),-abs(rg_min)),-N),])



#----
# Create tables
#----

table <- rg1[,.(trait1=trait2,trait2=trait1,rg,se,p,q,p_het,h2,h2_se,pmid)]

t1 <- dcast(table, trait2 + h2 + h2_se + p_het + pmid ~ trait1, value.var=c("rg","se","p","q"))
t2 <- t1[,.(trait=trait2, pmid, rg_healthspan, se_healthspan, p_healthspan, q_healthspan, rg_lifespan, se_lifespan, p_lifespan, q_lifespan, rg_longevity, se_longevity, p_longevity, q_longevity, p_het)]
t2 <- t2[order(p_het)]

fwrite(t2, "t008_rg_results_table.csv", sep=",", quote=F, na="NA")


tech_table <- rg1[,.(trait=trait1,trait2,pmid,h2,h2_se,n_snp)]
tech2 <- dcast(tech_table, trait + pmid + h2 + h2_se ~ trait2, value.var="n_snp",)
tech3 <- tech2[,.(trait, ancestry="European", pmid, `h2 (se)`=sprintf("%.3f (%.3f)",h2,h2_se), snp_health=healthspan, snp_life=lifespan, snp_long=longevity)]

fwrite(tech3, "t008_rg_technical_table.csv", sep=",", quote=F, na="NA")


#----
# Corrplot
#----

library(ggdendro)
library(corrplot)
library(corpcor)
library(RColorBrewer)


#----
# Full correlations
#----


q_thresh <- 0.01
sig_traits <- rg1[q < q_thresh,trait1]

rg2 <- dcast(rg[trait1 %in% sig_traits, trait2 %in% sig_traits,.(trait1,trait2,z=rg)], trait1 ~ trait2, value.var="z")
mat1 <- as.matrix(rg2[,-1])
rownames(mat1) <- rg2[,trait1]

clusters <- hclust(dist(mat1[,!colnames(mat1) %in% c("healthspan","lifespan","longevity")]))

dend_data <- dendro_data(clusters, type = "rectangle")
ord <- clusters$order

mat2 <- mat1[,c("healthspan", "lifespan", "longevity")]
mat2 <- mat2[ord,]


# Significance
p.mat <- as.matrix(t2[trait %in% sig_traits,.(healthspan=q_healthspan, lifespan=q_lifespan, longevity=q_longevity)])
rownames(p.mat) <- t2[trait %in% sig_traits,trait]
p.mat <- p.mat[match(rownames(mat2),rownames(p.mat)),]


# Pretty names
old <- c("age_at_menarche", "alzheimer", "bmi", "cancer_breast", "cancer_genital_male", "cancer_melanoma", "cancer_respiratory", "copd", "coronary_artery", "depression", "drinking", "hernia", "risk_taking", "smoking", "stroke", "type_2_diabetes", "years_of_schooling")
new <- c("Age at menarche", "Alzheimer's disease", "Body mass index", "Breast cancer", "Prostate/Testicular cancer", "Melanoma", "Lung cancers", "COPD", "Coronary artery disease", "Depression", "Alcohol intake", "Hernia", "Risk taking behaviour", "Ever smoked", "Stroke", "Type 2 diabetes", "Years of schooling")
rownames(mat2) <- new[match(rownames(mat2), old)]
colnames(mat2) <- c("Healthspan", "Lifespan", "Longevity")


# Bold
bold <- function(x) paste0(':bold("',x,'")')
het_traits <- rownames(mat2)[old[match(rownames(mat2), new)] %in% rg1[p_het < 0.05,unique(trait1)]] 
rownames(mat2)[rownames(mat2) %in% het_traits] <- bold(rownames(mat2)[rownames(mat2) %in% het_traits])



pdf("t008_rg_results_corrplot.pdf", width=12, height=8)
corrplot(mat2, method="circle", cl.ratio = 0.1, cl.pos="b", cl.length=3, cl.cex=1.1, tl.col = "black", p.mat = p.mat, insig="blank",
  sig.level = q_thresh, tl.srt=45, mar=unit(c(0,0,1,8), units="in"))
dev.off()


pdf("t008_rg_results_dendrogram.pdf", width=6, height=8)
ggdendrogram(dend_data, rotate=T)
dev.off()


#----
# Publication
#----

setEPS()
postscript("t008_rg_results_corrplot.eps", width=12, height=8)
corrplot(mat2, method="circle", cl.ratio = 0.1, cl.pos="b", cl.length=3, cl.cex=1.1, tl.col = "black", p.mat = p.mat, insig="blank",
  sig.level = q_thresh, tl.srt=45, mar=unit(c(0,0,1,8), units="in"))
dev.off()

postscript("t008_rg_results_dendrogram.eps", width=6, height=8)
ggdendrogram(dend_data, rotate=T)
dev.off()