#!/usr/bin/env Rscript

#----
# Setup environment
#----

set.seed(1)
options(width=200)
library(data.table)
library(xlsx)
library(pheatmap)


#----
# Get data
#----

df <- fread("st04_04_gws_hits_diseases_summary2.tsv")

table2 <-  data.table(read.xlsx("../q03_analyse/st03_06_triple_manova_table.xlsx", sheetIndex=4, startRow=3))[,1:2]
names(table2) <- c("gene","rsid")
table2 <- table2[!is.na(gene),]
table2[is.na(rsid),"gene"] <- ""


df <- df[match(table2$gene, gene)]

pdf <- as.matrix(df[,c(-1,-(ncol(df)-1),-ncol(df)),with=F])
rownames(pdf) <- df[,gene]

rownames(pdf)[is.na(rownames(pdf))][1] <- " "

pdf <- pdf[order(nrow(pdf):1),]
pdf[is.na(pdf)] <- -1

pdf1 <- pdf
pdf[pdf > 10] <- 11


annot <- data.frame(Phenotype = rev(c(rep("Novel",5)," ",rep("Known",5))))
rownames(annot) <- rownames(pdf)
annot_colors <- list(Phenotype = c(Novel="#404040", Known="#008837",` `="#FFFFFF"))

pheatmap(pdf, color=c("#FFFFFF",colorRampPalette(c("#f7fbff","#08306b"))(11),"#000038"),
 border_color="white", cluster_rows=FALSE, cluster_cols=FALSE, annotation_row = annot, annotation_colors = annot_colors, filename="st004_05_gws_hits_diseases.1.pdf", annotation_legend=FALSE, width=10, height=4)

