#!/usr/bin/env Rscript

#----
# Setup environment
#----

options(width=200)
library(data.table)
library(ComplexHeatmap)
library(UpSetR)

#----
# START
#----

table <- fread("st03_06_triple_manova_table.csv")

health_loci <- table[p_health < 0.05, gene]
life_loci <- table[p_life < 0.05, gene]
long_loci <- table[p_long < 0.05, gene]

listInput <- list(Healthspan = health_loci, Lifespan = life_loci, Longevity = long_loci)
matrixInput <- list_to_matrix(listInput)
m <- make_comb_mat(matrixInput)


pdf("st03_08_upset_plot.pdf", height=8, width=8)
UpSet(m, pt_size = unit(5, "mm"), lwd = 3, 
	top_annotation = HeatmapAnnotation("Number of shared loci" = anno_barplot(comb_size(m), 
            border = FALSE, 
            gp = gpar(fill = "black"), 
            height = unit(15, "cm")),
		annotation_name_side = "left", 
        annotation_name_rot = 90),
	right_annotation = rowAnnotation("Locus P < 0.05" = anno_barplot(set_size(m),
		border = FALSE, 
        gp = gpar(fill = "white", col="black"), 
        width = unit(5, "cm"),
        ylim=c(0,24))
        ))
dev.off()
