#!/usr/bin/env R

#----
# Setup environment
#----

options(width=200)
library(data.table)
library(ggplot2)

#----
# Load data
#----

phen_cors <- fread("st02_pheno_correlations.tsv")


cors <- fread("st02_genetic_correlations.tsv", select=1:6)
cors[,c("p1","p2"):=.(gsub(".sumstats.gz","",p1), gsub(".sumstats.gz","",p2))]
cors[,c("p1","p2"):=.(gsub("^.*_","",p1), gsub("^.*_","",p2))]
cors[,c("p1","p2"):=.(gsub("pct","longevity",p1), gsub("pct","longevity",p2))]
#cors[,c("p1","p2"):=.(gsub("lif","bothpl_life_",p1), gsub("lif","bothpl_life_",p2))]
cors[,c("p1","p2"):=.(factor(p1, levels=unique(p1)), factor(p2, levels=rev(unique(p1))))]

cors[,value:=0]
cors[p1 == "healthspan" | p2 == "healthspan", value:=value+100000]
cors[p1 == "lif4060" | p2 == "lif4060", value:=value+10000]
cors[p1 == "lif6080" | p2 == "lif6080", value:=value+1000]
cors[p1 == "lifegen" | p2 == "lifegen", value:=value+100]
cors[p1 == "lif80120" | p2 == "lif80120", value:=value+10]
cors[p1 == "longevity" | p2 == "longevity", value:=value+11]

cors1 <- cors[!duplicated(value),]
#cors1 <- rbind(cors1, data.table(p1=cors1[,unique(p1)],p2=cors1[,unique(p1)],rg=1,se=0,z=NA,p=NA,value=NA))

pdf("st03_genetic_correlations.1.pdf", height=6, width=8)
ggplot(cors1, aes(x=p1, y=p2, fill=rg)) +
 geom_tile(width=0.95, height=0.95) + 
 geom_text(aes(label=round(rg,2)), nudge_y= 0.118/2, col="white") +
 geom_text(aes(label=paste0("(",round(rg-qnorm(0.975)*se,2),"-",round(rg+qnorm(0.975)*se,2),")")), size=3, nudge_y=-0.118/2, col="white") +
 scale_fill_gradient(low="#6baed6", high="#08306b") +
 scale_x_discrete(position="top") +
 coord_equal() + labs(x="",y="") +
 theme_minimal() + theme(panel.grid=element_blank())
dev.off()

cors2 <- cors[p2 %in% c("healthspan","longevity"),]
cors2 <- rbind(cors2, data.table(p1=c("healthspan","longevity"), p2=c("healthspan","longevity"), rg=1, se=0, z=NA, p=NA, value=NA))
cors2[,c("p1","p2"):=.(factor(p1, levels=unique(p1)), factor(p2, levels=unique(p1)))]

pdf("st03_genetic_correlations.2.pdf", height=4, width=6)
ggplot(cors2, aes(x=p1, y=rg, fill=p1)) +
 geom_errorbar(aes(ymin=rg-qnorm(0.975) * se, ymax=rg+qnorm(0.975)*se), width=0.5) +
 geom_bar(stat="identity", col="black") + facet_wrap(~p2) +
 scale_fill_brewer(type = "seq", palette="Greys") + guides(fill=FALSE) +
 labs(x="",y="Genetic Correlation") + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
dev.off()

