#!/usr/bin/env R

#----
# Setup environment
#----

options(width=200)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(survival)

#----
# Basic three-way correlation
#----

cors <- fread("st005_05_genetic_correlations.tsv", select=1:6)
cors[,c("p1","p2"):=.(gsub(".sumstats.gz","",p1), gsub(".sumstats.gz","",p2))]
cors[,c("p1","p2"):=.(gsub("^.*_","",p1), gsub("^.*_","",p2))]
cors[,c("p1","p2"):=.(gsub("pct","longevity",p1), gsub("pct","longevity",p2))]
cors[,c("p1","p2"):=.(factor(p1, levels=unique(p1)), factor(p2, levels=rev(unique(p1))))]

cors[,value:=0]
cors[p1 == "healthspan" | p2 == "healthspan", value:=value+100]
cors[p1 == "lifegen" | p2 == "lifegen", value:=value+10]
cors[p1 == "longevity" | p2 == "longevity", value:=value+1]

cors1 <- cors[!duplicated(value),]

cors2 <- cors1
cors2[p1=="lifegen",p1:="Parental Lifespan"]
cors2[p2=="lifegen",p2:="Parental Lifespan"]
cors2[,p1:=paste0(toupper(substr(p1,0,1)),substr(p1,2,99))]
cors2[,p2:=paste0(toupper(substr(p2,0,1)),substr(p2,2,99))]
cors2[,label:=paste0(p1, "\nx\n", p2)]
cors2[,label:=factor(label, unique(label))]


lims <- data.table(x=cors2[,label], y=c(0,0.5,1))

pdf("st005_09_genetic_correlations.1.pdf", height=3, width=4)
g1 <- ggplot(cors2, aes(x=label, y=rg)) +
geom_hline(yintercept=seq(0, 1, 0.2), colour="grey90") +
geom_errorbar(aes(ymin=rg-qnorm(0.975) * se, ymax=rg+qnorm(0.975)*se), width=0.33) +
geom_bar(stat="identity", col="black", fill="grey", width=0.66) +
geom_rangeframe(data=lims, aes(x=x,y=y), colour="black") +
scale_y_continuous(breaks=seq(0,1,0.2), expand=c(0,0.15)) +
labs(x="",y="Genetic Correlation") + theme_bw() +
theme(legend.position="right", panel.grid.minor=element_blank(), panel.spacing = unit(1.5, "lines"), panel.border=element_blank(),panel.grid=element_blank(),axis.title.y=element_text(size=10))
print(g1)
dev.off()





#----
# Age-stratified correlation
#----

cors <- fread("st005_08_age_genetic_correlations.tsv", select=1:6)
cors[,c("p1","p2"):=.(gsub(".sumstats.gz","",p1), gsub(".sumstats.gz","",p2))]
cors[,c("p1","p2"):=.(gsub("^.*_","",p1), gsub("^.*_","",p2))]
cors[,c("p1","p2"):=.(gsub("pct","longevity",p1), gsub("pct","longevity",p2))]
cors[,c("p1","p2"):=.(factor(p1, levels=unique(p1)), factor(p2, levels=rev(unique(p1))))]

cors[,value:=0]
cors[p1 == "healthspan" | p2 == "healthspan", value:=value+100000]
cors[p1 == "lif4060" | p2 == "lif4060", value:=value+10000]
cors[p1 == "lif6080" | p2 == "lif6080", value:=value+1000]
cors[p1 == "lif80120" | p2 == "lif80120", value:=value+100]
cors[p1 == "longevity" | p2 == "longevity", value:=value+1]

cors1 <- cors[!duplicated(value),]

cors2 <- cors[!p1 %in% c("healthspan","longevity") & p2 %in% c("healthspan","longevity"),]
cors2[p1=="lif4060",p1:="40-60"]
cors2[p1=="lif6080",p1:="60-80\n\nParental Lifespan (UKBB)"]
cors2[p1=="lif80120",p1:="80+"]
cors2[,p1:=paste0(toupper(substr(p1,0,1)),substr(p1,2,99))]
cors2[,p2:=paste0(toupper(substr(p2,0,1)),substr(p2,2,99))]
cors2[,c("p1","p2"):=.(factor(p1, levels=unique(p1)), factor(p2, levels=unique(p2)))]


panels <- c("Healthspan","Longevity")
lims <- data.table(x=rep(c(1,length(cors2[,unique(p1)])),length(panels)), y=c(0,1, rep(NA,(length(panels)-1)*2)), p1=panels, p2=factor(rep(panels,each=2), levels=panels))


pdf("st005_09_genetic_correlations.2.pdf", height=3, width=5)
g2 <- ggplot(cors2, aes(x=p1, y=rg, fill=p2)) +
 geom_hline(yintercept=seq(0, 1, 0.2), colour="grey90") +
 geom_linerange(aes(x=p1, ymin=0, ymax=1.2), colour="grey90") +
 geom_errorbar(aes(ymin=rg-qnorm(0.975) * se, ymax=rg+qnorm(0.975)*se), width=0.5) +
 geom_bar(stat="identity", col="black") + facet_wrap(~p2) +
 geom_rangeframe(data=lims, aes(x=x,y=y), colour="black") +
 scale_fill_manual(values=c("white", "black")) +
 scale_y_continuous(breaks=seq(0,1,0.2), expand=c(0,0.05)) +
 guides(fill=FALSE) +
 labs(x="",y="Genetic Correlation") + theme_bw() +
 theme(legend.position="right", panel.grid.minor=element_blank(), panel.spacing = unit(1.5, "lines"), panel.border=element_blank(),panel.grid=element_blank(),axis.text.x=element_text(size=10), axis.title.y=element_text(size=10))
print(g2)
dev.off()


pdf("st005_09_genetic_correlations.pdf", height=6, width=5)
ggarrange(g1, g2, ncol=1, nrow=2, labels=c("(a)","(b)"))
dev.off()
