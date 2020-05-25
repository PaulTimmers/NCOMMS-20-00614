#!/usr/bin/env R

#----
# Setup environment
#----

library(ggplot2)
library(ggpubr)
library("yaml")
library(metafor)
library(data.table)
source("/opt/working/wilson/apps_by_us/gen_useful.R")
#R --args args ../impv1.sample /opt/working/wilson/projects/prj_022_longev/ukbb5/p21_form_res/st20_01_long_par_plink.phe t4 "fath_life_res" "fath_life_res_rn" "moth_life_res" "moth_life_res_rn"

if(!dir.exists("all")) dir.create("all")

options(width=180)

#----
# START
#----


table2_file<- "../../p010_table2/st010_01_table2.txt"
table2<-read_pj(table2_file)
#table2<-table2[,c(1:8,10,16,22,28,ncol(table2))]
table2 <- table2[,names(table2)[!grepl("bothor",names(table2))]]
names(table2)[names(table2)=="source"] <- "cohort"
names(table2) <- gsub("_alldr_bothpl","",names(table2))

#cand_snps <- read.table("../../r001_ukbb_pilling/p99_actual_sumstats/st99_03_candidates_dataframe.tsv",sep="\t",header=T)
#cand_snps$cohort <- "other"
#cand_snps$snpid <- paste(cand_snps$chr,cand_snps$pos,sep="_")
#names(cand_snps) <- sub("_replication","",names(cand_snps))
#cand_snps <- cand_snps[,names(cand_snps)%in%names(table2)]

#table2 <- rbind(table2,cand_snps)
table2<-table2[!duplicated(table2$gene) & !is.na(table2$beta1),]


phwas_snps <- read_pj("../../p020_extract_cands/st020_04_phewas_candidates.tsv")
#phwas_snps <- phwas_snps[phwas_snps$source!="other",]


fdr_disease_results <- data.frame(matrix(nrow=0,ncol=10))
colnames(fdr_disease_results) <- c('anal_id','score','n','trait','cases','beta','se','z','p','q')

results_age_bands <- NULL
res_age_bands_regress <- NULL

res<-read_pj("../p82_comb_assoc/st82_01_assoc_res.txt")
res<- res[!duplicated(res$formula),]
res$rsid<-gsub("_.","",res$score)
res <- res[res$rsid%in%table2$rsid,]
res$phewas <- ifelse(res$rsid %in% phwas_snps$rsid,T,F)
res_phewas <- res[!grepl("lif|extll|height|bmi",res$anal_id) & res$phewas,c("formula","p")]
res_phewas$q <- p.adjust(res_phewas$p,method="fdr")
res_phewas$p <- NULL
res <- merge_pj_many(df=res,new_data=res_phewas,id_var="formula")
res$phewas <- NULL
row.names(res) <- 1:nrow(res)




g2 <- list()

for (snp in unique(res$rsid)){
#   for (snp in "rs429358"){


heading(snp)


heading("Convert mortality hazards to longevity & set parental effects * 2")
hit_result<-table2[table2$rsid==snp,]
if(nrow(hit_result)<1) {cat("No data for",snp); next}
print(hit_result,row.names=F)
res1<-res[res$rsid==snp ,]


pdf(paste0("./all/st021_82_sex_n_age_",make.names(hit_result$gene),".pdf"),width=7,height=5)

score_name_res1<-res1[1,2]
effect_allele_res1<-gsub(".*_","",score_name_res1)
hit_result_effect_allele<-hit_result[1,"a1"]

replace_X_Y <- function(x,y,vector){

    vector <- gsub(paste("_",x,sep=""),paste("_",y,sep=""),vector)
    return(vector)
}

if(effect_allele_res1 != hit_result_effect_allele) {
    res1$beta<- -res1$beta
    res1$score<-replace_X_Y(effect_allele_res1,hit_result_effect_allele,res1$score)
} 

#convert to benefical - care this code does not generalise
res1$beta <- res1$beta * ifelse(grepl("age",res1$trait),-1,1)
#convert effect sizes to in self ie. parent *2
res1$beta <- res1$beta * ifelse(grepl("fath|moth",res1$trait),2,1)

res1 <- res1[order(res1$p),]

snp_results <- res1[!is.na(res1$q)&res1$q<0.05,c(1:9,16)]
if(snp %in% phwas_snps$rsid) fdr_disease_results <- rbind(fdr_disease_results, snp_results)

print(snp_results,row.names=F)
res1$formula <-  gsub("\\+ pc[0-9]* *","",res1$formula)


res1_age_bands <- res[res$rsid==snp & grepl("_age_",res$trait),]
res1_age_bands$sex <- ifelse(grepl("fath",res1_age_bands$trait),"M","F")
res1_age_bands$age_range <- gsub(".*_age_","",res1_age_bands$trait)

effect_allele_res1_age_bands<-gsub(".*_","",res1_age_bands$score[1])
if(effect_allele_res1_age_bands != hit_result_effect_allele) res1_age_bands$beta<- -res1_age_bands$beta
res1_age_bands$beta<-res1_age_bands$beta*-2
res1_age_bands$se<-res1_age_bands$se*2
#res1_age_bands$age<-as.numeric(gsub("_.*","",res1_age_bands$age_range))+5
res1_age_bands$age<-res1_age_bands$trait_mean

res1_age_bands$gene<-hit_result[1,1]
results_age_bands <- rbind(results_age_bands,res1_age_bands)


heading("Age and sex specific?")
 m1<-lm(beta~sex+age,data=res1_age_bands,weights=1/se^2)
#print( summary(m1))
 #p_values <-  summary(m1)$coeff[c(2,3),4]
 m2<-rma(beta,se^2,data=res1_age_bands,mods=~sex+age,method="FE")
 print(summary(m2))
 p_values<-summary(m2)$pval[c(2,3)]
p_sex_text<-paste("sex p:",ifelse(p_values[1]>0.001, sprintf("%0.3f", p_values[1]), sprintf("%0.0e", p_values[1])))
p_age_text<-paste("age p:",ifelse(p_values[2]>0.001, sprintf("%0.3f", p_values[2]), sprintf("%0.0e", p_values[2])))

res_age_bands_regress1<- t( c(summary(m2)$beta[c(2,3)],summary(m2)$se[c(2,3)],summary(m2)$pval[c(2,3)])[c(1,3,5,2,4,6)])
colnames(res_age_bands_regress1) <- c("sexM_beta","sexM_se","sexM_p","age_beta","age_se","age_p")
#add row descriptor
res_age_bands_regress1 <- cbind(res1_age_bands[1,c(20,2,4,3)],res_age_bands_regress1)
res_age_bands_regress <- rbind(res_age_bands_regress,res_age_bands_regress1)

t1<- paste(names(hit_result),collapse="\t")
t2<-  apply(hit_result,1,paste,collapse="\t")
#titl <- paste(t1,t2,sep="\n")
#title_size <- 10
titl <- paste0(hit_result[1,1]," (",snp,"_",hit_result[1,6],")")
title_size <- 12


pd <- position_dodge(0.1)
y_plot_lim <- with(res1_age_bands, max(abs(2*se+c(beta,-beta))))
g1<-ggplot(res1_age_bands, aes(x=age_range, y=beta, colour=sex,group=sex)) + 
    geom_errorbar(aes(ymin=beta-qnorm(0.975)*se, ymax=beta+qnorm(0.975)*se), width=.1,position=pd) +
    geom_hline(yintercept=0) +
    geom_line(position=pd) +
    geom_point(position=pd) + theme_light() +
   ggtitle(titl) + 
    annotate("text", label = p_sex_text, x = 0.5, y = y_plot_lim, color = "black", hjust = 0) +
    annotate("text", label = p_age_text, x = 0.5, y = y_plot_lim/1.1, color = "black", hjust = 0) +
    scale_color_manual(values=c("#F8766D","#00BFC4"),name="Sex",labels=c("Female","Male"))

g2[[hit_result[1,1]]] <- g1 +
 scale_x_discrete(name=NULL, labels=c("40-50","50-60","60-70","70-80","80-90","90-120")) +
 scale_y_continuous(name=NULL,breaks=seq(-2,2,ifelse(max(res1_age_bands$beta+qnorm(0.975)*res1_age_bands$se)>0.25,0.2,0.1)), minor_breaks=seq(-2,2,0.1), limits= c(-y_plot_lim,y_plot_lim)) +
 theme(title=element_text(size=title_size-2))

print(
    g1 + scale_x_discrete(name="Age range", labels=c("40-50","50-60","60-70","70-80","80-90","90-120")) + 
    scale_y_continuous(name="Beta",limits= c(-y_plot_lim,y_plot_lim)) + theme(title=element_text(size=title_size))
    )

dev.off()
}


fdr_disease_results$a1 <- substr(fdr_disease_results$score,nchar(fdr_disease_results$score),nchar(fdr_disease_results$score))
fdr_disease_results$rsid <- sub("_.$","",fdr_disease_results$score)

for (snp in unique(fdr_disease_results$score)) {
fdr_disease_results[fdr_disease_results$score==snp,"gene"] <- table2[paste0(table2$rsid,"_",table2$a1)==snp,"gene"]
fdr_disease_results[fdr_disease_results$score==snp,"cohort"] <- table2[paste0(table2$rsid,"_",table2$a1)==snp,"cohort"]
}

fdr_disease_results <- fdr_disease_results[order(match(fdr_disease_results$gene,table2$gene)),]



heading("Disease associations 5% FDR")
print(fdr_disease_results,row.names=F)

write.table(fdr_disease_results,"st021_83_phewas_results.tsv",sep="\t",quote=F,col.names=T,row.names=F)



results_age_bands_nice<-data.table(results_age_bands[,c(20,2,4,3,{5:9})])
results_age_bands_nice<- data.frame(results_age_bands_nice[order(match(gene,table2$gene),trait),])

heading("Age and sex specific HRRs (doubled)")
results_age_bands_nice


write.table(results_age_bands_nice,"st021_83_02_age_bands_results.tsv",sep="\t",quote=F,col.names=T,row.names=F)

res_age_bands_regress$sexM_q <- p.adjust(res_age_bands_regress$sexM_p,method="hochberg")
res_age_bands_regress$age_q <- p.adjust(res_age_bands_regress$age_p,method="hochberg")
res_age_bands_regress$weight <- log10(res_age_bands_regress$sexM_p)^3+log10(res_age_bands_regress$age_p)^3
res_age_bands_regress$weight2 <- scale(res_age_bands_regress$sexM_se,T,T)^3+scale(res_age_bands_regress$age_se,T,T)^3
res_age_bands_regress <- res_age_bands_regress[order(res_age_bands_regress$weight,-res_age_bands_regress$weight2),!grepl("weight",colnames(res_age_bands_regress))]

heading("Age and sex specificity")
print(res_age_bands_regress,row.names=F)

write.table(res_age_bands_regress,"st021_83_03_age_bands_results_q.tsv",sep="\t",quote=F,col.names=T,row.names=F)

heading("Age and sex specificity pass FDR 5%")
res_fdr <- res_age_bands_regress[res_age_bands_regress$sexM_q<.05|res_age_bands_regress$age_q<.05,]
print(res_fdr,row.names=F)




g2 <- g2[names(g2) %in% res_fdr$gene]


#----
# Arrangement #1
#----

#top <- g2[1:ceiling(length(g2)/2)]
#bottom <- g2[(1+ceiling(length(g2)/2)):length(g2)]

#png(paste0("st021_82_sex_n_age.png"),width=1200,height=675)
#top <- ggarrange(plotlist=top,ncol=length(top),nrow=1,common.legend=T,align="hv",legend="right")
#top <- annotate_figure(top, left="Beta (95% CI)", fig.lab.size=2)
#bottom <- ggarrange(plotlist=bottom,ncol=length(bottom),nrow=1,common.legend=T,align="hv",legend="right")
#bottom <- annotate_figure(bottom, left="Beta (95% CI)", bottom="Age range", fig.lab.size=2)
#ggarrange(plotlist=list(top,NULL,bottom), nrow=3, heights=c(25,1,25))
#dev.off()


#----
# Arrangement #2
#----

tl <- g2[res_fdr[res_fdr$age_q<0.05,]$gene[1:2]]
tr <- g2[res_fdr[res_fdr$sexM_q<0.05,]$gene[1]]

bl <- g2[res_fdr[res_fdr$age_q<0.05,]$gene[3:4]]
br <- g2[res_fdr[res_fdr$sexM_q<0.05,]$gene[2]]

pdf("st021_82_sex_n_age.pdf", width=12, height=8)
tl <- ggarrange(plotlist=tl,ncol=length(tl),nrow=1,common.legend=T,align="hv",legend="none")
tl <- annotate_figure(tl, left="Beta (95% CI)", fig.lab.size=2)
tr <- ggarrange(plotlist=tr,ncol=length(tr),nrow=1,common.legend=T,align="hv",legend="right")
tr <- annotate_figure(tr, left="Beta (95% CI)", fig.lab.size=2)
bl <- ggarrange(plotlist=bl,ncol=length(bl),nrow=1,common.legend=T,align="hv",legend="none")
bl <- annotate_figure(bl, left="Beta (95% CI)", bottom="Age range", fig.lab.size=2)
br <- ggarrange(plotlist=br,ncol=length(br),nrow=1,common.legend=T,align="hv",legend="right")
br <- annotate_figure(br, left="Beta (95% CI)", bottom="Age range", fig.lab.size=2)

left <- ggarrange(plotlist=list(tl,NULL,bl),ncol=1,nrow=3,common.legend=T,align="hv",legend="none", heights=c(30,1,30))
right <- ggarrange(plotlist=list(tr,NULL,br),ncol=1,nrow=3,align="hv", heights=c(30,1,30))
ggarrange(plotlist=list(left,NULL,right), ncol=3, widths=c(16,1,10),labels=c("A","","B"))

dev.off()

heading("done")
