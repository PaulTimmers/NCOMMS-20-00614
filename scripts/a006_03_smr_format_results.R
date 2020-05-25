#!/usr/bin/env Rscript

#----
# Setup environment
#----

options(width=200)
library(data.table)
flank <- 500000
smr_thresh <- 0.05
heidi_thresh <- 0.01

#----
# START
#----

# Read in data
res <- fread("st006_02_smr_results.txt")
res1 <- res[trait=="manova"]
#res1[,eqtls:=gsub("-smr-multi","",eqtls)]

# Make tissue and eqtl columns
res1[,tissue:=ifelse(grepl("gtex",eqtls), gsub("gtex-","",eqtls), "whole-blood")]
res1[,eqtls:=gsub("gtex.*$","gtex",eqtls)]


# Merge with data on loci of interest
loc <- fread("../q03_analyse/st03_03_triple_manova_res.tsv", select=c(1:3,5), col.names=c("rsid","chr","pos","locus"))
res2 <- loc[res1,,on="rsid"][order(chr,pos,p_SMR)]


# Subset to loci of interest FDR < 0.05
res3 <- res2[locus %in% c("SLC4A7", "LINC02513","FOXO3","ZW10","FGD6","LPA","CDKN2B-AS1","TOX3","LDLR","APOE")]
res3[,c("q_SMR"):=.(p.adjust(p_SMR, "fdr")),by=c("eqtls")]
res3 <- res3[q_SMR < smr_thresh & (p_HEIDI > heidi_thresh | is.na(p_HEIDI)), .(locus, chr, pos, gene=Gene, z_manova=b_GWAS/se_GWAS, z_eqtl=b_eQTL/se_eQTL, z_smr=b_SMR/se_SMR, p_smr=p_SMR, q_smr=q_SMR, p_heidi=p_HEIDI, n_heidi=nsnp_HEIDI, tissue, source=eqtls)]


# Set GWAS effect to positive, add direction column 
res3[z_manova < 0, c("z_manova","z_eqtl"):=.(-z_manova, -z_eqtl)]
res3[,dir:=ifelse(sign(z_smr)>0,"+","-")]


# Make tables with cis and trans genes
res3_cis <- res3[!grepl("trans",source),]
res3_trans <- res3[grepl("trans",source),]


# Output descriptives for manuscript
cat("\n\n====================== Cis Genes q_smr <",smr_thresh,"& p_heidi >",heidi_thresh,"======================\n\n")

tissue <- res3_cis[p_heidi > heidi_thresh,.(N_genes=length(unique(gene)), loci=paste0(unique(locus), collapse=", "), N_loci=length(unique(locus))),by="tissue"][order(-N_genes)]
tissue <- rbind(tissue, data.table(tissue=paste0("< ALL ",tissue[,.N]," tissues > "),N_genes=res3_cis[p_heidi > heidi_thresh,length(unique(gene))], loci=paste0(res3_cis[p_heidi > 0.05,unique(locus)], collapse=", "), N_loci=res3_cis[p_heidi > 0.05,length(unique(locus))]))
print(tissue)

cat("\nNumber of unique genes across all tissues:",res3_cis[,length(unique(gene))],"\n")


cat("\n\n====================== Cis Blood p_heidi >",heidi_thresh,"======================\n\n")

blood <- res3_cis[p_heidi > heidi_thresh & tissue == "whole-blood", .(locus, gene, tissue, dir)][!duplicated(paste0(locus,gene,tissue,dir))]
print(blood)

cat("\n\n====================== Cis Nerve Tibial p_heidi >",heidi_thresh,"======================\n\n")

tibial <- res3_cis[p_heidi > heidi_thresh & tissue == "nerve-tibial", .(locus, gene, tissue, dir)][!duplicated(paste0(locus,gene,tissue,dir))]
print(tibial)


cat("\n\n====================== Cis Top Loci ======================\n\n")

most_genes <- res3_cis[p_heidi > heidi_thresh,.(N_genes=length(unique(gene))),by="locus"][order(-N_genes)]
most_tissues <- res3_cis[p_heidi > heidi_thresh,.(N_tissues=length(unique(tissue))),by="locus"][order(-N_tissues)]

print(most_genes)
cat("\n\n")

print(most_tissues)



cat("\n\n====================== Trans Genes q_smr <",smr_thresh,"& p_heidi >",heidi_thresh,"======================\n\n")

trans <- res3_trans[p_heidi > heidi_thresh,.(locus, gene, tissue, dir)][!duplicated(paste0(locus,gene,tissue,dir))]
print(trans)


cat("\n\n====================== Cis Genes q_smr <",smr_thresh,"& is.na(p_heidi) ======================\n\n")

cis_strong <- res3_cis[p_heidi > heidi_thresh,.(gene, p_heidi, tissue, dir, .N),by="locus"][!duplicated(paste0(locus,gene,tissue,dir))]
cis_putative <- res3_cis[is.na(p_heidi),.(gene, p_heidi, tissue, dir, .N),by="locus"][!duplicated(paste0(locus,gene,tissue,dir))]
cis_add <- cis_putative[! gene %in% cis_strong$gene] 

print(cis_add)

cat("\nNumber of putative cis genes not previously reported in any tissue:",cis_add[,length(unique(gene))], "\n")

cis_add_tissue <- cis_add[,.(locus,tissue=paste0(tissue, collapse=", ")),by="gene"][!duplicated(paste0(locus,gene,tissue))]
print(cis_add_tissue)


cat("\n\n====================== Trans Genes q_smr <",smr_thresh,"& is.na(p_heidi) ======================\n\n")

trans_putative <- res3_trans[is.na(p_heidi),.(gene, p_heidi, tissue, dir, .N),by="locus"][!duplicated(paste0(locus,gene,tissue,dir))]
print(trans_putative)


cat("\n\n====================== Total unique genes q_smr <",smr_thresh,"& p_heidi >",heidi_thresh,"or NA ======================\n\n")

tot <- res3[!duplicated(gene),.(gene,trans=grepl("trans",source)),by=c("locus")]
print(tot)

trans_vs_cis <- res3[is.na(p_heidi) | p_heidi > heidi_thresh,.(locus, gene, trans=grepl("trans",source))][!duplicated(paste0(locus,gene,trans))][,.N,by=c("trans")]
print(trans_vs_cis)


cat("\n\n")

#----
# Format tables
#----


t1 <- res3[!grepl("trans",source),][order(chr,pos,p_heidi < 0.05, q_smr),][,.(cis=paste(unique(paste0(gene,dir)), collapse=", ")),by=c("locus","chr","pos")]
t2 <- res3[grepl("trans",source),][order(chr,pos,p_heidi < 0.05, q_smr),][,.(trans=paste(unique(paste0(gene,dir)), collapse=", ")),by=c("locus","chr","pos")]

t3 <- merge(t1, t2, all=T)

# Write tables to files

fwrite(format(res2[order(chr,pos,p_SMR)],digits=6), "st006_03_smr_loci_results.tsv", sep="\t", quote=F, na="NA")
fwrite(format(res3_cis[order(chr,pos,q_smr)],digits=6), "st006_04_smr_cis_5pct.tsv", sep="\t", quote=F, na="NA")
fwrite(format(res3_trans[order(chr,pos,q_smr)],digits=6), "st006_04_smr_trans_5pct.tsv", sep="\t", quote=F, na="NA")
fwrite(t3[order(chr,pos),], "st006_05_smr_genes_only.tsv", sep="\t", quote=F, na="NA")



#----
# CVD check
#----

#res_cvd <- res[trait=="cvd"][order(p_SMR, -p_HEIDI)]
#res_cvd1 <- res_cvd[p_SMR < 0.05 & p_HEIDI > 0.05,]
#res_cvd2 <- loc[res_cvd1,,on="rsid"][order(chr,pos,p_SMR)]