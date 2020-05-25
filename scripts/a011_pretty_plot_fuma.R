#!/usr/bin/env Rscript

#----
# Setup environment
#----

set.seed(1)
library(data.table)
options(width=200)
library(ggplot2)

# Colorblind-friendly palette 
okabe_ito_pal <- c(rgb(0,0,0), rgb(0.9,0.6,0), rgb(0.35, 0.7, 0.9), rgb(0,0.6,0.5), rgb(0.95, 0.9, 0.25), rgb(0,0.45,0.7), rgb(0.8,0.4,0), rgb(0.8,0.6,0.7), "grey")


#----
# Load data
#----

loci <- fread("../q06_smr_heidi/st006_04_smr_cis_5pct.tsv")
loci <- rbind(loci, fread("../q06_smr_heidi/st006_04_smr_trans_5pct.tsv"))
loci <- loci[!duplicated(gene),]

all_loci <- loci[order(chr,pos),unique(locus)]



#----
# Plot Hallmark
#----

collection <- "hallmark"

stat_df <- fread(paste0("st011_01_",collection,"_stats.tsv"))[q < 0.05, ][order(-n, p)]
gene_df <- fread(paste0("st011_01_",collection,"_genes.tsv"))[gene_set %in% stat_df$gene_set]
plot_df <- stat_df[gene_df,,on="gene_set"]

plot_df[,gene_set:=gsub("HALLMARK_","",gene_set)]


# Order by locus, then gene_set

plot_df[,gene_set:=factor(gene_set, levels=plot_df[order(q),unique(gene_set)])]
plot_df[,locus:=factor(locus, levels=plot_df[,.(.N,q=unique(q)),by=c("gene_set","locus")][order(-N,q),unique(locus)])]


plot_df <- plot_df[order(locus, gene_set, -n)]

plot_df[,c("gene_set","gene"):=.(factor(gene_set, rev(unique(gene_set))), factor(gene, levels=unique(gene)))]

# Plot

pdf(paste0("st011_02_",collection,"_plot.pdf"), width=8, height=3)

print(ggplot(plot_df, aes(x=gene, y=gene_set, fill=locus)) + 
 geom_tile(colour="white") + 
  scale_fill_manual(values=okabe_ito_pal) +
   coord_equal() + 
    theme_classic() + 
     guides(fill = guide_legend(nrow = 1, title="Genomic locus containing gene expression signal", title.position="bottom", title.hjust=0.5, title.theme=element_text(face="bold"))) +
      theme(axis.text.x=element_text(angle=45, hjust=1), axis.title.x=element_blank(), 
      	axis.title.y=element_blank(), legend.position="bottom"))

dev.off()

# Table

table1 <- plot_df[,.(n,n_set,p,q,locus,gene=paste0(gene,collapse=", ")),by=c("gene_set","locus")][!duplicated(paste0(gene_set,n,n_set,p,q,locus))][order(gene_set,nchar(gene),decreasing=T)][,.(gene_set, n, n_set, p, q, locus, gene)]
fwrite(table1, paste0("st011_03_",collection,"_table.csv"), sep=",", quote=T, na="NA")

#----
# Plot GO Biological Processes
#----

collection <- "go_bioproc"

stat_df <- fread(paste0("st011_01_",collection,"_stats.tsv"))[q < 0.05, ][order(-n, p)]
gene_df <- fread(paste0("st011_01_",collection,"_genes.tsv"))[gene_set %in% stat_df$gene_set]
plot_df <- stat_df[gene_df,,on="gene_set"]


# Cluster based on genes

heat_df <- fread(cmd=paste0("zcat datasets/",collection,".entrez.heat.gz"), header=TRUE)
heat_df <- heat_df[gene_set %in% plot_df$gene_set,]
heat_df <- heat_df[,unlist(lapply(heat_df, function(x) {!all(x == 0 | is.na(x))})), with=F]
heat_df$tot <- apply(heat_mat[,-1],1,sum)

gene2entrez <- fread("datasets/gene_to_entrez.tsv.gz")
plot_df$entrez <- gene2entrez[match(plot_df$gene,gene),entrez]


heat_mat <- as.matrix(heat_df[,-1])
rownames(heat_mat) <- heat_df[,gene_set]
heat_mat_subset <- heat_mat[,colnames(heat_mat) %in% plot_df$entrez]


# Determine number of clusters

wss <- (nrow(heat_mat_subset)-1)*sum(apply(heat_mat_subset,2,var))
for (i in 2:20) wss[i] <- sum(kmeans(heat_mat_subset,
   centers=i)$withinss)
plot(1:20, wss, type="b", xlab="Number of Clusters",
  ylab="Within groups sum of squares")
abline(v=seq(0,25,5), lty=2)

groups <- 8
k <- kmeans(dist(heat_mat_subset), centers=groups)
grouping <- data.table(gene_set=names(k$cluster), group=k$cluster)[order(group)]

print(grouping)

for (set in grouping$gene_set) {
  d <- system(paste0("curl -silent http://software.broadinstitute.org/gsea/msigdb/cards/",set," | grep -A 1 '<th>Brief description</th>' | sed '/Brief description/d; s/^.*<td>//; s=</td>==; s/ \\[.*$//'"), intern=T)
  cat(paste0(set,"(",grouping[gene_set==set,group],"):\n",d,"\n\n"))
  grouping[gene_set==set, description:=d]
}

plot_df <- grouping[plot_df,,on="gene_set"]


#----
# Create sensible names for clusters
#----

### THIS IS MANUAL STUFF ###
### RE-DO IF UPDATING CLUSTERS ###


plot_df1 <- plot_df[!duplicated(paste0(group,gene,locus)),]

cat(paste0(plot_df[group==1, sort(unique(gene_set))],collapse="\n"),"\n")
plot_df1[group==1, gene_set:="Breakdown of H2O2, antibiotics or cofactors"]

cat(paste0(plot_df[group==2, sort(unique(gene_set))],collapse="\n"),"\n")
plot_df1[group==2, gene_set:="Erythrocyte and myeloid cell homeostasis"]

cat(paste0(plot_df[group==3, sort(unique(gene_set))],collapse="\n"),"\n")
plot_df1[group==3, gene_set:="Catecholamine (dopamine) transport"]

cat(paste0(plot_df[group==4, sort(unique(gene_set))],collapse="\n"),"\n")
plot_df1[group==4, gene_set:="Apoptosis or chemical homeostasis"]

cat(paste0(plot_df[group==5, sort(unique(gene_set))],collapse="\n"),"\n")
plot_df1[group==5, gene_set:="Spleen development and cytoskeletal organisation"]

cat(paste0(plot_df[group==6, sort(unique(gene_set))],collapse="\n"),"\n")
plot_df1[group==6, gene_set:="Response to organic cyclic (steroid) compounds"]

cat(paste0(plot_df[group==7, sort(unique(gene_set))],collapse="\n"),"\n")
plot_df1[group==7, gene_set:="Organic cation (ammonium) transport"]

cat(paste0(plot_df[group==8, sort(unique(gene_set))],collapse="\n"),"\n")
plot_df1[group==8, gene_set:="Apoptotic signaling pathway by p53 class mediator"]




# Order by locus, then gene_set

plot_df1[,gene_set:=factor(gene_set, levels=plot_df1[order(q),unique(gene_set)])]
plot_df1[,locus:=factor(locus, levels=plot_df1[,.N,by=c("gene_set","locus")][order(-N),unique(locus)])]
plot_df1 <- plot_df1[order(locus, gene_set, -n)]

plot_df1[,c("gene_set","gene"):=.(factor(gene_set, rev(unique(gene_set))), factor(gene, levels=unique(gene)))]



# Plot

pdf(paste0("st011_02_",collection,"_plot.pdf"), width=12, height=8)

print(ggplot(plot_df1, aes(x=gene, y=gene_set, fill=locus)) + 
 geom_tile(colour="white") + 
  scale_fill_manual(values=okabe_ito_pal) +
   coord_equal() +
     theme_classic() + 
      guides(fill = guide_legend(nrow = 1, title="Genomic locus containing gene expression signal", title.position="bottom", title.hjust=0.5, title.theme=element_text(face="bold"))) +
      theme(axis.text.x=element_text(angle=45, hjust=1), axis.title.x=element_blank(), 
      	axis.title.y=element_blank(), legend.position="bottom"))


dev.off()

# Table

plot_df[,gene_set:=factor(gene_set, levels=plot_df[order(group,q,-n),unique(gene_set)])]
table2 <- plot_df1[,.(group_name=unique(gene_set)),by= "group"][plot_df,,on="group"]
table2 <- table2[,.(n,n_set,p,q,locus,gene=paste0(gene,collapse=", ")),by=c("group","group_name","gene_set","locus")][!duplicated(paste0(gene_set,n,n_set,p,q,locus,group))][order(gene_set,-nchar(gene),decreasing=F)][,.(group,group_name,gene_set,n,n_set,p,q,locus,gene)]
fwrite(table2, paste0("st011_03_",collection,"_table.csv"), sep=",", quote=T, na="NA")