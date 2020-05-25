#!/usr/bin/env Rscript

#----
# Setup environment
#----

options(width=200)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(mygene))

#----
# Convert Gene to Entrez ID
#----

# Load SMR highlighted genes

loci <- fread("../q06_smr_heidi/st006_04_smr_cis_5pct.tsv")
loci <- rbind(loci, fread("../q06_smr_heidi/st006_04_smr_trans_5pct.tsv"))
loci <- loci[!duplicated(gene),]

# Convert Names to Entrez ID -- method #1

loci[,entrez:=mapIds(org.Hs.eg.db, gene, 'ENTREZID', 'SYMBOL')]
loci[is.na(entrez), entrez:=mapIds(org.Hs.eg.db, gene, 'ENTREZID', 'ALIAS')]


# Convert Names to Entrez ID -- method #2

for (g in loci[is.na(entrez), gene]) {
	d <- try(data.frame(query(g)), silent=T)
	if(class(d) == "try-error") next
	if(is.null(d$hits.entrezgene)) next
	loci[gene==g, entrez:=d$hits.entrezgene[1]]
}

k <- nrow(loci[!duplicated(gene) & !is.na(entrez),])

cat("\n\n================================\n List of our eQTLs with Entrez IDs \n================================\n\n")

print(loci[!duplicated(gene),][order(as.numeric(entrez))][!is.na(entrez)])

cat("\n\nGenes missing entrez IDs:",paste0(loci[!duplicated(gene) & is.na(entrez),gene],collapse=" "),"\n")

# Load all possible SMR genes

if( file.exists("datasets/gene_to_entrez.tsv.gz")) {

	all_eqtls <- fread("datasets/gene_to_entrez.tsv.gz")


} else {

	all_eqtls <- data.table(gene=system("cat ../q06_smr_heidi/*.epi | tr -s ' ' '\t' | cut -f5 | sort -u", intern=T))
	all_eqtls <- all_eqtls[!duplicated(gene),]


	# Convert Names to Entrez ID -- method #1

	all_eqtls[,entrez:=mapIds(org.Hs.eg.db, gene, 'ENTREZID', 'SYMBOL')]
	all_eqtls[is.na(entrez),entrez:=mapIds(org.Hs.eg.db, gene, 'ENTREZID', 'ALIAS')]

	# Convert Names to Entrez ID -- method #2

	lookup <- all_eqtls[is.na(entrez), sort(unique(gene))]
	lookup_res <- data.table(data.frame(queryMany(lookup, scopes=c("symbol", "reporter","accession"), fields=c("entrezgene"), species="human")))
	lookup_res <- lookup_res[is.na(notfound),]
	lookup_res <- lookup_res[all_eqtls,,on=c(query="gene")]
	lookup_res[is.na(entrez), entrez:=entrezgene]

	all_eqtls <- lookup_res[!duplicated(entrez),.(gene=query, entrez)]
	#all_eqtls <- 19960 # https://www.eqtlgen.org/index.html

	fwrite(all_eqtls[!is.na(entrez)],"datasets/gene_to_entrez.tsv", sep="\t", quote=F, na="NA")
	system("gzip -9f datasets/gene_to_entrez.tsv")
}

cat("\n\n================================\n List of eQTLs with Entrez IDs \n================================\n\n")
print(all_eqtls[order(as.numeric(entrez))][!is.na(entrez)])


#----
# Process gene sets
#----

min_genes <- 3
for (collection in  c("hallmark","go_bioproc")) { # "go_bioproc" "canonical"

cat("gene_set", "n", "n_set", "p", "\n", file=paste0("st011_01_",collection,"_stats.tsv"), sep="\t") # Print header
cat("gene_set","locus","gene","\n", file=paste0("st011_01_",collection,"_genes.tsv"), sep="\t") # Print header

n_gene_sets <- system(paste0("wc -l < datasets/",collection,".entrez.gmt"), intern=T) # Find number of gene sets in collection

# Loop through gene sets

for (i in 1:n_gene_sets) {
	cat(collection, paste0(i,"/",n_gene_sets), "\r")
	gene_set <- fread(cmd=paste0("awk 'FNR == ",i," {print; exit}' datasets/",collection,".entrez.gmt | tr -s '\t' '\n' | sed '/http/d' "), data.table=FALSE)
	
	q <- sum(loci$entrez %in% gene_set[,1])
	
	if (q < min_genes) next
	m <- all_eqtls[entrez %in% gene_set[,1], .N]
	n <- all_eqtls[,.N] - m
	p <- phyper(q, m, n, k, lower.tail=FALSE) # Test for enrichment

	cat(names(gene_set), q, m, p, "\n", file=paste0("st011_01_",collection,"_stats.tsv"), sep="\t", append=TRUE) # Write results to file
	fwrite(data.table(gene_set=names(gene_set), loci[entrez %in% gene_set[,1], .(locus,gene)]), paste0("st011_01_",collection,"_genes.tsv"), sep="\t", append=TRUE, quote=FALSE) # Write results to file
}


# Sort and calculate FDR

collection_stats <- fread(paste0("st011_01_",collection,"_stats.tsv"), col.names=c("gene_set","n","n_set","p","q"))[order(p)]
collection_stats[,q:=p.adjust(p, "bonferroni")]; fwrite(collection_stats, paste0("st011_01_",collection,"_stats.tsv"), sep="\t", quote=F)

collection_genes <- fread(paste0("st011_01_",collection,"_genes.tsv"), select=1:3, fill=T, header=T)
collection_genes[,gene_set:=factor(gene_set, levels=collection_stats[,gene_set])]
collection_genes <- collection_genes[order(gene_set)]; fwrite(collection_genes, paste0("st011_01_",collection,"_genes.tsv"), sep="\t", quote=F)

cat("\n")
}
