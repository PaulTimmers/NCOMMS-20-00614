#!/usr/env/R

#----
# Setup environment 
#----

options(width=200)
library(data.table)
library(phenoscanner)

#----
# Load data 
#----

files <- list.files("proxy_info")
files <- paste0("proxy_info/",files[grepl("st04_00_rsids",files)])

results <- data.table()
for (file in files) {
rsids <- fread(file, data.table=FALSE)[,]
rsid_info <- phenoscanner(snpquery = rsids, catalogue = "GWAS", pvalue = 1, proxies = "None", r2 = 0.8, build = 37)
res1 <- data.table(rsid_info$results)
results <- rbind(results,res1)
}

fwrite(results, "st04_02_phenoscanner.tsv", na="NA", quote=FALSE, sep="\t")
system("gzip -9 -f st04_02_phenoscanner.tsv")