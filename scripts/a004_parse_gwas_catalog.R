#!/usr/bin/env R

#----
# Setup environment
#----

options(width=200)
library(data.table)
library(RSQLite)


# Define functions

shorten_text <- function(text, length=50) {
	words <- strsplit(text, split=" ")
	for (i in 1:length(words)) {
		word <- words[[i]]
		if(sum(nchar(word)) <= length | any(is.na(nchar(word)))) {words[[i]] <- paste(word, collapse=" "); next}
		repeat {
			if(sum(nchar(word)) <= length | any(is.na(nchar(word))) | nchar(word)[1] == sum(nchar(word))) break
			word <- word[-length(word)/2]
		}
		words[[i]] <- paste(c(word[1:(length(word)/2)],"...",word[((length(word)/2)+1):length(word)]),collapse=" ")
	}
	return(unlist(words))
}


print.top.bottom <- function(df,top=5,bottom=2,text.width=50) {
	min <- top+bottom
	for (j in colnames(df)[unlist(lapply(df,is.numeric))]) set(df, j = j, value = ifelse(abs(df[[j]])>0.5, round(df[[j]],4), ifelse(abs(df[[j]])<=1e-3, sprintf("%.2e",df[[j]]), sprintf("%.4f",df[[j]])) ))
	for (j in colnames(df)[unlist(lapply(df,is.character))]) set(df, j = j, value = shorten_text(df[[j]], text.width))

	if(nrow(df)<=min) return(print(df, nrow=nrow(df)))

	out <- df[c(1:top,(nrow(df)-bottom):nrow(df)),]
	out[, ] <- lapply(out[, ], as.character)
	out[,` `:=paste0(c(1:top,(nrow(df)-bottom):nrow(df)),":")]
	out[top+1,] <- ""
	out[top+1,ncol(out)] <- "---"
	return(print(out[, c(ncol(out),1:(ncol(out)-1)),with=F],nrow=nrow(df),row.names=F))
}






#----
# Create catalog
#----

if (!file.exists("st04_01_gwas_catalog.tsv.gz")) {

# Get association data

gwas_cat_association <- fread("curl -s ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/2020/03/09/gwas-catalog-associations.tsv",na.strings="NR",colClasses="character",select=c(2,8,9,10,21,27,28,29,31,32),col.names=c("pmid","trait","disco","repli","snp_a1","freq1","p","log_p","beta_or","ci"))

gwas_cat_ancestry <- fread("curl -s ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/2019/10/14/gwas-catalog-ancestry.tsv",na.strings="NR",colClasses="character",sep="\t",select=c(2,9),col.names=c("pmid","ancestry"),fill=T)
gwas_cat_ancestry <- gwas_cat_ancestry[!duplicated(pmid),]

# Load known hg19 chomosome and positions

snp_pos <- fread(paste0(root,"locuszoom_v1.4/data/gwas_catalog/gwas_catalog_hg19.txt"), select=c("snp","chr","pos"),sep="\t")
snp_pos2 <- fread("zcat ../q01_sumstats/st01_01_healthspan.tsv.gz", select=c("rsid","chr","pos"), col.names=c("snp","chr","pos"))
snp_pos <- rbind(snp_pos, snp_pos2)
snp_pos <- snp_pos[!duplicated(snp),]

# Subset to European ancestry
gwas_cat <- gwas_cat_association[gwas_cat_ancestry,on="pmid"]


# Reformat
gwas_cat$p <- 10^(-as.numeric(gwas_cat$log_p))
gwas_cat$a1 <- gsub("^.*-","",gwas_cat$snp_a1)
gwas_cat$snp <- gsub("-.*$","",gwas_cat$snp_a1)
gwas_cat$snp <- gsub(" ","",gwas_cat$snp)
gwas_cat[grep("[^rs|0-9]",snp),snp:=gsub("[^rs|0-9]","",snp)]

gwas_cat <- gwas_cat[grep("rs",snp),]
gwas_cat <- snp_pos[gwas_cat,on="snp"]
gwas_cat$beta1 <- with(gwas_cat, round(ifelse(grepl("crease",ci), as.numeric(beta_or), log(as.numeric(beta_or))),6))
gwas_cat$beta1 <- with(gwas_cat, ifelse(grepl("decrease",ci), -beta1, beta1))
gwas_cat$freq1 <- with(gwas_cat, as.numeric(gsub("[^0-9\\.]","",freq1)))
gwas_cat$units <- with(gwas_cat, gsub("\\[.*\\][ ]?| ..crease","",ci))
gwas_cat$units <- with(gwas_cat, gsub("NR ","",units))

# Subset to genome-wide significant
gwas_cat <- gwas_cat[nchar(beta_or)>0 & p<=5e-8,]

# Remove duplicates
gwas_cat <- gwas_cat[order(pmid,snp,trait,p),]
writeLines(paste0("GWAS CATALOG: Removing ",gwas_cat[duplicated(paste0(pmid,snp,trait)),.N], " snps with the same trait but higher p-value in same study"))
gwas_cat <- gwas_cat[!duplicated(paste0(pmid,snp,trait)),]

missing <- gwas_cat[is.na(pos),.(snp)]
if (nrow(missing) > 0) {
fwrite(missing, "st04_01_gwas_catalog.missing.tsv",sep="\t", na="NA", quote=F)
found <- fread("grep -wFf st04_01_gwas_catalog.missing.tsv ${ROOT}/data/reference_data/1000genomes/phase3/plink/1000genomes_chr*.bim ")
found <- found[,.(snp=V2, chr=as.numeric(gsub("^.*:","",V1)),pos=V4)]
gwas_cat <- found[gwas_cat,,on="snp"]
gwas_cat[,c("chr","pos","i.chr","i.pos"):=.(ifelse(is.na(chr),i.chr,chr), ifelse(is.na(pos),i.pos,pos), NULL, NULL)]
system("rm st04_01_gwas_catalog.missing.tsv")
}

# Get unknown chromosome and position data (slow)
SQLite()
con = dbConnect(drv=SQLite(), dbname=paste0(root,'locuszoom_v1.4/data/database/locuszoom_hg19.db'))
dbExecute(con,"CREATE TEMP VIEW snp_pos_trans AS SELECT rs_orig as snp,chr,pos FROM snp_pos p INNER JOIN refsnp_trans t ON (t.rs_current = p.snp)" )


i <- 0
rsids <- unique(gwas_cat[is.na(pos),snp])
writeLines(paste0("Looking up chromosome and position information for ",length(rsids)," SNPs..."))

for (rsid in sort(rsids)) {
i <- i + 1
if (i %% 20 == 1) cat(sprintf("%.2f%%\r",100*i/length(rsids)))
p1 = dbGetQuery(con,paste0("SELECT snp,chr,pos FROM snp_pos_trans WHERE snp='",rsid,"'"))
if(nrow(p1)==1) gwas_cat[snp==rsid,c("chr","pos"):=.(p1$chr,p1$pos)]
}
cat("100.00%\n")
rm(rsid)

# Tidy up
catalog <- gwas_cat[,.(rsid=snp,snpid=paste(chr,pos,sep="_"),a1,freq1=as.numeric(sprintf("%.2f", freq1)),trait_beta1=as.numeric(sprintf("%.4g",beta1)),trait_se=as.numeric(sprintf("%.4g",-beta1/qnorm(p/2))),trait_p=as.numeric(sprintf("%.4g",p)),trait_name=trait,units,pmid)]


fwrite(catalog, "st04_01_gwas_catalog.tsv", sep="\t", na="NA")
system("gzip -9 -f st04_01_gwas_catalog.tsv")

}



#----
# Load catalog
#----

catalog <- fread("zcat st04_01_gwas_catalog.tsv.gz", colClasses=c(rep("character",3), rep("numeric",3), rep("character",4)))
catalog[, trait_p := as.numeric(trait_p)]




#----
# Get longevity data
#----


# Set proxy data to match longevity allele

gws_hits <- fread("zcat proxy_info/st04_01_*_info.tsv.gz | awk '$1 != \"rsid\"'", header=FALSE, col.names=c("rsid","a1","a0","proxy","p.a1","p.a0","r2"))
life_df <- fread("zcat ../q03_analyse/st03_05_manova_lifespan_stats.tsv.gz") 
life_df[beta1 < 0,c("a1","a0","freq1","beta1"):=.(a0,a1,1-freq1,-beta1)] # Set longevity alleles
genes <- fread("../q03_analyse/st03_02_triple_manova_genes.tsv", select=c("summary","rsid"), col.names=c("gene","rsid"))
life_df <- genes[life_df,,on="rsid"][match(genes$gene,gene)]



#----
# Align LDLink
#----

# Align LDproxy alleles to longevity allele

for (rs in unique(life_df$rsid)) {
	life_allele <- life_df[rsid==rs, a1]
	gws_hits[rsid==rs, gene:=life_df[rsid==rs, gene]]
	gws_hits[rsid==rs & a1 != life_allele, c("a1","a0","p.a1","p.a0"):=.(a0,a1,p.a0,p.a1)]
}



#----
# Get catalog hits
#----


cat <- catalog[,.(proxy=rsid, snpid, trait_a1=a1, trait_freq1=freq1, trait_beta1, trait_se, trait_p, trait_name, units, pmid)]


# merge catalog hits with SNPs and close proxies (r2 > 0.6)

diseases <- cat[gws_hits, on="proxy"]
diseases[,trait_freq1:=as.numeric(trait_freq1)]



# If disease allele is known, set to longevity allele

diseases[trait_a1==p.a0, c("trait_a1","trait_freq1","trait_beta1"):=.(p.a1,1-trait_freq1,-trait_beta1)]



# Format nicely

diseases <- diseases[,.(gene,rsid,a1,a0,proxy,p.a1,p.a0,r2,trait_a1,trait_freq1,trait_beta1,trait_se,trait_p,trait_name,units,pmid)]



#----
# Add phenoscanner
#----

phenos <- fread("zcat st04_02_phenoscanner.tsv.gz", colClasses=c(rep("character",11), rep("numeric",4),rep("character",7)))
#phenos <- phenos[ancestry=="European",]
phenos <- phenos[p < 5e-8,]
phenos[,c("rsid","proxy"):=.(NULL, snp)]

phenos <- gws_hits[phenos,, on="proxy"]

# Align to longevity allele
phenos[i.a1==p.a0, c("i.a1","a2","beta","direction"):=.(p.a1,i.a1,-beta,ifelse(direction=="-","+","-"))]
phenos <- phenos[,.(gene,rsid, a1, a0, proxy, p.a1, p.a0, r2, trait_a1=i.a1, trait_freq1=NA, trait_beta1=beta, trait_se=se, trait_p=as.numeric(p), trait_name=trait, units=unit, pmid)]


# Correct unit decreases

phenos[grepl("decrease",units),c("trait_beta1","units"):=.(-trait_beta1,gsub("decrease","",units))]
phenos[grepl("crease",units),units:=gsub(" ..crease","",units)]
phenos[units=="-",units:=""]

# Other corrections


phenos[trait_name=="Cholesterol",trait_name:="Total cholesterol"]





#----
# Perform QC
#----


disease_df <- rbind(diseases, phenos)
disease_df <- disease_df[!is.na(pmid),]
disease_df <- disease_df[!duplicated(paste0(rsid,pmid,trait_p,trait_beta1,trait_se))]
disease_df[grepl("^unit",units), units:=""]
disease_df <- disease_df[!grepl("lifespan|longevity|healthspan|age at death|mortality|still alive|death in close genetic family",trait_name,ignore.case=TRUE),] # Don't care about traits in our own study
disease_df <- disease_df[!grepl("expression|transcript",trait_name,ignore.case=TRUE),] # Exclude expression levels
disease_df <- disease_df[!grepl("medication|treatment|therapy",trait_name,ignore.case=TRUE),] # Exclude medications and treatments
disease_df <- disease_df[!grepl("none of the above",trait_name,ignore.case=TRUE),] # Exclude useless traits
disease_df <- disease_df[order(trait_p),]


# Remove associations without effect sizes
no_effect <- disease_df[is.na(trait_beta1),]
writeLines(paste0("Removing ",no_effect[,.N], " disease associations without effect sizes")) 
print.top.bottom(no_effect[order(trait_p),], 10, 2)
writeLines("")
disease_df <- disease_df[!is.na(trait_beta1),]



# Remove tri-allelic

triallelic <- disease_df[trait_a1!="?" & trait_a1!=p.a1 & trait_a1!=p.a0,]
triallelic <- disease_df[proxy%in%triallelic$proxy,]

if(nrow(triallelic)>0) {
writeLines(paste0("Removing ",triallelic[,.N], " tri-allelic disease associations")) 
print.top.bottom(triallelic[order(proxy,trait_a1==p.a1),], 10, 2)
writeLines("")
disease_df <- disease_df[!proxy%in%triallelic$proxy,]
}




#----
# Duplicates
#----

#### MANUAL CURATION

disease_df[grepl("coronary artery calcification",trait_name,ignore.case=TRUE), trait_name:="coronary artery calcification"] # format for detecting duplicates
disease_df[,trait_name:=gsub("II","2",trait_name)] # format for detecting duplicates
disease_df[,trait_name:=gsub(" levels","",trait_name)] # format for detecting duplicates
disease_df[grepl("Lpa|Lp a|Lp \\(a\\)",trait_name),trait_name:="Lipoprotein a"] # format for detecting duplicates
disease_df[grepl("cholesterol",trait_name,ignore.case=TRUE) & grepl("total",trait_name,ignore.case=TRUE) & !grepl("HDL|LDL",trait_name),trait_name:="Total cholesterol"] # format for detecting duplicates
disease_df[grepl("cholesterol ldl",trait_name,ignore.case=TRUE),trait_name:="LDL cholesterol"] # format for detecting duplicates
disease_df[grepl("lipoproteins ldl",trait_name,ignore.case=TRUE),trait_name:="LDL cholesterol"] # format for detecting duplicates
disease_df[grepl("breast",trait_name,ignore.case=TRUE) & grepl("neoplasm", trait_name, ignore.case=T),trait_name:="Breast cancer"] # format for detecting duplicates
disease_df[,trait_name:=gsub("adenocarcinoma|carcinoma","cancer",trait_name)] # format for detecting duplicates
disease_df[grepl("^alzheimer",trait_name,ignore.case=TRUE) & grepl("disease$",trait_name,ignore.case=TRUE),trait_name:="Alzheimer's disease"] # format for detecting duplicates
disease_df[grepl("^Cerebrospinal",trait_name,ignore.case=TRUE) & grepl("tau$",trait_name,ignore.case=TRUE), trait_name:="Cerebrospinal fluid tau"] # format for detecting duplicates
disease_df[,trait_name:=gsub("Systolic |Diastolic ","",trait_name)] # format for detecting duplicates
disease_df[,trait_name:=gsub("percentage of white cells","count",trait_name)] # format for detecting duplicates
disease_df[,trait_name:=gsub("Coronary heart disease","Coronary artery disease",trait_name)] # format for detecting duplicates
disease_df[,trait_name:=gsub("Blue vs green eyes","Eye color",trait_name)] # format for detecting duplicates
disease_df[,trait_name:=gsub("Illnesses of father:","Father",trait_name)] # format for detecting duplicates
disease_df[,trait_name:=gsub("Pulse pressure","Blood pressure",trait_name)] # format for detecting duplicates
disease_df[grepl("pulse rate",trait_name,ignore.case=TRUE),trait_name:="Blood pressure"] # format for detecting duplicates
disease_df[grepl("skin cancer|malignant neoplasms of skin",trait_name,ignore.case=TRUE),trait_name:="Melanoma"] # format for detecting duplicates
disease_df[,trait_name:=gsub("^blood","Blood",trait_name)] # format for detecting duplicates
disease_df[grepl("Vitiligo",trait_name, ignore.case=TRUE),trait_name:="Vitiligo"] # format for detecting duplicates
disease_df[grepl(" tan ",trait_name, ignore.case=TRUE),trait_name:="Tanning"] # format for detecting duplicates
disease_df[grepl("sunburns|suntan",trait_name, ignore.case=TRUE),trait_name:="Sunburns"] # format for detecting duplicates
disease_df[grepl("freckl",trait_name, ignore.case=TRUE),trait_name:="Freckling"] # format for detecting duplicates
disease_df[grepl("blue vs. green eyes",trait_name, ignore.case=TRUE),trait_name:="Eye color"] # format for detecting duplicates
disease_df[grepl("hair color",trait_name, ignore.case=TRUE),trait_name:="Hair color"] # format for detecting duplicates
disease_df[grepl("logical memory",trait_name, ignore.case=TRUE),trait_name:="Logical memory"] # format for detecting duplicates
disease_df[grepl("balding pattern",trait_name,ignore.case=TRUE),trait_name:="Male-pattern baldness"]
disease_df[trait_name=="Diabetes mellitus type 1",trait_name:="Type 1 diabetes"]
disease_df[trait_name=="Diabetes mellitus type 2",trait_name:="Type 2 diabetes"]
disease_df[grepl("glucose",trait_name,ignore.case=TRUE),trait_name:=gsub("2 hour glucose", "2 hour fasting glucose", trait_name)]
disease_df[grepl("Vascular or heart problems diagnosed by doctor: ",trait_name), trait_name:=gsub("Vascular or heart problems diagnosed by doctor: ","",trait_name)]
disease_df[grepl("Pain type experienced in last month: ",trait_name), trait_name:=gsub("Pain type experienced in last month: ","",trait_name)]
disease_df[grepl("low density lipoprotein", trait_name, ignore.case=T), trait_name:=gsub("low density lipoprotein", "LDL cholesterol", trait_name, ignore.case=T)]
disease_df[grepl("high density lipoprotein", trait_name, ignore.case=T), trait_name:=gsub("high density lipoprotein", "HDL cholesterol", trait_name, ignore.case=T)]
disease_df[,trait_name:=gsub("Intracranial abdominal aortic or thoracic aortic aneurysm pleiotropy","Intracranial aneurysm",trait_name)]
disease_df[,trait_name:=gsub("Self-reported ", "", trait_name, ignore.case=T)]
disease_df[,trait_name:=gsub("Alzheimer's disease in hypertension.*","Alzheimer's disease",trait_name)]
disease_df[,trait_name:=gsub("^.*alzheimers disease or dementia","Alzheimer's disease",trait_name)]
disease_df[trait_name=="Supranuclear palsy progressive",trait_name:="Progressive supranuclear palsy"]
disease_df[grepl("qualifications", trait_name, ignore.case=T),trait_name:=gsub("Qualifications:.*","Qualifications",trait_name)]
disease_df[,trait_name:=gsub("BMI","Body mass index",trait_name)]
disease_df[grepl(" fat mass",trait_name,ignore.case=T),trait_name:="Fat mass"]
disease_df[grepl(" predicted mass",trait_name,ignore.case=T),trait_name:="Fat mass"]
disease_df[grepl(" fat percentage",trait_name,ignore.case=T),trait_name:="Fat percentage"]
disease_df[grepl("impedance of ",trait_name,ignore.case=T),trait_name:="Impedance"]


# Don't care about directionality
disease_df[grepl("left$",trait_name, ignore.case=T), trait_name:=gsub("left","",trait_name, ignore.case=T)]
disease_df[grepl("right$",trait_name, ignore.case=T), trait_name:=gsub("right","",trait_name, ignore.case=T)]

# Combine kin
disease_df[,trait_name:=gsub("Illnesses of mother: ","",trait_name)]
disease_df[,trait_name:=gsub("Illnesses of father: ","",trait_name)]
disease_df[,trait_name:=gsub("Illnesses of siblings: ","",trait_name)]
disease_df[grepl("father",trait_name,ignore.case=T), trait_name:=gsub("father","",trait_name,ignore.case=T)]
disease_df[grepl("mother",trait_name,ignore.case=T), trait_name:=gsub("mother","",trait_name,ignore.case=T)]

# Capitalise
disease_df[substr(trait_name,1,1)==" ", trait_name:=substr(trait_name,2,999)]
disease_df[,trait_name:=(paste0(substr(toupper(trait_name),1,1),substr(trait_name,2,999)))]


# Identify duplicates

writeLines("\n\nIdentify similar traits\n=======================\n")

disease_df$group <- 0
i <- 1

for (t in unique(disease_df[order(nchar(trait_name)),trait_name])) {
	if(disease_df[trait_name==t, any(group>0)]) next
	
	matches <- unique(agrep(paste0("^",t), disease_df[group==0,trait_name], value=T))
	
	if(grepl("HDL",t)) matches <- unique(matches[!grepl("LDL",matches)], unique(grep(paste0("high density lipo"), disease_df[group==0,trait_name], ignore.case=T, value=T)))
	if(grepl("LDL",t)) matches <- unique(matches[!grepl("HDL",matches) & !grepl("VLDL",matches)], unique(grep("low density lip", disease_df[group==0,trait_name], ignore.case=T, value=T)))
	if(grepl("Type 2",t)) matches <- matches[!grepl("Type 1",matches)]
	if(grepl("Type 1",t)) matches <- matches[!grepl("Type 2",matches)]
	if(grepl("Red",t)) matches <- matches[!grepl("White",matches)]
	if(grepl("White",t)) matches <- matches[!grepl("Red",matches)]
	if(grepl("APOE",t)) matches <- matches[!grepl("APO[BA]",matches)]
	if(grepl("APOB",t)) matches <- matches[!grepl("APO[EA]",matches)]
	if(grepl("APOA",t)) matches <- matches[!grepl("APO[BE]",matches)]
	if(grepl("Lipoprotein",t)) matches <- matches[!grepl("phospholipase",matches)]
	if(grepl("Myocardial",t)) matches <- matches[!grepl("Coronary",matches)]
 	if(grepl("Asthma",t)) matches <- c(matches, "Self-reported asthma")
	
	if(grepl("Hemoglobin",t)) matches <- matches[!grepl("Glycated",matches)]

	disease_df[trait_name%in%matches,group:=i]

	if(length(matches)>1) {
		
		display <- data.table(data.frame(disease_df[group==i]))
		print(display[, .(gene,traits=shorten_text(trait_name,130),pmid,group=i)][order(gene,nchar(traits)),][!duplicated(paste0(gene,traits,pmid,group)),])
		disease_df[trait_name%in%matches,trait_name:=t]
		cat("\n")	
	}

	i <- i+1
}




# Remove duplicates

disease_df <- disease_df[order(is.na(trait_beta1),trait_p),]

dup <- disease_df[duplicated(paste0(rsid,group)),]
dup <- rbind(dup, disease_df[!duplicated(paste0(rsid,group)), ][duplicated(paste0(rsid,trait_name)), ])

writeLines(paste0("\nRemoving ",dup[,.N], " non-unique disease associations per locus")) 
print.top.bottom(dup[order(group,trait_p),], 10, 2)
writeLines("")
disease_df <- disease_df[!duplicated(paste0(rsid,group)),]
disease_df <- disease_df[!duplicated(paste0(rsid,trait_name)),]
disease_df[,group:=NULL]




# Report number of ambigious direction of effect SNPs

ambiguous <- disease_df[(trait_a1=="?" | is.na(trait_a1)) & is.na(trait_freq1), ]
writeLines(paste0(ambiguous[,.N], " SNPs have ambiguous trait direction\n"))




#----
# Predict missing alleles and align
#----

writeLines("\n\nAlign alleles\n=============\n")

disease_df_aligned <- data.table(disease_df[0,])

for (rs in sort(unique(disease_df$rsid))) {
	writeLines(paste0("\n\n===== ", rs, " ====="))
	
	
	
	df <- disease_df[rsid==rs, ]
	
	freq1 <- mean(df[trait_a1 == p.a1, .SD[which.max(r2) & !is.na(trait_freq1),trait_freq1] ])

	mismatched_freq <- df[trait_a1!="?" & (0.5-trait_freq1)*(0.5-freq1)<0 & abs(trait_freq1-freq1)>0.3,]
	
	if(mismatched_freq[,.N]>0) {
		writeLines(paste0("Removing ", mismatched_freq[,.N], " SNPs with wildly different allele frequency (ref = ",round(freq1,4),")"))
		print.top.bottom(mismatched_freq[order(trait_p),], 5, 2)
		df <- df[!paste0(proxy,trait_name,pmid) %in% mismatched_freq[,paste0(proxy,trait_name,pmid)],]
		writeLines("")
	}
	
	freq1 <- mean(df[trait_a1 == p.a1, .SD[which.max(r2) & !is.na(trait_freq1),trait_freq1] ])

	df[trait_a1=="?" & (0.5-trait_freq1)*(0.5-freq1)>0 & abs(trait_freq1-freq1)<0.3, trait_a1 := a1]
	df[trait_a1=="?" & (0.5-trait_freq1)*(0.5-freq1)<0 & abs(trait_freq1-freq1)>0.3, c("trait_freq1","trait_beta1","trait_a1") := .(1-trait_freq1,-trait_beta1, a1)]
	
	df[trait_a1 == p.a1 & is.na(trait_freq1), trait_freq1:=freq1]
	
	df <- df[order(r2, rank(-log10(trait_p)), decreasing=TRUE),]

	disease_df_aligned <- rbind(disease_df_aligned, df)
	print.top.bottom(df,20,2,50)
}


#----
# Assign categories
#----

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5300907/pdf/main.pdf

writeLines("\n\nAssign categories\n================\n")

disease_df_aligned[,disease:="Miscellaneous"]
disease_df_aligned[grepl("ankle-brachial|heart|cholesterol|lipid|CVD|Carotid|LDL|HDL|IDL|DBP|SBP|myocardial|cardiovascular|coronary|hypertension|pressure|statin|triglyc|aort|Ankle brachial index|artery disease|thrombo|GpTot|lipoprotein|phospholipase|aneurysm|stroke|angina|vascular|wheeze or whistling in the chest",trait_name,ignore.case=T),"disease"] <- "Cardiovascular"
disease_df_aligned[grepl("homa b|diab|glucos|HbA1C|insulin|glycated|impedance",trait_name,ignore.case=T),"disease"] <- "Metabolic"
disease_df_aligned[grepl("neuro|cognitive impairment|alzh|parkin|hunting|synap|progranulin|depress|schizo|attention|amyloid|lewy body|epilepsy|cross disorder|dementia|insomnia|restless leg|palsy",trait_name,ignore.case=T),"disease"] <- "Neuropsychiatric"
disease_df_aligned[grepl("asthma|basophil|monocyte|eosinophil|lymphocyte|neutrophil|white blood|immun|rheumat|inflamm|crohns|ulcerative colitis|celiac|coeliac|multiple sclerosis|arthritis|irritible bowel|psoriasis|lupus|graves|myasthenia|primary biliary cirrhosis|thyroid peroxidase|allerg|viti[l]+igo|cholangitis|thyroidi|crohn|c-reactive|sarcoidosis",trait_name,ignore.case=T),"disease"] <- "Immune-related"
disease_df_aligned[grepl("waist|hip |fat | mass|^[^adjusted for]*bmi*$|weight|obesity|childhood bmi|BMI in|body mass|BMI tails|metabolic syndrome|basal metabolic| predicted mass|adipose tissue",trait_name,ignore.case=T),"disease"] <- "Metabolic"
disease_df_aligned[grepl("neoplasm|cancer|noma|glioma",trait_name,ignore.case=T),"disease"] <- "Cancer"
disease_df_aligned[grepl("copd|pulmonary|lung|airflow obstruction|respiratory|FEV|emphysema|exhale|smok|cigar|nicotine|tobacco",trait_name,ignore.case=T),"disease"] <- "Smoking-related"
disease_df_aligned[grepl("allergic|type 1",trait_name,ignore.case=TRUE),disease:="Immune-related"]
disease_df_aligned[grepl("cerebrospinal",trait_name,ignore.case=TRUE),disease:="Miscellaneous"]

# Ageing traits
disease_df_aligned[grepl("age[- ]related| ageing| aging|relative age|gr[ae]y|baldness|hearing|macular|age at menarche|infirm",trait_name,ignore.case=TRUE),disease:="Age-related"]

# Fix specific cases
disease_df_aligned[grepl("Number of non-cancer illnesses",trait_name,ignore.case=T),disease:="Age-related"]
disease_df_aligned[grepl("No blood clot, bronchitis, emphysema, asthma, rhinitis, eczema or allergy diagnosed by doctor",trait_name,ignore.case=T),disease:="Miscellaneous"]
disease_df_aligned[grepl("Lipoprotein-associated phospholipase A2",trait_name, ignore.case=T),disease:="Cardiovascular"]


writeLines("\n# Unassigned\n")
print(disease_df_aligned[disease=="Miscellaneous" & !duplicated(trait_name),.(gene,rsid,trait_p,trait_name)][order(nchar(trait_name))], 9999)


writeLines("\n# Assigned\n")
d <- data.table(data.frame(disease_df_aligned[,.(gene,rsid,trait_beta1,trait_p,trait_name,disease)]))
print(d[,.(gene,rsid,trait_beta1,trait_p,trait=shorten_text(trait_name,100),disease)][order(gene,-abs(trait_beta1), trait_p),], 9999)


#----
# Summary
#----

writeLines("\n\nSummary1\n========\n")

# Where longevity allele is protective (i.e. negative beta on disease), note which categories have at least 1 association

smry <- disease_df_aligned[,.(rsid, .N),by="gene"][!duplicated(gene),]#[!duplicated(gene),][order(gene),]
main_disease <- disease_df_aligned[disease!="Miscellaneous", .(disease=paste(names(which(table(disease) >= .N/length(table(disease)))),collapse="/")),by="gene"]
smry <- main_disease[smry,.(gene,rsid,N,disease),on="gene"]

smry <- smry[match(life_df$gene,gene),]
smry[,c("gene","rsid"):=.(life_df$gene,life_df$rsid)]
smry[is.na(smry)] <- 0
smry[gene=="HYKK",gene:="CHRNA3/5"]
smry[gene=="SLC22A2/SLC22A3", gene:="SLC22A2/3"]
smry[gene=="HLA-DRB1/HLA-DQA1", gene:="HLA-DRB1/DQA1"]
smry[gene=="NOL4L/NOL4L-DT", gene:="NOL4L/-DT"]
smry[disease=="0",disease:="-"]
print(smry)

fwrite(smry, "st04_04_gws_hits_diseases_summary1.tsv", sep="\t", na="NA")

cat("")

writeLines("\n\nSummary2\n========\n")

smry2 <- disease_df_aligned[,.N, by=c("gene","disease")][order(gene,-N),]
disease_levels <- c(smry2[disease!="Miscellaneous",.(N=sum(N)),by="disease"][order(-N),disease],"Miscellaneous")
smry2[,disease:=factor(disease,levels=disease_levels)]
smry2 <- dcast(smry2, gene ~ disease, value.var="N")


smry2 <- smry2[match(substr(life_df$gene,0,5),substr(gene,0,5)),]
smry2[,c("gene"):=.(life_df$gene)]
smry2[is.na(smry2)] <- 0
smry2$Total <- smry2[,apply(.SD,1,sum),.SDcols=-1]
smry2 <- rbind(smry2, data.table(gene="ALL",smry2[,lapply(.SD,sum),.SDcols=-1]))
smry2[gene=="HYKK",gene:="CHRNA3/5"]
smry2[gene=="SLC22A2/SLC22A3", gene:="SLC22A2/3"]
smry2[gene=="HLA-DRB1/HLA-DQA1", gene:="HLA-DRB1/DQA1"]
smry2[gene=="NOL4L/NOL4L-DT", gene:="NOL4L"]
print(smry2)

fwrite(smry2, "st04_04_gws_hits_diseases_summary2.tsv", sep="\t", na="NA")

#----
# Export
#----

disease_df_aligned[units=="",units:=NA]
disease_df_aligned[gene=="HYKK",gene:="CHRNA3/5"]
disease_df_aligned[gene=="SLC22A2/SLC22A3", gene:="SLC22A2/3"]
disease_df_aligned[gene=="HLA-DRB1/HLA-DQA1", gene:="HLA-DRB1/DQA1"]
disease_df_aligned[gene=="NOL4L/NOL4L-DT", gene:="NOL4L/-DT"]
fwrite(disease_df_aligned, "st04_04_gws_hits_diseases_full.tsv", sep="\t", na="NA")


classification_df <- disease_df_aligned[order(gene),.(locus=paste0(gene,collapse=", "), .N),by=c("disease","trait_name","pmid")][order(disease,-N, -nchar(locus))]
fwrite(classification_df, "st04_05_gws_hits_diseases_classification.tsv", sep="\t", na="NA")
