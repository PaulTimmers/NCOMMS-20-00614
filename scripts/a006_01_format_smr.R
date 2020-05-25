#!/usr/bin/env Rscript

#----
# Setup environment
#----

options(width=200)
library(data.table)
library(ggplot2)
library(ggpubr)

files <- paste0("manova/",list.files("manova")[grepl("manova.tsv$",list.files("manova"))])

# Functions

cmyk = function(c, m, y, k=0, alpha, maxColorValue=1){
  if(maxColorValue != 1) { c <- c/maxColorValue; m <- m/maxColorValue; y <- y/maxColorValue }
  c <- (c * (1-k) + k)
  m <- (m * (1-k) + k)  
  y <- (y * (1-k) + k)

  c <- 1-c; m <- 1-m; y <- 1-y
  hex <- function(v) substring(rgb(v,0,0),2,3)
  if(!missing(alpha)) alpha <- hex(alpha) else alpha <- ''
  paste0('#',hex(c), hex(m), hex(y), alpha)
}


map_interval <- function(v, interval=c(0,1)) {
	if(length(interval)!=2) stop("Interval must contain exactly two values")
	v <- v - min(v, na.rm=T)
	intercept <- interval[1]
	slope <- (interval[1]-interval[2])/(min(v, na.rm=T)-max(v, na.rm=T))
	return(intercept + v*slope)
}


#----
# START
#----

plots <- list()

for (file in files) {
	
	# Get info
	rsid <- gsub("_.*$","",gsub("^.*rs","rs",file))
	gene <- fread(cmd=paste0("snp2gene.sh -s ",rsid))$summary
	if(gene=="HYKK") gene <- "CHRNA3/5"

	cat(paste0(gene," (",rsid,")"),"\n")

	# Read file
	df <- fread(file)

	df[health_a1 == A2, c("health_a1", "health_beta1"):=.(A1, -health_beta1)]
	df[life_a1 == A2, c("life_a1", "life_beta1"):=.(A1, -life_beta1)]
	df[long_a1 == A2, c("long_a1", "long_beta1"):=.(A1, -long_beta1)]

	df[,c("health_z1","life_z1","long_z1"):=.(health_beta1/health_se, life_beta1/life_se, long_beta1/long_se)]
	df[,z_add := health_z1 + life_z1 + long_z1]
	df[,z1 := sign(z_add) * abs(qnorm(p/2))]


	# Create plot
	df[,c("health","life","long") := .(abs(health_z1 * z_add), abs(life_z1 * z_add), abs(long_z1 * z_add))]
	df[,colour:=cmyk(health, life, long, k=0.3, maxColorValue=max(health, life, long))]
	ylim <- df[,max(-log10(p))]	

	p1 <- ggplot(df[order(abs(health_z1+life_z1+long_z1)),], aes(x=POS/1000000, y=-log10(p), colour=colour)) +
	       geom_hline(yintercept=-log10(5e-8), linetype=2) +
	       geom_point() +
	        geom_point(data=df[SNP==rsid,], size=4, shape=18) +
	         coord_cartesian(ylim=c(0, ylim + ylim/2.5), expand=0.1) +
	         geom_text(data=df[SNP==rsid,], label=rsid, size=4, nudge_y=ylim/5 ) +
	           scale_colour_identity() + labs(y="", x="") +
	            ggtitle(gene) + theme_bw()

	plots[[gene]] <- p1


	# Effective sample sizes based on events
	df[,n_effective := round(84949 + 609139 + 2/(1/11262 + 1/25483))]

	# Simple total
	df[,n_simple := 300447 + 1012240 + (11262 + 25483)]


	# Calculate betas and SE based on sample size and allele frequency
	df[, se_effective :=  1 / sqrt( 2 * freq * (1 - freq) * (n_effective + z1^2)) ]
	df[, beta1_effective := z1 * se_effective]

	df[, se_simple :=  1 / sqrt( 2 * freq * (1 - freq) * (n_simple + z1^2)) ]
	df[, beta1_simple := z1 * se_simple]

	out1 <- df[,.(SNP,A1,A2,freq,b=beta1_effective,se=se_effective,p,n=n_effective)]
	fwrite(out1, file, quote=F, sep="\t", na="NA")

	out2 <- df[,.(SNP,A1,A2,freq,b=beta1_simple,se=se_simple,p,n=n_simple)]
	fwrite(out2, paste0(file,".simple"), quote=F, sep="\t", na="NA")	
}

bitmap("manova/st006_01_regional_plots.png", width=12, height=8, res=600)
gg <- ggarrange(plotlist=plots)
annotate_figure(gg, left="-log10(P)", bottom="GRCh37 position (Mb)")
dev.off()


bitmap("manova/st006_01_regional_plots_of_interest.png", width=12, height=4, res=600)
gg <- ggarrange(plotlist=plots[c("LPA","FOXO3", "FGD6", "LINC02513", "SLC4A7", "APOE", "TOX3", "ZW10", "LDLR", "CDKN2B-AS1")], nrow=2, ncol=5)
annotate_figure(gg, left="-log10(P)", bottom="GRCh37 position (Mb)")
dev.off()