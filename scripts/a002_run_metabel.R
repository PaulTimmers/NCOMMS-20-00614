#/usr/bin/env R

#----
# Setup environment
#----

options(width=200)
library(data.table)
library(MultiABEL)
sink("/dev/null");source(paste0(root,"apps_by_us/gen_useful.R"));sink()


args <- commandArgs(T)

health_file <-args[1]
long_file <-args[2]
life_file <- args[3]
out_name <- args[4]

indep_snp_file <- "../data/d001_indep.snps.RData"

sessionInfo()

#========== FIX PROBLEMATIC FUNCTION ============================================================

environment <- environment("MultiABEL")
MultiSummary_edit <- function (x, index = NULL, type = "direct", vars = NULL) {
    if (class(x) != "multi.summary") {
        stop("Wrong data type!")
    }
    cat("Multi-trait genome scan ... ")
    if (!is.null(index) & length(index) < nrow(x$cor.pheno)) {
        m <- nrow(x$cor.pheno)
        x$cor.pheno <- x$cor.pheno[index, index]
        x$var.pheno <- x$var.pheno[index]
        idx <- which(!(1:m %in% index))
        x$gwa <- x$gwa[, -c(idx * 2 - 1, idx * 2)]
    }
    m <- nrow(x$cor.pheno)
    f <- x$gwa[, 2 * m + 1]
    n <- x$gwa[, 2 * m + 2]
    k <- nrow(x$gwa)
    betamat <- x$gwa[, seq(1, 2 * m, 2)]
    semat <- x$gwa[, seq(2, 2 * m, 2)]
    tmat <- betamat/semat
    bad <- which(rowSums(is.na(tmat)) > 0)
    if (length(bad) > 0) {
        betamat <- betamat[-bad, ]
        semat <- semat[-bad, ]
        tmat <- tmat[-bad, ]
        x$gwa <- x$gwa[-bad, ]
        f <- f[-bad] # CHANGED
        n <- n[-bad] # CHANGED
        k <- nrow(x$gwa)
    }
    dimnames(betamat) <- list(NULL, NULL)
    res <- data.frame(marker = rownames(x$gwa), freq = f, n = n, stringsAsFactors = FALSE)
    rownames(res) <- rownames(x$gwa)
    if (type != "direct") {
        bnames <- paste0("beta.cond.", rownames(x$cor.pheno))
        snames <- paste0("se.cond.", rownames(x$cor.pheno))
        pnames <- paste0("p.cond.", rownames(x$cor.pheno))
        cns <- rep(NA, m * 3)
        cns[seq(1, m * 3, 3)] <- bnames
        cns[seq(2, m * 3, 3)] <- snames
        cns[seq(3, m * 3, 3)] <- pnames
        res3 <- matrix(NA, k, m * 3)
        colnames(res3) <- cns
        bnames <- paste0("est.coef.", rownames(x$cor.pheno))
        snames <- paste0("se.coef.", rownames(x$cor.pheno))
        pnames <- paste0("p.coef.", rownames(x$cor.pheno))
        cns[seq(1, m * 3, 3)] <- bnames
        cns[seq(2, m * 3, 3)] <- snames
        cns[seq(3, m * 3, 3)] <- pnames
        resc <- matrix(NA, k, m * 3)
        colnames(resc) <- cns
    }
    if (type == "outbred") {
        scan1 <- .Fortran("MultiSummaryLoopDirect", k = as.integer(k), 
            m = as.integer(m), nn = as.numeric(n), tmat = tmat, 
            invR = solve(x$cor.pheno), pil = numeric(k), fstat = numeric(k), 
            PACKAGE = "MultiABEL")
        scan2 <- .Fortran("MultiSummaryLoop", k = as.integer(k), 
            m = as.integer(m), nn = as.numeric(n), f = as.numeric(f), 
            betamat = betamat, R = x$cor.pheno, invR = solve(x$cor.pheno), 
            D = matrix(0, m, m), sdY = diag(sqrt(x$var.pheno)), 
            invsdY = diag(1/sqrt(x$var.pheno)), sY = sqrt(x$var.pheno), 
            b = betamat, s = betamat, pil = numeric(k), coef = betamat, 
            ss = betamat, PACKAGE = "MultiABEL")
        res2 <- data.frame(p = pf(scan1$fstat, m, n - m - 1, 
            lower.tail = FALSE), beta.score = scan2$pil)
        res2$se.score <- res2$beta.score/sqrt(qchisq(res2$p, 
            1, lower.tail = FALSE))
        res3[, seq(1, m * 3, 3)] <- scan2$b
        res3[, seq(2, m * 3, 3)] <- scan2$s
        res3[, seq(3, m * 3, 3)] <- pchisq(scan2$b^2/scan2$s^2, 
            1, lower.tail = FALSE)
        res <- cbind(res, res2, res3)
        resc[, seq(1, m * 3, 3)] <- scan2$coef
        resc[, seq(2, m * 3, 3)] <- sqrt(scan2$ss)
        resc[, seq(3, m * 3, 3)] <- pchisq(scan2$coef^2/scan2$ss, 
            1, lower.tail = FALSE)
    }
    else if (type == "inbred") {
        scan1 <- .Fortran("MultiSummaryLoopDirect", k = as.integer(k), 
            m = as.integer(m), nn = as.numeric(n), tmat = tmat, 
            invR = solve(x$cor.pheno), pil = numeric(k), fstat = numeric(k), 
            PACKAGE = "MultiABEL")
        scan2 <- .Fortran("MultiSummaryLoopInbred", k = as.integer(k), 
            m = as.integer(m), nn = as.numeric(n), f = as.numeric(f), 
            betamat = betamat, R = x$cor.pheno, invR = solve(x$cor.pheno), 
            D = matrix(0, m, m), sdY = diag(sqrt(x$var.pheno)), 
            invsdY = diag(1/sqrt(x$var.pheno)), sY = sqrt(x$var.pheno), 
            b = betamat, s = betamat, pil = numeric(k), coef = betamat, 
            ss = betamat, PACKAGE = "MultiABEL")
        res2 <- data.frame(p = pf(scan1$fstat, m, n - m - 1, 
            lower.tail = FALSE), beta.score = scan2$pil)
        res2$se.score <- res2$beta.score/sqrt(qchisq(res2$p, 
            1, lower.tail = FALSE))
        res3[, seq(1, m * 3, 3)] <- scan2$b
        res3[, seq(2, m * 3, 3)] <- scan2$s
        res3[, seq(3, m * 3, 3)] <- pchisq(scan2$b^2/scan2$s^2, 
            1, lower.tail = FALSE)
        res <- cbind(res, res2, res3)
        resc[, seq(1, m * 3, 3)] <- scan2$coef
        resc[, seq(2, m * 3, 3)] <- sqrt(scan2$ss)
        resc[, seq(3, m * 3, 3)] <- pchisq(scan2$coef^2/scan2$ss, 
            1, lower.tail = FALSE)
    }
    else if (type == "precise" & !is.null(vars)) {
        scan1 <- .Fortran("MultiSummaryLoopDirect", k = as.integer(k), 
            m = as.integer(m), nn = as.numeric(n), tmat = tmat, 
            invR = solve(x$cor.pheno), pil = numeric(k), fstat = numeric(k), 
            PACKAGE = "MultiABEL")
        scan2 <- .Fortran("MultiSummaryLoopPrecise", k = as.integer(k), 
            m = as.integer(m), nn = as.numeric(n), varg = vars, 
            betamat = betamat, R = x$cor.pheno, invR = solve(x$cor.pheno), 
            D = matrix(0, m, m), sdY = diag(sqrt(x$var.pheno)), 
            invsdY = diag(1/sqrt(x$var.pheno)), sY = sqrt(x$var.pheno), 
            b = betamat, s = betamat, pil = numeric(k), coef = betamat, 
            ss = betamat, PACKAGE = "MultiABEL")
        res2 <- data.frame(p = pf(scan1$fstat, m, n - m - 1, 
            lower.tail = FALSE), beta.score = scan2$pil)
        res2$se.score <- res2$beta.score/sqrt(qchisq(res2$p, 
            1, lower.tail = FALSE))
        res3[, seq(1, m * 3, 3)] <- scan2$b
        res3[, seq(2, m * 3, 3)] <- scan2$s
        res3[, seq(3, m * 3, 3)] <- pchisq(scan2$b^2/scan2$s^2, 
            1, lower.tail = FALSE)
        res <- cbind(res, res2, res3)
        resc[, seq(1, m * 3, 3)] <- scan2$coef
        resc[, seq(2, m * 3, 3)] <- sqrt(scan2$ss)
        resc[, seq(3, m * 3, 3)] <- pchisq(scan2$coef^2/scan2$ss, 
            1, lower.tail = FALSE)
    }
    else if (type == "direct") {
        scan <- .Fortran("MultiSummaryLoopDirect", k = as.integer(k), 
            m = as.integer(m), nn = as.numeric(n), tmat = tmat, 
            invR = solve(x$cor.pheno), pil = numeric(k), fstat = numeric(k), 
            PACKAGE = "MultiABEL")
        res$p <- pf(scan$fstat, m, n - m - 1, lower.tail = FALSE)
    }
    else {
        stop("Wrong type of analysis!")
    }
    cat("Done.\n")
    if (type %in% c("outbred", "inbred", "precise")) {
        colnames(scan2$coef) <- rownames(x$cor.pheno)
        rownames(scan2$coef) <- rownames(res)
        return(list(scan = res, coef = resc))
    }
    else {
        return(list(scan = res, coef = NULL))
    }
}

#========== END PROBLEMATIC FUNCTION ============================================================


#----
# Calculate correlation
#----


# Load data

heading(paste("loading",health_file,long_file, life_file))

health <- gzfread(health_file, select=c("rsid","beta1","a1","a0"))
names(health)[2:4] <- paste0(names(health)[2:4], "_health")

long <- gzfread(long_file, select=c("rsid","beta1","a1","a0"))
names(long)[2:4] <- paste0(names(long)[2:4], "_long")

life <- gzfread(life_file,select=c("rsid","beta1","a1","a0"))
names(life)[2:4] <- paste0(names(life)[2:4], "_life")


# Select independent SNPs

load(indep_snp_file)
print(data.table(indep.snps))

health <- health[rsid%in%indep.snps,]
health <- health[!duplicated(rsid),]

long <- long[rsid%in%indep.snps,]
long <- long[!duplicated(rsid),]

life <- life[rsid%in%indep.snps,]
life <- life[!duplicated(rsid),]

# Align alleles

df1 <- merge_pt_data_table(health,life,id_var="rsid")

df <- merge_pt_data_table(df1,long,id_var="rsid")
df <- df[complete.cases(df),]

cat("\n\n")
print(df)
cat("\n\n")

cat("N flipped alleles (healthspan):",df[a1_health == a0_life,.N],"\n")
df[a1_health == a0_life,c("a1_health","a0_health","beta1_health"):=.(a0_health,a1_health,-beta1_health)]
cat("N mismatched alleles (healthspan):",nrow(df) - df[a1_health == a1_life & a0_health == a0_life,.N], "\n")
df <- df[a1_health == a1_life & a0_health == a0_life, ]


cat("N flipped alleles (longevity):",df[a1_long == a0_life,.N],"\n")
df[a1_long == a0_life,c("a1_long","a0_long","beta1_long"):=.(a0_long,a1_long,-beta1_long)]
cat("N mismatched alleles (longevity):",nrow(df) - df[a1_long == a1_life & a0_long == a0_life,.N], "\n")
df <- df[a1_long == a1_life & a0_long == a0_life, ]

cat("\n\n")
print(df)
cat("\n\n")


# Calculate correlation


cormat <- cor(df[,.(beta1_health, beta1_long, beta1_life)])


#----
# Prepare for analysis
#----

manova_data <- load.summary(c(health_file,long_file,life_file), est.var=TRUE, cor.pheno = cormat, columnNames = c("rsid", "a1", "a0", "freq1", "beta1", "se", "n"), indep.snps = indep_snp_file)

print(cormat)
print(data.table(study=colnames(cormat),trait_var=manova_data$var.pheno))
print(data.table(manova_data$alleles))


#----
# Run MANOVA
#----

heading(paste("doing manova"))
manova_result_out <- MultiSummary_edit(manova_data, type = 'outbred')



#----
# Reformat results
#----

gwas_metad<-manova_result_out$scan
names(gwas_metad)[c(1,2)]<-c("rsid","freq1")
names(gwas_metad) <- gsub(".cond.*",".cond", names(gwas_metad))
gwas_metad <- data.frame(gwas_metad)

gwas_alleles <- manova_data$alleles
gwas_alleles$rsid <- rownames(gwas_alleles)
names(gwas_alleles) <- c("a1", "a0", "rsid")

gwas_alleles <- data.table(gwas_alleles, key="rsid")
gwas_metad <- data.table(gwas_metad, key="rsid")

gwas_metad1 <- gwas_alleles[gwas_metad,,on="rsid"]


heading(paste("loading snp info data from",life_file))
gwas_info <- gzfread(life_file,select=c(1:6,12))

gwas_info <- gwas_info[!duplicated(rsid),]
gwas_metad1 <- data.table(gwas_metad1)[!duplicated(rsid),]


gwas_refmt <- gwas_info[gwas_metad1,.(rsid,snpid,chr,pos,a1=i.a1,a0=i.a0,n,freq1,beta1=sprintf("%.6g",beta.score),se=sprintf("%.6g",se.score),p=sprintf("%.6g",p),info,
    beta.cond.health=sprintf("%.4g",beta.cond),se.cond.health=sprintf("%.4g",se.cond),beta.cond.long=sprintf("%.4g",beta.cond.1),se.cond.long=sprintf("%.4g",se.cond.1),
    beta.cond.life=sprintf("%.4g",beta.cond.2),se.cond.life=sprintf("%.4g",se.cond.2)),on="rsid"]



cat("N invalid statistics:",gwas_refmt[se==0,.N],"\n")
print(gwas_refmt[se==0,])
gwas_refmt <- gwas_refmt[se>0,]



#----
# Export results
#----

gwas_refmt <- gwas_refmt[order(chr,pos),]
print(gwas_refmt)



rare <- gwas_refmt[freq1 < 0.05 | freq1 > 0.95,.N]
common <- gwas_refmt[freq1 >= 0.05 & freq1 <= 0.95,.N]
maf_range <- gwas_refmt[,range(freq1)]

print(data.table(n_snps=gwas_refmt[,.N],rare,common,maf_range=paste0(maf_range,collapse="-")))

heading(paste("writing ",out_name))
fwrite(gwas_refmt, file=out_name, na="NA", sep="\t", quote=FALSE)

heading("DONE!")



