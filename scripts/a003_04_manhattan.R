#!/usr/bin/Rscript

#----
# Setup environment
#----

options(width=180)
library(data.table)


#----
# Define miami function
#----

miami <- function (x, y=NULL, chr = "chr", bp = "pos", p = "p", snp = "rsid", col1 = c("#e34a33", "#fdbb84"), col2 = c("#2b8cbe","#a6bddb"), chrlabs = NULL, suggestiveline = -log10(1e-05),
          genomewideline = -log10(5e-08), ymin = -15, ymax = 15, x.name = "x", y.name = "y", highlight = NULL, highlight_col="red", label = NULL, logp = TRUE, ymid = 0, ...)
  {
  CHR = BP = P = index = NULL
  ymid <- abs(ymid)
  ymin = ymin + ymid
  ymax = ymax - ymid

  #----
  # Check variables
  #----

  column_check <- function(column, df) {
    if (!(column %in% names(df)))
      stop(paste0("Column ", column, " not found in ",deparse(substitute(df)),"!"))
    if (!(column == snp) & !is.numeric(df[[column]]))
      stop(paste(column, "column in",deparse(substitute(df)),"should be numeric."))
  }

  for (column in c(bp,p,chr,snp)) column_check(column, x)
  if(!is.null(y)) for (column in c(bp,p,chr,snp)) column_check(column, y)

  #----
  # Create data frame
  #----
  print("##LOAD D1")
  
  d1 = data.table(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
  if (!is.null(x[[snp]])) d1$SNP <- x[[snp]]
  if (!is.null(highlight)) d1$highlight <- x[[highlight]]
  d1 <- subset(d1, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d1 <- d1[order(d1$CHR, d1$BP), ]
  if (logp) {
    d1$logp <- -log10(d1$P) - ymid
    d1 <- d1[d1$logp > 0,]
  } else {
    d1$logp <- d1$P - ymid
    d1 <- d1[d1$logp > 0,]
  }
  
  
  if(!is.null(y)) {
    print("##LOAD D2")
    d2 = data.table(CHR = y[[chr]], BP = y[[bp]], P = y[[p]])
    if (!is.null(y[[snp]])) d2$SNP <- y[[snp]]
    if (!is.null(highlight)) d2$highlight <- y[[highlight]]
    d2 <- subset(d2, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
    d2 <- d2[order(d2$CHR, d2$BP), ]
    if (logp) {
      d2$logp <- log10(d2$P) + ymid
      d2 <- d2[d2$logp < 0,]
    } else {
      d2$logp <- -d2$P + ymid
      d2 <- d2[d2$logp < 0,]
    }

  d <- data.frame(rbind(d1,d2))
  
  } else {
    d <- data.frame(d1)
  }

  #d <- d[abs(d$logp) > ymid,]
  
  #----
  # Define plotting parameters
  #----

  d$pos = NA
  d$index = NA
  ind = 0
  for (i in unique(d$CHR)) {
    ind = ind + 1
    d[d$CHR == i, ]$index = ind
  }
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    d$pos = d$BP
    ticks = floor(length(d$pos))/2 + 1
    xlabel = paste("Chromosome", unique(d$CHR), "position")
    labs = ticks
  } else {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + tail(subset(d, index ==
                                            i - 1)$BP, 1)
        d[d$index == i, ]$pos = d[d$index == i, ]$BP +
          lastbase
      }
      ticks = c(ticks, (min(d[d$index == i, ]$pos) + max(d[d$index ==
                                                             i, ]$pos))/2 + 1)
    }
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)

  print("##CALL PLOT")
  
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", yaxt = "n", cex.axis = 2, cex.lab = 2,
                   las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(ymin-1,ymax+1), xlab = xlabel, ylab = expression(-log[10](italic(p))),mar=5*c(5.1,6.1,4.1,2.1))
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))

  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  yticks <- seq(-5*ceiling(abs(ymin/5)), 5*ceiling(abs(ymax/5)), by = 5)
  yticks <- sort(c(yticks+sign(yticks)*(5-ymid),5-ymid,-5+ymid))
  if(is.null(y)) yticks <- yticks[yticks>=0]
  ylabels <- abs(yticks)+ymid
  
  if ((abs(suggestiveline) - ymid) > 0) {
    abline(h = suggestiveline - ymid, col = "blue")
    if(!is.null(y)) abline(h = -suggestiveline + ymid, col = "blue")
    }
  if ((abs(genomewideline) - ymid) > 0) {
    abline(h = genomewideline - ymid, col = "red")
    if(!is.null(y)) abline(h = -genomewideline + ymid, col = "red")
    }
  if(!is.null(y)) abline(h = 0, col = "black")
  #if(!is.null(y)) legend(0, y = ymax - ymid - 5*strheight(1), x.name, fill = col1[1], bty="n", cex = 3, title=NULL, yjust=0.5)
  #if(!is.null(y)) legend(0, y = -ymax + ymid + 5*strheight(1), y.name, fill = col2[1], bty="n", cex = 3, title=NULL, yjust=0.5)
  if(!is.null(y)) legend(0, y = -ymax + ymid + 5*strheight(1), c(x.name,y.name), fill = c(col1[1],col2[1]), bty="n", cex = 3, title=NULL, yjust=0.5)

  if (nchr == 1) {
    axis(1, cex.axis = 1.5, cex.lab = 1.5, las=1, ...)
    axis(2, at = yticks, labels = ylabels, las=1, cex.axis=1.5, cex.lab = 1.5, ...)
  } else {
    axis(1, at = ticks, labels = labs, las=1, cex.axis = 1.5, cex.lab = 1.5, ...)
    axis(2, at = yticks, labels = ylabels, las=1, cex.axis = 1.5, cex.lab = 1.5, ...)
  }
  col1 = rep(col1, max(d$CHR))
  col2 = rep(col2, max(d$CHR))
  
  d1 <- d[d$logp>0,]
  d2 <- d[d$logp<0,]

  if (nchr == 1) {
    with(d1[d1$logp<ymax,], points(pos, logp, pch = 20, col = col1[1], ...))
    with(d1[d1$logp>ymax,], points(pos, rep(ymax+0.5, length(pos)), pch = 2, col = col1[1], ...))
    with(d2[d$logp>ymin,], points(pos, logp, pch = 20, col = col2[1], ...))
    with(d2[d2$logp<ymin,], points(pos, rep(ymin-0.5, length(pos)), pch = 6, col = col2[1], ...))
  } else {
    print("##PLOT X")
    icol = 1
    for (i in unique(d1$index)) {
      
      if(!is.null(highlight)) {
        with(d1[d1$index == unique(d1$index)[i] & d1$logp<ymax & d1$highlight == 0, ], points(pos, logp, col = col1[icol], pch = 20, ...))  
        with(d1[d1$index == unique(d1$index)[i] & d1$logp<ymax & d1$highlight == 1, ], points(pos, logp, col = highlight_col, pch = 4, cex= 1.5, lwd= 3, ...))
        with(d1[d1$index == unique(d1$index)[i] & d1$logp<ymax & d1$highlight == 2, ], points(pos, logp, col = highlight_col, bg=highlight_col, pch = 24, cex= 1.5, lwd= 2, ...))
        with(d1[d1$index == unique(d1$index)[i] & d1$logp>ymax,], points(pos, rep(ymax+0.5, length(pos)), col = highlight_col, bg=highlight_col, pch = 24, cex=1.5, lwd = 2, ...))
      } else {
        with(d1[d1$index == unique(d1$index)[i] & d1$logp<ymax, ], points(pos, logp, col = col1[icol], pch = 20, ...))  
        with(d1[d1$index == unique(d1$index)[i] & d1$logp>ymax,], points(pos, rep(ymax+0.5, length(pos)), pch = 24, col = col1[icol], bg=col1[icol], cex=1.5, lwd = 2, ...))
      }
      
      
      icol = icol + 1
    }
    
    print("##PLOT Y")
    if(!is.null(y)) {

      icol = 1
      for (i in unique(d2$index)) {
      
      if(!is.null(highlight)) {
        with(d2[d2$index == unique(d2$index)[i] & d2$logp>ymin & d2$highlight == 0, ], points(pos, logp, col = col2[icol], pch = 20, ...))
        with(d2[d2$index == unique(d2$index)[i] & d2$logp<ymax & d2$highlight == 1, ], points(pos, logp, col = highlight_col, pch = 4, cex = 1.5, lwd=3, ...))
        with(d2[d2$index == unique(d2$index)[i] & d2$logp<ymax & d2$highlight == 2, ], points(pos, logp, col = highlight_col, bg=highlight_col, pch = 25, cex = 1.5, lwd=2, ...))
        with(d2[d2$index == unique(d2$index)[i] & d2$logp<ymin,], points(pos, rep(ymin-0.5, length(pos)), col = highlight_col, bg=highlight_col, pch = 25, cex=1.5, lwd=2, ...))
      } else {
        with(d2[d2$index == unique(d2$index)[i] & d2$logp>ymin, ], points(pos, logp, col = col2[icol], pch = 20, ...))
        with(d2[d2$index == unique(d2$index)[i] & d2$logp<ymin,], points(pos, rep(ymin-0.5, length(pos)), pch = 25, col = col2[icol], bg=col2[icol], cex=1.5, lwd=2, ...))
      }
      
      
      icol = icol + 1
      }
    }
  }



  #----
  # Plot names
  #----
  
  if (!is.null(label)) {
    required_cols <- c("LABEL","CHR","POS","P")
    if (is.null(colnames(label))) warning(paste("label must be a data.frame with col.names:",paste0(required_cols,collapse=" ")))
    label <- data.frame(label)
    
    if (all(required_cols %in% colnames(label))) {  
      label$pos <- label$POS
      for (i in unique(label$CHR)[!unique(label$CHR)==1]) {
        lastbase <- tail(d[d$index==(i-1),"pos"],n=1)
        label[label$CHR==i,"pos"] <- label[label$CHR==i,"pos"]+lastbase
      }
      
      #with(label, points(pos, P, col = "#974d26", pch = 4, cex=2, ...))

      label$P <- ifelse(label$P>0, label$P-ymid+2*strheight(label$LABEL), label$P+ymid-2*strheight(label$LABEL))
      label$P <- ifelse(label$P>ymax,ymax-strheight(label$LABEL)/2,label$P)
      label$P <- ifelse(label$P<ymin,ymin+strheight(label$LABEL)/2,label$P)
      label <- label[order(-sign(label$P),abs(label$P),label$pos),]
      row.names(label) <- 1:nrow(label)

      print("##POSITION LABELS")
      
      x.buffer <- strwidth("A")
      y.buffer <- strheight("A")

      label$YMIN <- label$P - strheight(label$LABEL)/2
      label$YMAX <- label$P + strheight(label$LABEL)/2

      label$XMIN <- label$pos - strwidth(label$LABEL)/2 - x.buffer
      label$XMAX <- label$pos + strwidth(label$LABEL)/2 + x.buffer
      
      forbidden <- data.frame(data.table(d)[highlight==TRUE,.(pos,P=logp,XMIN=pos-strwidth("A")/2, XMAX=pos+strwidth("A")/2, YMIN=ifelse(logp<ymin,ymin,ifelse(logp>ymax,ymax,logp - sign(logp) * strheight("A")/2)), YMAX=ifelse(logp>ymax, ymax, ifelse(logp<ymin,ymin,logp + sign(logp) * strheight("A")/2)))])

      #ggplot(label, aes(x=pos, xmin=XMIN, xmax=XMAX, y=(YMIN+YMAX)/2, ymin=YMIN, ymax=YMAX)) + geom_point(data=forbidden, aes(x=pos, y=YMAX), colour="red") + geom_rect() + geom_text(aes(label=LABEL), size=4, color="white") 

      rect_overlap <- function(df1, df2) {
        return(!(df2$XMAX <= df1$XMIN | df2$XMIN >= df1$XMAX | df2$YMAX <= df1$YMIN | df2$YMIN >= df1$YMAX))
      }


      for (precision in 1:10) {
      
      for (i in 1:nrow(label)) {
      

      n <- 0
      repeat{
      overlap_highlight <- forbidden[rect_overlap(df1=label[i,], df2=forbidden),]
      
      if(nrow(overlap_highlight) > 0) {
        if(all(overlap_highlight$P > 0)) {
          label[i,"YMAX"] <- max(overlap_highlight$YMAX) + strheight("A")/2*3 + y.buffer
          label[i,"YMIN"] <- max(overlap_highlight$YMAX) + strheight("A")/2 + y.buffer
        } else {
          label[i,"YMIN"] <- min(overlap_highlight$YMIN) - strheight("A")/2*3 - y.buffer
          label[i,"YMAX"] <- min(overlap_highlight$YMIN) - strheight("A")/2 - y.buffer
        }
      }

      overlap <- rbind(label[i,], label[-i,][rect_overlap(df1=label[i,], df2=label[-i,]),])
      

      if(nrow(overlap) == 1 ) break

      n <- n+1
      if(n == 1) cat("\n",label[i,"LABEL"],"")
      cat(n,"")
      

      if(all(overlap$P > 0)) {
        overlap[1,"YMAX"] <- max(overlap[-1,"YMAX"]) + strheight("A")/2*3 + y.buffer
        overlap[1,"YMIN"] <- max(overlap[-1,"YMAX"]) + strheight("A")/2 + y.buffer
      } else {
        overlap[1,"YMIN"] <- min(overlap[-1,"YMIN"]) - strheight("A")/2*3 - y.buffer
        overlap[1,"YMAX"] <- min(overlap[-1,"YMIN"]) - strheight("A")/2 - y.buffer
      }

      label[match(paste0(overlap$LABEL,overlap$P),paste0(label$LABEL,label$P)),] <- overlap
      }
      }
    }
    cat("\n")

      label$P <- ifelse(label$P>ymax,ymax-strheight(label$LABEL)/6,label$P)
      label$P <- ifelse(label$P<ymin,ymin+strheight(label$LABEL)/6,label$P)
      moved <- data.table(label)[P!=(YMIN+YMAX)/2, .(x=pos, y1=P, y2=(YMIN+YMAX)/2)]

      for(i in 1:nrow(moved)) lines(x=moved[i,x] + c(0,0), y=c(moved[i,y1], moved[i,y2]) - sign(moved[i,y1]) * strheight("A"))


      shadowtext <- function(x, y=NULL, labels, col='white', bg='black', 
                       theta= seq(0, 2*pi, length.out=50), r=0.1, ... ) {
        xy <- xy.coords(x,y)
        xo <- r*strwidth('A')
        yo <- r*strheight('A')
      # draw background text with small shift in x and y in background colour
      for (i in theta) {
        text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
      }
        # draw actual text in exact xy position in foreground colour
        text(xy$x, xy$y, labels, col=col, ... )
      }
      
      print("##DRAW LABELS")
      with(label, shadowtext(x=pos,y=(YMIN+YMAX)/2,labels=LABEL,cex=1.25,col="black",bg="white",font=4,r=0.2))

     } else {
      warning(paste("label data.frame is missing columns:",paste0(required_cols[!required_cols%in%colnames(label)],collapse=" ")))
     }
   }

  par(xpd = FALSE)
}



#----
# Load data
#----

writeLines("Loading data...")
gwas_data <- list()
hits <- list()

gwas_data[[1]] <- fread("zcat ../q02_manova/st02_01_triple_manova.tsv.gz | awk 'NR == 1 || $11 < 0.1'", select=c("rsid","chr","pos","p"))
hits[[1]] <- fread('st03_03_triple_manova_res.tsv')[order(chr,p)]

gwas_data[[2]] <- fread("zcat ../q01_sumstats/st01_01_healthspan.tsv.gz | awk 'NR==1 || $11 < 0.1'")
hits[[2]] <- fread('st03_03_healthspan_res.tsv')[order(chr,p)]

gwas_data[[3]] <- fread("cat ../../../prj_065b_lifegen_phase2/p007_gather_res/bothpl_alldr.tsv| awk 'NR==1 || $11 < 0.1'")
hits[[3]] <- fread('st03_03_lifespan_res.tsv')[order(chr,p)]

gwas_data[[4]] <- fread("zcat ../q01_sumstats/st01_01_longevity_90_pct.tsv.gz | awk 'NR==1 || $11 < 0.1'")
hits[[4]] <- fread("st03_03_longevity_res.tsv")[order(chr,p)]


#----
# Add highlights
#----

writeLines("Adding highlights...")

labels <- list()

for (n in 1:4) {
  for (i in 2:nrow(hits[[n]])) if(abs(hits[[n]][i,pos]-hits[[n]][i-1,pos]) < 250000) hits[[n]][i,p:=NA]
  hits[[n]] <- hits[[n]][!is.na(p),]
  labels[[n]] <- hits[[n]][p<5e-8,.(summary,chr,pos,p)]
  

  colnames(labels[[n]]) <- c("LABEL","CHR","POS","P")
  labels[[n]]$P <- -log10(labels[[n]]$P)
  labels[[n]]$LABEL <- gsub("/.*$","",labels[[n]]$LABEL)
  labels[[n]]$LABEL <- gsub("HYKK","CHRNA3/5",labels[[n]]$LABEL)
  
  hit_rsids <- hits[[n]]$rsid
  gwas_data[[n]][,highlight:=ifelse(rsid%in%hit_rsids, TRUE, FALSE)]
  labels[[n]] <- labels[[n]][!duplicated(LABEL),]
}


green <-c("#006d2c", "#238b45", "#41ab5d")
blue <- c("#2166ac","#2171b5","#4393c3")
purple <- c("#542788","#6a51a3","#8073ac")
black <- c("#252525", "#525252", "#737373")


#----
# Create plots
#----


writeLines("Plotting...")

png("st03_04_mh_triple_manova.png", width = 1900, height = 720)
par(mar=c(5.1,5.1,4.1,2.1))

miami(gwas_data[[1]],
  chr = "chr", bp = "pos", snp = "rsid", p = "p", 
  genomewideline = -log10(5e-8), ymin = 0, ymax = 29, 
  col1 = black, 
  suggestiveline = FALSE, ymid=0, 
  label=labels[[1]], highlight="highlight", 
  highlight_col="#cb181d")
dev.off()



png("st03_04_mh_healthspan_manova.png", width = 1900, height = 540)
par(mar=c(5.1,5.1,4.1,2.1))


png("st03_04_mh_healthspan.png", width = 1900, height = 540)
miami(gwas_data[[2]],
  chr = "chr", bp = "pos", snp = "rsid", p = "p", 
  genomewideline = -log10(5e-8), ymin = 0, ymax = 29, 
  col1 = green, #col1 = c("#9ecae1","#c6dbef"),
  suggestiveline = FALSE, ymid=0, 
  label=labels[[2]], highlight="highlight", 
  highlight_col="#cb181d")
dev.off()


png("st03_04_mh_lifespan.png", width = 1900, height = 540)
miami(gwas_data[[3]],
  chr = "chr", bp = "pos", snp = "rsid", p = "p", 
  genomewideline = -log10(5e-8), ymin = 0, ymax = 29, 
  col1 = blue, #col1 = c("#9ecae1","#c6dbef"),
  suggestiveline = FALSE, ymid=0, 
  label=labels[[3]], highlight="highlight", 
  highlight_col="#cb181d")
dev.off()

png("st03_04_mh_longevity.png", width = 1900, height = 540)
miami(gwas_data[[4]],
  chr = "chr", bp = "pos", snp = "rsid", p = "p", 
  genomewideline = -log10(5e-8), ymin = 0, ymax = 29, 
  col1 = purple, #col1 = c("#9ecae1","#c6dbef"),
  suggestiveline = FALSE, ymid=0, 
  label=labels[[4]], highlight="highlight", 
  highlight_col="#cb181d")
dev.off()





#----
# Publication
#----

setEPS()
postscript("st03_04_mh_triple_manova_publication.eps", width=19, height=7.2)
par(mar=c(5.1,5.1,4.1,2.1))
miami(gwas_data[[1]],
  chr = "chr", bp = "pos", snp = "rsid", p = "p", 
  genomewideline = -log10(5e-8), ymin = 0, ymax = 29, 
  col1 = black,
  suggestiveline = FALSE, ymid=0, 
  label=labels[[1]], highlight="highlight", 
  highlight_col="#cb181d")
dev.off()

# Open in Preview and export to .PNG at 330 dpi


