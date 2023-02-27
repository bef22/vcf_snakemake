### assess depth and missing individual data
#
## imiss
# - stripchart of missing per chromosome/scaffold
# - barplot of overall missing per individual
# - outlier summary, outliers per chromosome which are outside 1.5 times the interquartile range above the upper quartile
# - summary matrix of all missing proportions per individual and chromosome/scaffold
#
## depth (=DP)
# - histogram of DP distribution
# - depth summary
#
# Bettina Fischer, 20230126
###############################################################################

#Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
filePath <- args[1]
namePrefix <- args[2]
chrListFile <- args[3]
scaffoldListFile <- args[4]
type <- args[5]

if(type == "raw_imiss") {
  
  outlierFile <- paste0(filePath, namePrefix, "_raw_missing_individual_outlier_report.txt")
  summaryFile <- paste0(filePath, namePrefix, "_raw_missing_individual_summary_report.txt")
  plotFile <- paste0(filePath, namePrefix, "_raw_missing_individual.pdf")
  barplotFile <- paste0(filePath, namePrefix, "_raw_missing_individual_barplot.pdf")
  
  # read the main chromosome names
  chrList <- read.delim(chrListFile, sep="\t", header=F, stringsAsFactors = F)[,1]
  
  mList <- list()
  x.outliers <- c()
  miss.df <- c()
  nmiss.df <- c()
  nChrom <- c()
  
  # read each file
  for(i in 1:length(chrList)) {
    
    mFile <- paste0(filePath, namePrefix, ".raw.", chrList[i], ".imiss")
    x <- read.delim(mFile, sep="\t")
    
    # get outlier
    bx <- boxplot(x$F_MISS, plot = F)
    x.upper <- bx$stats[5,1]
    
    # also report any above 0.9
    if(x.upper > 0.9) {
      x.upper <- 0.9
    }
    
    x.out <- x[x$F_MISS > x.upper,]
    if(nrow(x.out) > 0) {
      x.out$CHROM <- chrList[i]
      x.outliers <- rbind(x.outliers, x.out)
    }
    
    if(i == 1) {
      nChrom <- x$N_DATA[1]
      mm <- x
      miss.df <- x[,c(1,5)]
      nmiss.df <- x[,c(1,4)]
    }
    if(i > 1) {
      nChrom <- c(nChrom, x$N_DATA[1])
      mm <- rbind(mm, x)
      miss.df$tmp <- x[,5]
      nmiss.df$tmp <- x[,4]
    }
    colnames(miss.df)[i+1] <- chrList[i]
    colnames(nmiss.df)[i+1] <- chrList[i]
    mList[[i]] <- x$F_MISS
    names(mList)[i] <- chrList[i]
  }
  
  
  
  # multipage pdf with max 25 per row
  pdf(plotFile, height=8.27, width=11.69, onefile=TRUE)
  max.n <- 25
  max.plotrows <- ceiling(length(mList)/max.n)
  par(mar=c(15,4,2,1))
  if(length(mList) > max.n) {
    for(j in 1:max.plotrows) {
      mListj <- mList[((max.n*j)-max.n+1):(max.n*j)]
      stripchart(mListj, pch=16, cex=0.5, method = "jitter", jitter = 0.3, vertical = T, ylim=c(0,1), las=2,
                 ylab="F_MISS", main="missing individual")
    }
  }
  if(length(mList) <= max.n) {
    stripchart(mList, pch=16, cex=0.5, method = "jitter", jitter = 0.3, vertical = T, ylim=c(0,1), las=2,
               ylab="F_MISS", main="missing individual")
  }
  invisible(dev.off())
  

  write.table(x.outliers, file=outlierFile, sep="\t", row.names=F, quote=F)
  write.table(miss.df, file=summaryFile, sep="\t", row.names=F, quote=F)
  
  
  # plot the overall missing per individual proportion across all chr/scaffolds
  mt <- nmiss.df
  rownames(mt) <- nmiss.df[,1]
  mt[,1] <- NULL
  mt <- t(mt)
  cs <- colSums(mt)
  
  cs <- cs/sum(nChrom)
  max.missing <- max(cs) + 0.01
  

  # multipage pdf with max 50 per row
  pdf(barplotFile, height=8.27, width=11.69, onefile=TRUE)
  max.n <- 55
  max.plotrows <- ceiling(length(mList)/max.n)
  par(mar=c(15,4,2,1))
  
  if(length(cs) > max.n) {
    for(j in 1:max.plotrows) {
      csj <- cs[((max.n*j)-max.n+1):(max.n*j)]
      barplot(csj, ylim=c(0,max.missing), las=2, col="grey", ylab="overall missing")
    }
  }
  if(length(cs) <= max.n) {
    barplot(cs, ylim=c(0,max.missing),las=2, col="grey", ylab="overall missing")
  }
  
  invisible(dev.off())
  
}


if(type == "imiss") {
  
  outlierFile <- paste0(filePath, namePrefix, "_missing_individual_outlier_report.txt")
  summaryFile <- paste0(filePath, namePrefix, "_missing_individual_summary_report.txt")
  plotFile <- paste0(filePath, namePrefix, "_missing_individual.pdf")
  barplotFile <- paste0(filePath, namePrefix, "_missing_individual_barplot.pdf")
  
  # read the main chromosome names
  chrList <- read.delim(chrListFile, sep="\t", header=F, stringsAsFactors = F)[,1]
  
  if(scaffoldListFile != "no") {
    chrList <- c(chrList, "Scaffolds")
  }
  
  mList <- list()
  x.outliers <- c()
  miss.df <- c()
  nmiss.df <- c()
  nChrom <- c()
  
  # read each file
  for(i in 1:length(chrList)) {
    
    mFile <- paste0(filePath, namePrefix, ".", chrList[i], ".imiss")
    x <- read.delim(mFile, sep="\t")
    
    # get outlier
    bx <- boxplot(x$F_MISS, plot = F)
    x.upper <- bx$stats[5,1]
    
    # also report any above 0.9
    if(x.upper > 0.9) {
      x.upper <- 0.9
    }
    
    x.out <- x[x$F_MISS > x.upper,]
    if(nrow(x.out) > 0) {
      x.out$CHROM <- chrList[i]
      x.outliers <- rbind(x.outliers, x.out)
    }
    
    if(i == 1) {
      nChrom <- x$N_DATA[1]
      mm <- x
      miss.df <- x[,c(1,5)]
      nmiss.df <- x[,c(1,4)]
    }
    if(i > 1) {
      nChrom <- c(nChrom, x$N_DATA[1])
      mm <- rbind(mm, x)
      miss.df$tmp <- x[,5]
      nmiss.df$tmp <- x[,4]
    }
    colnames(miss.df)[i+1] <- chrList[i]
    colnames(nmiss.df)[i+1] <- chrList[i]
    mList[[i]] <- x$F_MISS
    names(mList)[i] <- chrList[i]
  }
  
  
  
  # multipage pdf with max 25 per row
  pdf(plotFile, height=8.27, width=11.69, onefile=TRUE)
  max.n <- 25
  max.plotrows <- ceiling(length(mList)/max.n)
  par(mar=c(15,4,2,1))
  if(length(mList) > max.n) {
    for(j in 1:max.plotrows) {
      mListj <- mList[((max.n*j)-max.n+1):(max.n*j)]
      stripchart(mListj, pch=16, cex=0.5, method = "jitter", jitter = 0.3, vertical = T, ylim=c(0,1), las=2,
                 ylab="F_MISS", main="missing individual")
    }
  }
  if(length(mList) <= max.n) {
    stripchart(mList, pch=16, cex=0.5, method = "jitter", jitter = 0.3, vertical = T, ylim=c(0,1), las=2,
               ylab="F_MISS", main="missing individual")
  }
  invisible(dev.off())
  
  
  write.table(x.outliers, file=outlierFile, sep="\t", row.names=F, quote=F)
  write.table(miss.df, file=summaryFile, sep="\t", row.names=F, quote=F)
  
  
  # plot the overall missing per individual proportion across all chr/scaffolds
  mt <- nmiss.df
  rownames(mt) <- nmiss.df[,1]
  mt[,1] <- NULL
  mt <- t(mt)
  cs <- colSums(mt)
  
  cs <- cs/sum(nChrom)
  max.missing <- max(cs) + 0.01
  
  
  # multipage pdf with max 50 per row
  pdf(barplotFile, height=8.27, width=11.69, onefile=TRUE)
  max.n <- 55
  max.plotrows <- ceiling(length(mList)/max.n)
  par(mar=c(15,4,2,1))
  
  if(length(cs) > max.n) {
    for(j in 1:max.plotrows) {
      csj <- cs[((max.n*j)-max.n+1):(max.n*j)]
      barplot(csj, ylim=c(0,max.missing), las=2, col="grey", ylab="overall missing")
    }
  }
  if(length(cs) <= max.n) {
    barplot(cs, ylim=c(0,max.missing),las=2, col="grey", ylab="overall missing")
  }
  
  invisible(dev.off())
  
}




if(type == "depth") {
  
  # read the main chromosome names
  chrList <- read.delim(chrListFile, sep="\t", header=F, stringsAsFactors = F)
  chrs <- chrList[,1]
  
  # get the unique depth values for all chromosomes
  for(i in 1:length(chrs)) {
    mFile <- paste0(filePath, namePrefix, ".", chrs[i], "_freq.txt.gz")
    x <- read.delim(mFile, sep="", header=F)
    x2 <- data.frame(Val = as.integer(x[,2]),
                     N = as.integer(x[,1]))
    if(i == 1) {
      dpvals <- x2$Val
    }
    if(i > 1) {
      dpvals <- unique(c(dpvals, x2$Val))
    }
    rm(x, x2)
  }
  
  dpvals <- sort(dpvals)
  cdf <- data.frame(Val = dpvals)
  rm(dpvals)
  
  # sum the counts for each dpval
  for(i in 1:length(chrs)) {
    mFile <- paste0(filePath, namePrefix, ".", chrs[i], "_freq.txt.gz")
    cat("Reading", mFile, "\n")
    
    # read the frequency tables and merge across all chromosomes
    x <- read.delim(mFile, sep="", header=F)
    x2 <- data.frame(Val = as.integer(x[,2]),
                     N = as.integer(x[,1]))
    
    if(i == 1) {
      cdf <- merge(cdf, x2, by="Val", all.x=T)
      colnames(cdf)[2] <- "Sum"
      cdf[is.na(cdf[,2]),2] <- 0
    }
    
    if(i > 1) {
      tmp <- merge(cdf, x2, by="Val", all.x=T)
      tmp$N[is.na(tmp$N)] <- 0
      tmp$Sum <- tmp$Sum + tmp$N
      cdf <- tmp[,c(1:2)]
      rm(tmp)
    }
    rm(x, x2)
  }
  
  # remove the 0-counts
  cdf <- cdf[cdf$Val > 0,]
  
  # create vector of depths (this is memory intensive)
  depths <- rep(cdf$Val, cdf$Sum)
  
  
  # median, mean, mode across sites
  median.depth <- median(depths)
  mean.depth <- mean(depths)
  
  # calculate the mode, but exclude very low values
  getMode <- function(x, min.cutoff=10) {
    ux <- unique(x[x > min.cutoff])
    ux[which.max(tabulate(match(x, ux)))]
  }
  Mode <- getMode(depths)
  
  summary.depth <- as.matrix(summary(depths))
  summary.depth <- rbind(summary.depth, Mode)
  
  # set nice break values for histogram
  max.depth <- max(depths)
  breaks <- as.integer(max.depth/50)
  
  # plot
  plotName <- paste0(filePath, "depth_histogram.png")
  png(plotName, width=1200)
  hist(depths, breaks=breaks, xlim=c(0, 2*summary.depth[5,1]), main="Depth", xlab="depth")
  abline(v=median.depth, lwd=2, lty=2, col="red")
  abline(v=mean.depth, lwd=2, lty=2, col="blue")
  abline(v=Mode, lwd=2, lty=2, col="forestgreen")
  abline(v=Mode-(Mode*0.25), lwd=1, lty=2)
  abline(v=Mode+(Mode*0.25), lwd=1, lty=2)
  abline(v=Mode-(Mode*0.35), lwd=1, lty=3)
  abline(v=Mode+(Mode*0.35), lwd=1, lty=3)
  
  legend("topright", legend=c("median", "mean", "mode", "+/-25% mode", "+/-35% mode"), lwd=2, lty=c(1,1,1,2,3), 
         col=c("red", "blue", "forestgreen", "black", "black"), bty="n")
  invisible(dev.off())
  
  outFile <- paste0(filePath, "depth_summary.txt")
  write.table(summary.depth, file=outFile, sep="\t", col.names=F,  quote=F)
  
  
}




