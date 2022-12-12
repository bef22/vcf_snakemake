### report missing individuals
# the script reports missing individual outliers which are outside 1.5 times the interquartile range above the upper quartile
# 
#
# 20221129, Bettina Fischer
# 20221212, add missing barplot per individual, changed the number of plotrows for scaffolds
###############################################################################

#Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
filePath <- args[1]
namePrefix <- args[2]
chrListFile <- args[3]
scaffoldListFile <- args[4]
type <- args[5]

### missing individual
if(type == "individual") {

  outlierFile <- paste0(filePath, namePrefix, "_missing_individual_outlier_report.txt")
  summaryFile <- paste0(filePath, namePrefix, "_missing_individual_summary_report.txt")
  plotMainFile <- paste0(filePath, namePrefix, "_mainChrom_missing_individual.png")
  barplotFile <- paste0(filePath, namePrefix, "_sum_missing_individual_barplot.png")
  
  # read the main chromosome names
  chrList <- read.delim(chrListFile, sep="\t", header=F, stringsAsFactors = F)
  
  # create the list of files in chromosome order
  mFiles <- paste0(filePath, namePrefix, ".", chrList[,1], ".imiss")
  
  mList <- list()
  x.outliers <- c()
  
  # rad each file
  for(i in 1:length(mFiles)) {
    x <- read.delim(mFiles[i], sep="\t")
    
    # get outliers
    bx <- boxplot(x$F_MISS, plot = F)
    x.upper <- bx$stats[5,1]
    
    x.out <- x[x$F_MISS > x.upper,]
    if(nrow(x.out) > 0) {
      x.out$CHROM <- chrList[i,1]
      x.outliers <- rbind(x.outliers, x.out)
    }
    
    if(i == 1) {
      mm <- x
    }
    if(i > 1) {
      mm <- rbind(mm, x)
    }
    mList[[i]] <- x$F_MISS
    names(mList)[i] <- chrList[i,1]
  }
  
  png(plotMainFile, width=1000, height=400)
  stripchart(mList, pch=16, cex=0.5, method = "jitter", jitter = 0.3, vertical = T, ylim=c(0,1), las=2,
             ylab="F_MISS", main="missing individual")
  invisible(dev.off())
  
  
  if(scaffoldListFile != "no") {
    
    plotScaffoldFile <- paste0(filePath, namePrefix, "_Scaffolds_missing_individual.png")
    
    # read the main chromosome names
    scaffoldList <- read.delim(scaffoldListFile, sep="\t", header=F, stringsAsFactors = F)
    
    # create the list of files in chromosome order
    sFiles <- paste0(filePath, namePrefix, ".", scaffoldList[,1], ".imiss")
    
    sList <- list()
    s.outliers <- c()
    
    # rad each file
    for(i in 1:length(sFiles)) {
      x <- read.delim(sFiles[i], sep="\t")
      
      # get outliers
      bx <- boxplot(x$F_MISS, plot = F)
      x.upper <- bx$stats[5,1]
      
      x.out <- x[x$F_MISS > x.upper,]
      if(nrow(x.out) > 0) {
        x.out$CHROM <- scaffoldList[i,1]
        s.outliers <- rbind(s.outliers, x.out)
      }
      
      if(i == 1) {
        mm <- x
      }
      if(i > 1) {
        mm <- rbind(mm, x)
      }
      mList[[i]] <- x$F_MISS
      names(mList)[i] <- scaffoldList[i,1]
    }
    
    # plot max 100 per row
    max.n <- 100
    max.plotrows <- ceiling(length(mList)/max.n)
    plot.width <- 1200
    plot.height <- 300 * max.plotrows
    
    png(plotScaffoldFile, width=plot.width, height=plot.height)
    par(mfrow=c(max.plotrows, 1), mar=c(10,4,2,1))
    for(j in 1:max.plotrows) {
      mListj <- mList[((max.n*j)-max.n+1):(max.n*j)]
      stripchart(mListj, pch=16, cex=0.5, method = "jitter", jitter = 0.3, vertical = T, ylim=c(0,1), las=2,
                 ylab="F_MISS")
    }
    invisible(dev.off())
    
    x.outliers <- rbind(x.outliers, s.outliers)
    
  }

  write.table(x.outliers, file=outlierFile, sep="\t", row.names=F, quote=F)

  # plot the sum of missing per individual across all chr/scaffolds
  mt <- miss.df
  rownames(mt) <- miss.df[,1]
  mt[,1] <- NULL
  mt <- t(mt)
  cs <- colSums(mt)
  
  # plot max 50 per row
  max.n <- 50
  max.plotrows <- ceiling(length(cs)/max.n)
  plot.width <- 1200
  plot.height <- 300 * max.plotrows
  
  png(barplotFile, width=plot.width, height=plot.height)
  par(mfrow=c(max.plotrows, 1), mar=c(10,4,2,1))
  if(length(cs) > max.n) {
    for(j in 1:max.plotrows) {
      csj <- cs[((max.n*j)-max.n+1):(max.n*j)]
      barplot(csj, las=2, col="grey", ylab="sum missing")
    }
  }
  if(length(cs) <= max.n) {
    barplot(cs, las=2, col="grey", ylab="sum missing")
  }
  invisible(dev.off())

}

