### report missing individuals
# the script reports missing individual outliers which are outside 1.5 times the interquartile range above the upper quartile
# 
#
# 20221129, Bettina Fischer
# needs more work for plotting scaffolds... and a usage
###############################################################################

#Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
filePath <- args[1]
namePrefix <- args[2]
chrListFile <- args[3]
scaffoldListFile <- args[4]
type <- args[5]

#filePath <- "../bcf_qc/missing_individual/"
#type <- "individual"              # individual or site
#namePrefix <- "zebra_callainos_pearly.MZeb"
#chrListFile <- "../chromList.txt"
#scaffoldListFile <- "../scaffoldList.txt"



### missing individual
if(type == "individual") {

  outlierFile <- paste0(filePath, namePrefix, "_missing_individual_report.txt")
  plotMainFile <- paste0(filePath, namePrefix, "_mainChrom_missing_individual.png")
  
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
    
    png(plotScaffoldFile, width=6000, height=400)
    stripchart(mList, pch=16, cex=0.5, method = "jitter", jitter = 0.3, vertical = T, ylim=c(0,1), las=2,
               ylab="F_MISS", main="missing individual")
    invisible(dev.off())
    
    x.outliers <- rbind(x.outliers, s.outliers)
    
  }

  write.table(x.outliers, file=outlierFile, sep="\t", row.names=F, quote=F)

}

