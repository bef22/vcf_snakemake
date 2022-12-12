# calculate the genome-wide median depth
# 20221212 - now using frequency table
###############################################################################


#Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
filePath <- args[1]
namePrefix <- args[2]
chrListFile <- args[3]

# read the main chromosome names
chrList <- read.delim(chrListFile, sep="\t", header=F, stringsAsFactors = F)
chrs <- chrList[,1]

for(i in 1:length(chrs)) {
  mFile <- paste0(filePath, namePrefix, ".", chrs[i], ".ldepth.gz")
  cat("Reading", mFile, "\n")
  
  
  # read the frequency tables and merge acrcross all chromosomes
  x <- read.delim(mFile, sep="", header=F)
  x2 <- data.frame(Val = x[,2],
                   N = x[,1])
  
  if(i == 1) {
    cdf <- data.frame(Val = c(0:max(x2[,1])))
    cdf <- merge(cdf, x2, by="Val", all.x=T)
    colnames(cdf)[2] <- "Sum"
    cdf[is.na(cdf[,2]),2] <- 0
  }
  
  if(i > 1) {
    if(max(x2[,1]) > nrow(cdf)) {
      tmp2 <- data.frame(Val = c(0:max(x2[,1])))
      tmp2 <- merge(tmp2, cdf, by="Val", all.x=T)
      tmp2$Sum[is.na(tmp2$Sum)] <- 0
      tmp <- merge(tmp2, x2, by="Val", all.x=T)
      tmp$N[is.na(tmp$N)] <- 0
      tmp$Sum <- tmp$Sum + tmp$N
      cdf <- tmp[,c(1:2)]
      rm(tmp2)
      rm(tmp)
    }
    
    if(max(x2[,1]) <= nrow(cdf)) {
      tmp <- merge(cdf, x2, by="Val", all.x=T)
      tmp$N[is.na(tmp$N)] <- 0
      tmp$Sum <- tmp$Sum + tmp$N
      cdf <- tmp[,c(1:2)]
      rm(tmp)
    }
  }
  rm(x2)
  rm(x)
}


# remove the 0-counts
cdf <- cdf[cdf$Sum > 0,]
cdf <- cdf[cdf$Val > 0,]

# create vector of depths
depths <- rep(cdf$Val, cdf$Sum)

# median and mean across sites
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

# plot
plotName <- paste0(filePath, "site_depth_histogram.png")
png(plotName)
hist(depths, breaks=5000, xlim=c(0, 2*summary.depth[5,1]), main="Sum depth at sites", xlab="depth")
abline(v=median.depth, lwd=2, col="red")
abline(v=mean.depth, lwd=2, col="blue")
abline(v=median.depth-(median.depth*0.25), lwd=1, lty=2)
abline(v=median.depth+(median.depth*0.25), lwd=1, lty=2)
abline(v=median.depth-(median.depth*0.35), lwd=1, lty=3)
abline(v=median.depth+(median.depth*0.35), lwd=1, lty=3)

legend("topright", legend=c("median", "mean", "+/-25%", "+/-35%"), lwd=2, lty=c(1,1,2,3), col=c("red", "blue", "black", "black"), bty="n")
dev.off()

outFile <- paste0(filePath, "summary_site_depth.txt")
write.table(summary.depth, file=outFile, sep="\t", col.names=F,  quote=F)

