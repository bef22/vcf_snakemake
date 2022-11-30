# calculate the genome-wide median depth
#
###############################################################################


#Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
filePath <- args[1]
namePrefix <- args[2]
chrListFile <- args[3]


# read the main chromosome names
chrList <- read.delim(chrListFile, sep="\t", header=F, stringsAsFactors = F)
chrs <- chrList[,1]


depths <- c()

for(i in 1:length(chrs)) {
  mFile <- paste0(filePath, namePrefix, ".", chrs[i], ".ldepth.gz")
  cat("Reading", mFile, "\n")
  x <- read.delim(mFile, sep="\t", header=T)
  depths <- c(depths, x$SUM_DEPTH)
}


# median across sites
median.depth <- median(depths)
mean.depth <- mean(depths)

summary.depth <- as.matrix(summary(depths))

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

