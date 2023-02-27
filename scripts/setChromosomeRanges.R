# split chromosomes (and scaffolds) into ranges for vcf calling
#
# this also creates the list files for concatenating these back to chromosome level
# scaffolds will be merged into one file
###############################################################################

#Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
namePrefix <- args[1]
chrListFile <- args[2]
scaffoldListFile <- args[3]
faiFile <- args[4]
chunkSize <- as.integer(args[5])

#namePrefix <- "test"
#chrListFile <- "../chromList.txt"
#scaffoldListFile <- "../scaffoldList.txt"
#faiFile <- "../../../gtf_genomes_bams/fAstCal1.2_chromnames/GCA_900246225.3_fAstCal1.2_genomic_chromnames_mt.fa.fai"
#chunkSize <- 10000000




rawDir <- "bcf_raw"
callDir <- "bcf_call"


###############################################################################

if(!dir.exists(rawDir)) {
  dir.create(rawDir)
}
if(!dir.exists(callDir)) {
  dir.create(callDir)
}



### read fai
fai <- read.delim(faiFile, sep="\t", header=F, stringsAsFactors=F)
rownames(fai) <- fai[,1]


### read chromosomes (and scaffold) names to include
chrs <- read.delim(chrListFile, sep="\t", header=F, stringsAsFactors=F)[,1]


### per chromosome ranges
regions <- c()
for(i in 1:length(chrs)) {
  starts <- ends <- chunks <- NULL
  len <- fai[chrs[i], 2]
  starts <- as.integer(seq(1, len, by=chunkSize))
  ends <- starts[-1] -1
  ends <- as.integer(c(ends, len))
  chunks <- paste0(chrs[i], "_", starts, "-", ends)
  regions <- c(regions, chunks)
  
  ### write per chromosome bcf lists for concat step
  outName <- paste0(callDir, "/", chrs[i], "_list.txt")
  bcfRegionList <- paste0(rawDir, "/", namePrefix, ".raw.", chunks, ".bcf.gz")
  write.table(bcfRegionList, file=outName, sep="\t", row.names=F, col.names=F, quote=F)
  
}


### include scaffold ranges
if(scaffoldListFile != "no") {
  scaffolds <- read.delim(scaffoldListFile, sep="\t", header=F, stringsAsFactors=F)[,1]
  scaffoldRegions <- c()
  
  for(i in 1:length(scaffolds)) {
    starts <- ends <- NULL
    len <- fai[scaffolds[i], 2]
    starts <- as.integer(seq(1, len, by=chunkSize))
    ends <- starts[-1] -1
    ends <- as.integer(c(ends, len))
    chunks <- paste0(scaffolds[i], "_", starts, "-", ends)
    scaffoldRegions <- c(scaffoldRegions, chunks)
  }
  
  outName <- paste0(callDir, "/Scaffolds_list.txt")
  bcfRegionList <- paste0(rawDir, "/", namePrefix, ".raw.", scaffoldRegions, ".bcf.gz")
  write.table(bcfRegionList, file=outName, sep="\t", row.names=F, col.names=F, quote=F)
  
  regions <- c(regions, scaffoldRegions)
  
}


### write all regions to file
outName <- paste0(rawDir, "/regions.txt")
write.table(regions, file=outName, sep="\t", row.names=F, col.names=F, quote=F)

