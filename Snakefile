### Snakefile to create bcf files using bcftools, splitting chromosome into chunks and merging Scaffolds

"""

REQUIRED FILES AND UPDATE PARAMETRS IN SCRIPT:
- user = crsid
- cramListFile = file with the full paths to each cram (if >200 samples merge 50 samples into new sorted bam/cram files!)
- mainChromFile = file with main chromosome names per row
- scaffoldFile = file with scaffold names per row, or set to "no" if scaffolds should not be included
- speciesTableFile = tab delimited file with 2 columns: 1st is the sampleID, 2nd is the species
- reference = fasta file of reference genome
- fai = fai index of the reference genome


UPDATE IN SCRIPT:
- namePrefix = name of the sample set
- mutation_rate = 0.001                   # chichlids = 0.001, moths = 0.01
- chromChunkSize = 1000000                # default is 1000000, for <50 samples could increase this to 5000000
- biallelicType = "snps,indels"           # "snps,indels" or "snps", it is more efficient to include the indels at this stage if they might be required!

- passQC1 = "no"                          # review missing QC for chromosome chunks, change to yes after QC passed
- passQC2 = "no"                          # review missing and depth QC after normalisation, change to yes after QC passed


# after both QC steps have passed update the following
- medianDP = update this with the median or mode depth value from file bcf_qc/depth/depth_summary.txt, '<INT>'
- percentDP = medianDP +/- percentDP %, '<INT>', a stringent value would be 25, less stringent 50
- minQ = minimum QUAL vlaue, '<INT>', less stringet use 20 (in GATK forum it is suggested not to use this for filtering!)
- minGQ = minimum mead GQ value, '<INT>', suggest to use 30
- maxMissing = maximum number of missing samples, '<INT>', less stringent use 75


RUN
login to icelake node using login-icelake.hpc.cam.ac.uk
start a screen terminal (screen -RD)
make sure to use snakemake version 7.8.5

check how many jobs are required AT EACH STAGE, this command also shows the shell commands it would use:
- snakemake -n --printshellcmds -p

to run use:
- snakemake --profile .full_path_to/profile/ --latency 20 --jobs n
where n is the number of jobs to be allowed to start in parallel


NOTES
This script has to be run 4 times!
1. creates the regions for subsetting chromosomes into chunks (very fast ~1-2 min)
2. performs mpileup and call and does the missing QC per chunk (this takes the longest and depends on n samples and chunk size)
     after reviewing the QC plots change passQC1 = "yes" (but keep passQC2 = "no")
3. merges to chromosome level, performs normalisation step and does the missing QC and depth QC
     after reviewing the QC plots change passQC2 = "yes" an update the medianDP (and other filter parameters)
4. sets filter to PASS and fail terms, extracts biallelic 


AUTHORS and UPDATES
Bettina Fischer, 20230209
20230706, removed tabix_biallelic from localrules
20240504, updated to icelake-himem and hard coded the path to softwares in software_RHEL8/bin

"""

###############################################################################


### set global parameters
namePrefix = "example"
mainChromFile = "chromList.txt"            # file with chromosome names per row
scaffoldFile = "scaffoldList.txt"           # set this to "no" if scaffolds should not be included, otherwise provide the scaffold file
cramListFile = "cramList.txt"                # file with cramPaths
speciesTableFile = "speciesTable.txt"        # tab delimited file with 2 columns: 1st is the sampleID, 2nd is the species

reference = "/rds/project/rd109/rds-rd109-durbin-group/ref/fish/Astatotilapia_calliptera/fAstCal1.2/GCA_900246225.3_fAstCal1.2_genomic_chromnames_mt.fa"
fai = "/rds/project/rd109/rds-rd109-durbin-group/ref/fish/Astatotilapia_calliptera/fAstCal1.2/GCA_900246225.3_fAstCal1.2_genomic_chromnames_mt.fa.fai"

mutation_rate = 0.001                   # chichlids = 0.001, moths = 0.01
chromChunkSize = 1000000                # default is 1000000, for <50 samples could increase this to 5000000


### paths to dependant scripts
R_setChromosomeRanges_script = "./scripts/setChromosomeRanges.R"
R_QC_script = "./scripts/vcf_QC.R"
perl_filter_script = "./scripts/setPassFilter.pl"


### biallelic
# choose either "snps" or "snps,indels"
# if there is a chance that you might want to use the indels select them now and subset 
# snps at GWAS stage when creating input files
#biallelicType = "snps"                 # snps only
biallelicType = "snps,indels"           # snps and indels


### QC stops
# set this to "no" if the script hasn't been run before, update to "yes" after running next stage
passQC1 = "no"                          # review missing QC for chromosome chunks
passQC2 = "no"                          # review missing and depth QC after normalisation


### filter parameters for PASS sites - update these after first pass of script
# the result file name string will be named as _DP<percentDP>_Q<minQ>_GQ<minGQ>_MM<maxMissing>
# to switch off a parameter set it to 0
medianDP = 0                          # update this with the median depth value from file bcf_qc/depth/summary_site_depth.txt, '<INT>'
percentDP = 35                          # medianDP +/- percentDP %, '<INT>'
minQ = 20                               # minimum QUAL vlaue, '<INT>'
minGQ = 30                              # minimum mead GQ value, '<INT>'
maxMissing = 10                         # % max missing per site, '<INT>'


###############################################################################
### the following shouldn't be changed to ensure working software versions on icelake are used
bcfPath = f'/home/{user}/rds/rds-durbin-group-8b3VcZwY7rY/software_RHEL8/bin/bcftools'
vcfPath = f'/home/{user}/rds/rds-durbin-group-8b3VcZwY7rY/software_RHEL8/bin/vcftools'
tabixPath = f'/home/{user}/rds/rds-durbin-group-8b3VcZwY7rY/software_RHEL8/bin/tabix'

###############################################################################
# the following shouldn't require changes

import pdb
import subprocess
import os
import re


# create the .slurm directory
if not os.path.exists('.slurm'):
    os.makedirs('.slurm')


# need to constrain the namePrefix, otherwsie unexpected behavior in creating file names
wildcard_constraints:
    namePrefix=namePrefix


filterString="_DP"+str(percentDP)+"_Q"+str(minQ)+"_GQ"+str(minGQ)+"_MM"+str(maxMissing)
print("filterString: " + str(filterString) + '\n')



### Setting up the chrs wildcard
chromosomes=[]
with open(mainChromFile) as f:
    chromosomes=f.readlines()
main_chrs=[]
chrs=[]
chrsPlus=[]

for string in chromosomes:
    chromosomes_e = string.replace("\n","")
    main_chrs.append(chromosomes_e)
    chrs.append(chromosomes_e)
    chrsPlus.append(chromosomes_e)
#print("chromosome names used: " + str(main_chrs) + '\n')


### scaffolds wildcards
scaffolds=[]
if scaffoldFile != "no":
    # add scaffold name 
    chrsPlus.append("Scaffolds")
    scaff=[]
    with open(scaffoldFile) as f:
        scaff=f.readlines()
    scaffolds=[]

    for string in scaff:
        scaff_e = string.replace("\n","")
        scaffolds.append(scaff_e)
        chrs.append(scaff_e)
    #print("scaffold names used: " + str(scaffolds) + '\n')

#print("all names used: " + str(chrs) + '\n')


### Setting up the regions wildcard
regionsExists = os.path.isfile("bcf_raw/regions.txt")
regions=[]
#test=[]
if (regionsExists):
    subregions=[]
    with open("bcf_raw/regions.txt") as f:
        subregions=f.readlines()
    for string in subregions:
        regions_e = string.replace("\n","")
#        s = 'mississipi'
#        old = '_'
#       new = ':'
#        maxreplace = 1
#        test_e = new.join(regions_e.rsplit(old, maxreplace))
        regions.append(regions_e)
#        test.append(test_e)
    print("regions: " + str(regions) + '\n')
#    print("test: " + str(test) + '\n')


### replace the last _ with : in regions
def rreplace(s, old, new, occurrence):
 li = s.rsplit(old, occurrence)
 return new.join(li)


# rules not passed to SLURM
#localrules: all, set_chrom_ranges, tabix_biallelic
localrules: all, set_chrom_ranges


# create list of output files for first stage
myoutput = list()
myoutput.append("bcf_raw/regions.txt")
if (regionsExists):
    myoutput.append(expand("bcf_raw/{namePrefix}.raw.{regions}.bcf.gz", regions=regions, namePrefix=namePrefix))
    myoutput.append(expand("bcf_qc/raw_missing_individual/{namePrefix}.raw.{regions}.imiss", regions=regions, namePrefix=namePrefix))
    myoutput.append(expand("bcf_qc/raw_missing_individual/{namePrefix}_raw_missing_individual_summary_report.txt", namePrefix=namePrefix)) 



# steps after reviewing the 1st QC files
    if passQC1 == "yes":       
        myoutput.append(expand("bcf_call/{namePrefix}.call.{chrsPlus}.bcf.gz", chrsPlus=chrsPlus, namePrefix=namePrefix))
        myoutput.append(expand("bcf_norm/{namePrefix}.norm.{chrsPlus}.bcf.gz", chrsPlus=chrsPlus, namePrefix=namePrefix))
        myoutput.append(expand("bcf_qc/missing_individual/{namePrefix}.{chrsPlus}.imiss", chrsPlus=chrsPlus, namePrefix=namePrefix))
        myoutput.append(expand("bcf_qc/missing_individual/{namePrefix}_missing_individual_summary_report.txt", namePrefix=namePrefix)) 
        myoutput.append(expand("bcf_qc/depth/{namePrefix}.{main_chrs}_freq.txt.gz", main_chrs=main_chrs, namePrefix=namePrefix))
        myoutput.append("bcf_qc/depth/depth_histogram.png")
        myoutput.append("bcf_qc/depth/depth_summary.txt")


# steps after reviewing the 2nd QC files
    if passQC2 == "yes":
            myoutput.append(expand("bcf_normf/{namePrefix}{filterString}.normf.{chrsPlus}.bcf.gz.csi", chrsPlus=chrsPlus, namePrefix=namePrefix, filterString=filterString))
            myoutput.append(expand("bcf_biallelic_tmp/{namePrefix}{filterString}.biallelic.{chrsPlus}.bcf.gz", chrsPlus=chrsPlus, namePrefix=namePrefix, filterString=filterString))
            myoutput.append(expand("bcf_biallelic/{namePrefix}{filterString}.biallelic.bcf.gz", namePrefix=namePrefix, filterString=filterString))
            myoutput.append(expand("bcf_biallelic/{namePrefix}{filterString}.biallelic.bcf.gz.csi", namePrefix=namePrefix, filterString=filterString))

	

rule all:
    input:
        myoutput


### create chromosome regions
# R script
# this creates the regions for mpileup and call step
# it also creates the bcf list files for concat step, all scaffolds are merged into one file
rule set_chrom_ranges:
    input:
        {mainChromFile}, {fai}
#        {mainChromFile}, {scaffoldFile}, {fai}
    output:
        "bcf_raw/regions.txt"
    shell:
        "Rscript {R_setChromosomeRanges_script} {namePrefix} {mainChromFile} {scaffoldFile} {fai} {chromChunkSize}"




### run bcftools mpileup for each chr/LG (and optional the scaffolds)
# Note: bcftools uses threads only at the compression step so is not worthwhile using!!!
# bcftools mpileup:
# -f 		{reference} 
# -b 		{input.cramList} 
# -r 		<chrXX:START-END>: Genotype by region
# -C 0 		Coefficient for downgrading mapping quality for reads containing excessive mismatches. A zero value (the default) disables this functionality.
# -d 250000 Basically disables read limits per file. Want to take all reads for every sample.
# -L 250000 Basically disables read limits per file for INDEL calling.
# -q 20		skip alignments with mapQ smaller than INT
# -Q 13		skip bases with baseQ/BAQ smaller than INT 
# --ff UNMAP,SECONDARY,QCFAIL,DUP 	filter flags: skip reads with mask bits
# -a FORMAT/AD,FORMAT/DP,QS,SP,FORMAT/SCR,INFO/AD,INFO/SCR Additional annotations to add
# -p 		Apply -m and -F thresholds per sample to increase sensitivity of calling
# -O u		Output compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v)
# bcftools call:
# --ploidy 2 	predefined ploidy
# -a PV4,GQ,GP Comma-separated list of FORMAT and INFO tags to output
# -m 		alternative model for multiallelic and rare-variant calling designed to overcome known limitations in -c calling model (conflicts with -c)	
# -P 0.001	Mutation rate FLOAT (use bigger for greater sensitivity), use with -m.
# -G {input.speciesTable}  		by default, all samples are assumed to come from a single population. This option allows to group samples into populations and apply the HWE assumption within but not across the populations. FILE is a tab-delimited text file with sample names in the first column and group names in the second column. The -G option requires the presence of per-sample FORMAT/QS or FORMAT/AD tag generated with bcftools mpileup -a QS (or -a AD).
# bcftools +fill-tags
# -O b 		output bcf
# -t 'AF,ExcHet,NS'	which tags to update
#
rule mpileup_call_filltags:
    input:
        cramList = {cramListFile},	
        speciesTable = {speciesTableFile}
    output:
        "bcf_raw/{namePrefix}.raw.{regions}.bcf.gz"
    params:
        chrRegions = lambda wildcards: rreplace(wildcards.regions, '_', ':', 1)

    shell:
        "{bcfPath} mpileup -f {reference} -b {input.cramList} -r {params.chrRegions} -C 0 -d 250000 -L 250000 -q 20 -Q 13 --ff UNMAP,SECONDARY,QCFAIL,DUP -a FORMAT/AD,FORMAT/DP,QS,SP,FORMAT/SCR,INFO/AD,INFO/SCR -p -O u | {bcfPath} call --ploidy 2 -a PV4,GQ,GP -m -P {mutation_rate} -O u -G {input.speciesTable} | {bcfPath} +fill-tags -O b -o {output} -- -t 'AF,ExcHet,NS'"


### check for missing individual in raw files
rule raw_missing_individual:
    input:
        "bcf_raw/{namePrefix}.raw.{regions}.bcf.gz"
    output:
        "bcf_qc/raw_missing_individual/{namePrefix}.raw.{regions}.imiss"
    shell:
        "{vcfPath} --bcf {input} --missing-indv --stdout > {output}"


### plot and report missing individual R in raw files
# creates a stripchart for chromosomes and scaffolds and a report file which contains individuals which are 
# above 1.5 times the interquartile range above the upper quartile (e.g. whiskers in boxplot) or have
# > 0.9 missing
rule raw_report_missing_individuals:
    input:
        expand("bcf_qc/raw_missing_individual/{namePrefix}.raw.{regions}.imiss", regions=regions, namePrefix=namePrefix)
    params:
        filePath = "bcf_qc/raw_missing_individual/"
    output:
        "bcf_qc/raw_missing_individual/{namePrefix}_raw_missing_individual_summary_report.txt"
    shell:
        "Rscript {R_QC_script} {params.filePath} {namePrefix} bcf_raw/regions.txt no raw_imiss"


#***** check-point, script stops here and <passQC1> must be updated to "yes" to continue *****#


### concat to chromosome level and all scaffolds
rule concat:
    input:
        "bcf_call/{chrsPlus}_list.txt"
    output:
        "bcf_call/{namePrefix}.call.{chrsPlus}.bcf.gz"
    shell:
        "{bcfPath} concat -f {input} -O b -o {output}"


### normalise and merge multiallelic sites
# bcftools norm:
# -m +any	-|+[snps|indels|both|any]     split multiallelic sites into biallelic records (-) or join biallelic sites into multiallelic records (+). An optional type string can follow which controls variant types which should be split or merged together: If only SNP records should be split or merged, specify snps; if both SNPs and indels should be merged separately into two records, specify both; if SNPs and indels should be merged into a single record, specify any.
#
rule norm_merge:
    input:
        "bcf_call/{namePrefix}.call.{chrsPlus}.bcf.gz"
    output:
        "bcf_norm/{namePrefix}.norm.{chrsPlus}.bcf.gz"
    shell:
        "{bcfPath} norm -f {reference} -m +any -O u {input} | {bcfPath} +fill-tags -O b -o {output} -- -t 'AF,AC,AN,ExcHet,NS'"


### check for missing individual
rule missing_individual:
    input:
        "bcf_norm/{namePrefix}.norm.{chrsPlus}.bcf.gz"
    output:
        "bcf_qc/missing_individual/{namePrefix}.{chrsPlus}.imiss"
    shell:
        "{vcfPath} --bcf {input} --missing-indv --stdout > {output}"


### plot and report missing individual R
# creates a stripchart for chromosomes and scaffolds and a report file which contains individuals which are 
# above 1.5 times the interquartile range above the upper quartile (e.g. whiskers in boxplot) or have
# > 0.9 missing
rule report_missing_individuals:
    input:
        expand("bcf_qc/missing_individual/{namePrefix}.{chrsPlus}.imiss", chrsPlus=chrsPlus, namePrefix=namePrefix)
    params:
        filePath = "bcf_qc/missing_individual/"
    output:
        "bcf_qc/missing_individual/{namePrefix}_missing_individual_summary_report.txt"
    shell:
        "Rscript {R_QC_script} {params.filePath} {namePrefix} {mainChromFile} {scaffoldFile} imiss"


### get site depths across main chromosomes
# this extracts the depth in the vcf files at all positions usig the AD field and creates a frequency table
rule site_depth:
    input:
        "bcf_norm/{namePrefix}.norm.{main_chrs}.bcf.gz"
    output:
        temp("bcf_qc/depth/{namePrefix}.{main_chrs}_freq.txt")
    shell:
        "{bcfPath} view {input} | {vcfPath} --vcf - --site-depth --stdout | grep -v 'SUM_DEPTH' | cut -f 3 | sort -n | uniq -c > {output}"


### gzip depth files
rule gzip_site_depth:
    input:
        "bcf_qc/depth/{namePrefix}.{main_chrs}_freq.txt"
    output:
        "bcf_qc/depth/{namePrefix}.{main_chrs}_freq.txt.gz"
    shell:
        "gzip {input}"


### calculate genome-wide median depth 
# get the genome wise median and plot a histogram of the site depths
rule median_DP:
    input:
        expand("bcf_qc/depth/{namePrefix}.{main_chrs}_freq.txt.gz", main_chrs=main_chrs, namePrefix=namePrefix)
    params:
        filePath = "bcf_qc/depth/"
    output:
        "bcf_qc/depth/depth_histogram.png", "bcf_qc/depth/depth_summary.txt"
    shell:
        "Rscript {R_QC_script} {params.filePath} {namePrefix} {mainChromFile} no depth"


#***** check-point, script stops here and <passQC2> must be updated to "yes" to continue *****#


### set filter to PASS and lowDP|highDP|lowQ|lowGQ|highMiss flags
# this only sets the filter column flag but doesn't remove any sites
rule set_vcf_filter:
    input:
        "bcf_norm/{namePrefix}.norm.{chrsPlus}.bcf.gz"
    output:
        bcf="bcf_normf/{namePrefix}_DP{percentDP}_Q{minQ}_GQ{minGQ}_MM{maxMissing}.normf.{chrsPlus}.bcf.gz",
        csi="bcf_normf/{namePrefix}_DP{percentDP}_Q{minQ}_GQ{minGQ}_MM{maxMissing}.normf.{chrsPlus}.bcf.gz.csi"
    run:
        shell("{bcfPath} view {input} | {perl_filter_script} --medianDP {medianDP} --percentDP {percentDP} --minQ {minQ} --minGQ {minGQ} --maxMissing {maxMissing} /dev/stdin | {bcfPath} view -O b -o {output.bcf}")
        shell("{tabixPath} -p bcf {output.bcf}")


### biallelic per chromosome - temporary file
rule extract_biallelic:
    input:
        "bcf_normf/{namePrefix}{filterString}.normf.{chrsPlus}.bcf.gz"
    output:
        temp("bcf_biallelic_tmp/{namePrefix}{filterString}.biallelic.{chrsPlus}.bcf.gz")
    shell:
        "{bcfPath} view -m2 -M2 -v {biallelicType} -O b -o {output} {input}"


### merge biallelic
rule merge_biallelic_bcfs:
    input:
        expand("bcf_biallelic_tmp/{namePrefix}{filterString}.biallelic.{chrsPlus}.bcf.gz", chrsPlus=chrsPlus, namePrefix=namePrefix, filterString=filterString)
    output:
        "bcf_biallelic/{namePrefix}{filterString}.biallelic.bcf.gz"
    shell:
        "{bcfPath} concat -O b -o {output} {input}"


### tabix
rule tabix_biallelic:
    input:
        "bcf_biallelic/{namePrefix}{filterString}.biallelic.bcf.gz"
    output:
        "bcf_biallelic/{namePrefix}{filterString}.biallelic.bcf.gz.csi"
    shell:
        "{tabixPath} -p bcf {input}"


