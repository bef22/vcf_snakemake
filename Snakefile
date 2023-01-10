### Snakefile to create bcf files using bcftools

"""

REQUIRED FILES:
- mainChromFile = file with main chromosome names per row
- scaffoldFile = file with scaffold names per row, or set to "no" if scaffolds should not be included
- cramListFile = file with the full paths to each cram
- speciesTableFile = tab delimited file with 2 columns: 1st is the sampleID, 2nd is the species

UPDATE
the following parameters in the script below:
namePrefix = this is the name of the resulting vcf
reference = path to the genome reference fasta
mainChromFile, scaffoldFile, cramListFile, speciesTableFile - see above

R_QC_script, R_depth_script, perl_filter_script paths

missingQCpass = "no" this should be "no" the first time the script is run, then it will stop for QC review

# after QC files have been created update the following
medianDP = update this with the median depth value from file bcf_qc/depth/summary_site_depth.txt, '<INT>'
percentDP = medianDP +/- percentDP %, '<INT>', a stringent value would be 25, less stringent 50
minQ = minimum QUAL vlaue, '<INT>', less stringet use 20
minGQ = minimum mead GQ value, '<INT>', suggest to use 30
maxMissing = maximum number of missing samples, '<INT>', less stringent use 75


RUN
in a screen terminal and load the following modules
- module load python-3.9.6-gcc-5.4.0-sbr552h

check how many jobs are required, this command also shows the shell commands it would use:
- snakemake -n --printshellcmds -p

to run use:
- snakemake --profile path2profile/ --latency 20 --jobs n




NOTES
- tried very hard to better deal with the scaffolds but failed to get snakemake to deal with this


AUTHORS and UPDATES
20221129 Written by Bettina Fischer (bef22) - version 1.0
20221212 updated by Bettina Fischer - version 1.1

"""

###############################################################################


### set global parameters
namePrefix = "zebra_callainos_pearly_MZeb"

reference = "/rds/project/rd109/rds-rd109-durbin-group/ref/fish/Maylandia_zebra/M_zebra_UMD2a/bwa-mem2_index/GCF_000238955.4_M_zebra_UMD2a_genomic_LG.fa"

mutation_rate = 0.001                                # chichlids = 0.001

mainChromFile = "example/chromList.txt"              # file with chromosome names per row
scaffoldFile = "example/scaffoldList.txt"            # set this to "no" if scaffolds should not be included, otherwise provide the scaffold file
cramListFile = "example/cramList.txt"                # file with cramPaths
speciesTableFile = "example/speciesTable.txt"        # tab delimited file with 2 columns: 1st is the sampleID, 2nd is the species

# set this to "no" if the script hasn't been run before
missingQCpass = "no"                        # update this to "yes" after reviewing the QC files, update the filter parameters and start snakemake again


### filter parameters for PASS sites - update these after first pass of script
# the result file name string will be named as _DP<percentDP>_Q<minQ>_GQ<minGQ>_MM<maxMissing>
# to switch off a parameter set it to 0
medianDP = 0                                 # update this with the median depth value from file bcf_qc/depth/summary_site_depth.txt, '<INT>'
percentDP = 25                               # medianDP +/- percentDP %, '<INT>'
minQ = 20                                    # minimum QUAL vlaue, '<INT>'
minGQ = 30                                   # minimum mead GQ value, '<INT>'
maxMissing = 50                              # % max missing per site, '<INT>'


### paths to dependant scripts
R_QC_script = "scripts/QC.R"
R_depth_script = "scripts/depth.R"
perl_filter_script = "scripts/setPassFilter.pl"


###############################################################################
# the following shouldn't require changes

import pdb
import subprocess
import os


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

for string in chromosomes:
    chromosomes_e = string.replace("\n","")
    main_chrs.append(chromosomes_e)
    chrs.append(chromosomes_e)
#print("chromosome names used: " + str(main_chrs) + '\n')


### scaffolds wildcards
scaffolds=[]
#if callScaffolds == "yes":
if scaffoldFile != "no":
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



# rules not passed to SLURM
localrules: all, tabix_biallelic


# create list of output files for first stage
myoutput = list()
myoutput.append(expand("bcf_call/{namePrefix}.call.{chrs}.bcf.gz", chrs=chrs, namePrefix=namePrefix))
myoutput.append(expand("bcf_norm/{namePrefix}.norm.{chrs}.bcf.gz", chrs=chrs, namePrefix=namePrefix))
myoutput.append(expand("bcf_qc/missing_individual/{namePrefix}.{chrs}.imiss", chrs=chrs, namePrefix=namePrefix))
myoutput.append(expand("bcf_qc/missing_individual/{namePrefix}_missing_individual_summary_report.txt", namePrefix=namePrefix)) 
myoutput.append(expand("bcf_qc/depth/{namePrefix}.{main_chrs}_freq.txt.gz", main_chrs=main_chrs, namePrefix=namePrefix))
myoutput.append("bcf_qc/depth/site_depth_histogram.png")
myoutput.append("bcf_qc/depth/summary_site_depth.txt")


# steps after reviewing the QC and depth files
if missingQCpass == "yes":
    myoutput.append(expand("bcf_normf/{namePrefix}{filterString}.normf.{chrs}.bcf.gz.csi", chrs=chrs, namePrefix=namePrefix, filterString=filterString))
    myoutput.append(expand("bcf_biallelic_tmp/{namePrefix}{filterString}.biallelic.{chrs}.bcf.gz", chrs=chrs, namePrefix=namePrefix, filterString=filterString))
    myoutput.append(expand("bcf_biallelic/{namePrefix}{filterString}.biallelic.bcf.gz", namePrefix=namePrefix, filterString=filterString))
    myoutput.append(expand("bcf_biallelic/{namePrefix}{filterString}.biallelic.bcf.gz.csi", namePrefix=namePrefix, filterString=filterString))
	

rule all:
    input:
        myoutput


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
        "bcf_call/{namePrefix}.call.{chrs}.bcf.gz"
    params:
        chromosome = lambda wildcards: wildcards.chrs
    shell:
        "bcftools mpileup -f {reference} -b {input.cramList} -r {params.chromosome} -C 0 -d 250000 -L 250000 -q 20 -Q 13 --ff UNMAP,SECONDARY,QCFAIL,DUP -a FORMAT/AD,FORMAT/DP,QS,SP,FORMAT/SCR,INFO/AD,INFO/SCR -p -O u | bcftools call --ploidy 2 -a PV4,GQ,GP -m -P {mutation_rate} -O u -G {input.speciesTable} | bcftools +fill-tags -O b -o {output} -- -t 'AF,ExcHet,NS'"


### normalise and merge multiallelic sites
# bcftools norm:
# -m +any	-|+[snps|indels|both|any]     split multiallelic sites into biallelic records (-) or join biallelic sites into multiallelic records (+). An optional type string can follow which controls variant types which should be split or merged together: If only SNP records should be split or merged, specify snps; if both SNPs and indels should be merged separately into two records, specify both; if SNPs and indels should be merged into a single record, specify any.
#
rule norm_merge:
    input:
        "bcf_call/{namePrefix}.call.{chrs}.bcf.gz"
    output:
        "bcf_norm/{namePrefix}.norm.{chrs}.bcf.gz"
    shell:
        "bcftools norm -f {reference} -m +any -O u {input} | bcftools +fill-tags -O b -o {output} -- -t 'AF,AC,AN,ExcHet,NS'"


### check for missing individual
rule missing_individual:
    input:
        "bcf_norm/{namePrefix}.norm.{chrs}.bcf.gz"
    output:
        "bcf_qc/missing_individual/{namePrefix}.{chrs}.imiss"
    params:
        chromosome = lambda wildcards: wildcards.chrs
    shell:
        "vcftools --bcf {input} --chr {params.chromosome} --missing-indv --stdout > {output}"


### plot and report missing individual R
# creates a stripchart for chromosomes and scaffolds and a report file which contains individuals which are above 1.5 times the interquartile range above the upper quartile (e.g. whiskers in boxplot)
rule report_missing_individuals:
    input:
        expand("bcf_qc/missing_individual/{namePrefix}.{chrs}.imiss", chrs=chrs, namePrefix=namePrefix)
    params:
        filePath = "bcf_qc/missing_individual/"
    output:
        "bcf_qc/missing_individual/{namePrefix}_missing_individual_summary_report.txt"
    shell:
        "Rscript {R_QC_script} {params.filePath} {namePrefix} {mainChromFile} {scaffoldFile} individual"



### get site depths across main chromosomes
# this extracts the depth in the vcf files at all positions usig the AD field and creates a frequency table
rule site_depth:
    input:
        "bcf_norm/{namePrefix}.norm.{main_chrs}.bcf.gz"
    output:
        temp("bcf_qc/depth/{namePrefix}.{main_chrs}_freq.txt")
    shell:
        "bcftools view {input} | vcftools --vcf - --site-depth --stdout | grep -v 'SUM_DEPTH' | cut -f 3 | sort -n | uniq -c > {output}"
	#"bcftools view --min-ac 1 {input} | vcftools --vcf - --site-depth --stdout > {output}"


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
        "bcf_qc/depth/site_depth_histogram.png", "bcf_qc/depth/summary_site_depth.txt"
    shell:
        "Rscript {R_depth_script} {params.filePath} {namePrefix} {mainChromFile}"



#***** check-point, script stops here and <missingQCpass> must be updated to "yes" to continue *****#



### set filter to PASS and lowDP|highDP|lowQ|lowGQ|highMiss flags
# this only sets the filter column flag but doesn't remove any sites
rule set_vcf_filter:
    input:
        "bcf_norm/{namePrefix}.norm.{chrs}.bcf.gz"
    output:
        bcf="bcf_normf/{namePrefix}_DP{percentDP}_Q{minQ}_GQ{minGQ}_MM{maxMissing}.normf.{chrs}.bcf.gz",
        csi="bcf_normf/{namePrefix}_DP{percentDP}_Q{minQ}_GQ{minGQ}_MM{maxMissing}.normf.{chrs}.bcf.gz.csi"
    run:
        shell("bcftools view {input} | {perl_filter_script} --medianDP {medianDP} --percentDP {percentDP} --minQ {minQ} --minGQ {minGQ} --maxMissing {maxMissing} /dev/stdin | bcftools view -O b -o {output.bcf}")
        shell("tabix -p bcf {output.bcf}")


### biallelic per chromosome - temporary file
rule extract_biallelic:
    input:
        "bcf_normf/{namePrefix}{filterString}.normf.{chrs}.bcf.gz"
    output:
        temp("bcf_biallelic_tmp/{namePrefix}{filterString}.biallelic.{chrs}.bcf.gz")
    shell:
        "bcftools view -m2 -M2 -v snps -O b -o {output} {input}"


### merge biallelic
rule merge_biallelic_bcfs:
    input:
        expand("bcf_biallelic_tmp/{namePrefix}{filterString}.biallelic.{chrs}.bcf.gz", chrs=chrs, namePrefix=namePrefix, filterString=filterString)
    output:
        "bcf_biallelic/{namePrefix}{filterString}.biallelic.bcf.gz"
    shell:
        "bcftools concat -O b -o {output} {input}"


### tabix
rule tabix_biallelic:
    input:
        "bcf_biallelic/{namePrefix}{filterString}.biallelic.bcf.gz"
    output:
        "bcf_biallelic/{namePrefix}{filterString}.biallelic.bcf.gz.csi"
    shell:
        "tabix -p bcf {input}"




