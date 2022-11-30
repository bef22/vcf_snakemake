# vcf_snakemake

Pipeline for creating vcf files from cram files on the Cambridge HPC cluster using Snakemake

The pipeline runs in 2 stages:
1. The first stage performs the mpileup+call, norm, missing individual QC and calculating the genome wide median depth
2. The second stage sets the FILTER column to PASS and fail and extracts biallelic sites


### Dependencies
- snakemake v7.8.5
the following must be in your path (for Durbin group you can set your PATH to include rds/rds-durbin-group-8b3VcZwY7rY/software/bin/)
- bcftools v1.15.1
- vcftools v0.1.17

### Required Files:
see example folder
- mainChromFile = file with main chromosome names per row
- scaffoldFile = file with scaffold names per row (optional)
- cramListFile = file with the full paths to each cram file per row
- speciesTableFile = tab delimited file with two columns: 1st is the sampleID, 2nd is the species

### profile
update the cluster.yaml file with the path to the location of the **profile** directory

### Update the Snakemake file:
set the following parameters
- namePrefix = provide a name for the vcf
- reference = path to the genome reference fasta
- mainChromFile = path to the mainChromFile (see above)
- scaffoldFile = path to scaffoldFile (optional), or set to "no" if the scaffolds should not be included
- cramListFile = path to cramListFile (see above)
- speciesTableFile = path to speciesTableFile (see above)
paths to the scripts
- R_QC_script
- R_depth_script
- perl_filter_script
- missingQCpass = "no" this should be "no" the first time the script is run, then it will stop for QC review


#### after QC files have been created update the following:
- medianDP = INT, median depth from file bcf_qc/depth/summary_site_depth.txt
- percentDP = INT, +/- percentDP % of medianDP, a nore stringent value would be 25, less stringent 50
- minQ = INT, minimum QUAL vlaue, for less stringet use 20
- minGQ = INT, minimum mead GQ value, suggest to use 30
- maxMissing = INT, maximum number of missing samples, less stringent use 75

### Run:
in a screen terminal and load the following module
```
module load python-3.9.6-gcc-5.4.0-sbr552h
```

To check how many jobs are required and display the shell commands use:
 ```
snakemake -n --printshellcmds -p
```

To run use:
```
snakemake --profile path2profile/ --latency 20 --jobs n
```

