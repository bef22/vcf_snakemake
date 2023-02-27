# vcf_snakemake

Pipeline for creating vcf files from cram files on the Cambridge HPC cluster using Snakemake. The pipeline can also be found on rds in /rds/rds-durbin-group-8b3VcZwY7rY/projects/cichlid/vcfs/vcf_call_snakemake_resource

The pipeline runs in 4 stages:
1. The first stage creates regions which splits chromosomes into user specified sized chunks (very fast)
2. The second run performs the mpileup+call on each of the regions and creates missing individual QC plots
3. The third merges the regions back to chromosome level and creates missing individual QC plots and genome wide depth assessment plots
6. The last stage sets the FILTER column to PASS and fail and extracts biallelic sites


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
update the cluster.yaml file with the path to the location of the **profile** directory or use relative paths

### Update the Snakemake file:
set the following parameters
- namePrefix = provide a name for the vcf
- reference = path to the genome reference fasta
- fai = fai index of the reference
- mainChromFile = path to the mainChromFile (see above)
- scaffoldFile = path to scaffoldFile (optional), or set to "no" if the scaffolds should not be included
- cramListFile = path to cramListFile (see above)
- speciesTableFile = path to speciesTableFile (see above)
- mutation_rate = mutation rate for this species
- chromChunkSize = size into which chromsomes are split during mpileup+call step, this will depend on how many samples are processed at the same time
- biallelicType = choose either "snps,indels" or "snps"

paths to the scripts
- R_setChromosomeRanges_script
- R_QC_script
- perl_filter_script
- passQC1 = "no" this should be "no" the first time the script is run, then it will stop for 1st QC review, update to "yes" to continue
- passQC2 = "no" this should be "no" the first time the script is run, then it will stop for 2nd QC review, update to "yes" to continue

#### after both QC steps have been complated update the following:
- medianDP = INT, median or mode depth from file bcf_qc/depth/summary_site_depth.txt
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
snakemake --profile ./profile/ --latency 20 --jobs n
```

### Troubleshooting:
You will receive an email for each failed submission, you can check the log files per rule and file/submission in .slurm directory.


#### To overwrite filter parameters
and add maf 10% and PASS filtering (as example using medianDPValue = 2040, increase the percentDP to +/- 50%, mores stringent on the maxMissing 10%)

```
bcftools view bcf_file | scripts/setPassFilter.pl --medianDP 2040 --percentDP 50 --minQ 20 --minGQ 30 --maxMissing 10 /dev/stdin | bcftools view -f PASS -q 0.1:minor -O b -o new_bcf_file
```

