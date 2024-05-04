# vcf_snakemake

Pipeline for creating vcf files from cram files on the Cambridge HPC cluster login-icelake.hpc.cam.ac.uk node using Snakemake.

The pipeline runs in 4 stages:
1. The first stage creates regions which splits chromosomes into user specified sized chunks (very fast)
2. The second run performs the mpileup+call on each of the regions and creates missing individual QC plots. After reviewing the raw QC plots update passQC1 = "yes" (but leave passQC2 = "no")
3. The third merges the regions back to chromosome level and creates missing individual QC plots and genome wide depth assessment plots. After reviewing the QC plots update passQC2 = "yes" and update the filter parameters
6. The last stage sets the FILTER column to PASS and fail and extracts biallelic sites


### Dependencies
- snakemake v7.8.5  (It only seems to work for snakemake v7.8.5 and not any later versions!)

The following software versions have been hard coded in the Snakefile to use the compiled versions to run on icelake RHLE8 (rds/rds-durbin-group-8b3VcZwY7rY/software_RHEL8/bin/)
- bcftools v1.20
- vcftools v0.1.17
- tabix v1.20

### profile
Update the cluster.yaml file with the full path to the location of the **profile** directory (I found that relative paths do not always work!). Update the config.yaml account to your groups account.

### scripts
Please make sure you run the following commands in the scripts folder
```
chmod ug+rwx *
dos2unix *                          (otherwise you can get a /usr/bin/perl^M: bad interpreter error)
```

### Performance considerations
If you have 200 or more cram files you should merge 50 files into position sorted cram files and use these in the cramListFile instead of the 200 or more individual cram files to limit IO overload. If you have less than 50 cram files you can select use large chromosome chunks.


### Required Files:
see example folder which is for Maylandia genome aligned data
- mainChromFile = file with main chromosome names per row (make sure you use the correct chromosome names for your reference genome - see below)
- scaffoldFile = file with scaffold names per row (optional)
- cramListFile = file with the full paths to each cram file per row
- speciesTableFile = tab delimited file with two columns: 1st is the sampleID, 2nd is the species

#### Main chromosome names for cichlid refrerence assemblies
- fAstCal1.2: chr1, chr2, chr3 ... chr20, chr22, chr23
- fAulStu2:  chr1, chr2, chr3 ... chr20, chr22, chr23
- M_zebra_UMD2: LG1, LG2, LG3 ... LG20, LG21, LG23

#### Reference fasta locations on rds
- /rds/project/rd109/rds-rd109-durbin-group/ref/fish/Astatotilapia_calliptera/fAstCal1.2/GCA_900246225.3_fAstCal1.2_genomic_chromnames_mt.fa
- /rds/project/rd109/rds-rd109-durbin-grou/ref/fish/Aulonocara_stuartgranti/fAulStu2/fAulStu2.mm.asm.FINAL.2.mt.fa
- /rds/project/rd109/rds-rd109-durbin-group/ref/fish/Maylandia_zebra/M_zebra_UMD2a/bwa-mem2_index/GCF_000238955.4_M_zebra_UMD2a_genomic_LG.fa

### Update the Snakemake file:
set the following parameters
- user = crsid
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
open a screen terminal (take a not of the login node so you can attach back to the same later)
```
screen -RD
```

To check how many jobs are required and display the shell commands use:
 ```
snakemake -n --printshellcmds -p
```
This will print a DAG table in yellow which states the number of total jobs this run/stage will require.


To run use the full path to the profile folder, and where n is the maximum number of jobs to be submitted in parallel:
```
snakemake --profile full_path_to/profile/ --latency 20 --jobs n
```

The pipeline runs in 4 stages:
1. The first stage creates regions which splits chromosomes into user specified sized chunks (very fast)
2. The second run performs the mpileup+call on each of the regions and creates missing individual QC plots. After reviewing the raw QC plots update passQC1 = "yes" (but leave passQC2 = "no")
3. The third merges the regions back to chromosome level and creates missing individual QC plots and genome wide depth assessment plots. After reviewing the QC plots update passQC2 = "yes" and update the filter parameters
6. The last stage sets the FILTER column to PASS and fail and extracts biallelic sites


### Troubleshooting:
You will receive an email for each failed submission, you can check the log files per rule and file/submission in .slurm directory.


#### To overwrite filter parameters
and add maf 10% and PASS filtering (as example using medianDPValue = 2040, increase the percentDP to +/- 50%, mores stringent on the maxMissing 10%)

```
bcftools view bcf_file | scripts/setPassFilter.pl --medianDP 2040 --percentDP 50 --minQ 20 --minGQ 30 --maxMissing 10 /dev/stdin | bcftools view -f PASS -q 0.1:minor -O b -o new_bcf_file
```

