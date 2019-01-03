# *Snakemake workflow for TITAN analysis of 10X Genomics WGS*

# Description
This workflow will run the TITAN copy number analysis for set of tumour-normal pairs, starting from the BAM files aligned using [Long Ranger](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/what-is-long-ranger) software. The analysis includes haplotype-based copy number prediction and post-processing of results. It will also perform model selection at the end of the workflow to choose the optimal ploidy and clonal cluster solutions.  
[Viswanathan SR*, Ha G*, Hoff A*, et al. Structural Alterations Driving Castration-Resistant Prostate Cancer Revealed by Linked-Read Genome Sequencing. *Cell* 174, 433–447.e19 (2018).](https://www.cell.com/cell/abstract/S0092-8674(18)30648-2)

# Contact
Gavin Ha  
Fred Hutchinson Cancer Research Center  
contact: <gavinha@gmail.com> or <gha@fredhutch.org>  
Date: August 7, 2018  
Website: [GavinHaLab.org](https://gavinhalab.org/)

# Requirements
## Software packages or libraries
 - R-3.4
   - [TitanCNA](https://github.com/gavinha/TitanCNA) (v1.15.0) or higher
   		- TitanCNA imports: GenomicRanges, GenomeInfoDb, VariantAnnotation, dplyr, data.table, foreach
   - [ichorCNA](https://github.com/broadinstitute/ichorCNA) (v0.1.0) 
   - HMMcopy
   - optparse
   - stringr
   - SNPchip
   - doMC
 - Python 3.5
   - snakemake-5.2.0
   - PySAM-0.11.2.1
   - PyYAML-3.12
 - [bxtools](https://github.com/walaj/bxtools)

# Files in the workflow
### Scripts used by the workflow
The following scripts are used by this snakemake workflow:
 - [getMoleculeCoverage.R](code/getMoleculeCoverage.R) Normalizing/correcting molecule-level coverage
 - [getPhasedHETSitesFromLLRVCF.R](code/getPhasedHETSitesFromLLRVCF.R) - Extracts phased germline heterozygous SNP sites from the Long Ranger analysis of the normal sample
 - [getTumourAlleleCountsAtHETSites.py](code/getTumourAlleleCountsAtHETSites.py) - Extracts allelic counts from the tumor sample at the germline heterozygous SNP sites
 - [titanCNA_v1.15.0_TenX.R](code/titanCNA_v1.15.0_TenX.R) - Main R script to run TitanCNA
 - [selectSolution.R](code/selectSolution.R) - R script to select optimal solution for each sample
 - [combineTITAN-ichor.R](code/combineTITAN-ichor.R) - R script to merge autosomes and chrX results, plus post-processing steps including adjusting max copy values.

### Tumour-Normal sample list [config/samples.yaml](config/samples.yaml)
The list of tumour-normal paired samples should be defined in a YAML file. In particular, the [Long Ranger](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/what-is-long-ranger) (v2.2.2) analysis directory is listed under samples.  See [config/samples.yaml](config/samples.yaml) for an example.  Both fields `samples` and `pairings` must to be provided.  `pairings` key must match the tumour sample while the value must match the normal sample.
```
samples:
  tumor_sample_1:  /path/to/tumor/longranger/dir
  normal_sample_1:  /path/to/normal/longranger/dir


pairings:
  tumor_sample_1:  normal_sample_1
```

### snakefiles
1. [moleculeCoverage.snakefile](moleculeCoverage.snakefile)
2. [getPhasedAlleleCounts.snakefile](getPhasedAlleleCounts.snakefile)
3. [TitanCNA.snakefile](TitanCNA.snakefile)

### [config.yaml](config/config.yaml)
See below for details about [config/config.yaml](config/config.yaml)

# Run the analysis

## 1. Invoking the full snakemake workflow for TITAN on a local machine
This will also run both [moleculeCoverage.snakefile](moleculeCoverage.snakefile) and [getPhasedAlleleCounts.snakefile](getPhasedAlleleCounts.snakefile) which generate the necessary inputs for [TitanCNA.snakefile](TitanCNA.snakefile).
```
# show commands and workflow
snakemake -s TitanCNA.snakefile -np
# run the workflow locally using 5 cores
snakemake -s TitanCNA.snakefile --cores 5
```
## 2. Invoking the TITAN snakemake workflow on a cluster
Here are instructions for running workflow on a cluster using specific resource settings for memory and runtime limits, and parallel environments.  
There are two cluster configurations provided: `qsub` and `slurm`

#### a. `qsub`
There are 2 separate files in use for `qsub`, which are provided as a template:
	`config/cluster_qsub.sh` - This file contains other `qsub` parameters. *Note that these settings are used for the Broad's UGER cluster so users will need to modify this for their own clusters.*  
	`config/cluster_qsub.yaml` - This file contains the memory, runtime, and number of cores for certain tasks.  

To invoke the snakemake pipeline for `qsub`:
```
snakemake -s  TitanCNA.snakefile --jobscript config/cluster_qsub.sh --cluster-config config/cluster_qsub.yaml --cluster-sync "qsub -l h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -pe {cluster.pe} -binding {cluster.binding}" -j 100
```
Here, the `h_vmem` (max memory), `h_rt` (max runtime) are used. For `runTitanCNA` task, the default setting is to use 1 core but additional number of cpus (per task) can help to speed up the analysis. This can be set with `-pe` and `-binding`. Your SGE settings may be different and users should adjust accordingly.

#### b. `slurm`
There is only one file in use for `slurm`:
	`config/cluster_slurm.yaml` - This file contains the memory, runtime, and number of cores for certain tasks. 
To invoke the snakemake pipeline for `qsub`:
```
snakemake -s  TitanCNA.snakefile --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output}" -j 50
```


## 3. Invoking individual steps in the workflow
Users can run the snakemake files individually. This can be helpful for testing each step or if you only wish to generate results for a particular step. The snakefiles need to be run in this same order since input files are generated by the previous steps.
  ### a. [moleculeCoverage.snakefile](moleculeCoverage.snakefile)
  i.   Run [bxtools](https://github.com/walaj/bxtools) to compute counts of unique molecules in each window.  
  ii.  Perform GC-content bias correction for barcode counts.  
  iii. Perform ichorCNA analysis to generate initial molecule coverage-based copy number. For male samples, chrX results will be used from this step.  
  ```
  snakemake -s moleculeCoverage.snakefile -np
  snakemake -s moleculeCoverage.snakefile --cores 5
  # OR
  snakemake -s  moleculeCoverage.snakefile --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output}" -j 50
  # OR
  snakemake -s moleculeCoverage.snakefile --jobscript config/cluster_qsub.sh --cluster-config config/cluster_qsub.yaml --cluster-sync "qsub -l h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -pe {cluster.pe} -binding {cluster.binding}" -j 100
  ```
  
  ### b. [getPhasedAlleleCounts.snakefile](getPhasedAlleleCounts.snakefile) 
  i.   Read the Long Ranger output file `phased_variants.vcf.gz` and extract heterozygous SNP sites (that overlap a SNP database, e.g. [hapmap_3.3.hg38.vcf.gz](https://storage.cloud.google.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz?_ga=2.110868357.-1633399588.1531762721)).  
  ii.   Extract the allelic read counts from the Long Ranger tumor bam file `phased_possorted_bam.bam` for each chromosome.  
  iii. Cat the allelic read counts from each chromosome file into a single counts file.
  ```
  snakemake -s getPhasedAlleleCounts.snakefile -np
  snakemake -s getPhasedAlleleCounts.snakefile --cores 5
  # OR
   snakemake -s getPhasedAlleleCounts.snakefile --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output}" -j 50
   # OR
   snakemake -s getPhasedAlleleCounts.snakefile --jobscript config/cluster_qsub.sh --cluster-config config/cluster_qsub.yaml --cluster-sync "qsub -l h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -pe {cluster.pe} -binding {cluster.binding}" -j 100
  ``` 
  ### c. [TitanCNA.snakefile](TitanCNA.snakefile)
  i.   Run the [TitanCNA](https://github.com/gavinha/TitanCNA) analysis and generates solutions for different ploidy initializations and each clonal cluster.  
  ii.  Merge results with ichorCNA output generate by [moleculeCoverage.snakefile](moleculeCoverage.snakefile) and post-processes copy number results.  
  iii. Select optimal solution for each samples and copies these to a new folder. The parameters are compiled in a text file.  
  
# Configuration and settings
All settings for the workflow are contained in [config/config.yaml](config/config.yaml). The settings are organized by paths to scripts and reference files and then by each step in the workflow.

### 1. Path to tools
These tools are used by various `snakefiles` and are required.
```
samTools:  /path/to/samtools ## need to specify
bxTools:  /path/to/bxtools ## need to specify
```

### 2. Path to scripts
These are provided in this repo under [code/](code/). The prefix to the name of the setting indicates which `snakefile` it is used in. 
```
molCov_script:  code/getMoleculeCoverage.R
phaseCounts_hetSites_script:  code/getPhasedHETSitesFromLLRVCF.R
phaseCounts_counts_script:  code/getTumourAlleleCountsAtHETSites.py
TitanCNA_rscript: code/titanCNA_v1.15.0_TenX.R
TitanCNA_selectSolutionRscript: code/selectSolution.R
TitanCNA_combineTitanIchorCNA:  code/combineTITAN-ichor.R
```

### 3. Path to R package files
Specify the directory in which [TitanCNA](https://github.com/gavinha/TitanCNA) and [ichorCNA](https://github.com/broadinstitute/ichorCNA) are installed.  
*Set these if the R files in these libraries have been modified or updated but not yet installed or updated in R*.
```
TitanCNA_libdir:  /path/to/TitanCNA/ ## optional
ichorCNA_libdir:  /path/to/ichorCNA/ ## optional
```

### 4. Reference files and settings
Global reference files used by many of the `snakefiles` and scripts.  
- `snpVCF` you can download the HapMap file (used for filtering heterozygous SNPs) here: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0  
- `genomeStyle` specifies the chromosome naming convention to used for **output** files. Input files can be any convention as long as it is the same genome build. Only use `UCSC` (e.g. chr1) or `NCBI` (e.g. 1). 
- `sex` set to `male` or `female`, otherwise `None` if both females and males are in sample set.
```
genomeBuild: hg38
genomeStyle:  UCSC
snpVCF:  /path/to/hapmap_3.3.hg38.vcf.gz ## optional
cytobandFile:  data/cytoBand_hg38.txt # only need if hg38
centromere:  data/GRCh38.GCA_000001405.2_centromere_acen.txt
sex:  male   # use None if both females and males are in sample set
```

### 5. Long Ranger filenames
Set this to the filenames that are used for the BAM and variant files generated by Long Ranger. The current filenames are ones generated by [Long Ranger](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/what-is-long-ranger) v2.2.2
```
bamFileName:  phased_possorted_bam.bam
phaseVariantFileName:  phased_variants.vcf.gz
```

### 6. bxtools settings
```
bx_mapQual:  60  # mapping quality threshold
bx_bedFileRoot:  data/10kb_hg38/10kb_hg38  # bed files to specify intervals for analysis
``` 

### 7. [moleculeCoverage.snakefile](moleculeCoverage.snakefile) settings
Settings for the analysis of molecule coverage.  
- `molCov_minReadsPerBX` specify the minimum number of reads required for a barcode to be counted in the coverage.  
- `molCov_chrs` specifies the chromosomes to analyze; users do not need to be concerned about chromosome naming convention here as the code will handle it based on the `genomeStyle` set in the reference settings above.  
- The GC and Map wig files must have bin sizes that match the `bx_bedFileRoot` bed files. At the moment, only 10kb is supported.
```
molCov_minReadsPerBX:  2
molCov_chrs:  c(1:22, \"X\")
molCov_gcWig: data/gc_hg38_10kb.wig
molCov_mapWig:  data/map_hg38_10kb.wig
molCov_maxCN:  8
```

### 8. [getPhasedAlleleCounts.snakefile](getPhasedAlleleCounts.snakefile) settings: Heterozygous SNP 
Minimum thresholds used when determining heterozygous SNP sites from the Long Ranger `phased_variants.vcf.gz` file for the matched normal sample.
```
het_minVCFQuality:  100
het_minDepth:  10
het_minVAF:  0.25
```

### 9. [getPhasedAlleleCounts.snakefile](getPhasedAlleleCounts.snakefile) settings: Tumor allelic counts 
Minimum thresholds to use for extracting allelic read counts from the tumor sample.
```
het_minBaseQuality:  10
het_minMapQuality:  20
```

### 10. [TitanCNA.snakefile](TitanCNA.snakefile) settings
Most settings can be left as default.  
- `TitanCNA_maxNumClonalClusters` specifies the maximum number of clonal clusters to consider. For example, if set to 5, then 5 solutions are generated, each one considering a different number of cluster(s).  
- `TitanCNA_maxPloidy` specifies the maximum ploidy to initialize. This be set to either `2` (only considers diploid solutions), `3` (considers diploid and triploid, and usually accounts for tetraploid), or `4` (for diploid, triploid, tetraploid or higher ploidies). Usually, `3` is suitable for most samples unless you know that your samples are tetraploid or even higher. For example, if set to `3`, then solutions for diploid and triploid will be generated. [code/selectSolution.R](code/selectSolution.R) will try to select the optimal solution; however, users should inspect to make sure results are accurate.  
- `TitanCNA_numCores` specifies the number of cores to use on a single machine. `TitanCNA_pe` should also be set as to be consistent.
```
TitanCNA_maxNumClonalClusters: 2
TitanCNA_chrs:  c(1:22, \"X\")
TitanCNA_normalInit:  0.5
TitanCNA_maxPloidy:  3
TitanCNA_haplotypeBinSize: 1e5
TitanCNA_estimateNormal:  map
TitanCNA_estimatePloidy:  TRUE
TitanCNA_estimateClonality: TRUE
TitanCNA_alphaK:  10000
TitanCNA_alphaR:  5000
TitanCNA_txnExpLen: 1e15
TitanCNA_plotYlim:  c(-2,4)
TitanCNA_solutionThreshold: 0.05
TitanCNA_numCores:  1  #must match the settings for number of cpus for cluster settings
```

  
