#' getPhasedHETsitesFromLRVCF.R
#' author: Gavin Ha 
#' institution: Fred Hutchinson Cancer Research Center
#' contact: <gha@fredhutch.org>
#' date: July 19, 2018

#' requires R-3.3+
#' @import VariantAnnotation
#' @import TitanCNA

library(optparse)

option_list <- list(
  make_option(c("-i", "--inVCF"), type = "character", help = "LongRanger 2.1 phased variant result VCF file; typically has filename suffix \"*phased_variants.vcf.gz\". [Required]"),
  make_option(c("-q", "--minQuality"), type = "numeric", default = 100, help = "Heterozygous variants with QUAL greater than or equal to this value are considered. [Default: %default]"),
  make_option(c("-d", "--minDepth"), type = "numeric", default=10, help = "Heterozygous variants with read depth greater than or equal to this value are considered. [Default: %default]"),
  make_option(c("-v", "--minVAF"), type = "numeric", default=0.25, help = "Heterozygous variants with variant allele fraction or reference allele fraction greater than this value are considered. [Default: %default]"),
  make_option(c("--genomeBuild"), type = "character", default="hg19", help = "Genome build (e.g. hg19 or hg38"),
  make_option(c("--chrs"), type = "character", default="c(1:22,\"X\")", help = "Chromosomes to analyze. [Default: %default]"),
  make_option(c("--genomeStyle"), type = "character", default="NCBI", help = "Chr naming convention. NCBI (e.g. 1) or UCSC (e.g. chr1). Default: [%default]"),
  make_option(c("-s","--snpDB"), type="character", default=NULL, help= "VCF file of list of known SNPs (e.g. HapMap, 1000 Genomes, dbSNP)"),
  make_option(c("--altCountField"), type="character", default="AD", help="Alternate allele count field name in genotype. [Default: %default]"),
  make_option(c("-o","--outVCF"), type = "character", help = "Output VCF file with suffix \"_phasedHets.vcf\" [Required]"),
  make_option(c("--libdir"), type="character", default=NULL, help="Path to TitanCNA package directory. If provided, will source scripts after loading installed TitanCNA package.")
)

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
print(opt)

library(TitanCNA)
library(GenomeInfoDb)
library(VariantAnnotation)

vcfFile <- opt$inVCF
snpFile <- opt$snpDB
genomeBuild <- opt$genomeBuild
genomeStyle <- opt$genomeStyle
chrs <- eval(parse(text = opt$chrs))
minQUAL <- opt$minQuality
minDepth <- opt$minDepth
minVAF <- opt$minVAF
outVCF <- opt$outVCF
outVCFgz <- paste0(outVCF, ".gz")
libdir <- opt$libdir
altCountField <- opt$altCountField

if (!is.null(libdir) && libdir != "None"){
  source(paste0(libdir, "/R/haplotype.R"))
}

if (is.null(snpFile) || snpFile == "None"){
	snpFile <- NULL
} 

## set genome style for chromosome names
seqlevelsStyle(chrs) <- genomeStyle

filterFlags <- c("PASS", "10X_RESCUED_MOLECULE_HIGH_DIVERSITY")
#minQUAL <- 100
keepGenotypes <- c("1|0", "0|1", "0/1")



#vcfFiles <- list.files(LRdir, pattern = "_phased_variants.vcf.gz$", full.name = TRUE)

#for (i in 1:length(vcfFiles)){
#id <- gsub("_phased_variants.vcf.gz", "", basename(vcfFile))
hap <- getHaplotypesFromVCF(vcfFile, 
			chrs = chrs, build = genomeBuild, genomeStyle = genomeStyle,
			filterFlags = filterFlags, minQUAL = minQUAL, minDepth = minDepth,
			minVAF = minVAF, keepGenotypes = keepGenotypes, 
			altCountField = altCountField, snpDB = snpFile)

#outFile <- paste0(outDir, "/", id, "_phasedHets.vcf")
## remove BX genotype field to make things faster
geno(hap$vcf)$BX <- NULL

message("Writing to file: ", outVCF)
writeVcf(hap$vcf, filename = outVCF)
bgzipFile <- bgzip(outVCF, dest = outVCFgz, overwrite = TRUE)
indexTabix(bgzipFile, format = "vcf")

#	outFile <- paste0(outDir, "/", id, "_phasedGR.rds")
#	saveRDS(hap$geno, file = outFile)
								
#}