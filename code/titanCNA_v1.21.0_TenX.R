#' titanCNA_v1.15.0_TenX.R
#' author: Gavin Ha 
#' Dana-Farber Cancer Institute
#' Broad Institute
#' contact: <gavinha@gmail.com> or <gavinha@broadinstitute.org>
#' date:	  May 24, 2017
#' Notes: This script was tested for TitanCNA v1.15.0 and higher

library(optparse)

option_list <- list(
	make_option(c("--id"), type = "character", help = "Sample ID"),
	make_option(c("--hetFile"), type = "character", help = "File containing allelic read counts at HET sites. [Required]"),
	make_option(c("--cnFile"), type = "character", help = "File containing normalized coverage as log2 ratios. [Required]"),
	make_option(c("--outDir"), type = "character", help = "Output directory to output the results. [Required]"),
	make_option(c("--ichorParamFile"), type = "character", default = NULL, help = "ichorCNA parameter file to use to determine sex if not provided. [Default: %default]"),
	make_option(c("--numClusters"), type = "integer", default = 1, help = "Number of clonal clusters. [Default: 1]"),
	make_option(c("--numCores"), type = "integer", default = 1, help = "Number of cores to use. [Default: %default]"),
	make_option(c("--ploidy_0"), type = "numeric", default = 2, help = "Initial ploidy value; float [Default: %default]"),
	make_option(c("--estimatePloidy"), type = "logical", default = TRUE, help = "Estimate ploidy; TRUE or FALSE [Default: %default]"),
	make_option(c("--normal_0"), type = "numeric", default = 0.5, help = "Initial normal contamination (1-purity); float [Default: %default]"),
	make_option(c("--estimateNormal"), type = "character", default = "map", help = "Estimate normal contamination method; string {'map', 'fixed'} [Default: %default]"),
	make_option(c("--estimateClonality"), type="logical", default=TRUE, help="Estimate cellular prevalence. [Default: %default]"),
	make_option(c("--maxCN"), type = "integer", default = 8, help = "Maximum number of copies to model; integer [Default: %default]"),
	make_option(c("--alphaK"), type = "numeric", default = NA, help = "Hyperparameter of the prior on Gaussian variance for log ratio data; for WES, use 2500; for WGS, use 10000; float [Default: %default]"),
	make_option(c("--diploidStrength"), type = "numeric", default = 0, help = "Multiplier for likelihood model to favor diploid; 0 - no effect, 3 - strongly favor diploid. [Default: %default]"),
	make_option(c("--alleleModel"), type = "character", default = "binomial", help = "Emission density to use for allelic input data (binomial or Gaussian). [Default: %default]"),
	make_option(c("--alphaR"), type = "numeric", default = NA, help = "Hyperparaemter on the Gaussian variance for allelic fraction data; used if --alleleModel=\"Gaussian\". [Defaule: %default]"),
	make_option(c("--txnExpLen"), type = "numeric", default = 1e20, help = "Expected length of segments; higher leads to longer (less sensitive) segments; float [Default: %default]"),
	make_option(c("--txnZStrength"), type = "numeric", default = 1, help = "Expected length of clonal cluster segmentation (factor of txnExpLen); float [Default: %default]"),
	make_option(c("--minDepth"), type = "integer", default = 10, help = "Minimum read depth of a HET site to include in analysis; integer [Default: %default]"),
	make_option(c("--maxDepth"), type = "integer", default = 10000, help = "Maximum read depth of a HET site to include in analysis; integer [Default: %default]"),
	make_option(c("--skew"), type = "numeric", default=0, help = "Allelic reference skew for all states except heterozygous states (e.g. 1:1, 2:2, 3:3). Value is additive to baseline allelic ratios. float [Default: %default]"),
	make_option(c("--hetBaselineSkew"), type="numeric", default=NULL, help="Allelic reference skew for heterozygous states (e.g. 1:1, 2:2, 3:3). Value is the additive to baseline allelic ratios. float [Default: %default]"), 
	make_option(c("--genomeStyle"), type = "character", default = "NCBI", help = "NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have \"chr\" string. [Default: %default]"),
	make_option(c("--genomeBuild"), type = "character", default = "hg38", help="Genome build to use; will load Seqinfo from GenomeInfoDb."),
	make_option(c("--chrs"), type = "character", default = "c(1:22, 'X')", help = "Chromosomes to analyze; string [Default: %default"),
	make_option(c("--sex"), type = "character", default = "None", help = "User specified sex: male or female or None [Default: %default]"),
	make_option(c("--cytobandFile"), type = "character", default = "", help = "Cytoband file should be provided only if reference genome is hg38."),
	make_option(c("--mapWig"), type = "character", default = NULL, help = "Mappability score file for bin sizes matching cnfile. [Default: %default]"),
	make_option(c("--mapThres"), type = "numeric", default = 0.9, help = "Minimum mappability score threshold to use; float [Default: %default]"),
	make_option(c("--centromere"), type = "character", default=NULL, help = "Centromere gap file. [Default: %default]"),
	make_option(c("--libdir"), type = "character", default=NULL, help = "Directory containing source code. Specify if changes have been made to source code and want to over-ride package code. [Default: %default]"),
	make_option(c("--outFile"), type = "character", default = NULL, help = "Output file to write position-level file. (default uses extension: *.titan.txt]"),
	make_option(c("--outSeg"), type = "character", default = NULL, help = "Output file to write detailed segments. (default uses extension: *.segs.txt]"),
	make_option(c("--outIGV"), type = "character", default = NULL, help = "Output file to write segments for loading into IGV. (default uses extension: *.seg]"),
	make_option(c("--outParam"), type = "character", default = NULL, help = "Output file to write parameters. [Default: %default]"),
	make_option(c("--outPlotDir"), type = "character", default = NULL, help = "Output directory to save plots. [Default: %default]"),
	make_option(c("--plotYlim"), type = "character", default = "c(-2,4)", help = "The Y-axis limits to use for plotting log ratio coverage results. [Default: %default]"),
	make_option(c("--haplotypeBinSize"), type = "integer", default = 1e5, help = "Bin size for summarizing phased haplotypes. [Default: %default]"),
	make_option(c("--phaseSummarizeFun"), type = "character", default = "sum", help = "Strategy to summarize specified haplotype bins: mean, sum, SNP. [Default: %default]")
)

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
print(opt)

library(TitanCNA)
library(data.table)
library(GenomicRanges)
#library(GenomeInfoDb)
library(dplyr)
library(doMC)
library(SNPchip)
sessionInfo()
options(bitmapType='cairo', scipen=0)

libdir <- opt$libdir
if (!is.null(libdir) && libdir != "None"){
	source(paste0(libdir, "/R/plotting.R"))
	source(paste0(libdir, "/R/utils.R"))
	source(paste0(libdir, "/R/hmmClonal.R"))
	source(paste0(libdir, "/R/paramEstimation.R"))
	source(paste0(libdir, "/R/correction.R"))
	source(paste0(libdir, "/R/haplotype.R"))
}

id <- opt$id
hetfile <- opt$hetFile
cnfile <- opt$cnFile
ichorParamFile <- opt$ichorParamFile
numClusters <- opt$numClusters
numCores <- opt$numCores
ploidy_0 <- opt$ploidy_0
boolEstPloidy <- opt$estimatePloidy
norm_0 <- opt$normal_0
normEstMeth <- opt$estimateNormal
estimateS <- opt$estimateClonality
maxCN <- opt$maxCN
alphaK <- opt$alphaK
diploidStrength <- opt$diploidStrength
alleleEmissionModel <- opt$alleleModel
alphaR <- opt$alphaR
txn_exp_len <- opt$txnExpLen
txn_z_strength <- opt$txnZStrength
mapThres <- opt$mapThres
minDepth <- opt$minDepth
maxDepth <- opt$maxDepth
skew <- opt$skew
hetBaselineSkew <- opt$hetBaselineSkew
chrs <- eval(parse(text = opt$chrs))
genomeStyle <- opt$genomeStyle
genomeBuild <- opt$genomeBuild
cytobandFile <- opt$cytobandFile
sex <- opt$sex
mapWig <- opt$mapWig
centromere <- opt$centromere
haplotypeBinSize <- opt$haplotypeBinSize
phaseSummarizeFun <- opt$phaseSummarizeFun
outdir <- opt$outDir
outfile <- opt$outFile
outparam <- opt$outParam
outseg <- opt$outSeg
outigv <- opt$outIGV
outplot <- opt$outPlotDir
plotYlim <- eval(parse(text = opt$plotYlim))

## check arguments ##
if (!normEstMeth %in% c("map", "fixed")){
	stop("--estimateNormal must be \"map\" or \"fixed\"")
}

### SETUP OUTPUT FILE NAMES ###
numClustersStr <- as.character(numClusters)
if (numClusters < 10) { 
	numClustersStr <- paste0("0", numClusters)
}
if (is.null(outfile)){
	outfile <- paste0(outdir, "/", id, "_cluster", numClustersStr, ".titan.txt")
}
if (is.null(outparam)){
	outparam <- gsub(".titan.txt", ".params.txt", outfile)
}
if (is.null(outseg)){
	outseg <- gsub(".titan.txt", ".segs.txt", outfile)
}
if (is.null(outigv)){
	outigv <- gsub(".titan.txt", ".seg", outfile)
}
if (is.null(outplot)){
	outplot <- paste0(outdir, "/", id, "_cluster", numClustersStr, "/")
	dir.create(outplot)
}
outImage <- gsub(".titan.txt", ".RData", outfile)

## set up chromosome naming convention ##
bsg <- paste0("BSgenome.Hsapiens.UCSC.", genomeBuild)
if (!require(bsg, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE)) {
	seqinfo <- Seqinfo(genome=genomeBuild)
} else {
	seqinfo <- seqinfo(get(bsg))
}
seqlevelsStyle(seqinfo) <- genomeStyle
seqlevelsStyle(chrs) <- genomeStyle

## set up sex for analysis
if (sex == "None" || is.null(sex)){
	if (file.exists(ichorParamFile)){
		ichorParams <- read.delim(ichorParamFile, header=T, as.is=T)
		sex <- as.character(ichorParams[3, 2])
	}else{
		stop("TitanCNA: one of sex (male/female) or ichorParamFile must be provided.")
	}
	if (sex == "unknown" || sex == "Unknown"){
		sex <- "male"
	} 
}
message("Analyzing sample as sex: ", sex)
## exclude chrX if sex==male ##
if (sex == "male" || sex == "Male" || sex == "MALE" || sex == "None" || is.null(sex)){
	chrs <- chrs[!grepl("X", chrs)]
}

#pseudo_counts <- 1e-300
centromereFlank <- 100000
maxI <- 50

message('Running TITAN...')
save.image(file=outImage)
#### LOAD DATA ####
data <- loadHaplotypeAlleleCounts(hetfile, cnfile, chrs = chrs,
	seqinfo = seqinfo, genomeStyle = genomeStyle,
	haplotypeBinSize = haplotypeBinSize, fun=phaseSummarizeFun,
	map=mapWig, mapThres=mapThres, centromere=centromere, minDepth=minDepth, maxDepth=maxDepth) 
data <- data$haplotypeData

## reassign chromosomes ##
chrs <- unique(data$chr)


#### LOAD PARAMETERS ####
message('titan: Loading default parameters')
params <- loadDefaultParameters(copyNumber=maxCN,numberClonalClusters=numClusters, 
								skew=skew, hetBaselineSkew=hetBaselineSkew, 
								alleleEmissionModel = alleleEmissionModel, data=data)

#### MODEL SELECTION USING EM (FWD-BACK) TO SELECT NUMBER OF CLUSTERS ####
registerDoMC()
options(cores=numCores)
message("Using ",getDoParWorkers()," cores.")
K <- length(params$genotypeParams$rt)
N <- nrow(data)
params$ploidyParams$phi_0 <- ploidy_0
params$normalParams$n_0 <- norm_0
#params$ploidyParams$alphaPHyper <- params$ploidyParams$alphaPHyper
#params$ploidyParams$betaPHyper <- params$ploidyParams$betaPHyper
params$genotypeParams$var_0 <- rep(var(data$logR), K)
params$genotypeParams$varR_0 <- rep(var(data$ref / data$tumDepth, na.rm = TRUE), K)
params$genotypeParams$betaKHyper <- rep(2, K)
params$genotypeParams$betaRHyper <- rep(2, K)
params$genotypeParams$alphaKHyper <- params$genotypeParams$betaKHyper / params$genotypeParams$var_0
params$genotypeParams$alphaRHyper <- params$genotypeParams$betaRHyper / params$genotypeParams$varR_0
if (diploidStrength > 0){
	params$genotypeParams$alphaKHyper[c(1,5:K)] <- params$genotypeParams$alphaKHyper[c(1,5:K)] * sqrt(N) / 10 * diploidStrength
	params$genotypeParams$betaKHyper[c(1,5:K)] <- params$genotypeParams$betaKHyper[c(1,5:K)] * sqrt(N) / 2 * diploidStrength
}
convergeParams <- runEMclonalCN(data, params=params,
                                maxiter=maxI,maxiterUpdate=15,
                                txnExpLen=txn_exp_len,txnZstrength=txn_z_strength,
                                useOutlierState=FALSE,
                                normalEstimateMethod=normEstMeth,estimateS=estimateS,
                                estimatePloidy=boolEstPloidy, pseudoCounts=0)
    
#### COMPUTE OPTIMAL STATE PATH USING VITERBI ####
message("Using ",getDoParWorkers()," cores.")
optimalPath <- viterbiClonalCN(data,convergeParams)
save.image(file=outImage)
#### PRINT RESULTS TO FILES ####
results <- outputTitanResults(data,convergeParams,optimalPath,is.haplotypeData=TRUE,
			filename=NULL,posteriorProbs=FALSE,subcloneProfiles=TRUE, recomputeLogLik = FALSE,
			proportionThreshold = 0.05, proportionThresholdClonal = 0.05, verbose=FALSE)
convergeParams <- results$convergeParams
results <- results$corrResults
norm <- tail(convergeParams$n,1)
ploidy <- tail(convergeParams$phi,1)

# save specific objects to a file
convergeParams$rhoG <- NULL; convergeParams$rhoZ <- NULL
save.image(file=outImage)

#### OUTPUT SEGMENTS ####
segs <- outputTitanSegments(results, id, convergeParams, filename = NULL, igvfilename = outigv)
corrIntCN.results <- correctIntegerCN(results, segs, 1 - norm, ploidy, maxCNtoCorrect.autosomes = maxCN, 
			maxCNtoCorrect.X = NULL, correctHOMD = TRUE, minPurityToCorrect = 0.2, gender = sex, chrs = chrs)
results <- corrIntCN.results$cn
segs <- corrIntCN.results$segs
message("Writing results to ", outfile, ",\n\t", outseg, ",\n\t", outparam)
write.table(results, file = outfile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(segs, file = outseg, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
outputModelParameters(convergeParams, results, outparam)

save.image(file=outImage)

#### PLOT RESULTS ####
dir.create(outplot)
if (genomeBuild == "hg38" && file.exists(cytobandFile)){
	cytoband <- fread(cytobandFile)
	names(cytoband) <- c("chrom", "start", "end", "name", "gieStain")
	cytoband$chrom <- setGenomeStyle(cytoband$chrom, genomeStyle = genomeStyle, filterExtraChr = FALSE)
	cytoband <- cytoband[chrom %in% chrs]
	cytoband <- as.data.frame(cytoband)
}

if (sex == "male" || sex == "Male" || sex == "MALE" || sex == "None" || is.null(sex)){
	chrsToPlot <- chrs[!grepl("X", chrs)]
}else{
	chrsToPlot <- chrs
}
for (chr in chrsToPlot){
	if (grepl("chr", chr)){ 
		chrStr <- chr 
	}else{
		chrStr <- paste0("chr", chr)
	}	
	outfig <- paste0(outplot, "/", id, "_cluster", numClustersStr, "_", chrStr, ".png")
	png(outfig,width=1200,height=1000,res=100)
	par(mfrow=c(5,1))  

	if (genomeStyle == "UCSC"){
		titleStr <- chr
	}else{
		titleStr <- paste("Chr ",chr,sep="")
	}
	plotCNlogRByChr(results, chr=chr, segs = NULL, ploidy=ploidy, normal = norm, geneAnnot=NULL,  cex.axis=1.5, 
					cex.lab=1.5, ylim=plotYlim, cex=0.5, xlab="", main=titleStr)
	plotAllelicRatio(results, chr=chr, geneAnnot=NULL, spacing=4, cex.axis=1.5, cex.lab=1.5, 
					ylim=c(0,1), xlab="", cex=0.5)
	#plotHaplotypeFraction(data, type = "AllelicRatio", colType = "haplotypes", xlab="", chr="6", cex=0.25)
	#plotHaplotypeFraction(results, chr=chr, type = "AllelicRatio", colType = "haplotype", 
	#  xlab="", cex=0.5, cex.axis=1.5, cex.lab=1.5)
	plotHaplotypeFraction(results, chr=chr, resultType = "HaplotypeRatio", colType = "Haplotypes", 
	  xlab="", ylim=c(0,1), cex=0.5, cex.axis=1.5, cex.lab=1.5)
	#plotHaplotypeFraction(results, chr, type = "HaplotypeRatio", colType = "titan", 
	#  xlab="", cex=0.5, cex.axis=1.5, cex.lab=1.5)
	#plotHaplotypeFraction(data, chr=chr, type = "HaplotypeRatio", colType = "haplotype", xlab="", ylim=c(0,1), cex=0.5, cex.axis=1.5, cex.lab=1.5)
	maxCorCN <- segs[Chromosome==chr, max(Corrected_Copy_Number, na.rm = TRUE)]
	plotSegmentMedians(segs, chr=chr, resultType = "LogRatio", plotType = "CopyNumber", 
				plot.new=TRUE, ylim=c(0,maxCorCN), xlab="", cex.axis=1.5, cex.lab=1.5, spacing=4)
	plotClonalFrequency(results, chr, normal=norm, geneAnnot=NULL, spacing=4, 
					cex.axis=1.5, ylim=c(0,1), xlab="", cex=0.5, cex.axis=1.5, cex.lab=1.5,
					main=paste("Chr ",chr,sep=""))
  
  	par(xpd = NA)
  	if (genomeBuild == "hg38" && file.exists(cytobandFile)){
  		sl <- seqlengths(seqinfo[chr])
  		pI <- plotIdiogram.hg38(chr, cytoband=cytoband, seqinfo=seqinfo, xlim=c(0, max(sl)), unit="bp", label.y=-0.35, new=FALSE, ylim=c(-0.2,-0.1))	
  	}else{
  		pI <- plotIdiogram(chr, build="hg19", unit="bp", label.y=-0.35, new=FALSE, ylim=c(-0.2,-0.1))	
  	}
	dev.off()

}

################################################
############## GENOME WIDE PLOTS ###############
################################################
outFile <- paste0(outplot, "/", id, "_cluster", numClustersStr, "_CNA.png")
png(outFile,width=1200,height=400,res=100)
#pdf(outFile,width=20,height=6)
plotCNlogRByChr(dataIn=results, chr=chrs, segs = segs, ploidy=ploidy,  normal = norm, geneAnnot=genes, spacing=4, main=id, xlab="", ylim=plotYlim, cex=0.5, cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
dev.off()

outFile <- paste0(outplot, "/", id, "_cluster", numClustersStr, "_LOH.png")
png(outFile,width=1200,height=400,res=100)
#pdf(outFile,width=20,height=6)
plotAllelicRatio(dataIn=results, chr=chrs, geneAnnot=genes, spacing=4, main=id, xlab="", ylim=c(0,1), cex=0.5, cex.axis=1.5, cex.lab=1.5, cex.main=1.5)	
dev.off()

outFile <- paste0(outplot, "/", id, "_cluster", numClustersStr, "_HAP-PHASE.png")
png(outFile,width=1200,height=400,res=100)
#pdf(outFile,width=20,height=6)
par(mfrow=c(1,1))
plotHaplotypeFraction(results, chr=chrs, resultType = "HaplotypeRatio", colType = "Haplotypes", xlab="", ylim=c(0,1), cex=0.5, cex.axis=1.5, cex.lab=1.5)
dev.off()

outFile <- paste0(outplot, "/", id, "_cluster", numClustersStr, "_HAP-FRAC.png")
png(outFile,width=1200,height=400,res=100)
#pdf(outFile,width=20,height=6)
par(mfrow=c(1,1))
plotHaplotypeFraction(results, chr=chrs, resultType = "HaplotypeRatio", colType = "CopyNumber", xlab="", ylim=c(0,1), cex=0.5, cex.axis=1.5, cex.lab=1.5)
dev.off()

outFile <- paste0(outplot, "/", id, "_cluster", numClustersStr, "_CP.png")
png(outFile,width=1200,height=400,res=100)
#pdf(outFile,width=20,height=6)
plotClonalFrequency(dataIn=results, chr=chrs, norm, geneAnnot=genes, spacing=4, main=id, xlab="", ylim=c(0,1), cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
dev.off()

outFile <- paste0(outplot, "/", id, "_cluster", numClustersStr, "_LOH-SEG.png")
maxCorCN <- segs[Chromosome %in% chrs, max(Corrected_Copy_Number, na.rm = TRUE)]
png(outFile,width=1200,height=400,res=100)
#pdf(outFile, width=20, height=6)
plotSegmentMedians(dataIn=segs, chr=chrs, resultType = "AllelicRatio", plotType = "CopyNumber", plot.new=T, ylim=c(0, maxCorCN), cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
dev.off()


if (as.numeric(numClusters) <= 2){
	outFile <- paste0(outplot, "/", id, "_cluster", numClustersStr, "_subclone.png")
	png(outFile,width=1200,height=400,res=100)
	#pdf(outFile,width=20,height=6)
	plotSubcloneProfiles(dataIn=results, chr=chrs, cex = 0.5, spacing=4, main=id, cex.axis=1.5, xlab="")
	dev.off()
}

