#' combineTITAN-ichor.R
#' author: Gavin Ha 
#' institution: Fred Hutchinson Cancer Research Center
#' contact: <gha@fredhutch.org>
#' date: July 23, 2018

#' requires R-3.3+
#' @import data.table
#' @import GenomicRanges
#' @import stringr
#' @import optparse

library(optparse)

option_list <- list(
	make_option(c("--titanSeg"), type="character", help="TitanCNA segs.txt file. Required."),
	make_option(c("--titanBin"), type="character", help="TitanCNA titan.txt file. Required."),
	make_option(c("--titanParams"), type="character", help="TitanCNA params.txt file. Required."),
  	make_option(c("--ichorSeg"), type="character", help="ichorCNA segs.txt file. Required."),
	make_option(c("--ichorBin"), type="character", help="ichorCNA cna.seg file. Required."),
	make_option(c("--ichorParams"), type="character", help="ichorCNA params.txt file. Required."),
	make_option(c("--sex"), type="character", default="female", help="female or male. Default [%default]."),
	make_option(c("--libdir"), type="character", help="TitanCNA directory path to source R files if custom changes made."),
  	make_option(c("--outSegFile"), type="character", help="New combined segment file. Required"),
  	make_option(c("--outBinFile"), type="character", help="New combined bin-level file. Required"),
  	make_option(c("--centromere"), type="character", default=NULL, help="Centromere table.")  
  )

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
print(opt)

titanSeg <- opt$titanSeg
titanBin <- opt$titanBin
titanParams <- opt$titanParams
ichorSeg <- opt$ichorSeg
ichorBin <- opt$ichorBin
ichorParams <- opt$ichorParams
gender <- opt$sex
outSegFile <- opt$outSegFile
outBinFile <- opt$outBinFile
centromere <- opt$centromere
libdir <- opt$libdir

library(TitanCNA)
library(stringr)
library(data.table)
library(GenomicRanges)

if (!is.null(libdir) && libdir != "None"){
	source(paste0(libdir, "/R/utils.R"))
}

options(stringsAsFactors=F, width=150, scipen=999)

## copy number state mappings ##
ichorCNmap <- list("0"="HOMD", "1"="DLOH", "2"="NEUT", "3"="GAIN", "4"="AMP", "5"="AMP")
maxichorcn <- 5

## load segments 
titan <- fread(titanSeg)
ichor <- fread(ichorSeg)
setnames(ichor, c("ID", "chrom", "start", "end", "num.mark", "seg.median.logR", "copy.number", "call"), 
		c("Sample", "Chromosome", "Start_Position.bp.", "End_Position.bp.", 
		  "Length.snp.", "Median_logR", "Copy_Number", "TITAN_call"))

## load data points ##
titan.cn <- fread(titanBin)
titan.cn <- cbind(Sample=titan[1,Sample], titan.cn)
#titan.cn[, chr := as.character(Chr)]
id <- titan[1, Sample]
ichor.cn <- fread(ichorBin)
ichor.cn <- cbind(Sample = id, ichor.cn)
#ichor.cn[, CopyNumber := state - 1]
ichor.cn[, Position := start]
setnames(ichor.cn, c("chr", "start", paste0(id,".copy.number"), paste0(id,".event"), paste0(id,".logR"), "end"), 
		c("Chr", "Start", "Copy_Number", "TITANcall", "LogRatio", "End"))


## get chromosome style
genomeStyle <- seqlevelsStyle(titan$Chr)
chrs <- c(1:22, "X")
seqlevelsStyle(chrs) <- genomeStyle

## load parameters ##
params <- read.delim(titanParams, header=F, as.is=T)
purity <- 1 - as.numeric(params[1,2])
ploidyT <- as.numeric(params[2,2])
ploidy <- purity * ploidyT + (1-purity) * 2
params.ichor <- read.delim(ichorParams, header=T, as.is=T)

## get gender
if (is.null(gender) || gender == "None"){
	gender <- params.ichor[3, 2]
}

## get bin overlap with SNPs - include ichor bins even if no SNPs overlap 
titan.gr <- titan.cn[, .(Chr, Position)]
titan.gr[, Start := Position]; titan.gr[, End := Position]
titan.gr <- as(titan.gr, "GRanges")
ichor.gr <- as(ichor.cn, "GRanges")
hits <- findOverlaps(query = titan.gr, subject = ichor.gr)
titan.cn[queryHits(hits), Start := ichor.cn[subjectHits(hits), Start]]
titan.cn[queryHits(hits), End := ichor.cn[subjectHits(hits), End]]
titan.ichor.cn <- merge(titan.cn, ichor.cn, by=c("Sample", "Chr", "Start", "End"), all=T, suffix=c("",".ichor"))
titan.ichor.cn[is.na(LogRatio), LogRatio := LogRatio.ichor] # assign ichor log ratio to missing titan SNPs
titan.ichor.cn <- titan.ichor.cn[, -c(grep("ichor", colnames(titan.ichor.cn),value=T)), with=F] 

## combine TITAN (chr1-22) and ichorCNA (chrX) segments and bin/SNP level data ##
## if male only ##
if (gender == "male"){
	cn <- rbind(titan.ichor.cn[Chr %in% chrs[1:22]], ichor.cn[Chr == chrs[grep("X", chrs)]], fill = TRUE)
	segs <- rbind(titan[Chromosome %in% chrs[1:22]], ichor[Chromosome == chrs[grep("X", chrs)]], fill = TRUE)
}else{
	cn <- titan.ichor.cn
	segs <- titan
}

## sort column order
setnames(segs, c("Start_Position.bp.", "End_Position.bp."), c("Start", "End"))
cols <- c("Sample", "Chr", "Position", "Start", "End")
setcolorder(cn, c(cols, colnames(cn)[!colnames(cn) %in% cols]))

## get major/minor CN from segs and place in SNP/level data ##
cn.gr <- cn[, .(Chr, Start, End)]
cn.gr <- as(na.omit(cn.gr), "GRanges")
segs.gr <- as(segs, "GRanges")
hits <- findOverlaps(query = cn.gr, subject = segs.gr)
cn[queryHits(hits), MajorCN := segs[subjectHits(hits), MajorCN]]
cn[queryHits(hits), MinorCN := segs[subjectHits(hits), MinorCN]]
cn[is.na(CopyNumber), CopyNumber := MajorCN + MinorCN]

## correct copy number beyond maximum CN state based on purity and logR
correctCN <- correctIntegerCN(cn, segs, purity, ploidyT, maxCNtoCorrect.autosomes = NULL, 
		maxCNtoCorrect.X = NULL, minPurityToCorrect = 0.2, gender = gender, chrs = chrs)
segs <- correctCN$segs
cn <- correctCN$cn
## extend segments to remove gaps
centromeres <- fread(centromere)
segs <- extendSegments(segs, removeCentromeres = TRUE, centromeres = centromeres, extendToTelomeres = FALSE,
	chrs = chrs, genomeStyle = genomeStyle)

## write segments to file ##
write.table(segs, file = outSegFile, col.names=T, row.names=F, quote=F, sep="\t")
write.table(cn, file = outBinFile, col.names=T, row.names=F, quote=F, sep="\t")
## write segments without germline SNPs
outSegNoSNPFile <- gsub(".txt", ".noSNPs.txt", outSegFile)
write.table(segs[, -c("Start.snp", "End.snp")], file = outSegNoSNPFile, col.names=T, row.names=F, quote=F, sep="\t")

outImageFile <- gsub(".segs.txt", ".RData", outSegFile)
save.image(outImageFile)
