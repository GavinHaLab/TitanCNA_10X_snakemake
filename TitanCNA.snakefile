configfile: "config/config.yaml"
configfile: "config/samples.yaml"

include: "moleculeCoverage.snakefile"
include: "getPhasedAlleleCounts.snakefile"
import os.path

#CHRS = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X']
CHRS = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']
CLUST = {1:[1], 2:[1,2], 3:[1,2,3], 4:[1,2,3,4], 5:[1,2,3,4,5], 6:[1,2,3,4,5,6], 7:[1,2,3,4,5,6,7], 8:[1,2,3,4,5,6,7,8], 9:[1,2,3,4,5,6,7,8,9], 10:[1,2,3,4,5,6,7,8,9,10]}
PLOIDY = {2:[2], 3:[2,3], 4:[2,3,4]}


rule all:
	input: 
		expand("results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt", tumor=config["pairings"], clustNum=CLUST[config["TitanCNA_maxNumClonalClusters"]], ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		#expand("results/titan/hmm/titanCNA_ploidy{ploidy}/", ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		expand("results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.ichor.seg.txt", tumor=config["pairings"], clustNum=CLUST[config["TitanCNA_maxNumClonalClusters"]], ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		expand("results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.ichor.cna.txt", tumor=config["pairings"], clustNum=CLUST[config["TitanCNA_maxNumClonalClusters"]], ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		"results/titan/hmm/optimalClusterSolution.txt",
		#"results/titan/hmm/optimalClusterSolution/"
		
rule makeOutDir:
	output:
		"results/titan/hmm/titanCNA_ploidy{ploidy}/"
	shell:
		"mkdir -p {output}"
		
rule runTitanCNA:
	input:
		alleleCounts="results/phasedCounts/tumCounts/{tumor}.tumCounts.txt",
		corrDepth="results/moleculeCoverage/{tumor}/{tumor}.BXcounts.txt"		
	output:
		outRoot="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}/",
		titan="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt",
		param="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.params.txt",
		segTxt="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.segs.txt",
		seg="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.seg"
	params:
		titanRscript=config["TitanCNA_rscript"],
		libdir=config["TitanCNA_libdir"],
		genomeBuild=config["genomeBuild"],
		genomeStyle=config["genomeStyle"],
		numCores=config["TitanCNA_numCores"],
		normal=config["TitanCNA_normalInit"],
		chrs=config["TitanCNA_chrs"],
		gender=config["gender"],
		estimatePloidy=config["TitanCNA_estimatePloidy"],
		estimateClonality=config["TitanCNA_estimateClonality"],
		estimateNormal=config["TitanCNA_estimateNormal"],
		centromere=config["centromere"],
		haplotypeBinSize=config["TitanCNA_haplotypeBinSize"],
		alphaK=config["TitanCNA_alphaK"],
		alphaR=config["TitanCNA_alphaR"],
		#alleleModel=config["TitanCNA_alleleModel"],
		txnExpLen=config["TitanCNA_txnExpLen"],
		plotYlim=config["TitanCNA_plotYlim"]
	log:
		"logs/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.log"
	shell:
		"Rscript {params.titanRscript} --id {wildcards.tumor} --hetFile {input.alleleCounts} --cnFile {input.corrDepth} --numClusters {wildcards.clustNum} --numCores {params.numCores} --normal_0 {params.normal} --ploidy_0 {wildcards.ploidy} --chrs \"{params.chrs}\" --gender {params.gender} --haplotypeBinSize {params.haplotypeBinSize} --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateClonality {params.estimateClonality}  --centromere {params.centromere} --genomeBuild {params.genomeBuild} --genomeStyle {params.genomeStyle} --libdir {params.libdir} --alphaK {params.alphaK} --alphaR {params.alphaR} --alleleModel Gaussian --txnExpLen {params.txnExpLen} --plotYlim \"{params.plotYlim}\" --outFile {output.titan} --outSeg {output.segTxt} --outParam {output.param} --outIGV {output.seg} --outPlotDir {output.outRoot} > {log} 2> {log}"
	
				
rule combineTitanAndIchorCNA:
	input:
		titanSeg="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.segs.txt", 
		titanBin="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt",
		titanParam="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.params.txt",
		ichorSeg="results/moleculeCoverage/{tumor}/{tumor}.seg.txt",
		ichorBin="results/moleculeCoverage/{tumor}/{tumor}.cna.seg",
		ichorParam="results/moleculeCoverage/{tumor}/{tumor}.params.txt"
	output:
		segFile="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.ichor.seg.txt",
		binFile="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.ichor.cna.txt",
	params:
		combineScript=config["TitanCNA_combineTitanIchorCNA"],
		libdir=config["TitanCNA_libdir"],
		centromere=config["centromere"],
		gender=config["gender"]
	log:
		"logs/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.combineTitanIchorCNA.log"
	shell:
		"Rscript {params.combineScript} --libdir {params.libdir} --titanSeg {input.titanSeg} --titanBin {input.titanBin} --titanParam {input.titanParam} --ichorSeg {input.ichorSeg} --ichorBin {input.ichorBin} --ichorParam {input.ichorParam} --gender {params.gender} --outSegFile {output.segFile} --outBinFile {output.binFile} --centromere {params.centromere} > {log} 2> {log}"	


rule selectSolution:
	input:
		ploidyDirs=expand("results/titan/hmm/titanCNA_ploidy{ploidy}/", ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		resultFiles=expand("results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt", tumor=config["pairings"], clustNum=CLUST[config["TitanCNA_maxNumClonalClusters"]], ploidy=PLOIDY[config["TitanCNA_maxPloidy"]])
	output:
		solutionsTxt="results/titan/hmm/optimalClusterSolution.txt",
		solutionsDir="results/titan/hmm/optimalClusterSolution/"
	params:
		solutionRscript=config["TitanCNA_selectSolutionRscript"],
		threshold=config["TitanCNA_solutionThreshold"]
	log:
		"logs/titan/selectSolution.log"
	shell:
		"""
		if [ -d results/titan/hmm/titanCNA_ploidy3/ ]; then
			ploidyRun3=results/titan/hmm/titanCNA_ploidy3/
		else
			ploidyRun3=NULL
		fi
		if [ -d results/titan/hmm/titanCNA_ploidy4/ ]; then
			ploidyRun4=results/titan/hmm/titanCNA_ploidy4/
		else
			ploidyRun4=NULL
		fi
		Rscript {params.solutionRscript} --ploidyRun2 {input.ploidyDirs[0]} --ploidyRun3 $ploidyRun3 --ploidyRun4 $ploidyRun4 --threshold {params.threshold} --outFile {output.solutionsTxt} --outDir {output.solutionsDir} > {log} 2> {log}
		"""
	
	
		