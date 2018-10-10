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
		expand("results/titan/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt", tumor=config["pairings"], clustNum=CLUST[config["TitanCNA_maxNumClonalClusters"]], ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		#expand("results/titan/titanCNA_ploidy{ploidy}/", ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		expand("results/titan/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.ichor.seg.txt", tumor=config["pairings"], clustNum=CLUST[config["TitanCNA_maxNumClonalClusters"]], ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		expand("results/titan/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.ichor.cna.txt", tumor=config["pairings"], clustNum=CLUST[config["TitanCNA_maxNumClonalClusters"]], ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		"results/titan/optimalClusterSolution.txt",
		"results/titan/optimalClusterSolution/"
		
rule makeOutDir:
	params:
		mem=config["std_mem"],
		runtime=config["std_runtime"],
		pe=config["std_numCores"]
	output:
		"results/titan/titanCNA_ploidy{ploidy}/"
	shell:
		"mkdir -p {output}"
		
rule runTitanCNA:
	input:
		alleleCounts="results/phasedCounts/tumCounts/{tumor}.tumCounts.txt",
		corrDepth="results/moleculeCoverage/{tumor}/{tumor}.BXcounts.txt",
		ichorParam="results/moleculeCoverage/{tumor}/{tumor}.params.txt"	
	output:
		outRoot="results/titan/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}/",
		titan="results/titan/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt",
		param="results/titan/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.params.txt",
		segTxt="results/titan/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.segs.txt",
		seg="results/titan/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.seg"
	params:
		titanRscript=config["TitanCNA_rscript"],
		libdir=config["TitanCNA_libdir"],
		genomeBuild=config["genomeBuild"],
		genomeStyle=config["genomeStyle"],
		cytobandFile=config["cytobandFile"],
		numCores=config["TitanCNA_numCores"],
		normal=config["TitanCNA_normalInit"],
		chrs=config["TitanCNA_chrs"],
		sex=config["sex"],
		maxCN=config["TitanCNA_maxCN"],
		estimatePloidy=config["TitanCNA_estimatePloidy"],
		estimateClonality=config["TitanCNA_estimateClonality"],
		estimateNormal=config["TitanCNA_estimateNormal"],
		centromere=config["centromere"],
		haplotypeBinSize=config["TitanCNA_haplotypeBinSize"],
		alphaK=config["TitanCNA_alphaK"],
		alphaR=config["TitanCNA_alphaR"],
		#alleleModel=config["TitanCNA_alleleModel"],
		txnExpLen=config["TitanCNA_txnExpLen"],
		plotYlim=config["TitanCNA_plotYlim"],
		mem=config["TitanCNA_mem"],
		runtime=config["TitanCNA_runtime"],
		pe=config["TitanCNA_pe"]
	log:
		"logs/titan/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.log"
	shell:
		"Rscript {params.titanRscript} --id {wildcards.tumor} --hetFile {input.alleleCounts} --cnFile {input.corrDepth} --ichorParam {input.ichorParam} --numClusters {wildcards.clustNum} --numCores {params.numCores} --normal_0 {params.normal} --ploidy_0 {wildcards.ploidy} --chrs \"{params.chrs}\" --maxCN {params.maxCN} --sex {params.sex} --haplotypeBinSize {params.haplotypeBinSize} --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateClonality {params.estimateClonality}  --centromere {params.centromere} --genomeBuild {params.genomeBuild} --genomeStyle {params.genomeStyle} --libdir {params.libdir} --alphaK {params.alphaK} --alphaR {params.alphaR} --alleleModel Gaussian --txnExpLen {params.txnExpLen} --plotYlim \"{params.plotYlim}\" --cytobandFile {params.cytobandFile} --outFile {output.titan} --outSeg {output.segTxt} --outParam {output.param} --outIGV {output.seg} --outPlotDir {output.outRoot} > {log} 2> {log}"
	
				
rule combineTitanAndIchorCNA:
	input:
		titanSeg="results/titan/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.segs.txt", 
		titanBin="results/titan/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt",
		titanParam="results/titan/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.params.txt",
		ichorSeg="results/moleculeCoverage/{tumor}/{tumor}.seg.txt",
		ichorBin="results/moleculeCoverage/{tumor}/{tumor}.cna.seg",
		ichorParam="results/moleculeCoverage/{tumor}/{tumor}.params.txt"
	output:
		segFile="results/titan/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.ichor.seg.txt",
		binFile="results/titan/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.ichor.cna.txt",
	params:
		combineScript=config["TitanCNA_combineTitanIchorCNA"],
		libdir=config["TitanCNA_libdir"],
		centromere=config["centromere"],
		sex=config["sex"],
		mem=config["std_mem"],
		runtime=config["std_runtime"],
		pe=config["std_numCores"]
	log:
		"logs/titan/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.combineTitanIchorCNA.log"
	shell:
		"Rscript {params.combineScript} --libdir {params.libdir} --titanSeg {input.titanSeg} --titanBin {input.titanBin} --titanParam {input.titanParam} --ichorSeg {input.ichorSeg} --ichorBin {input.ichorBin} --ichorParam {input.ichorParam} --sex {params.sex} --outSegFile {output.segFile} --outBinFile {output.binFile} --centromere {params.centromere} > {log} 2> {log}"	


rule selectSolution:
	input:
		ploidyDirs=expand("results/titan/titanCNA_ploidy{ploidy}/", ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		resultFiles=expand("results/titan/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.ichor.seg.txt", tumor=config["pairings"], clustNum=CLUST[config["TitanCNA_maxNumClonalClusters"]], ploidy=PLOIDY[config["TitanCNA_maxPloidy"]])
	output:
		solutionsTxt="results/titan/optimalClusterSolution.txt",
	params:
		solutionRscript=config["TitanCNA_selectSolutionRscript"],
		threshold=config["TitanCNA_solutionThreshold"],
		mem=config["std_mem"],
		runtime=config["std_runtime"],
		pe=config["std_numCores"]
	log:
		"logs/titan/optSolution/selectSolution.log"
	shell:
		"""
		if [ -d results/titan/titanCNA_ploidy3/ ]; then
			ploidyRun3=results/titan/titanCNA_ploidy3/
		else
			ploidyRun3=NULL
		fi
		if [ -d results/titan/titanCNA_ploidy4/ ]; then
			ploidyRun4=results/titan/titanCNA_ploidy4/
		else
			ploidyRun4=NULL
		fi
		Rscript {params.solutionRscript} --ploidyRun2 {input.ploidyDirs[0]} --ploidyRun3 $ploidyRun3 --ploidyRun4 $ploidyRun4 --threshold {params.threshold} --outFile {output.solutionsTxt} > {log} 2> {log}
		"""
	
rule copyOptSolution:
	input:
		"results/titan/optimalClusterSolution.txt"
	output:
		"results/titan/optimalClusterSolution/"
	params:
		mem=config["std_mem"],
		runtime=config["std_runtime"],
		pe=config["std_numCores"]
	log:
		"logs/titan/optSolution/copyOptSolution.log"
	shell:
		"""
		curDir=`pwd`
		for i in `cut -f11 {input} | grep -v "path"`;
		do
			echo -e "Creating sym links for $curDir/${{i}} to {output}"
			ln -s ${{curDir}}/${{i}}* {output}
		done		
		"""
	
	
	
	
		