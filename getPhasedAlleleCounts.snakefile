configfile: "config/config.yaml"
configfile: "config/samples.yaml"

CHRS = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']
#CHRS = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X']

import glob
def getLRFullPath(base, filename):
  return glob.glob(''.join([base, "/*/outs/", filename]))

rule phasedCounts:
	input: 
		expand("results/phasedCounts/tumCounts/{tumor}.tumCounts.txt", tumor=config["pairings"])
		
rule getHETsites:
	input:
		lambda wildcards: getLRFullPath(config["samples"][config["pairings"][wildcards.tumor]], config["phaseVariantFileName"])
	output:
		"results/phasedCounts/hetPosns/{tumor}.phasedHetsFromNormal.vcf"
	params:
		getHETsitesScript=config["phaseCounts_hetSites_script"],
		snpDB=config["snpVCF"],
		genomeBuild=config["genomeBuild"],
		genomeStyle=config["genomeStyle"],
		minQual=config["het_minVCFQuality"],
		minDepth=config["het_minDepth"],
		minVAF=config["het_minVAF"],
		libdir=config["TitanCNA_libdir"],
		mem=config["het_mem"],
		runtime=config["het_runtime"],
		pe=config["std_numCores"]
	log:
		"logs/phasedCounts/hetPosns/{tumor}.phasedHETsites.log"
	shell:
		"Rscript {params.getHETsitesScript} --inVCF {input} --genomeBuild {params.genomeBuild} --genomeStyle {params.genomeStyle} --snpDB {params.snpDB} --minQuality {params.minQual} --minDepth {params.minDepth} --minVAF {params.minVAF} --altCountField AD --libdir {params.libdir} --outVCF {output} > {log} 2> {log}"


rule getAlleleCountsByChr:
	input:
		hetSites="results/phasedCounts/hetPosns/{tumor}.phasedHetsFromNormal.vcf",
		tumBam=lambda wildcards: getLRFullPath(config["samples"][wildcards.tumor], config["bamFileName"])
	output:
		"results/phasedCounts/tumCounts/{tumor}/{tumor}.tumCounts.{chr}.txt"
	params:
		countScript=config["phaseCounts_counts_script"],
		mapQ=config["het_minMapQuality"],
		baseQ=config["het_minBaseQuality"],
		vcfQ=config["het_minVCFQuality"],
		mem=config["std_mem"],
		runtime=config["std_runtime"],
		pe=config["std_numCores"]
	log:
		"logs/phasedCounts/tumCounts/{tumor}/{tumor}.tumCounts.{chr}.log"
	shell:
		"python {params.countScript} {wildcards.chr} {input.hetSites} {input.tumBam} {params.baseQ} {params.mapQ} {params.vcfQ} > {output} 2> {log}"

rule catAlleleCountFiles:
	input:
		expand("results/phasedCounts/tumCounts/{{tumor}}/{{tumor}}.tumCounts.{chr}.txt", chr=CHRS)
	output:
		"results/phasedCounts/tumCounts/{tumor}.tumCounts.txt"
	log:
		"logs/phasedCounts/tumCounts/{tumor}/{tumor}.cat.log"
	shell:
		"cat {input} | grep -v Chr > {output} 2> {log}"






