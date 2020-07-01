#!/usr/bin/env nextflow

BED_GLOB = params.bed_glob

// Take human peaks
// Lift human to rat

process lift {

	publishDir "${params.results}"

	input:
	file(peaks) from Channel.fromPath(BED_GLOB)

	output:
	file("snp-to-rn6.bed")

	"""
	cat $peaks | cut -f1-3 > peaks.bed
	cat peaks.bed | perl -pe 's/\\t/:/g' > peaknames.txt
	paste peaks.bed peaknames.txt > peaks.tmp
	mv peaks.tmp peaks.bed
	bnMapper.py peaks.bed -o rn6-to-peak.bed /lab/work/porchard/data/chain/hg19ToRn6.over.chain.gz
	cut -f4 rn6-to-peak.bed | perl -pe 's/:/\t/g' > hg19.bed
	cut -f1-3 rn6-to-peak.bed | perl -pe 's/\t/:/' > rn6.txt
	paste hg19.bed rn6.txt > snp-to-rn6.bed
	"""

}
