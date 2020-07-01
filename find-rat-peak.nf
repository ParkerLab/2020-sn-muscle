#!/usr/bin/env nextflow

PEAK_GLOB = params.peak_glob

// Take human peaks
// Lift human to rat

process lift {

	publishDir "${params.results}"

	input:
	file(peaks) from Channel.fromPath(PEAK_GLOB)

	output:
	file("${cluster}.peak-to-rn6.bed")

	script:
	cluster = 'cluster_' + peaks.getName().tokenize('-')[0]

	"""
	cat $peaks | cut -f1-3 > peaks.bed
	cat peaks.bed | perl -pe 's/\\t/:/g' > peaknames.txt
	paste peaks.bed peaknames.txt > peaks.tmp
	mv peaks.tmp peaks.bed
	bnMapper.py peaks.bed -o ${cluster}.rn6-to-peak.bed /lab/work/porchard/data/chain/hg19ToRn6.over.chain.gz
	cut -f4 ${cluster}.rn6-to-peak.bed | perl -pe 's/:/\t/g' > hg19.bed
	cut -f1-3 ${cluster}.rn6-to-peak.bed | perl -pe 's/\t/:/' > rn6.txt
	paste hg19.bed rn6.txt > ${cluster}.peak-to-rn6.bed
	"""

}
