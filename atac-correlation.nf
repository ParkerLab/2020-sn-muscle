#!/usr/bin/env nextflow

BAMS = Channel.fromPath(params.bam_glob)
PEAKS = Channel.fromPath(params.peak_glob)

process make_master_peaks {
	
	executor 'local'

	input:
	file(beds) from PEAKS.toSortedList()

	output:
	file('master-peaks.bed') into master_peaks

	"""
	cat ${beds.join(' ')} | sort -k1,1 -k2n,2 -T . | bedtools merge -i stdin > master-peaks.bed
	"""
}


process get_counts {

	maxForks 5
	publishDir "${params.results}/counts"

	input:
	file(peaks) from master_peaks
	each file(bam) from BAMS

	output:
	file("${library}.counts.bed")

	script:
	library = bam.getName().replaceAll('.bam', '')

	"""
	sort-bed-by-bam.py $peaks $bam > peaks.txt
	coverageBed -counts -sorted -a peaks.txt -b $bam | sort -k1,1 -k2n,2 > ${library}.counts.bed
	"""

}
