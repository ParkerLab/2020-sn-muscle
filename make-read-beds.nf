#!/usr/bin/env nextflow

ATAC_BAM_GLOB = params.atac_bam_glob
RNA_BAM_GLOB = params.rna_bam_glob

process make_bed_atac {

	publishDir "${params.results}/beds/atac"

	input:
	file(bam) from Channel.fromPath(ATAC_BAM_GLOB)

	output:
	file("${library}.pruned.bed.gz")

	script:
	library = bam.getName().replaceAll('.pruned.bam', '')

	"""
	bedtools bamtobed -i $bam | awk '\$5="."' | perl -pe 's/ /\\t/g' | gzip -c > ${library}.pruned.bed.gz
	"""

}

process add_barcodes_to_readname {

	input:
	file(bam) from Channel.fromPath(RNA_BAM_GLOB)

	output:
	file("${library}.with-barcodes.bam") into rna_in

	script:
	library = bam.getName().replaceAll('.before-dedup.bam', '')

	"""
	add-barcode-and-umi-to-readname.py $bam ${library}.with-barcodes.bam
	"""

}

process make_bed_rna {

	publishDir "${params.results}/beds/rna"

	input:
	file(bam) from rna_in

	output:
	file("${library}.pruned-no-dedup.bed.gz")

	script:
	library = bam.getName().replaceAll('.with-barcodes.bam', '')

	"""
	bedtools bamtobed -i $bam | awk '\$5="."' | perl -pe 's/ /\\t/g' | gzip -c > ${library}.pruned-no-dedup.bed.gz
	"""

}
