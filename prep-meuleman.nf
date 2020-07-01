#!/usr/bin/env nextflow

INDEX = params.dhs_index
CHAIN = params.chain
METADATA = params.metadata
BINARY_DHS_MATRIX = params.dhs_matrix

process liftover {

	publishDir "${params.results}/common", pattern: "common_open_chromatin.bed"

	input:
	file(ix) from Channel.fromPath(INDEX)
	file(chain) from Channel.fromPath(CHAIN)

	output:
	file('hg19.bed') into lifted
	file('common_open_chromatin.bed') 

	"""
	zcat $ix | awk 'NR>1' > hg38.bed
	cut -f1-4 hg38.bed > liftover.bed
	liftOver -bedPlus=1 liftover.bed $chain hg19-lifted.bed not_lifted.txt
	cat hg19-lifted.bed | sort -k1,1 -k2n,2 | grep -v -e 'chrUn' -e 'chrX' -e 'chrY' -e 'random' > hg19.bed
	grep 'Tissue invariant' hg38.bed | awk '\$6>=500' | cut -f4 > tissue_invariant_ids.txt 
	filter.py hg19.bed tissue_invariant_ids.txt 4 | cut -f1-3 | sort -k1,1 -k2n,2 | bedtools merge -i stdin > common_open_chromatin.bed
	"""

}

process get_tissue_dhs_hg38 {

	input:
	file(mtx) from Channel.fromPath(BINARY_DHS_MATRIX)
	file(md) from Channel.fromPath(METADATA)

	output:
	file("*.hg38-peaks.txt") into tissue_dhs_hg38

	"""
	get-meuleman-tissue-peaks.R $mtx $md 	
	"""

}

process hg38_peak_to_index {

	input:
	file(ix) from Channel.fromPath(INDEX)
	each file(x) from tissue_dhs_hg38.flatten()

	output:
	set val(tissue), file('tissue_ids.txt') into tissue_indices

	script:
	tissue = x.getName().replaceAll('.hg38-peaks.txt', '')

	"""
	zcat $ix | awk 'NR>1' | cut -f1-4 | perl -pe 's/\\t/:/; s/\\t/-/' > hg38_peak_to_id.txt
	filter.py hg38_peak_to_id.txt $x 1 | cut -f2 > tissue_ids.txt
	"""

}

process index_to_hg19 {

	publishDir "${params.results}/tissues"

	input:
	set val(tissue), file(ids), file(hg19) from tissue_indices.combine(lifted)
	
	output:
	set val(tissue), file("${tissue}.bed")

	"""
	filter.py $hg19 $ids 4 | cut -f1-3 | sort -k1,1 -k2n,2 | bedtools merge -i stdin | sort -k1,1 -k2n,2 > ${tissue}.bed
	"""

}
