#!/usr/bin/env nextflow

SNP_BED_GLOB = params.snp_bed_glob
ROADMAP_BED_GLOB = params.roadmap_bed_glob

process sort_snps {

	input:
	file(x) from Channel.fromPath(SNP_BED_GLOB)

	output:
	file('snps.sorted.bed') into sort_snps_out

	"""
	cat $x | sort -k1,1 -k2n,2 > snps.sorted.bed
	"""

}


process add_cell_type_and_state {

	input:
	file(x) from Channel.fromPath(ROADMAP_BED_GLOB)

	output:
	set val(cell_type), file('tmp.bed') into concat_by_cell_type_in

	script:
	cell_type = x.getName().replaceAll('.bed', '').tokenize('.')[0]
	state = x.getName().replaceAll('.bed', '').tokenize('.')[1]

	"""
	cat $x | perl -pe 's/\$/\\t${cell_type}\\t${state}/' > tmp.bed
	"""

}

process concat_by_cell_type {

	input:
	set val(cell_type), file("x.*.bed") from concat_by_cell_type_in.groupTuple()

	output:
	file("${cell_type}.bed") into intersect_in

	"""
	cat x.*.bed | sort -k1,1 -k2n,2 > ${cell_type}.bed
	"""
}


process intersect {

	input:
	file(snps) from sort_snps_out
	each file(bed) from intersect_in

	output:
	file('intersected.bed') into intersect_out

	"""
	bedtools intersect -sorted -a $snps -b $bed -wa -wb | cut -f1-3,7,8 > intersected.bed
	"""

}


process concat {

	publishDir "${params.results}/intersected"

	input:
	file("intersected.*.bed") from intersect_out.toSortedList()

	output:
	file('snp-overlap-with-roadmap-states.bed') into concat_out

	"""
	cat intersected.*.bed | sort -k1,1 -k2n,2 > snp-overlap-with-roadmap-states.bed
	"""

}

process analyze {

	publishDir "${params.results}/analyze"

	input:
	file(x) from concat_out

	"""
	check-snp-chromatin-state.R $x /lab/work/porchard/data/roadmap/roadmap_cell_types.txt
	"""

}
