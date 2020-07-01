#!/usr/bin/env nextflow

ROOT = params.projroot

LIBRARIES = ['133155', '133156', '133157', '133158', '63_20_rna', '63_40_rna'].collect({it -> it + '-hg19'})

library_to_matrix = {
	library ->
	"${ROOT}/work/rnaseq/results/starsolo/" + library + "/Solo.out/GeneFull/raw"
}

uncorrected_in = Channel.from(LIBRARIES).map({it -> [it, file(library_to_matrix(it))]})

process uncorrected_counts {

	publishDir "${params.results}/uncorrected-counts"
	executor 'local'
	maxForks 1

	input:
	set val(lib), file(x) from uncorrected_in
	
	output:
	set val("uncorrected"), val(lib), file("${lib}.features.txt") into uncorrected_counts_in

	"""
	cat ${x}/barcodes.tsv > barcodes.txt
	cat ${x}/features.tsv | cut -f2 > features.txt
	cat ${x}/matrix.mtx > matrix.mtx
	starsolo-mtx-to-feature-file.py features.txt barcodes.txt matrix.mtx | perl -pe 's/^/${lib}\\t/' > ${lib}.features.txt
	"""

}

