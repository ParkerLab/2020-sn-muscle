#!/usr/bin/env nextflow

ROOT = params.projroot
ukbb_gwas = Channel.fromPath("${ROOT}/data/ukb-summary-stats/*.both_sexes.tsv.gz").into{list_snps_1; reformat_in}

process list_ukb_snps_1 {
	
	memory '15 GB'

	input:
	file(x) from list_snps_1.collate(20)

	output:
	file('snps.txt') into list_ukb_snps_2_in

	"""
	list-all-ukb-snps.py ${x.join(' ')} > snps.txt
	"""

}

process list_ukb_snps_2 {

	publishDir "${params.results}/convert-to-rsid"
	memory '30 GB'

	input:
	file("snps.*.txt") from list_ukb_snps_2_in.toSortedList()

	output:
	file('all-ukb-snps.txt') into list_ukb_snps_2_out

	"""
	sort_uniq.py snps*.txt > all-ukb-snps.txt
	"""

}

process make_conversions {

	publishDir "${params.results}/convert-to-rsid"
	memory '20 GB'

	input:
	file(x) from list_ukb_snps_2_out

	output:
	file("snp-conversions.txt") into make_conversions_out

	"""
	make-snp-id-conversions.py $x
	"""


}

process reformat_ukbb_for_ldsc {

	maxForks 20
	memory '15 GB'
	publishDir "${params.results}/results/ldsc"

	input:
	file(conversions) from make_conversions_out
	each file(x) from reformat_in

	output:
	set val(trt), file("${trt}.txt")

	script:
	trt = x.getName().replaceAll('.gwas.imputed_v3.both_sexes.tsv.gz', '')

	"""
	convert-ukbb-snp-id.py $x $conversions ${trt}.txt
	"""

}
