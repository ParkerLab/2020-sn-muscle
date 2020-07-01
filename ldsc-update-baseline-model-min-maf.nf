#!/usr/bin/env nextflow

ROOT = params.projroot
IONICE = 'ionice -c2 -n7'

process strip_annotations {

	container "${params.containers.general}"

	input:
	file(annot) from Channel.fromPath("${ROOT}/data/ldsc-data/baseline/*.annot.gz")

	output:
	set val(chrom), file("baselineLD.${chrom}.annot") into no_cell_type_annots

	script:
	chrom = annot.getName().replaceAll('baselineLD.', '').replaceAll('.annot.gz', '')

	"""
	remove-ldsc-baseline-cell-type-specific-annotations.py $annot > baselineLD.${chrom}.annot
	"""

}

process make_ldsc_ld_scores {

	container "${params.containers.ldsc}"
	publishDir "${params.results}/baseline-model"

	input:
	set val(chrom), file(annot) from no_cell_type_annots

	output:
	set file("baselineLD.${chrom}.l2.M"), file("baselineLD.${chrom}.l2.M_5_50"), file("baselineLD.${chrom}.l2.ldscore.gz"), file("baselineLD.${chrom}.annot.gz") 

	"""
	cat $annot | gzip -c > ${annot}.gz
	cut -f1 ${ROOT}/data/ldsc-data/w_hm3.snplist | grep -v -w 'SNP' | sort | uniq > print-these-snps.txt
	source activate ldsc && python ~/sw/ldsc/ldsc.py --l2 --maf 0.01 --bfile ${ROOT}/data/ldsc-data/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom} --ld-wind-cm 1 --print-snps print-these-snps.txt --annot ${annot}.gz --out baselineLD.${chrom}
	"""

}
