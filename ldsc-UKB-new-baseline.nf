#!/usr/bin/env nextflow

IONICE = 'ionice -c2 -n7'

CHROMS = (1..22)
ROOT = params.projroot
peaks = Channel.fromPath(params.peak_glob)
ldsc_annot_file_individual = Channel.from(CHROMS).combine(peaks)
OTHER_ANNOTATIONS = params.other_annotations
TRTS = Channel.fromPath(params.trait_list).splitText().map({it -> it.replaceAll('\n', '')})
CLUSTER_NAMES = Channel.fromPath(params.cluster_names)

process make_ldsc_annot_file {

	container "${params.containers.ldsc}"
	memory '30 GB'

	input:
	set val(chrom), file(bed) from ldsc_annot_file_individual.groupTuple()

	output:
	set val(chrom), file("tissues.${chrom}.annot.gz") into ldsc_annot_individual

	"""
	source activate ldsc && make_annot_custom.py --bimfile ${ROOT}/data/ldsc-data/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom}.bim --bed-files ${bed.join(' ')} $OTHER_ANNOTATIONS --annot-file tissues.${chrom}.annot.gz
	"""

}

process make_ldsc_ld_scores_individual {

	container "${params.containers.ldsc}"

	input:
	set val(chrom), file(annot) from ldsc_annot_individual

	output:
	set val("${cluster}"), file("${cluster}.${chrom}.l2.M"), file("${cluster}.${chrom}.l2.M_5_50"), file("${cluster}.${chrom}.l2.ldscore.gz"), file(annot) into ldsc_ld_scores_individual

	script:
	cluster = 'tissues'

	"""
	cat ${ROOT}/data/ldsc-data/w_hm3.snplist | grep -v "^SNP" | cut -f1 | sort | uniq > w_hm3.onlySNPid.snplist
	source activate ldsc && export MKL_NUM_THREADS=1 && export NUMEXPR_NUM_THREADS=1 && export OPENBLAS_NUM_THREADS=1 && export OMP_NUM_THREADS=1 && python /sw/ldsc/ldsc.py --l2 --maf 0.01 --thin-annot --bfile ${ROOT}/data/ldsc-data/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom} --ld-wind-cm 1 --print-snps w_hm3.onlySNPid.snplist --annot $annot --out ${cluster}.${chrom}
	"""

}

process ldsc_sumstats {

	container "${params.containers.ldsc}"
	publishDir "${params.results}/sumstats"
	memory '6 GB'
	maxForks 40

	input:
	val(trt) from TRTS

	output:
	set val(trt), file("${trt}.sumstats.gz") into sumstats_out 

	"""
	source activate ldsc && filter-sumstat-snps.py ${ROOT}/work/reformat-ukbb/results/ldsc/${trt}.txt ${trt}.sumstats-reformatted.txt ${ROOT}/data/ldsc-data/w_hm3.snplist
        source activate ldsc && python /sw/ldsc/munge_sumstats.py --sumstats ${trt}.sumstats-reformatted.txt --merge-alleles ${ROOT}/data/ldsc-data/w_hm3.snplist --out $trt
	"""

}


process run_ldsc {

	container "${params.containers.ldsc}"
	publishDir "${params.results}/partitioned-heritability"

	input:
	set val(cluster), file(M1), file(M2), file(ldscore), file(annot), val(trt), file(sumstats) from ldsc_ld_scores_individual.groupTuple(size: 22).combine(sumstats_out)

	output:
	file("*.log")
	set val(cluster), file("${cluster}.${trt}.results") into ldsc_results

	"""
	source activate ldsc && export MKL_NUM_THREADS=1 && export NUMEXPR_NUM_THREADS=1 && export OPENBLAS_NUM_THREADS=1 && export OMP_NUM_THREADS=1 && python /sw/ldsc/ldsc.py --h2 $sumstats --maf 0.05 --ref-ld-chr ${ROOT}/work/ldsc-baseline-min-maf/results/baseline-model/baselineLD.,${cluster}. --w-ld-chr ${ROOT}/data/ldsc-data/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. --overlap-annot --frqfile-chr ${ROOT}/data/ldsc-data/1000G_Phase3_frq/1000G.EUR.QC. --out ${cluster}.${trt} --print-coefficients
	"""

}


process plot_per_trait {

	publishDir "${params.results}/plot-per-trait"

	input:
	file(cluster_names) from CLUSTER_NAMES
	each file(x) from ldsc_results.map({x -> x[1]})

	output:
	file("*.pdf")

	"""
	# plot-ldsc-ukb.R ${x.join(' ')}
	plot-trait-ldsc-result.R $x $cluster_names
	"""

}
