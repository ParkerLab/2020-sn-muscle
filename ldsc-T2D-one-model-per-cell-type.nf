#!/usr/bin/env nextflow

ROOT = params.projroot
CLUSTER_NAMES = params.cluster_names
IONICE = 'ionice -c2 -n7'

CHROMS = (1..22)
peaks = Channel.fromPath(params.peak_glob)
ldsc_annot_file_individual = Channel.from(CHROMS).combine(peaks)
OTHER_ANNOTATIONS = params.other_annotations

process make_ldsc_annot_file_individual {

	container "${params.containers.ldsc}"
	maxForks 10

	input:
	set val(chrom), file(bed) from ldsc_annot_file_individual

	output:
	set val(chrom), val(cluster), file("${cluster}.${chrom}.annot.gz") into ldsc_annot_individual

	script:
	cluster = bed.getName().replaceAll('.bed', '')

	"""
	source activate ldsc && make_annot_custom.py --bimfile ${ROOT}/data/ldsc-data/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom}.bim --bed-files $bed $OTHER_ANNOTATIONS --annot-file ${cluster}.${chrom}.annot.gz
	"""

}

process make_ldsc_ld_scores_individual {

	container "${params.containers.ldsc}"

	input:
	set val(chrom), val(cluster), file(annot) from ldsc_annot_individual

	output:
	set val("${cluster}"), file("${cluster}.${chrom}.l2.M"), file("${cluster}.${chrom}.l2.M_5_50"), file("${cluster}.${chrom}.l2.ldscore.gz"), file(annot) into ldsc_ld_scores_individual

	"""
	cat ${ROOT}/data/ldsc-data/w_hm3.snplist | grep -v "^SNP" | cut -f1 | sort | uniq > w_hm3.onlySNPid.snplist
	source activate ldsc && export MKL_NUM_THREADS=1 && export NUMEXPR_NUM_THREADS=1 && export OPENBLAS_NUM_THREADS=1 && export OMP_NUM_THREADS=1 && python /sw/ldsc/ldsc.py --l2 --maf 0.01 --thin-annot --bfile ${ROOT}/data/ldsc-data/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom} --ld-wind-cm 1 --print-snps w_hm3.onlySNPid.snplist --annot $annot --out ${cluster}.${chrom}
	"""

}


process get_diamante_sumstats {

        container "${params.containers.general}"
        memory '10 GB'

        input:
        file(x) from Channel.fromPath(ROOT + '/data/diamante-summary-stats/diamante*.txt')

        output:
        set val(trt), file("${trt}.stats.txt") into get_diamante

        script:
        trt = x.getName().replaceAll('.txt', '')
        bmi_adjusted_or_unadjusted = trt.contains('bmiadj') ? 'adjusted' : 'unadjusted'

        """
        adjust-diamante-columns-use-neff.py $x > ${trt}.stats.txt
        """

}

process get_manning_sumstats {

	container "${params.containers.general}"
	executor 'local'

	input:
	file(x) from Channel.fromPath(ROOT + '/data/manning-summary-stats/manning*.txt')

	output:
	set val(trt), file("${trt}.stats.txt") into get_manning

	script:
	trt = x.getName().replaceAll('.txt', '')

	"""
	cat $x > ${trt}.stats.txt
	"""

}

process munge_sumstats {

        container "${params.containers.ldsc}"
        publishDir "${params.results}/sumstats"
        memory '6 GB'

        input:
        set val(trt), file(x) from get_diamante.mix(get_manning)

        output:
        set val(trt), file("${trt}.sumstats.gz") into sumstats_out

        """
        source activate ldsc && python /sw/ldsc/munge_sumstats.py --sumstats ${trt}.stats.txt --merge-alleles ${ROOT}/data/ldsc-data/w_hm3.snplist --out $trt
        """

}

process run_ldsc {

	container "${params.containers.ldsc}"
	publishDir "${params.results}/partitioned-heritability"

	input:
	set val(cluster), file(M1), file(M2), file(ldscore), file(annot), val(trt), file(sumstats) from ldsc_ld_scores_individual.groupTuple(size: 22).combine(sumstats_out)

	output:
	file("*.log")
	set val(cluster), file("${cluster}.${trt}.results") into ldsc_magic_results

	"""
	source activate ldsc && export MKL_NUM_THREADS=1 && export NUMEXPR_NUM_THREADS=1 && export OPENBLAS_NUM_THREADS=1 && export OMP_NUM_THREADS=1 && python /sw/ldsc/ldsc.py --h2 $sumstats --maf 0.05 --ref-ld-chr ${ROOT}/work/ldsc-baseline-min-maf/results/baseline-model/baselineLD.,${cluster}. --w-ld-chr ${ROOT}/data/ldsc-data/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. --overlap-annot --frqfile-chr ${ROOT}/data/ldsc-data/1000G_Phase3_frq/1000G.EUR.QC. --out ${cluster}.${trt} --print-coefficients
	"""

}

process plot {

	publishDir "${params.results}/plot"
	container "${params.containers.rplot}"

	input:
	file(x) from ldsc_magic_results.map({x -> x[1]}).toSortedList()

	output:
	file("*.pdf")

	"""
	plot-ldsc-T2D.R ${CLUSTER_NAMES} ${x.join(' ')}
	"""

}
