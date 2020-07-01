#!/usr/bin/env nextflow

ROOT = params.projroot
GENOME = params.genome
TSS = params.tss[GENOME]
CHROM_SIZES = params.chrom_sizes[GENOME]
DOWNSTREAM = params.downstream + '/results'

// Generic data
AUTOSOMAL_REFERENCES = ['hg19': (1..22).collect({it -> 'chr' + it}),
			'hg38': (1..22).collect({it -> 'chr' + it}),
			'rn5': (1..20).collect({it -> 'chr' + it}),
			'rn6': (1..20).collect({it -> 'chr' + it}),
			'mm9': (1..19).collect({it -> 'chr' + it}),
			'mm10': (1..19).collect({it -> 'chr' + it})
]

get_autosomes = {
	genome ->
	AUTOSOMAL_REFERENCES[genome]
}

// Name master peaks
// Filter master peaks to TSS-distal
// Get accessibility of each (peak, cluster)
// Lift to human coordinates, if necessary

process name_and_filter_peaks {

	executor 'local'
	container "${params.containers.general}"

	input:
	file(master_peaks) from Channel.fromPath("${DOWNSTREAM}/process-by-cluster-round-1/master-peaks/master-peaks.${GENOME}.bed")

	output:
	file('master-peaks.bed') into accessibility_peaks_in
	file('master-peaks.bed') into master_peaks_in

	"""
	bedtools slop -b 5000 -i $TSS -g $CHROM_SIZES | sort -k1,1 -k2n,2 > tss-proximal.bed
	cat $master_peaks | awk '{print(\$1, \$2, \$3, "peak_" NR)}' | perl -pe 's/ /\\t/g' | sort -T . -k1,1 -k2n,2 | bedtools intersect -a stdin -b tss-proximal.bed -v > master-peaks.bed
	"""	

}

liftover_in = Channel.create()
not_lifted = Channel.create()
master_peaks_in.choice(liftover_in, not_lifted) { a -> GENOME == 'rn6' ? 0 : 1}

process lift_master_peaks {

	executor 'local'
	container "${params.containers.bnmapper}"

	input:
	file(master_peaks) from liftover_in

	output:
	file('master-peaks.lifted-sorted.bed') into lifted

	"""
	bedtools slop -b 5000 -i ${params.tss['hg19']} -g ${params.chrom_sizes['hg19']} | sort -k1,1 -k2n,2 > tss-proximal.hg19.bed
	bnMapper.py -o master-peaks.lifted.bed $master_peaks ${params.chain}
	cat master-peaks.lifted.bed | grep -w -e ${get_autosomes('hg19').join(' -e ')}| sort -k1,1 -k2n,2 -T . | bedtools intersect -a stdin -b tss-proximal.hg19.bed -v > master-peaks.lifted-sorted.bed
	"""

}


process enhancer_regression_accessibility {
	
	container "${params.containers.general}"

	input:
	file(peaks) from Channel.fromPath("${DOWNSTREAM}/process-by-cluster-round-1/peaks/*-${GENOME}_peaks.broadPeak.noblacklist").toSortedList()
	file(master_peaks) from accessibility_peaks_in

	output:
	file("accessibility.bed") into enhancer_regression_accessibility_out

	"""
	cat ${peaks.join(' ')} | sort -T . -k1,1 -k2n,2 | perl -pe 's/\\t(\\d+)-${GENOME}_peak.*/\\tcluster_\$1/' > cluster-peaks.txt
	bedtools intersect -a $master_peaks -b cluster-peaks.txt -sorted -wa -wb | cut -f4,8 > accessibility.bed # peakname, cluster it's accessible in (one per line)
	"""

}

// Run the enhancer regression
process get_enhancer_posteriors {

	memory '10 GB'
	maxForks 50
	tag "${cell_type} ${chrom}"
	container "${params.containers.general}"

	input:
	file(master_peaks) from lifted.mix(not_lifted)
	each file(posterior_file) from Channel.fromPath("${ROOT}/data/roadmap-posteriors/*posterior.txt.gz")

	output:
	set val(cell_type), file("scores.bed") into enhancer_posteriors_out

	script:
	m = posterior_file.getName() =~ /(.*)_15_coreMarks_(chr.*)_posterior.txt.gz/
	cell_type = m[0][1]
	chrom = m[0][2]

	"""
	fetch-enhancer-posteriors-multiple-states.py $master_peaks $posterior_file E6 E7 E12 > scores.bed
	"""

}

process concat_enhancer_posteriors {

	tag "${cell_type}"
	container "${params.containers.general}"

	input:
	set val(cell_type), file("scores*.bed") from enhancer_posteriors_out.groupTuple()

	output:
	file("${cell_type}.bed") into enhancer_posteriors

	"""
	cat scores*.bed | grep -v '^name' | sort -T . > ${cell_type}.bed
	"""

}


process enhancer_regression {

	publishDir "${params.results}/enhancer-regression"
	memory '20 GB'

	input:
	file(accessibility) from enhancer_regression_accessibility_out
	file(posteriors) from enhancer_posteriors.toSortedList()

	output:
	file("model_results.txt") into enhancer_regression_out

	"""
	logistic-regression-on-enhancer-posteriors.py $accessibility ${posteriors.join(' ')} | grep cluster > model_results.txt
	"""

}
