#!/usr/bin/nextflow

ROOT = params.projroot
MAIN = ROOT + '/work/downstream-new-features/results'
CLUSTER_NAMES = ROOT + '/2019-01-04-cluster-names.txt'
FANS_AND_LOADING_POSITION = 'chr8:41,427,421-41,535,659'

get_tss = {
	genome ->
	params.tss[genome]
}

get_chrom_sizes = {
	genome ->
	params.chrom_sizes[genome]
}



process plot_qc_thresholds_atac {

	publishDir "${params.results}/fans-vs-no-fans", pattern: "*-fans-vs-no-fans.png"
	publishDir "${params.results}/20k-vs-40k", pattern: "*-20k-vs-40k.png"
	publishDir "${params.results}/feb", pattern: "*-feb.png"
	publishDir "${params.results}/qc-for-downstream-libraries", pattern: "*-used-downstream.png"
	publishDir "${params.results}/tables", pattern: "*.csv"
	executor 'local'
	container "${params.containers.rplot}"

	input:
	file(metrics) from Channel.fromPath(MAIN + '/nucleus-qc/metrics.txt')
	file(initial_thresholds) from Channel.fromPath("${ROOT}/initial-thresholds-atac.txt")
	file(library_labels) from Channel.fromPath("${ROOT}/library-labels.txt")
	file(upper_thresholds) from Channel.fromPath("${MAIN}/upper-atac-thresholds/atac-max-hqaa-thresholds.txt")

	output:
	file("*.png")
	file("atac-qc-thresholds.csv")

	"""
	call-nuclei-atac.R $metrics $initial_thresholds $upper_thresholds $library_labels
	"""

}

process plot_qc_thresholds_rna {

	publishDir "${params.results}/fans-vs-no-fans", pattern: "*-fans-vs-no-fans.png"
	publishDir "${params.results}/20k-vs-40k", pattern: "*-20k-vs-40k.png"
	publishDir "${params.results}/qc-for-downstream-libraries", pattern: "*-used-downstream.png"
	publishDir "${params.results}/tables", pattern: "*.csv"
	executor 'local'
	container "${params.containers.rplot}"
	
	input:
	file(qc) from Channel.fromPath("${MAIN}/nucleus-qc/*-hg19.txt").toSortedList()
	file(initial_thresholds) from Channel.fromPath("${ROOT}/initial-thresholds-rna.txt")
	file(library_labels) from Channel.fromPath("${ROOT}/library-labels.txt")

	output:
	file("*.png")
	file("rna-qc-thresholds.csv")

	"""
	call-nuclei-rna-with-mitochondrial.R $initial_thresholds $library_labels ${qc.join(' ')}
	"""

}

process read_and_nuclei_counts_and_summary_stats {

	publishDir "${params.results}/fans-vs-no-fans", pattern: "fans*"
	publishDir "${params.results}/20k-vs-40k", pattern: "loading*"
	publishDir "${params.results}/tables", pattern: "summary_stats.tsv"
	container "${params.containers.rplot}"
	executor 'local'
	
	input:
	file(library_labels) from Channel.fromPath("${ROOT}/library-labels.txt")

	output:
	file("*.png")

	"""
	plot-nuclei-counts-and-coverage.R --nuclei ${MAIN}/nucleus-qc/nuclei.txt --atac-metrics ${MAIN}/nucleus-qc/metrics.txt --rna-metrics ${MAIN}/nucleus-qc/hqaa.txt --library-labels $library_labels
	"""
}

process ataqv_metrics {

	publishDir "${params.results}/fans-vs-no-fans", pattern: "fans*"
	publishDir "${params.results}/20k-vs-40k", pattern: "loading*"
	executor 'local'

	input:
	file(ataqv) from Channel.fromPath("${ROOT}/work/process-as-bulk/results/ataqv/*.json.gz").mix(Channel.fromPath("${ROOT}/work/bulk-atacseq/results/ataqv/320-NM-*.ataqv.json.gz")).toSortedList()

	output:
	file("*.pdf")

	"""
	cp ${ROOT}/library-labels.txt .
	extractAtaqvMetric.py --files ${ataqv.join(' ')} --metrics total_peaks hqaa_overlapping_peaks_percent | perl -pe 's/.ataqv.json.gz//' > metrics.txt
	extractFragmentLengths.py --files ${ataqv.join(' ')} | perl -pe 's/.ataqv.json.gz//' > fragment-lengths.txt
	extractTssCoverage.py --files ${ataqv.join(' ')} | perl -pe 's/.ataqv.json.gz//' > tss-coverage.txt
	plot-ataqv-metrics.R
	"""

}

process chromatin_state_overlap {
	
	container "${params.containers.general}"
	executor 'local'

	input:
	file(peaks) from Channel.fromPath("${ROOT}/work/process-as-bulk/results/macs2/*.noblacklist")

	output:
	file("${library}.chromhmm_overlap.txt") into plot_aggregate_chromatin_state_overlap_in

	script:
	genome = 'hg19'
	library = peaks.getName().replaceAll('_peaks.broadPeak.noblacklist', '')

	"""
	chromhmm_overlap.py --tss-extension 5000 $peaks ${get_tss(genome)} ${get_chrom_sizes(genome)} ${params.chromatin_state_glob} > ${library}.chromhmm_overlap.txt
	"""

}


process plot_aggregate_chromatin_state_overlap {
	
	publishDir "${params.results}/fans-vs-no-fans", pattern: "fans*"
	publishDir "${params.results}/20k-vs-40k", pattern: "loading*"
	container "${params.containers.rplot}"
	executor 'local'

	input:
	file(coverage) from plot_aggregate_chromatin_state_overlap_in.toSortedList()

	output:
	file("*.pdf")

	"""
	plot_chromatin_state_overlap.R ${ROOT}/library-labels.txt ${coverage.join(' ')}
	"""

}


process fans_screenshot {

	publishDir "${params.results}/fans-vs-no-fans"
	executor 'local'
	cache false

	output:
	file("fans-gb.pdf")

	"""
	echo ${FANS_AND_LOADING_POSITION} | perl -pe 's/,//g; s/[:-]/\t/g' > locus.txt
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session 2020-sn-muscle-fans-no-nucleus-qc --bed locus.txt --prefix fans --viewer-width 700 --label-width 15 --extension 0 --text-size 24
	mv fans* fans-gb.pdf
	"""

}

process loading_screenshot {

	publishDir "${params.results}/20k-vs-40k"
	executor 'local'
	cache false

	output:
	file("loading-gb.pdf")

	"""
	echo ${FANS_AND_LOADING_POSITION} | perl -pe 's/,//g; s/[:-]/\t/g' > locus.txt
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session 2020-sn-muscle-loading-no-nucleus-qc --bed locus.txt --prefix loading --viewer-width 700 --label-width 22 --extension 0 --text-size 24
	mv loading* loading-gb.pdf
	"""

}

process per_nucleus_myh {
	
	executor 'local'
	publishDir "${params.results}/per-nucleus-plots"

	input:
	file(x) from Channel.fromPath("${MAIN}/features/same-space/*counts.txt").toSortedList()
	file(clusters) from Channel.fromPath("${MAIN}/liger/round-1/second-louvain-clusters.txt")
	file(library_labels) from Channel.fromPath("${ROOT}/library-labels.txt")

	output:
	file("*.png")

	"""
	cat ${x.join(' ')} | grep MYH > per-nucleus-gene-counts.txt
	plot-MYH-counts.R $clusters per-nucleus-gene-counts.txt $CLUSTER_NAMES $library_labels
	"""

}


process per_library_gene_counts {
	
	executor 'local'

	input:
	file(x) from Channel.fromPath("${MAIN}/features/raw/*counts.txt").toSortedList()

	output:
	file("per-library-gene-counts.txt") into per_library_gene_counts_out

	"""
	aggregate-counts-per-library.py ${x.join(' ')} > per-library-gene-counts.txt
	"""

}

process rnaseq_aggregate_correlation {
	
	executor 'local'
	
	publishDir "${params.results}/fans-vs-no-fans", pattern: "fans*"
	publishDir "${params.results}/20k-vs-40k", pattern: "loading*"

	input:
	file(x) from per_library_gene_counts_out

	output:
	file("*.png")

	"""
	cp ${ROOT}/library-labels.txt .
	rnaseq-aggregate-correlation.R $x
	"""

}

process atac_bulk_correlation_fans {

	publishDir "${params.results}/fans"
	container "${params.containers.rplot}"
	executor 'local'

	input:
	file(x) from Channel.fromPath("${ROOT}/work/atac-correlation-fans/results/counts/*.bed").toSortedList()

	output:
	file("fans-atac-correlation.png")

	"""
	atacseq-aggregate-correlation.R ${ROOT}/library-labels.txt ${x.join(' ')}
	ln -s all-atac-correlation.png fans-atac-correlation.png
	"""

}

process atac_bulk_correlation_loading {

	publishDir "${params.results}/loading"
	container "${params.containers.rplot}"
	executor 'local'

	input:
	file(x) from Channel.fromPath("${ROOT}/work/atac-correlation-loading/results/counts/*.bed").toSortedList()

	output:
	file("loading-atac-correlation.png")

	"""
	atacseq-aggregate-correlation.R ${ROOT}/library-labels.txt ${x.join(' ')}
	ln -s all-atac-correlation.png loading-atac-correlation.png
	"""

}


process nuclei_stats_table {

	publishDir "${params.results}/tables"
	executor 'local'

	input:
	file(library_labels) from Channel.fromPath("${ROOT}/library-labels.txt")
	file(nuclei_with_individual_assignments) from Channel.fromPath("${MAIN}/nucleus-qc/nuclei-with-individual-assignments.txt")
	file(hqaa) from Channel.fromPath("${MAIN}/nucleus-qc/hqaa.txt")
	file(ataqv_metrics) from Channel.fromPath("${MAIN}/nucleus-qc/metrics.txt")

	output:
	file('nuclei-summary-stats.csv')
	file('per-species-and-modality-fragment-stats.csv')

	"""
	make-table-of-nuclei-stats.R $library_labels $nuclei_with_individual_assignments $hqaa $ataqv_metrics
	make-table-of-per-species-per-modality-fragment-stats.R $library_labels $nuclei_with_individual_assignments $hqaa $ataqv_metrics
	"""

}

process cell_type_proportions {

	publishDir "${params.results}/cell-type-proportions"
	executor 'local'

	input:
	file(clusters) from Channel.fromPath("${MAIN}/liger/round-1/second-louvain-clusters.txt")
	file(cluster_names) from Channel.fromPath(CLUSTER_NAMES)

	output:
	file("*.pdf")

	"""
	plot-cell-type-proportions.R $cluster_names $clusters
	"""

}

process cell_type_counts_table {

	publishDir "${params.results}/tables"
	executor 'local'

	input:
	file(clusters) from Channel.fromPath("${MAIN}/liger/round-1/second-louvain-clusters.txt")
	file(cluster_names) from Channel.fromPath(CLUSTER_NAMES)
	file(library_labels) from Channel.fromPath("${ROOT}/library-labels.txt")

	output:
	file("*.tsv")

	"""
	make-cell-type-count-table.R $cluster_names $clusters $library_labels
	"""

}

process cell_type_counts_plots {

	publishDir "${params.results}/cell-type-counts"
	executor 'local'

	input:
	file(clusters) from Channel.fromPath("${MAIN}/liger/round-1/second-louvain-clusters.txt")
	file(cluster_names) from Channel.fromPath(CLUSTER_NAMES)
	file(library_labels) from Channel.fromPath("${ROOT}/library-labels.txt")

	output:
	file("*.pdf")

	"""
	plot-cell-type-numbers.R $cluster_names $clusters $library_labels
	"""

}


process cell_type_proportions_by_sample {

	publishDir "${params.results}/cell-type-proportions"
	executor 'local'

	input:
	file(clusters) from Channel.fromPath("${MAIN}/liger/round-1/second-louvain-clusters.txt")
	file(assignments) from Channel.fromPath("${MAIN}/nucleus-qc/nuclei-with-individual-assignments.txt")
	file(library_labels) from Channel.fromPath("${ROOT}/library-labels.txt")

	output:
	file("*.pdf")

	"""
	plot-fraction-of-library-by-metadata.R $clusters ${CLUSTER_NAMES} $assignments $library_labels
	"""

}

process fiber_type_rna_differences {

	publishDir "${params.results}/fiber-type-rna-diffs"
	container "${params.containers.rplot}"
	executor 'local'
	
	input:
	file(clusters) from Channel.fromPath("${MAIN}/liger/round-1/second-louvain-clusters.txt")
	file(counts_20k) from Channel.fromPath("${MAIN}/features/same-space/63_20_rna-hg19.counts.txt")
	file(counts_40k) from Channel.fromPath("${MAIN}/features/same-space/63_40_rna-hg19.counts.txt")
	file(counts_133155) from Channel.fromPath("${MAIN}/features/same-space/133155-hg19.counts.txt")
	file(counts_133156) from Channel.fromPath("${MAIN}/features/same-space/133156-hg19.counts.txt")
	file(counts_133157) from Channel.fromPath("${MAIN}/features/same-space/133157-hg19.counts.txt")
	file(counts_133158) from Channel.fromPath("${MAIN}/features/same-space/133158-hg19.counts.txt")
	
	output:
	file("*.pdf")

	"""
	compare-rubenstein.R $clusters ${ROOT}/data/rubenstein-fiber-type-differential-genes/fiber-type-differential-genes.txt rubenstein-vs-our-fiber-type-lfcs.pdf $counts_20k $counts_40k $counts_133155 $counts_133156 $counts_133157 $counts_133158
	"""

}


process logistic_regression_heatmap_human {

	publishDir "${params.results}/logistic-regression/hg19"
	executor 'local'

	input:
	file(x) from Channel.fromPath(ROOT + '/work/logistic-regression-multiple-states-hg19/results/enhancer-regression/model_results.txt')
	
	output:
	file("enhancer-similarity-heatmap.pdf")

	"""
	plot-logistic-regression-results.R $x $CLUSTER_NAMES ${ROOT}/data/roadmap-posteriors/roadmap_cell_types.txt
	"""

}

process logistic_regression_heatmap_rat {

	publishDir "${params.results}/logistic-regression/rn6"
	cache 'deep'
	executor 'local'

	input:
	file(x) from Channel.fromPath(ROOT + '/work/logistic-regression-multiple-states-rn6/results/enhancer-regression/model_results.txt')
	file(cluster_names) from Channel.fromPath(CLUSTER_NAMES)
	
	output:
	file("enhancer-similarity-heatmap.pdf")

	"""
	plot-logistic-regression-results.R $x $cluster_names ${ROOT}/data/roadmap-posteriors/roadmap_cell_types.txt
	"""

}

process logistic_regression_heatmap {

	publishDir "${params.results}/logistic-regression/both"
	cache 'deep'
	executor 'local'

	input:
	file('human.txt') from Channel.fromPath(ROOT + '/work/logistic-regression-multiple-states-hg19/results/enhancer-regression/model_results.txt')
	file('rat.txt') from Channel.fromPath(ROOT + '/work/logistic-regression-multiple-states-rn6/results/enhancer-regression/model_results.txt')
	file(cluster_names) from Channel.fromPath(CLUSTER_NAMES)
	
	output:
	file("enhancer-similarity-heatmap.pdf")

	"""
	plot-logistic-regression-results-both-species.R human.txt rat.txt $cluster_names ${ROOT}/data/roadmap-posteriors/roadmap_cell_types.txt
	"""

}

process plot_umap_by_modality_and_species {

	publishDir "${params.results}/umaps"
	container "${params.containers.rplot}"
	cache 'deep'
	executor 'local'

	input:
	file(umap) from Channel.fromPath("${MAIN}/liger/round-1/umap.txt")
	file(clusters) from Channel.fromPath("${MAIN}/liger/round-1/second-louvain-clusters.txt")

	output:
	file('full-umap.png')
	file('split-umap.png')

	"""
	umap-modality-and-species.R $umap $clusters
	"""

}

process marker_gene_screenshots_hg19 {

	publishDir "${params.results}/gb-screenshots/marker-genes"
	executor 'local'
	cache 'false'
	
	output:
	file("*.pdf")

	script:
	session = '2020-sn-muscle-cell-types-new-features-hg19'

	"""
	echo 'chr12:6,642,596-6,648,525' | perl -pe 's/,//g; s/[:-]/\t/g' > gapdh.bed
	echo 'chr17:10,414,886-10,423,629' | perl -pe 's/,//g; s/[:-]/\t/g' > myh1.bed
	echo 'chr14:23,878,228-23,912,613' | perl -pe 's/,//g; s/[:-]/\t/g' > myh7.bed
	echo 'chr4:55,091,467-55,101,711' | perl -pe 's/,//g; s/[:-]/\t/g' > pdgfra.bed
	echo 'chr12:6,221,410-6,238,772' | perl -pe 's/,//g; s/[:-]/\t/g' > vwf.bed
	echo 'chr16:15,948,630-15,953,189' | perl -pe 's/,//g; s/[:-]/\t/g' > myh11.bed
	echo 'chr12:7,650,408-7,657,741' | perl -pe 's/,//g; s/[:-]/\t/g' > cd163.bed
	echo 'chr1:18,956,440-18,961,112' | perl -pe 's/,//g; s/[:-]/\t/g' > pax7.bed
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session ${session} --bed gapdh.bed --prefix hg19-gapdh --viewer-width 500 --label-width 28 --extension 0 --text-size 18
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session ${session} --bed myh1.bed --prefix hg19-myh1 --viewer-width 500 --label-width 28 --extension 0 --text-size 18
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session ${session} --bed myh7.bed --prefix hg19-myh7 --viewer-width 500 --label-width 28 --extension 0 --text-size 18
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session ${session} --bed pdgfra.bed --prefix hg19-pdgfra --viewer-width 500 --label-width 28 --extension 0 --text-size 18
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session ${session} --bed vwf.bed --prefix hg19-vwf --viewer-width 500 --label-width 28 --extension 0 --text-size 18
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session ${session} --bed myh11.bed --prefix hg19-myh11 --viewer-width 500 --label-width 28 --extension 0 --text-size 18
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session ${session} --bed cd163.bed --prefix hg19-cd163 --viewer-width 500 --label-width 28 --extension 0 --text-size 18
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session ${session} --bed pax7.bed --prefix hg19-pax7 --viewer-width 500 --label-width 28 --extension 0 --text-size 18
	"""

}

process marker_gene_screenshots_rn6 {

	publishDir "${params.results}/gb-screenshots/marker-genes"
	executor 'local'
	
	output:
	file("*.pdf")

	"""
	echo "chr10:53,588,847-53,909,741" >> regions.bed # MYH loci
	echo "chr10:760,194-769,669" >> regions.bed # MYH11
	echo "chr14:35,552,738-35,605,843" >> regions.bed # PDGFRA
	echo "chr4:158,080,951-158,094,049" >> regions.bed # VWF
	echo "chr4:156,742,565-156,775,939" >> regions.bed # CD163
	echo "chr5:158,278,488-158,327,893" >> regions.bed # PAX7
	cat regions.bed | perl -pe 's/,//g; s/[:-]/\t/g' > regions.tmp
	mv regions.tmp regions.bed
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session 2020-sn-muscle-cell-types-rn6 --bed regions.bed --prefix rn6 --viewer-width 1200 --label-width 28 --extension 0 --text-size 18
	"""

}

process marker_gene_heatmaps {

	publishDir "${params.results}/marker-gene-heatmaps"
	executor 'local'

	input:
	file(clusters) from Channel.fromPath("${ROOT}/work/downstream-new-features/results/liger/round-1/second-louvain-clusters.txt")
	file(counts) from Channel.fromPath("${ROOT}/work/downstream-new-features/results/features/same-space/*.counts.txt").toSortedList()

	output:
	file("*.pdf")

	"""
	make-cluster-feature-counts.py $clusters ${ROOT}/sample_info/sample_info.txt ${counts.join(' ')} > per-cluster-counts.txt
	plot-marker-gene-heatmap-atac-and-rna.R per-cluster-counts.txt $clusters
	"""

}

process ukb_enrichments_hg19 {

	publishDir "${params.results}/UKB-GWAS-enrichments/hg19"
	executor 'local'

	input:
	file(cluster_names) from Channel.fromPath(CLUSTER_NAMES)
	file(phenotype_tsv) from Channel.fromPath("${ROOT}/data/ukb-summary-stats/ukb31063_h2_all.02Oct2019.tsv.gz")
	file(results) from Channel.fromPath("${ROOT}/work/ldsc/UKB/hg19/joint-model/results/partitioned-heritability/*.results").toSortedList()

	output:
	file("*.pdf")

	"""
	plot-ldsc-ukb.R $cluster_names $phenotype_tsv ${results.join(' ')}
	"""

}

process ukb_enrichments_rn6 {

	publishDir "${params.results}/UKB-GWAS-enrichments/rn6"
	executor 'local'

	input:
	file(cluster_names) from Channel.fromPath(CLUSTER_NAMES)
	file(phenotype_tsv) from Channel.fromPath("${ROOT}/data/ukb-summary-stats/ukb31063_h2_all.02Oct2019.tsv.gz")
	file(results) from Channel.fromPath("${ROOT}/work/ldsc/UKB/rn6/joint-model/results/partitioned-heritability/*.results").toSortedList()

	output:
	file("*.pdf")

	"""
	plot-ldsc-ukb.R $cluster_names $phenotype_tsv ${results.join(' ')}
	"""

}


process plot_hg19_creatintine_coefficients {
	
	publishDir "${params.results}/UKB-GWAS-enrichments/hg19"
	executor 'local'

	input:
	file(cluster_names) from Channel.fromPath(CLUSTER_NAMES)
	file(results) from Channel.fromPath("${ROOT}/work/ldsc/UKB/hg19/joint-model/results/partitioned-heritability/tissues.30700_irnt.results")

	output:
	file("creatinine_coefficients.pdf")

	"""
	plot-ldsc-open-chromatin.R $results $cluster_names creatinine_coefficients.pdf Creatinine
	"""

}


process T2D_gwas_enrichments_hg19 {

	publishDir "${params.results}/T2D-GWAS-enrichments/hg19"
	executor 'local'
	cache false
	
	input:
	file(single_human) from Channel.fromPath("${ROOT}/work/ldsc/T2D/hg19/one-model-per-cluster/results/partitioned-heritability/*.results").toSortedList()
	file(joint_human) from Channel.fromPath("${ROOT}/work/ldsc/T2D/hg19/joint-model/results/partitioned-heritability/*.results").toSortedList()

	output:
	file("*.pdf")

	"""
	plot-T2D-enrichment.R ${CLUSTER_NAMES} ${single_human.join(' ')} ${joint_human.join(' ')}
	"""

}

process T2D_gwas_enrichments_rn6 {

	publishDir "${params.results}/T2D-GWAS-enrichments/rn6"
	executor 'local'
	cache false
	
	input:
	file(single_human) from Channel.fromPath("${ROOT}/work/ldsc/T2D/rn6/one-model-per-cluster/results/partitioned-heritability/*.results").toSortedList()
	file(joint_human) from Channel.fromPath("${ROOT}/work/ldsc/T2D/rn6/joint-model/results/partitioned-heritability/*.results").toSortedList()

	output:
	file("*.pdf")

	"""
	plot-T2D-enrichment.R ${CLUSTER_NAMES} ${single_human.join(' ')} ${joint_human.join(' ')}
	"""

}

process T2D_gwas_enrichments_both_species {

	publishDir "${params.results}/T2D-GWAS-enrichments/both"
	executor 'local'
	cache false
	
	input:
	file("human/*") from Channel.fromPath("${ROOT}/work/ldsc/T2D/hg19/one-model-per-cluster/results/partitioned-heritability/*.results").toSortedList()
	file("rat/*") from Channel.fromPath("${ROOT}/work/ldsc/T2D/rn6/one-model-per-cluster/results/partitioned-heritability/*.results").toSortedList()

	output:
	file("*.pdf")

	"""
	plot-T2D-enrichment-both-species.R ${CLUSTER_NAMES} human/* rat/*
	"""

}

process T2D_gwas_enrichments_hg19_no_bmiadj {

	publishDir "${params.results}/T2D-GWAS-enrichments-no-bmiadj/hg19"
	executor 'local'
	cache false
	
	input:
	file(single_human) from Channel.fromPath("${ROOT}/work/ldsc/T2D/hg19/one-model-per-cluster/results/partitioned-heritability/*.results").toSortedList()
	file(joint_human) from Channel.fromPath("${ROOT}/work/ldsc/T2D/hg19/joint-model/results/partitioned-heritability/*.results").toSortedList()

	output:
	file("*.pdf")

	"""
	plot-T2D-enrichment-no-bmiadj.R ${CLUSTER_NAMES} ${single_human.join(' ')} ${joint_human.join(' ')}
	"""

}

process T2D_gwas_enrichments_rn6_no_bmiadj {

	publishDir "${params.results}/T2D-GWAS-enrichments-no-bmiadj/rn6"
	executor 'local'
	cache false
	
	input:
	file(single_human) from Channel.fromPath("${ROOT}/work/ldsc/T2D/rn6/one-model-per-cluster/results/partitioned-heritability/*.results").toSortedList()
	file(joint_human) from Channel.fromPath("${ROOT}/work/ldsc/T2D/rn6/joint-model/results/partitioned-heritability/*.results").toSortedList()

	output:
	file("*.pdf")

	"""
	plot-T2D-enrichment-no-bmiadj.R ${CLUSTER_NAMES} ${single_human.join(' ')} ${joint_human.join(' ')}
	"""

}

process T2D_gwas_enrichments_both_species_no_bmiadj {

	publishDir "${params.results}/T2D-GWAS-enrichments-no-bmiadj/both"
	executor 'local'
	cache false
	
	input:
	file("human/*") from Channel.fromPath("${ROOT}/work/ldsc/T2D/hg19/one-model-per-cluster/results/partitioned-heritability/*.results").toSortedList()
	file("rat/*") from Channel.fromPath("${ROOT}/work/ldsc/T2D/rn6/one-model-per-cluster/results/partitioned-heritability/*.results").toSortedList()

	output:
	file("*.pdf")

	"""
	plot-T2D-enrichment-both-species-T2D-no-bmiadj.R ${CLUSTER_NAMES} human/* rat/*
	"""

}


process ukb_enrichment_table {

	publishDir "${params.results}/tables"
	executor 'local'

	input:
	file(cluster_names) from Channel.fromPath(CLUSTER_NAMES)
	file(phenotype_tsv) from Channel.fromPath("${ROOT}/data/ukb-summary-stats/ukb31063_h2_all.02Oct2019.tsv.gz")
	file(results) from Channel.fromPath("${ROOT}/work/ldsc/UKB/hg19/joint-model/results/partitioned-heritability/*.results").toSortedList()

	output:
	file('LDSC-UKB-Z-scores.tsv')

	"""
	make-ldsc-ukb-table.R $cluster_names $phenotype_tsv ${results.join(' ')}
	"""

}

process diamante_overlap_table {

	publishDir "${params.results}/tables"
	executor 'local'

	input:
	file(cluster_names) from Channel.fromPath(CLUSTER_NAMES)
	file(beta_peaks) from Channel.fromPath("${ROOT}/data/other-annotations/beta_ATAC.bed")
	file(adipose_peaks) from Channel.fromPath("${ROOT}/data/other-annotations/Adipose*.broadPeak").toSortedList()
	file(islet_peaks) from Channel.fromPath("${ROOT}/data/other-annotations/HP*").mix(Channel.fromPath("${ROOT}/data/other-annotations/AB*")).mix(Channel.fromPath("${ROOT}/data/other-annotations/AC*")).toSortedList()

	output:
	file('diamante-overlap-summary.csv')

	"""
	cat ${ROOT}/data/diamante-credible-sets/genetic_credible_sets/*.txt | grep -v IndexSNP | awk '{print(\$2, \$3, \$1)}' | perl -pe 's/ /\\t/g; s/^/chr/' > diamante-loci.txt
	mkdir bed
	cp ${ROOT}/data/chromhmm/Islets.all_promoter.bed bed/Islets.TSS_state.bed
	cp ${ROOT}/data/chromhmm/Adipose.all_promoter.bed bed/Adipose.TSS_state.bed
	cp ${ROOT}/data/chromhmm/Liver.all_promoter.bed bed/Liver.TSS_state.bed
	cp ${ROOT}/data/chromhmm/Islets.all_enhancer.bed bed/Islets.enhancer_state.bed
	cp ${ROOT}/data/chromhmm/Adipose.all_enhancer.bed bed/Adipose.enhancer_state.bed
	cp ${ROOT}/data/chromhmm/Liver.all_enhancer.bed bed/Liver.enhancer_state.bed
	cp ${ROOT}/data/gencode-coding/coding.bed bed/
	cat ${adipose_peaks.join(' ')} | sort -k1,1 -k2n,2 | bedtools merge -i stdin > bed/adipose_ATAC.bed
	cat ${islet_peaks.join(' ')} | sort -k1,1 -k2n,2 | bedtools merge -i stdin > bed/islet_ATAC.bed
	cat $beta_peaks > bed/beta_ATAC.bed
	cp ${MAIN}/process-by-cluster-round-1/peaks/0-hg19_peaks.broadPeak.noblacklist bed/cluster_0_ATAC.bed
	cp ${MAIN}/process-by-cluster-round-1/peaks/1-hg19_peaks.broadPeak.noblacklist bed/cluster_1_ATAC.bed
	cp ${MAIN}/process-by-cluster-round-1/peaks/2-hg19_peaks.broadPeak.noblacklist bed/cluster_2_ATAC.bed
	cp ${MAIN}/process-by-cluster-round-1/peaks/3-hg19_peaks.broadPeak.noblacklist bed/cluster_3_ATAC.bed
	cp ${MAIN}/process-by-cluster-round-1/peaks/4-hg19_peaks.broadPeak.noblacklist bed/cluster_4_ATAC.bed
	cp ${MAIN}/process-by-cluster-round-1/peaks/5-hg19_peaks.broadPeak.noblacklist bed/cluster_5_ATAC.bed
	cp ${MAIN}/process-by-cluster-round-1/peaks/6-hg19_peaks.broadPeak.noblacklist bed/cluster_6_ATAC.bed
	annot-snp-bed.py --bed-files bed/* --pos-file diamante-loci.txt > annotated-loci.txt	
	make-diamante-table.R $cluster_names annotated-loci.txt
	"""

}

process NOG_screenshots {

	publishDir "${params.results}/gb-screenshots/NOG"
	executor 'local'
	cache false

	output:
	file("*.pdf")

	"""
	echo "chr17:54735000-54800000" > regions.bed	
	cat regions.bed | perl -pe 's/,//g; s/[:-]/\t/g' > regions.tmp
	mv regions.tmp regions.bed
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session 2020-sn-muscle-creatinine-GWAS-overlap --bed regions.bed --prefix NOG --viewer-width 1200 --label-width 24 --extension 0 --text-size 18 --highlight-color FF0000 --highlight-region chr17:54776935-54776975 --genome hg19
	"""

}


process ARL15_screenshots {

	publishDir "${params.results}/gb-screenshots/ARL15"
	executor 'local'
	cache false

	output:
	file("*.pdf")

	"""
	echo "chr5:53,269,334-53,278,795" >> regions.bed	
	cat regions.bed | perl -pe 's/,//g; s/[:-]/\t/g' > regions.tmp
	mv regions.tmp regions.bed
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session 2020-sn-muscle-GWAS-overlap-muscle-cell-types-only --bed regions.bed --prefix ARL15 --viewer-width 1200 --label-width 24 --extension 0 --text-size 18 --highlight-color FF0000 --highlight-region chr5:53271399-53271440 --genome hg19
	echo "chr5:53,269,196-53,277,711" > regions.bed	
	cat regions.bed | perl -pe 's/,//g; s/[:-]/\t/g' > regions.tmp
	mv regions.tmp regions.bed
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session 2020-sn-muscle-GWAS-overlap-ATAC --bed regions.bed --prefix ARL15-other-cell-types --viewer-width 1000 --label-width 24 --extension 0 --text-size 24 --highlight-color FF0000 --highlight-region chr5:53271399-53271440 --genome hg19
	echo "chr2:45977100-45977546" > regions.rn6.bed
	echo "chr2:45977438-45977469" >> regions.rn6.bed
	cat regions.rn6.bed | perl -pe 's/,//g; s/[:-]/\t/g' > regions.tmp
	mv regions.tmp regions.rn6.bed
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session 2020-sn-muscle-cell-types-new-features-rn6 --genome rn6 --bed regions.rn6.bed --prefix ARL15_rn6 --viewer-width 1200 --label-width 28 --extension 10000 --highlight --highlight-color FF0000 --text-size 18
	"""

}

process ITPR2_screenshots {

	publishDir "${params.results}/gb-screenshots/ITPR2"
	executor 'local'
	cache false

	output:
	file("*.pdf")

	"""
	echo "chr12:26,436,640-26,491,955" >> regions.bed
	echo "chr12:26,469,532-26,475,063" >> regions.bed
	cat regions.bed | perl -pe 's/,//g; s/[:-]/\t/g' > regions.tmp
	mv regions.tmp regions.bed
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session 2020-sn-muscle-GWAS-overlap-muscle-cell-types-only --bed regions.bed --prefix ITPR2 --viewer-width 1200 --label-width 24 --extension 0 --text-size 18 --highlight-color FF0000 --highlight-region chr12:26472541-26472582 --genome hg19
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session 2020-sn-muscle-GWAS-overlap-ATAC --bed regions.bed --prefix ITPR2-other-cell-types --viewer-width 1200 --label-width 24 --extension 0 --text-size 18 --highlight-color FF0000 --highlight-region chr12:26472541-26472582 --genome hg19
	echo "chr4:180399117-180400404" > regions.rn6.bed
	echo "chr4:180400064-180400103" >> regions.rn6.bed
	cat regions.rn6.bed | perl -pe 's/,//g; s/[:-]/\t/g' > regions.tmp
	mv regions.tmp regions.rn6.bed
	python3 /home/porchard/src/gb_screenshot_utility/gb_screenshot.py --user porchard --session 2020-sn-muscle-cell-types-new-features-rn6 --bed regions.rn6.bed --genome rn6 --prefix ITPR2_rn6 --viewer-width 1200 --label-width 28 --extension 10000 --highlight --highlight-color FF0000 --text-size 18
	"""

}


process chromatin_state_legend {

	publishDir "${params.results}/chromatin-state-legend"
	executor 'local'

	input:
	file(bed) from Channel.fromPath("${ROOT}/data/chromhmm/Adipose.dense.bed")

	output:
	file("*.pdf")

	"""
	make_chromatin_state_legend.R $bed
	"""	

}
