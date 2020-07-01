#!/usr/bin/env nextflow

IONICE = 'ionice -c2 -n7'

// Generic data
AUTOSOMAL_REFERENCES = ['hg19': (1..22).collect({it -> 'chr' + it}),
			'hg38': (1..22).collect({it -> 'chr' + it}),
			'rn5': (1..20).collect({it -> 'chr' + it}),
			'rn6': (1..20).collect({it -> 'chr' + it}),
			'mm9': (1..19).collect({it -> 'chr' + it}),
			'mm10': (1..19).collect({it -> 'chr' + it})
]

ORGANISMS = ['hg19': 'human', 
		'hg38': 'human',
		'rn5': 'rat',
		'rn6': 'rat',
		'mm9': 'mouse',
		'mm10': 'mouse']

MACS2_GENOME_SIZE = [
    'rn4': 'mm',
    'rn5': 'mm',
    'rn6': 'mm',
    'mm9': 'mm',
    'mm10': 'mm',
    'hg19': 'hs',
    'hg38': 'hs'
]

get_autosomes = {
	genome ->
	AUTOSOMAL_REFERENCES[genome]
}

get_gtf = {
	genome ->
	params.gtf[genome]
}

get_macs2_genome_size = {
	genome ->
	return MACS2_GENOME_SIZE[genome]
}

get_rnaseq_qc = {
	library ->
	return params.libraries[library].qc
}

get_ensembl = {
	genome ->
	return params.ensembl[genome]
}

get_primary_ataqv_json = {
	library ->
	return params.libraries[library].ataqv_json
}

get_primary_counts_matrix = {
	library ->
	return params.libraries[library].counts
}

make_excluded_regions_arg = {
	genome ->
	return params.blacklist[genome].collect({'--excluded-region-file ' + it}).join(' ')
}

get_genome_size = {
	genome ->
	MACS2_GENOME_SIZE[genome]
}

get_genome = {
	library ->
	params.libraries[library].genome
}

get_tss = {
	genome ->
	params.tss[genome]
}

get_organism = {
	genome ->
	ORGANISMS[genome]
}

get_chrom_sizes = {
	genome ->
	params.chrom_sizes[genome]
}

get_gene_bed = {
	genome ->
	params.gene_bed[genome]
}

get_samples = {
	library ->
	params.libraries[library].samples
}

get_pruned = {
	library ->
	params.libraries[library].pruned
}

get_modality = {
	library ->
	params['libraries'][library]['modality']
}

get_starsolo_counts = {
	library ->
	params['libraries'][library]['starsolo'] // features, barcodes, matrix.mtx
}

libraries = params.libraries.keySet()

ATAC_LIBRARIES = []
RNA_LIBRARIES = []

for(library in libraries) {
	if (get_modality(library) == 'ATAC') {
		ATAC_LIBRARIES << library
	}
	if (get_modality(library) == 'RNA') {
		RNA_LIBRARIES << library
	}
}


get_ataqv_metrics_in = []
call_nuclei_rna_in = []
counts_to_tpm_matrix_atac_in = []
make_rna_qc_in = []
doubletfinder_counts_in = []

for (library in ATAC_LIBRARIES) {
	get_ataqv_metrics_in << [library, file(get_primary_ataqv_json(library))]
}

for (library in RNA_LIBRARIES) {
	call_nuclei_rna_in << [library, file(get_primary_counts_matrix(library))]
	doubletfinder_counts_in << [library, file(get_primary_counts_matrix(library))]
	make_rna_qc_in << [library, file(get_rnaseq_qc(library))]
}

for (library in ATAC_LIBRARIES) {
	counts_to_tpm_matrix_atac_in << [library, file(get_primary_counts_matrix(library))]
}

process make_rna_qc {

	publishDir "${params.results}/nucleus-qc", mode: 'rellink'
	memory '20 GB'
	container "${params.containers.general}"
	cache 'lenient'
	maxForks 1
	executor 'local'

	input:
	set val(library), file(counts) from Channel.from(call_nuclei_rna_in)

	output:
	file("${library}.txt") into make_rna_qc_out
	set val(library), file("${library}.counts.txt") into counts_to_tpm_matrix_rna_in

	"""
	ln -s $counts ${library}.counts.txt
	count-matrix-and-gtf-to-rnaseq-qc-input.py $counts ${get_gtf(get_genome(library))} > ${library}.txt
	"""

}


process get_ataqv_metrics {
	
	memory '40 GB'
	container "${params.containers.general}"
	cache 'lenient'
	executor 'local'
	maxForks 1

	input:
	set val(library), file(ataqv_json) from Channel.from(get_ataqv_metrics_in)

	output:
	file("${library}.metrics.txt") into call_nuclei_atac_in

	"""
	extractAtaqvMetric.py --files $ataqv_json --metrics tss_enrichment percent_hqaa hqaa total_reads total_autosomal_reads percent_mitochondrial percent_autosomal_duplicate percent_duplicate max_fraction_reads_from_single_autosome | perl -pe 's@.*.ataqv.json.gz\t@${library}-@' > ${library}.metrics.txt
	"""

}

initial_atac_thresholds = Channel.fromPath(params.qc_thresholds['initial-atac'])
process call_nuclei_atac {

	publishDir "${params.results}/nucleus-qc", mode: 'rellink'
	publishDir "${params.results}/figures", mode: 'rellink'
	container "${params.containers.general}"
	executor 'local'

	input:
	file(ataqv_metrics) from call_nuclei_atac_in.toSortedList()
	file(thresholds) from initial_atac_thresholds

	output:
	set file('atac-nuclei.txt'), file('hqaa-vs-tss-enrichment.png'), file('hqaa-vs-max-fraction-reads-from-single-autosome.png'), file('hg19-rn6-ratio-threshold.png')
	file('atac-nuclei.txt') into call_nuclei_atac_out
	file('atac-nuclei.txt') into demuxlet_barcodes
	file('metrics.txt') into call_nuclei_atac_plot_stats
	file('metrics.txt') into set_upper_hqaa_threshold_atac_metrics
	
	"""
	cat ${ataqv_metrics.join(' ')} > metrics.txt
	call-nuclei-atac.R metrics.txt $thresholds
	"""

}

process call_nuclei_rna {

	publishDir "${params.results}/nucleus-qc", mode: 'rellink'
	publishDir "${params.results}/figures", mode: 'rellink'
	container "${params.containers.general}"

	input:
	file(x) from make_rna_qc_out.toSortedList()

	output:
	set file('rna-nuclei.txt'), file('hqaa-vs-mitochondrial.png')
	file('rna-nuclei.txt') into call_nuclei_rna_out
	file('rna-nuclei.txt') into souporcell_barcodes
	file('rna-nuclei.txt') into doubletfinder_barcodes
	file('hqaa.txt') into call_nuclei_rna_plot_stats

	"""
	call-nuclei-rna-with-mitochondrial.R ${params.qc_thresholds['initial-rna']} ${x.join(' ')}
	"""

}

SOUPORCELL_LIBRARIES = ['63_20_rna-hg19', '63_40_rna-hg19']
DEMUXLET_LIBRARIES = ['63_20-hg19', '63_40-hg19']
DOUBLETFINDER_LIBRARIES = RNA_LIBRARIES


process doubletfinder {

	publishDir "${params.results}/doubletfinder", mode: 'rellink'
	container "${params.containers.doubletfinder}"
	cache 'deep'
	memory '40 GB'
	maxForks 5

	input:
	set val(library), file(counts), file(barcodes) from Channel.from(doubletfinder_counts_in).combine(doubletfinder_barcodes)

	output:
	file("${library}.doubletfinder-assignments.txt") into doubletfinder_out
	
	"""
	filter-features.py --keep-nuclei $barcodes $counts > counts.txt
	mkdir -p $library
	doublet_finder.R --input counts.txt --max_dims 20 --outdir $library --sample $library
	ln -s assignments.txt ${library}.doubletfinder-assignments.txt
	"""
}


process souporcell {

	publishDir "${params.results}/souporcell", mode: 'rellink'
	container "${params.containers.souporcell}"
	cache 'deep'
	memory '75 GB'
	cpus 8

	input:
	file(barcodes) from souporcell_barcodes
	each name from SOUPORCELL_LIBRARIES
	
	output:
	set val(name), file("${name}.souporcell.txt"), file("${name}.souporcell.vcf") into souporcell_out
	
	script:
	bam = get_pruned(name)

	"""
	mkdir -p ${name}
        grep -w ${name} $barcodes | cut -f2 > barcodes.txt
        souporcell_pipeline.py -i $bam -b barcodes.txt -f ${params.fasta['hg19']} -t 8 -o $name -k 2
	cp ${name}/clusters.tsv ${name}.souporcell.txt
	cp ${name}/cluster_genotypes.vcf ${name}.souporcell.vcf
        """

}

process match_souporcell_to_ksm {

	publishDir "${params.results}/souporcell", mode: 'rellink'
	container "${params.containers.general}"
	cache 'deep'

	input:
	set val(library), file(clusters), file(vcf) from souporcell_out

	output:
	file("${library}.souporcell-ksm.txt") into souporcell_refactored

	"""
	match-souporcell-individuals-to-samples.py $vcf ${params.genotypes} $clusters $library ${library}.souporcell-ksm.txt
	"""
	
}


process process_souporcell_and_doubletfinder_out {

	publishDir "${params.results}/souporcell", mode: 'rellink'
	container "${params.containers.general}"
	cache 'deep'

	input:
	file(x) from souporcell_refactored.toSortedList()
	file(y) from doubletfinder_out.toSortedList()

	output:
	file("snRNA-doublets.pdf")
	file("rna-souporcell-assignments.txt") into souporcell_individual_assignments
	file('rna-nuclei-with-individuals.txt') into rna_nuclei_individuals

	"""
	cat ${y.join(' ')} > doubletfinder.txt
	cat ${x.join(' ')} > souporcell.txt
	process-souporcell-out.R souporcell.txt doubletfinder.txt
	grep -v -w doublet rna-souporcell-assignments.txt > rna-nuclei-with-individuals.txt
	"""
	
}


process demuxlet {

	publishDir "${params.results}/demuxlet", mode: 'rellink'
	container "${params.containers.demuxlet}"
	memory '50 GB'
	cache 'deep'
	
	input:
	file(barcodes) from demuxlet_barcodes
	each library from DEMUXLET_LIBRARIES

	output:
	file("${library}.best")
	file("${library}.single")
	file("${library}.sing2")
	file("${library}.best") into demuxlet_out

	script:
	bam = get_pruned(library)

	"""
        grep -w ${library} $barcodes | cut -f2 > barcodes.txt
	ln -s $bam ${library}.bam
	samtools index ${library}.bam
	demuxlet --sam ${library}.bam --vcf ${params.genotypes} --group-list barcodes.txt --field PL --out ${library}
	"""

}

process process_demuxlet_out {

	publishDir "${params.results}/demuxlet", mode: 'rellink'
	container "${params.containers.general}"
	cache 'deep'

	input:
	file(x) from demuxlet_out.toSortedList()

	output:
	file("snATAC-doublets.pdf")
	file("atac-demuxlet-assignments.txt") into demuxlet_individual_assignments

	"""
	process-demuxlet-out.R ${x.join(' ')}
	"""

}

process tighten_qc_thresholds_atac {

	publishDir "${params.results}/upper-atac-thresholds"
	container "${params.containers.rplot}"
	cache 'deep'

	input:
	file(metrics) from set_upper_hqaa_threshold_atac_metrics
	file(demuxlet_assignments) from demuxlet_individual_assignments
	file(atac_nuclei) from call_nuclei_atac_out

	output:
	file("*.pdf")
	file("atac-max-hqaa-thresholds.txt")
	file('atac-nuclei-with-individuals.txt') into atac_nuclei_individuals

	"""
	tighten-qc-thresholds-atac.R $demuxlet_assignments $atac_nuclei $metrics
	"""

}


process concat_nuclei {

	publishDir "${params.results}/nucleus-qc", mode: 'rellink'
	container "${params.containers.general}"
	cache 'deep'

	input:
	file(atac) from atac_nuclei_individuals
	file(rna) from rna_nuclei_individuals

	output:
	file('nuclei-with-individual-assignments.txt') into liger_individuals
	file('nuclei-with-individual-assignments.txt') into liger_individuals_2
	file('nuclei-with-individual-assignments.txt') into concat_nuclei_individuals
	file('nuclei.txt') into concat_nuclei_out
	file('nuclei.txt') into concat_nuclei_aggregate
	file('nuclei.txt') into concat_nuclei_plot_stats

	"""
	cat $atac $rna | sort | uniq > nuclei-with-individual-assignments.txt
	cut -f1-2 nuclei-with-individual-assignments.txt > nuclei.txt
	"""

}


process plot_nuclei_counts_and_coverage {

	publishDir "${params.results}/figures", mode: 'rellink'
	container "${params.containers.general}"

	input:
	file(nuclei) from concat_nuclei_plot_stats
	file(ataqv_metrics) from call_nuclei_atac_plot_stats
	file(hqaa) from call_nuclei_rna_plot_stats

	output:
	set file('nuclei-read-counts.pdf'), file('nuclei-per-library.pdf'), file('nuclei-per-modality.pdf')

	"""
	plot-nuclei-counts-and-coverage.R --nuclei $nuclei --atac-metrics $ataqv_metrics --rna-metrics $hqaa
	"""

}


filter_counts_matrix_in = Channel.from(counts_to_tpm_matrix_atac_in).mix(counts_to_tpm_matrix_rna_in).combine(concat_nuclei_out)
process filter_counts_matrix {

	publishDir "${params.results}/features/raw", mode: 'rellink'
	container "${params.containers.general}"
	
	input:
	set val(library), file("unfiltered.txt"), file(nuclei) from filter_counts_matrix_in

	output:
	set val(library), file("${library}.counts.txt") into counts_to_tpm_in
	set val(library), val('counts'), file("${library}.counts.txt") into filter_counts_matrix_same_space
	set val(library), val('counts'), file("${library}.counts.txt") into translate_rat_counts

	"""
	filter-features.py --keep-nuclei $nuclei unfiltered.txt | perl -pe 's/\\tchr.*\\d:/\\t/' | sort -k1,1 -k2,2 -k3,3 | bedtools groupby -g 1,2,3 -c 4 -o sum | grep -w -v 'nan' > ${library}.counts.txt
	"""

}


process counts_to_tpm {
	
	publishDir "${params.results}/features/raw", mode: 'rellink'
	container "${params.containers.general}"

	input:
	set val(library), file("counts") from counts_to_tpm_in

	output:
	set val(library), file("${library}.tpm.txt") into tpm_out
	set val(library), val('tpm'), file("${library}.tpm.txt") into tpm_same_space
	set val(library), val('tpm'), file("${library}.tpm.txt") into translate_rat_tpm

	"""
	feature-counts-to-tpm.py counts > ${library}.tpm.txt
	"""

}


process translate_gene_names_rat {

	container "${params.containers.general}"

	input:
	set val(library), val(category), file("features.txt") from translate_rat_tpm.mix(translate_rat_counts)

	output:
	set val(library), val(category), file("${library}.${category}.txt") into translate_gene_names_rat_out

	when:
	get_genome(library) == 'rn6'

	"""
	translate-feature-file-feature-names.py --from rn6 --to hg19 --drop-if-missing features.txt ${params.orthologues['gene-name']} > ${library}.${category}.txt
	"""	

}


same_space_in = filter_counts_matrix_same_space.filter({get_genome(it[0].toString()) == 'hg19'}).mix(tpm_same_space.filter({get_genome(it[0].toString()) == 'hg19'})).mix(translate_gene_names_rat_out)

process same_space {

	cache 'deep'
	publishDir "${params.results}/features/same-space", mode: 'rellink'
	container "${params.containers.general}"

	input:
	set val(library), val(category), file(x) from same_space_in

	output:
	set val(library), val(category), file(x) into seurat_in
	set val(library), val(category), file(x) into liger_round_1_in
	set val(library), val(category), file(x) into liger_round_1_in_2
	set val(library), val(category), file(x) into interesting_genes_tpm_in
	set val(library), val(category), file(x) into per_library_per_cluster_counts_features
	set val(library), val(category), file(x) into scina_in

	"""
	echo 'pass'
	"""	

}

// Process the aggregate of quality nuclei for each library
aggregate_rna_in = []
aggregate_atac_in = []
aggregate_atac_by_individual_in = []

for (library in RNA_LIBRARIES) {
	aggregate_rna_in << [library, file(get_pruned(library))]
}
for (library in ATAC_LIBRARIES) {
	aggregate_atac_in << [library, file(get_pruned(library))]
	for (samp in get_samples(library)) {
		aggregate_atac_by_individual_in << [library, samp, file(get_pruned(library))]
	}
}

process aggregate_bam {

	publishDir "${params.results}/aggregate-of-quality-nuclei/bam", mode: 'rellink'
	container "${params.containers.general}"

	input:
	set val(library), file(bam), file(nuclei) from Channel.from(aggregate_rna_in).mix(Channel.from(aggregate_atac_in)).combine(concat_nuclei_aggregate)

	output:
	set val(library), file("${library}.bam") into per_library_per_cluster_aggregate_bam_in
	set val(library), file("${library}.bam"), file("${library}.bam.bai") into aggregate_macs2_in
	set val(library), file("${library}.bam"), file("${library}.bam.bai") into ataqv_bam

	"""
	grep -w $library $nuclei | cut -f2 > keep-nuclei.txt
	filter-bam-by-barcode.py $bam ${library}.bam keep-nuclei.txt
	samtools index ${library}.bam
	rm keep-nuclei.txt
	"""

}


process aggregate_macs2 {

	publishDir "${params.results}/aggregate-of-quality-nuclei/peaks", mode: 'rellink'
	container "${params.containers.macs2}"

	input:
	set val(library), file(bam), file(index) from aggregate_macs2_in

	output:
	set val(library), file("${library}_peaks.broadPeak"), file("${library}_peaks.broadPeak.noblacklist"), file("${library}_treat_pileup.bdg") into macs2_out
	set val(library), file("${library}_treat_pileup.bdg") into aggregate_bigwig_in
	set val(library), file("${library}_peaks.broadPeak.noblacklist") into ataqv_peaks
	set val(library), file("${library}_peaks.broadPeak.noblacklist") into chromatin_state_overlap

	when:
	get_modality(library) == 'ATAC'

	"""
	bedtools bamtobed -i $bam > ${library}.bed
	macs2 callpeak -t ${library}.bed --outdir . -f BED -n ${library} -g ${get_macs2_genome_size(get_genome(library))} --nomodel --shift -100 --seed 762873 --extsize 200 -B --SPMR --broad --keep-dup all
	bedtools intersect -a ${library}_peaks.broadPeak -b ${params.blacklist[get_genome(library)].collect().join(' ')} -v > ${library}_peaks.broadPeak.noblacklist
	rm ${library}.bed
	"""

}


process aggregate_atac_bigwig {

	publishDir "${params.results}/aggregate-of-quality-nuclei/bigwig", mode: 'rellink'
	memory '40 GB'
	container "${params.containers.general}"

	input:
	set val(library), file(bedgraph) from aggregate_bigwig_in

	output:
	file("${library}.bw")

	"""
	bedSort $bedgraph sorted.bdg
	bedClip sorted.bdg ${get_chrom_sizes(get_genome(library))} trimmed.bdg
	bedGraphToBigWig trimmed.bdg ${get_chrom_sizes(get_genome(library))} ${library}.bw
	"""

}


process ataqv {

	publishDir "${params.results}/aggregate-of-quality-nuclei/ataqv", mode: 'rellink'
	container "${params.containers.snATAC}"

	input:
	set val(library), val(bam), file(index), file(peaks) from ataqv_bam.combine(ataqv_peaks, by: 0)

	output:
	set file("${library}.ataqv.json.gz"), file("${library}.ataqv.out")

	when:
	get_modality(library) == 'ATAC'

	"""
	ataqv --name $library --peak-file $peaks --metrics-file ${library}.ataqv.json.gz --tss-file ${get_tss(get_genome(library))} ${make_excluded_regions_arg(get_genome(library))} --ignore-read-groups ${get_organism(get_genome(library))} $bam > ${library}.ataqv.out
	"""

}


process chromatin_state_overlap {

        publishDir "${params.results}/aggregate-of-quality-nuclei/chromhmm-overlap", mode: 'rellink'
	container "${params.containers.general}"

        input:
        set val(library), file(peaks) from chromatin_state_overlap

        output:
        set file("${library}.chromhmm_overlap.txt"), file("${library}.chromhmm_overlap.pdf")

        when:
        get_genome(library) == 'hg19'
        
	"""
        chromhmm_overlap.py $peaks ${get_tss(get_genome(library))} ${get_chrom_sizes(get_genome(library))} ${params.chromatin_state_glob} > ${library}.chromhmm_overlap.txt
        plot_chromatin_state_overlap.R --overlap ${library}.chromhmm_overlap.txt --out ${library}.chromhmm_overlap.pdf
        """

}


// Filter TPM matrices down to just interesting genes
process interesting_genes_tpm {

	publishDir "${params.results}/features/interesting-genes", mode: 'rellink'
	container "${params.containers.general}"

	input:
	set val(library), val(category), file("features.txt") from interesting_genes_tpm_in

	output:
	file("${library}.${category}.txt") into interesting_genes_tpm_out
	set val(library), file("${library}.${category}.txt") into plot_gene_expression_on_umap

	when:
	category == 'tpm'
	
	"""
	echo ${params.interesting_genes.join(' ')} | perl -pe 's/ /\\n/g' > keep.txt
	filter-features.py --keep-features keep.txt features.txt > ${library}.${category}.txt
	"""

}


prep_counts_for_merge_in = liger_round_1_in_2.filter({x -> params.liger_exclude.indexOf(x[0].toString()) == -1 && x[1].toString() == 'counts'}).map({x -> x[2]}).toSortedList()
process prep_counts_for_liger_merge {

	publishDir "${params.results}/liger/round-1/input-merged", mode: 'rellink'
	container "${params.containers.general}"
	memory '50 GB'

	input:
	file(counts) from prep_counts_for_merge_in
	file(individuals) from liger_individuals_2

	output:
	set val("x"), file("*.liger.txt") into liger_round_1_in_filtered

	"""
	cat ${counts.join(' ')} > all_counts.txt
	make-liger-input-matrices-merge.py /lab/work/porchard/sn-muscle-project/library-labels.txt all_counts.txt $individuals
	rm all_counts.txt
	"""

}

process liger_factorize {

	publishDir "${params.results}/liger/round-1", mode: 'rellink', overwrite: true
	container "${params.containers.liger}"
	memory '75 GB'

	input:
	file(mats) from liger_round_1_in_filtered.map({x -> x[1]}).flatten().toSortedList()

	output:
	file("liger-factorized.Rda") into liger_normalize_visualize_in

	"""
	liger-factorize-with-individuals.R --mats ${mats.join(',')} --factorization_k 15 --factorization_lambda 5 --out liger-factorized.Rda
	"""

}

process liger_normalize_visualize {

	publishDir "${params.results}/liger/round-1", mode: 'rellink'
	container "${params.containers.liger}"
	memory '75 GB'

	input:
	file(rda) from liger_normalize_visualize_in

	output:
	file("shared-factor-markers.txt")
	file("umap.txt") into plot_clusters_dim_in
	file("liger-normalized.Rda") into liger_recluster_in

	"""
	liger-normalize-visualize.R --rda $rda --out liger-normalized.Rda --norm_resolution 0.05 --norm_knnk 10 --dim umap.txt
	"""

}

process liger_recluster {

	publishDir "${params.results}/liger/round-1", mode: 'rellink'
	container "${params.containers.liger}"
	memory '75 GB'

	input:
	file(rda) from liger_recluster_in

	output:
	file("liger-reclustered.Rda") into liger_reclustered_rda
	file("second-louvain-clusters.txt") into plot_clusters_clusters_in
	file("second-louvain-clusters.txt") into per_library_per_cluster_clusters_in
	file("second-louvain-clusters.txt") into per_library_per_cluster_counts_clusters_in
	file("second-louvain-clusters.txt") into clusters

	"""
	liger-recluster.R --rda $rda --out liger-reclustered.Rda --k 17 --resolution 0.05 --cluster_file second-louvain-clusters.txt
	"""

}

process liger_qc_plots {

	publishDir "${params.results}/liger/round-1/qc-plots", mode: 'rellink'
	container "${params.containers.liger}"
	memory '110 GB'

	input:
	file(rda) from liger_reclustered_rda

	output:
	file("*.pdf")

	"""
	liger-plot-factors.R --rda $rda --prefix liger-qc.
	"""

}

plot_clusters_gene_expression_in = plot_gene_expression_on_umap.filter({x -> get_modality(x[0].toString()) == 'RNA'}).map({x -> x[1]}).toSortedList()
process plot_clusters {

	publishDir "${params.results}/liger/round-1", mode: 'rellink'
	container "${params.containers.rplot}"
	executor 'local'

	input:
	file(clusters) from plot_clusters_clusters_in
	file(dim) from plot_clusters_dim_in
	file(tpms) from plot_clusters_gene_expression_in

	output:
	set file("liger.umap-by-modality.pdf"), file("liger.umap-by-library.pdf"), file("liger.number-nuclei-per-cluster-per-library.pdf"), file("liger.fraction-cluster-per-library.pdf"), file("liger.clusters.pdf"), file("liger.fraction-library-per-cluster.pdf"), file("liger.gene_expression.png"), file("liger.clusters.facetted.png"), file("liger.gene_expression*.png")

	"""
	plot-dim-by-library-facets.R $dim liger.
	plot-library-representation-per-cluster.R $clusters liger.
	plot-clusters.R $dim $clusters liger.clusters.pdf
	plot-clusters-facets.R $dim $clusters liger.clusters.facetted.png
	cat ${tpms.join(' ')} > tpms.txt
	overlay-tpm-on-umap.R liger.gene_expression.png $dim tpms.txt
	"""	

}



// get count/TPM matrices for each cluster
// TODO: for rat and RNA, should collect the counts from the untranslated count matrix -- otherwise, will overestimate TPM for all genes that were successfully translated...
process per_cluster_counts {

	publishDir "${params.results}/process-by-cluster-round-1/features", mode: 'rellink'
	container "${params.containers.general}"
	memory '100 GB'

	input:
	file(counts) from per_library_per_cluster_counts_features.filter({x -> x[1].toString() == 'counts'}).map({x -> x[2]})
	file(clusters) from per_library_per_cluster_counts_clusters_in
	
	output:
	file('per-cluster-gene-counts.txt')

	"""
	make-cluster-feature-counts.py $clusters /lab/work/porchard/sn-muscle-project/sample_info/sample_info.txt ${counts.join(' ')} > per-cluster-gene-counts.txt
	"""

}

// Process by cluster
process per_library_per_cluster_bam {

	publishDir "${params.results}/process-by-cluster-round-1/bam", mode: 'rellink'
	container "${params.containers.general}"
	maxForks 10

	input:
	set val(library), file("bam-to-split"), file(clusters) from per_library_per_cluster_aggregate_bam_in.combine(per_library_per_cluster_clusters_in)

	output:
	file("*.bam") into per_library_per_cluster_bam_out

	when:
	params.liger_exclude.indexOf(library) == -1

	script:
	genome = get_genome(library)
	modality = get_modality(library)
	
	"""
	grep -w $library $clusters | cut -f2,3 > barcode_to_cluster.txt
        split-bam-by-cluster.py bam-to-split ${library}. barcode_to_cluster.txt
	"""

}


process per_library_per_cluster_bam_redirect {

	executor 'local'

	input:
	file(bam) from per_library_per_cluster_bam_out.flatten()

	output:
	set val(cluster), val(genome), val(modality), file(bam) into merge_in
	set val(cluster), val(genome), val(library), file(bam) into peak_counts_bam_in

	script:
	library = bam.getName().split(/\./)[0]
	cluster = bam.getName().split(/\./)[1]
	genome = get_genome(library)
	modality = get_modality(library)

	"""
	echo $bam
	"""

}

process merge_cluster {

	publishDir "${params.results}/process-by-cluster-round-1/bam", mode: 'rellink'
	container "${params.containers.general}"
	maxForks 5

	input:
	set val(cluster), val(genome), val(modality), file(bams) from merge_in.groupTuple(by: [0, 1, 2])

	output:
	set val(cluster), val(genome), val(modality), file("${cluster}.${genome}.${modality}.bam") into bamtobed_in

	"""
	samtools merge ${cluster}.${genome}.${modality}.bam ${bams.join(' ')}
	"""

}


process per_cluster_bamtobed {

	container "${params.containers.general}"
	tag "${cluster}"
	maxForks 5

	input:
	set val(cluster), val(genome), val(modality), file(bam) from bamtobed_in

	output:
	set val(cluster), val(genome), file("reads.bed") into broadpeaks_in
	set val(cluster), val(genome), file("reads.bed") into narrowpeaks_in

	when:
	modality == 'ATAC'

	"""
	bedtools bamtobed -i $bam > reads.bed
	"""

}

process peaks {

	publishDir "${params.results}/process-by-cluster-round-1/peaks", mode: 'rellink'
	container "${params.containers.macs2}"
	tag "${cluster}"

	input:
	set val(cluster), val(genome), file(reads) from broadpeaks_in

	output:
	set val(cluster), val(genome), file("${cluster}-${genome}_treat_pileup.bdg") into cluster_bigwigs_in
	set val(genome), file("${cluster}-${genome}_peaks.broadPeak.noblacklist") into master_peaks_in
	set val(cluster), val(genome), file("${cluster}-${genome}_peaks.broadPeak.noblacklist") into peak_counts_peaks_in
	set val(cluster), val(genome), file("${cluster}-${genome}_peaks.broadPeak.noblacklist") into cluster_chromatin_state_overlap_hg19_in
	set val(cluster), val(genome), file("${cluster}-${genome}_peaks.broadPeak.noblacklist") into enhancer_regression_accessibility_in
	set val(cluster), val(genome), file("${cluster}-${genome}_peaks.broadPeak.noblacklist") into peak_liftover_in

	"""
        macs2 callpeak -t $reads --outdir . -f BED -n ${cluster}-${genome} --SPMR -g ${get_genome_size(genome)} --nomodel --shift -100 --seed 762873 --extsize 200 -B --broad --keep-dup all
        bedtools intersect -a ${cluster}-${genome}_peaks.broadPeak -b ${params.blacklist[genome].collect().join(' ')} -v | sort -k1,1 -k2n,2 > ${cluster}-${genome}_peaks.broadPeak.noblacklist
	"""

}

process narrow_peaks {

	publishDir "${params.results}/process-by-cluster-round-1/narrowpeaks", mode: 'rellink'
	container "${params.containers.macs2}"
	tag "${cluster}"

	input:
	set val(cluster), val(genome), file(reads) from narrowpeaks_in

	output:
	set val(genome), val(cluster), file("${cluster}-${genome}_peaks.narrowPeak.noblacklist") into enrichment_features_in

	"""
        macs2 callpeak -t $reads --outdir . -f BED -n ${cluster}-${genome} --SPMR -g ${get_genome_size(genome)} --nomodel --shift -37 --seed 762873 --extsize 73 -B --keep-dup all
        bedtools intersect -a ${cluster}-${genome}_peaks.narrowPeak -b ${params.blacklist[genome].collect().join(' ')} -v | sort -k1,1 -k2n,2 > ${cluster}-${genome}_peaks.narrowPeak.noblacklist
	"""

}

process rat_peak_liftover {

	publishDir "${params.results}/rat-peak-liftover", mode: 'rellink'
	container "${params.containers.bnmapper}"

	input:
	set val(cluster), val(genome), file(peaks) from peak_liftover_in

	output:
	set val(cluster), val(genome), file(outfile) into cluster_chromatin_state_overlap_rn6_in

	when:
	genome == "rn6"

	script:
	outfile = cluster + '-' + genome + '_peaks_in_hg19.bed'

	"""
	bnMapper.py -o ${outfile}.tmp $peaks /lab/work/porchard/data/chain/rn6ToHg19.over.chain.gz
	cat ${outfile}.tmp | sort -k1,1 -k2n,2 | bedtools merge -i stdin -d 50 > $outfile
	"""	

}


process cluster_chromatin_state_overlap {

        publishDir "${params.results}/process-by-cluster-round-1/chromhmm-overlap", mode: 'rellink'
	container "${params.containers.general}"

        input:
        set val(cluster), val(genome), file(peaks) from cluster_chromatin_state_overlap_hg19_in.filter({x -> x[1].toString() == 'hg19'}).mix(cluster_chromatin_state_overlap_rn6_in)

        output:
        set file("${cluster}-${genome}.chromhmm_overlap.txt"), file("${cluster}-${genome}.chromhmm_overlap.pdf")

	"""
        chromhmm_overlap.py $peaks ${get_tss(genome)} ${get_chrom_sizes(genome)} ${params.chromatin_state_glob} > ${cluster}-${genome}.chromhmm_overlap.txt
        plot_chromatin_state_overlap.R --overlap ${cluster}-${genome}.chromhmm_overlap.txt --out ${cluster}-${genome}.chromhmm_overlap.pdf
        """

}

peak_counts_in = peak_counts_bam_in.filter({x -> get_modality(x[2].toString()) == 'ATAC'}).combine(peak_counts_peaks_in, by: [0, 1])
process peak_counts {

	memory '10 GB'
	cpus 5
	container "${params.containers.general}"
	maxForks 10

	input:
	set val(cluster), val(genome), val(library), file(bam), file(peaks) from peak_counts_in

	output:
	set val(cluster), val(genome), file("${library}.${cluster}.counts.txt") into peak_counts_out

	"""
	sort-bed-by-bam.py $peaks $bam > peaks.sorted.bed
	bedtools intersect -wa -wb -bed -sorted -a $bam -b peaks.sorted.bed | cut -f4,13-15 | perl -pe 's/\\t/--/g; s/--/\\t/; s/--/:/g' | sort --parallel=5 -S 5G | perl -pe 's@.*_(.*)/\\d+\\t(.*)@\$1\\t\$2@' | perl -pe 's/:/_/g' | uniq -c > counts.intermediate.txt
        cat counts.intermediate.txt | perl -pe 's/^\\s+//; s/\\s+/\\t/' | awk '{print(\$2, \$3, \$1)}' | perl -pe 's/ /\\t/g; s/^/${library}\\t/' | sort --parallel=5 -S 5G | bedtools groupby -g 1,2,3 -c 4 -o sum > ${library}.${cluster}.counts.txt
	"""

}

process per_cluster_peak_counts {

	publishDir "${params.results}/peak-counts", mode: 'rellink'
	container "${params.containers.general}"
	memory '5 GB'

	input:
	set val(cluster), val(genome), file("counts*.txt") from peak_counts_out.groupTuple(by: [0, 1])

	output:
	set val(cluster), val(genome), file("${cluster}-${genome}.counts.txt") into per_cluster_peak_counts_out

	"""
	cat counts*.txt > ${cluster}-${genome}.counts.txt
	"""

}

process master_peaks {

	publishDir "${params.results}/process-by-cluster-round-1/master-peaks", mode: 'rellink'
	container "${params.containers.general}"

	input:
	set val(genome), file(peaks) from master_peaks_in.groupTuple()

	output:
	file("master-peaks.${genome}.bed") into enhancer_posteriors_master_peaks_in

	"""
	cat ${peaks.join(' ')} | sort -k1,1 -k2n,2 | bedtools merge -i stdin > master-peaks.${genome}.bed
	"""

}

process cluster_bigwigs {

	publishDir "${params.results}/process-by-cluster-round-1/atac-bigwigs", mode: 'rellink'
	container "${params.containers.general}"
	memory '50 GB'
	tag "${cluster}-${genome}"

	input:
	set val(cluster), val(genome), file(bedgraph) from cluster_bigwigs_in

	output:
	file("${cluster}-${genome}.bw") into bigwig_out

	"""
	bedClip $bedgraph ${get_chrom_sizes(genome)} bedgraph.clipped.bdg
        LC_COLLATE=C sort -k1,1 -k2n,2 bedgraph.clipped.bdg > bedgraph.sorted.bdg
        bedGraphToBigWig bedgraph.sorted.bdg ${get_chrom_sizes(genome)} ${cluster}-${genome}.bw
	"""
	
}
