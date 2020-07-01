#!/usr/bin/env nextflow

// Mapping has already been performed; just need the post-mapping merged files
MAPPED_BAM_GLOB = params.mapped_bam_glob 

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


make_excluded_regions_arg = {
	genome ->
	return params.blacklist[genome].collect({'--excluded-region-file ' + it}).join(' ')
}


get_blacklists = {
	genome ->
	params.blacklist[genome]
}


get_genome = {
	library ->
	return library.tokenize('-')[1]
}


get_tss = {
	genome ->
	params.tss[genome]
}


get_organism = {
	genome ->
	ORGANISMS[genome]
}

get_macs2_genome_size = {
	genome ->
	MACS2_GENOME_SIZE[genome]
}


get_chrom_sizes = {
	genome ->
	params.chrom_sizes[genome]
}


process mark_duplicates {

	publishDir "${params.results}/mark_duplicates", mode: 'rellink'
	errorStrategy 'retry'
	maxRetries 1
	time '5h'
	maxForks 15

	input:
	file(bam) from Channel.fromPath(MAPPED_BAM_GLOB)

	output:
	set val(library), file("${library}.md.bam"), file("${library}.md.bam.bai"), file("${library}.metrics") into md_bams
	set val(library), file("${library}.md.bam"), file("${library}.md.bam.bai") into prune_in
	set val(library), file("${library}.md.bam"), file("${library}.md.bam.bai") into ataqv_md_in

	script:
	library = bam.getName().replaceAll('.bam', '')

	"""
	java -Xmx4g -Xms4g -jar \$PICARD_JAR MarkDuplicates I=$bam O=${library}.md.bam ASSUME_SORTED=true METRICS_FILE=${library}.metrics VALIDATION_STRINGENCY=LENIENT
	samtools index ${library}.md.bam
	"""
}


process prune {

	memory '3 GB'
	time '5h'
	errorStrategy 'retry'
	maxRetries 3
	maxForks 10

	publishDir "${params.results}/prune", mode: 'rellink'

	input:
	set val(library), file(bam), file(bam_index) from prune_in

	output:
	set val(library), file("${library}.pruned.bam") into pruned_bams

	"""
	${IONICE} samtools view -h -b -f 3 -F 4 -F 8 -F 256 -F 1024 -F 2048 -q 30 $bam ${AUTOSOMAL_REFERENCES[get_genome(library)].join(' ')} > ${library}.pruned.bam 
	"""

}


process bamtobed {

	time '4h'
	maxForks 10

	input:
	set val(library), file(bam) from pruned_bams
	
	output:
	set val(library), file("${library}.bed") into beds

	"""
	bedtools bamtobed -i $bam > ${library}.bed
	"""

}


process macs2 {

	publishDir "${params.results}/macs2", mode: 'rellink'
	time '5h'

	input:
	set val(library), file(bed) from beds

	output:
	set val(library), file("${library}_peaks.broadPeak") into blacklist_in
	set val(library), file("${library}_peaks.broadPeak") into ataqv_macs2_in
	set val(library), file("${library}_treat_pileup.bdg") into bigwig_in

	"""
	macs2 callpeak -t $bed --outdir . --SPMR -f BED -n $library -g ${get_macs2_genome_size(get_genome(library))} --nomodel --shift -100 --seed 762873 --extsize 200 -B --broad --keep-dup all
	"""

}


process blacklist_filter_peaks {

	publishDir "${params.results}/macs2", mode: 'rellink'
	time '1h'

	input:
	set val(library), file(peaks) from blacklist_in

	output:
	file("${library}_peaks.broadPeak.noblacklist")

	"""
	bedtools intersect -a $peaks -b ${get_blacklists(get_genome(library)).join(' ')} -v > ${library}_peaks.broadPeak.noblacklist
	"""

}


process bigwig {

	time '5h'
	publishDir "${params.results}/bigwig", mode: 'rellink'

	input:
	set val(library), file(bedgraph) from bigwig_in

	output:
	file("${library}.bw")

	"""
	LC_COLLATE=C sort -k1,1 -k2n,2 $bedgraph > sorted.bedgraph
	bedClip sorted.bedgraph ${get_chrom_sizes(get_genome(library))} clipped.bedgraph
	bedGraphToBigWig clipped.bedgraph ${get_chrom_sizes(get_genome(library))} ${library}.bw
	rm sorted.bedgraph clipped.bedgraph
	"""	
	
}


ataqv_in = ataqv_md_in.combine(ataqv_macs2_in, by: 0)

process ataqv {
	
	publishDir "${params.results}/ataqv", mode: 'rellink'
	errorStrategy 'retry'
	maxRetries 1
	memory '5 GB'
	time '3h'
	
	input:
	set val(library), file(bam), file(bam_index), file(peaks) from ataqv_in
	
	output:
	set file("${library}.ataqv.json.gz"), file("${library}.ataqv.out") into ataqv_out

	"""
	${IONICE} ataqv --peak-file $peaks --name ${library} --metrics-file ${library}.ataqv.json.gz --tss-file ${get_tss(get_genome(library))} ${make_excluded_regions_arg(get_genome(library))} --ignore-read-groups ${get_organism(get_genome(library))} $bam > ${library}.ataqv.out
	"""	

}
