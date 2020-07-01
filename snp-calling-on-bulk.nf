#!/usr/bin/env nextflow

bams = Channel.fromPath(params.bam_glob)

process sort_bam {

	input:
	file(bam) from bams

	output:
	set val(library), file("${library}.sorted.bam") into mismatch_filter_in

	script:
	library = bam.getName().replaceAll('.pruned.bam', '')	

	"""
	samtools sort -m 3G -O bam -T ${library}.sort -o ${library}.sorted.bam $bam
	"""

}


process mismatch_filter {

	input:
	set val(library), file(bam) from mismatch_filter_in

	output:
	set val(individual), file("${library}.filtered.bam") into mismatch_filter_out

	script:
	m = library =~ /(KSM1|KSM2)_/
	individual = m[0][1]

	"""
	filter-bam-by-edit-distance.py --max-edit-distance 2 $bam ${library}.filtered.bam
	"""

}


process merge {

	input:
	set val(individual), file(bams) from mismatch_filter_out.groupTuple()

	output:
	set val('merged'), file("${individual}.bam"), file("${individual}.bam.bai") into pileup_in
	
	"""
	samtools merge ${individual}.bam ${bams.join(' ')}
	samtools index ${individual}.bam
	"""

}



process pileup_aggregate {

	time '48h'
	memory '100 GB'

	input:
	set val(throwaway), file(bams), file(indices) from pileup_in.groupTuple()
	
	output:
	file("pileup.txt") into pileup_out

	"""
	samtools mpileup -R -t DP,DV,DPR -Q 20 -f ${params.fasta} -I -d 10000 -E -u ${bams.join(' ')} > pileup.txt
	"""

}


process call_on_aggregate {

	publishDir "${params.results}/calls"
	time '48h'
	memory '100 GB'

	input:
	file(pileup) from pileup_out

	output:
	file("calls.vcf.gz") into calls_out
	
	"""
	bcftools call -cv -Oz -f GQ $pileup > calls.vcf.gz
	"""

}

process filter_to_quality_genotypes {

	publishDir "${params.results}/genotypes"
	time '48h'

	input:
	file(calls) from calls_out

	output:
	file('snp-calls.vcf.gz')

	"""
	bcftools filter -Oz -i '%MIN(GQ)>90' $calls > snp-calls.vcf.gz
	"""

}
