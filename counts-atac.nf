#!/usr/bin/env nextflow

ROOT = params.projroot
BAM_CHAN = Channel.fromPath(params.bam_glob)

process uncorrected_counts {

        memory '40 GB'
	container "${params.containers.snATAC}"
        maxRetries 1
        time '24h'
        errorStrategy 'retry'
        cache false

        publishDir "${params.results}/uncorrected-counts", mode: 'rellink'

        input:
	file(bam) from BAM_CHAN

        output:
        set val(library), val(genome), file("${library}-${genome}.counts.txt")

	script:
	m = bam.getName() =~ /(.*)-(.*)\.pruned\.bam/
	library = m[0][1]
	genome = m[0][2]

        """
	cat ${ROOT}/data/liger-features/${genome}.features.noblacklist | sort -k1V,1 -k2n,2 > genes.sorted.bed
        sort-bed-by-bam.py genes.sorted.bed $bam > genes.sorted_by_bam.bed
        bedtools intersect -wa -wb -bed -sorted -a $bam -b genes.sorted_by_bam.bed | cut -f4,16 | perl -pe 's@.*_(.*)/\\d+\\t(.*)@\$1\\t\$2@' > overlaps.txt 
	atac-overlaps-to-feature-counts.py overlaps.txt ${library}-${genome} > ${library}-${genome}.counts.txt
        """

}
