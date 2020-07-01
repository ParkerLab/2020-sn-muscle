#!/usr/bin/env nextflow

explain_channel = Channel.fromPath(params.snp_file)
FIMO_BACKGROUND = params.fimo_background
MEME_GLOB = params.meme_glob

process make_fastas {

	container "${params.containers.mkfasta}"

	input:
	file(snps) from explain_channel
	
	output:
	set file("ref.fa"), file("alt.fa") into fastas
	file("*.fa") into fastas_for_fimo

	"""
	make_ref_alt_flanking_fastas.py $snps ${params.fasta['hg19']} --flank_size 50
	"""

}

process fimo_scan {

	publishDir "${params.results}/fimo/scan"
	maxForks 50
	
	input:
	file(bg) from Channel.fromPath(FIMO_BACKGROUND)	
	each file(fasta) from fastas_for_fimo.flatten()
	each file(motif) from Channel.fromPath(MEME_GLOB)

	output:
        set val("${ref_or_alt}"), file("${ref_or_alt}.${motif_name}.fimo.txt") into concat_in

        script:
        ref_or_alt = fasta.getName().replaceAll('.fa', '')
        motif_name = motif.getName().replaceAll('.meme', '').replaceAll('::', '_')

        """
        fimo --text --bgfile $bg $motif $fasta > ${ref_or_alt}.${motif_name}.fimo.txt
	"""

}

process fimo_concat {

	publishDir "${params.results}/concat"
        executor 'local'

        input:
        set val(ref_or_alt), file(fimo) from concat_in.groupTuple()

        output:
        file("${ref_or_alt}.fimo.txt") into fimo_out

        """
        cat ${fimo.join(' ')} > ${ref_or_alt}.fimo.txt
        """

}
