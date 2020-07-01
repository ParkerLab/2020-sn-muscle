#!/usr/bin/env nextflow

explain_channel = Channel.fromPath(params.snp_file)
explain_channel_2 = Channel.fromPath(params.snp_file)
model_channel = Channel.fromPath(params.model_glob)
FIMO_BACKGROUND = params.fimo_background
MEME_GLOB = params.meme_glob
PLAIN_MOTIF_GLOB = params.plain_motif_glob

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

process explain {

	container "${params.containers.gkmsvm}"
	publishDir "${params.results}/explained"

	input:
	set file(model), file(fasta) from model_channel.combine(fastas.flatten())

	output:
	file("${cluster}.${ref_or_alt}.explained.txt") into explain_out
	
	script:
	cluster = model.getName().replaceAll('.model.txt', '')
	ref_or_alt = fasta.getName().replaceAll('.fa', '')

	"""
	gkmexplain $fasta $model ${cluster}.${ref_or_alt}.explained.txt
	"""

}

process reformat {
	
	publishDir "${params.results}/explained"
	container "${params.containers.general}"
	executor 'local'

	input:
	file(explained) from explain_out

	output:
	file(out) into reformat_out
	file(out) into reformat_out_2

	script:
	out = explained.getName().replaceAll('.explained.txt', '.explained.reformatted.txt')

	"""
	reformat-importance-scores.py $explained > $out
	"""

}

process plot {

	publishDir "${params.results}/plot"
	executor 'local'
	container "${params.containers.rplot}"

	input:
	file(x) from reformat_out.toSortedList()

	output:
	file("*.pdf")

	"""
	plot-gkmexplain.R ${x.join(' ')}
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


process plot_fimo {

        publishDir "${params.results}/fimo/plot"
        errorStrategy 'ignore'

        input:
        file(fimo) from fimo_out.toSortedList()
        file(gkm) from reformat_out_2.toSortedList()
	file(plain_motifs) from Channel.fromPath(PLAIN_MOTIF_GLOB).toSortedList()
        each snp from explain_channel_2.splitText().map({x -> x.tokenize(' ')[2]})

        output:
        file("*.pdf")

        when:
        snp != 'SNP'

        """
        plot-fimo-general.R . . $snp .
        """

}
