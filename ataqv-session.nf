#!/usr/bin/env nextflow

process ataqv_ignore_readgroups {

	memory '10 GB'
	clusterOptions '--constraint=wolverine'
	publishDir "${params.results}/ataqv"

	input:
	file(

}
