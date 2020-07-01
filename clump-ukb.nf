#!/usr/bin/env nextflow

TRAIT = params.trait

process prep_ukb {

	input:
	file(x) from Channel.fromPath("/lab/work/porchard/sn-muscle-project/work/reformat-ukbb/results/ldsc/${TRAIT}.txt")

	output:
	set val(trt), file("${trt}.txt") into clump_snps_in
	set val(trt), file("${trt}.txt") into subset_vcf_snps_in
	file("${trt}.txt") into ld_buddy_snps_in

	script:
	trt = x.getName().replaceAll('.txt', '')

	"""
	mv $x ${trt}.tmp
	prep-ukb-clump.py ${trt}.tmp > ${trt}.txt
	"""

}


process get_1000G_samples {

	output:
	file("samples.txt") into samples_out

	"""
	grep -w European /lab/data/genomes/human/hg19/1000GenomesDownloads/igsr_3990samples_allPopulations.tsv | cut -f1 > samples.txt
	"""

}

subset_vcf_in = Channel.fromPath('/lab/data/genomes/human/hg19/1000GenomesDownloads/ALL.*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz').combine(samples_out).combine(subset_vcf_snps_in)
process subset_vcf {

	input:
	set file(vcf), file(samples), val(trt), file(rsids) from subset_vcf_in

	output:
	set val(trt), file("subsetted.recode.vcf.gz") into subset_vcf_out

	"""
	vcftools --gzvcf $vcf --keep $samples --remove-indels --out subsetted-maybe-dups --recode
	# remove duplicated SNPs
	cat subsetted-maybe-dups.recode.vcf | grep -v "^#" | cut -f3 | sort | uniq -d > duplicated-snps.txt
	vcftools --vcf subsetted-maybe-dups.recode.vcf --exclude duplicated-snps.txt --out subsetted --recode
        bgzip subsetted.recode.vcf
        tabix subsetted.recode.vcf.gz
	"""

}

process make_plink_files {

	input:
	set val(trt), file("subsetted.*.vcf.gz") from subset_vcf_out.groupTuple()

	output:
	set val(trt), file("subset.map"), file("subset.ped") into make_plink_out
	set val(trt), file("subset.map"), file("subset.ped") into make_plink_out_2

	"""
	export PERL5LIB="/home/porchard/scripts/vcftools/src/perl/:\$PERL5LIB"
	vcf-concat subsetted.*.vcf.gz | grep -v "^#" | cut -f3 | sort | uniq -d > duplicated-snps.txt
	vcf-concat subsetted.*.vcf.gz | bgzip -c > subsetted.vcf.gz
        vcftools --gzvcf subsetted.vcf.gz --exclude duplicated-snps.txt --plink --out subset
	"""

}

process clump {

	publishDir "${params.results}/clumped"
	memory '50 GB'

	input:
	set val(trt), file(snps), file(mp), file(ped) from clump_snps_in.combine(make_plink_out, by: 0)

	output:
	file("${trt}.clumped")

	"""
	/lab/sw/modules/plink/1.9/bin/plink --file subset --clump $snps --clump-r2 0.8 --clump-p1 5e-8 --clump-p2 0.99999999 --out ${trt}
	"""
}

process ld_with_gws {	
	
	publishDir "${params.results}/ld-08"
	memory '50 GB'

	input:
	file(snps) from ld_buddy_snps_in
	set val(trt), file(mp), file(ped) from make_plink_out_2

	output:
	file("${trt}*")

	"""
	cat $snps | awk '\$2<=5e-8' | cut -f1 | grep -v SNP > snps.txt
	/lab/sw/modules/plink/1.9/bin/plink --file subset --r2 --ld-snp-list snps.txt --ld-window-kb 500 --ld-window 99999 --ld-window-r2 0.8 --out ${trt}
	"""

}
