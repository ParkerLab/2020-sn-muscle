###### CHANGE THESE PATHS
ROOT=/lab/work/porchard/sn-muscle-project
HG19_FASTA=/lab/data/reference/human/hg19/hg19.fa
HG19_CHROM_SIZES=/lab/data/reference/human/hg19/hg19.chrom_sizes
HG19_BLACKLISTS=/lab/data/reference/human/hg19/annot/wgEncodeDacMapabilityConsensusExcludable.bed.gz /lab/data/reference/human/hg19/annot/wgEncodeDukeMapabilityRegionsExcludable.bed.gz
RN6_FASTA=/lab/work/porchard/data/fasta/rn6.fa
RN6_CHROM_SIZES=/lab/work/porchard/data/chrom_sizes/rn6.chrom_sizes
RN6_BLACKLISTS=/lab/work/porchard/data/mappability.new/rn6.blacklist.1.bed.gz
HG38_TO_HG19_CHAIN=/lab/work/porchard/data/chain/hg38ToHg19.over.chain.gz

##### THERE SHOULD BE NO NEED TO TOUCH ANYTHING FROM HERE DOWN

WORK=$(ROOT)/work
DATA=$(ROOT)/data
FIGURES=$(ROOT)/figures
SRC=$(ROOT)/src
BIN=$(ROOT)/bin
CONTROL=$(ROOT)/control
SAMPLE_INFO=$(ROOT)/sample_info/sample_info.txt
NAME = sn-muscle-project
CLUSTER_NAMES=$(ROOT)/2019-01-04-cluster-names.txt
SHARED_OPEN_CHROMATIN=$(WORK)/meuleman/results/common/common_open_chromatin.bed
PATH_TO_BULK=$(WORK)/bulk-atacseq/results

.PHONY: data sample_info all

ANALYSIS = $(WORK)/$@
CONFIG = $(ANALYSIS)/config.json
PIPELINE = $(ANALYSIS)/pipeline
JOBNAME = $(NAME)-$@
DOWNSTREAM=$(WORK)/downstream-new-features

##### https://stackoverflow.com/questions/7039811/how-do-i-process-extremely-long-lists-of-files-in-a-make-recipe
define NL


endef
#####

### SET UP DATA FILES ###
data: repeats liger-features dbSNP-vcf orthologues roadmap-posteriors ukb-summary-stats chromhmm other-annotations reformat-ukbb ldsc-baseline diamante-summary-stats manning-summary-stats gencode-coding meuleman-data meuleman encode-cisregulatory-elements encode-cres 1000G-SNPs

repeats: ANALYSIS=$(DATA)/$@
repeats:
	mkdir -p $(ANALYSIS)
	mysql --host=genome-mysql.cse.ucsc.edu --user=genome -D hg19 -e "SELECT * FROM simpleRepeat" > $(ANALYSIS)/trf_table.txt
	cut -f2-4 $(ANALYSIS)/trf_table.txt | grep -v chromStart | sort -k1,1 -k2n,2 | bedtools merge -i stdin > $(ANALYSIS)/trf.bed


meuleman-data: ANALYSIS=$(DATA)/meuleman
meuleman-data:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && wget https://www.meuleman.org/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz
	cd $(ANALYSIS) && wget https://www.meuleman.org/DHS_Index_and_Vocabulary_metadata.tsv

meuleman:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --dhs_index $(DATA)/meuleman/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz --chain $(HG38_TO_HG19_CHAIN) --metadata $(DATA)/meuleman/DHS_Index_and_Vocabulary_metadata.tsv --dhs_matrix $(DATA)/meuleman/dat_bin_FDR01_hg38.RData $(ROOT)/prep-meuleman.nf &

encode-cisregulatory-elements: ANALYSIS=$(DATA)/$@
encode-cisregulatory-elements:
	cd $(ANALYSIS) && bash commands

encode-cres:
	mkdir -p $(ANALYSIS)
	echo "#!/bin/bash" > $(PIPELINE)
	echo "#SBATCH --mem=10G" >> $(PIPELINE)
	cd $(ANALYSIS) && echo "python $(ROOT)/bin/sort_uniq_gzip.py $(DATA)/encode-cisregulatory-elements/*.bed.gz | sort -k1,1 -k2n,2 > encode-cres.bed" >> $(PIPELINE) && sbatch $(PIPELINE)

1000G-SNPs: ANALYSIS=$(DATA)/$@
1000G-SNPs:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/*.vcf.gz
	cd $(ANALYSIS) && wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/*.vcf.gz.tbi
	cd $(ANALYSIS) && rm *chrMT* *chrX* *chrY* *wgs*
	bcftools concat $(ANALYSIS)/ALL*.vcf.gz | bcftools view -Oz -G --types snps -o $(ANALYSIS)/1000G-snps.vcf.gz

gencode-coding: ANALYSIS=$(DATA)/$@
gencode-coding:
	mkdir $(ANALYSIS) && cd $(ANALYSIS) && wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz && zcat gencode.v19.annotation.gtf.gz | awk '$$3=="CDS"' | cut -f1,4,5 | sort -k1,1 -k2n,2 -k3n,3 | bedtools merge -i stdin > coding.bed

liger-features: ANALYSIS=$(DATA)/$@
liger-features:
	# For each genome:
	# get the GTF
	# get the annotation report
	# convert the chromosome names to ucsc style
	# then, using python script:
	# keep curated transcripts
	# for each gene, merge all the transcripts, and keep 3 kb upstream (even if it overlaps with other stuff)
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Rattus_norvegicus/annotation_releases/current/GCF_000001895.5_Rnor_6.0/GCF_000001895.5_Rnor_6.0_assembly_report.txt && wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Rattus_norvegicus/annotation_releases/current/GCF_000001895.5_Rnor_6.0/GCF_000001895.5_Rnor_6.0_genomic.gtf.gz
	cd $(ANALYSIS) && wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/105.20190906/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz && wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/105.20190906/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt
	python $(ROOT)/bin/refseq-annotation-and-report-to-feature-file.py $(ANALYSIS)/GCF_000001895.5_Rnor_6.0_genomic.gtf.gz $(ANALYSIS)/GCF_000001895.5_Rnor_6.0_assembly_report.txt | grep -v -e 'random' -e 'chrUn' | sort -k1,1 -k2n,2 -k3n,3 > $(ANALYSIS)/rn6.genes.bed
	python $(ROOT)/bin/make-liger-features.py $(ANALYSIS)/rn6.genes.bed $(RN6_CHROM_SIZES) | sort -k1,1 -k2n,2 -k3n,3 > $(ANALYSIS)/rn6.features.bed
	python $(ROOT)/bin/refseq-annotation-and-report-to-feature-file.py $(ANALYSIS)/GCF_000001405.25_GRCh37.p13_genomic.gtf.gz $(ANALYSIS)/GCF_000001405.25_GRCh37.p13_assembly_report.txt | grep -v -e 'random' -e 'chrUn' | sort -k1,1 -k2n,2 -k3n,3 > $(ANALYSIS)/hg19.genes.bed
	python $(ROOT)/bin/make-liger-features.py $(ANALYSIS)/hg19.genes.bed $(HG19_CHROM_SIZES) | sort -k1,1 -k2n,2 -k3n,3 > $(ANALYSIS)/hg19.features.bed
	# blacklist filter
	bedtools intersect -a $(ANALYSIS)/rn6.features.bed -b $(RN6_BLACKLISTS) -v > $(ANALYSIS)/rn6.features.noblacklist 
	bedtools intersect -a $(ANALYSIS)/hg19.features.bed -b $(HG19_BLACKLISTS) -v > $(ANALYSIS)/hg19.features.noblacklist 

dbSNP-vcf: ANALYSIS=$(DATA)/$@
dbSNP-vcf:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz
	cd $(ANALYSIS) && wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz.tbi

diamante-summary-stats: ANALYSIS=$(DATA)/$@
diamante-summary-stats:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup python $(ROOT)/bin/prep-diamante-for-ldsc.py --vcf $(DATA)/dbSNP-vcf/All_20170710.vcf.gz $(ANALYSIS)/Mahajan.NatGenet2018b.T2D.European.txt diamante.T2D.European.ldsc.txt.gz &
	cd $(ANALYSIS) && nohup python $(ROOT)/bin/prep-diamante-for-ldsc.py --vcf $(DATA)/dbSNP-vcf/All_20170710.vcf.gz $(ANALYSIS)/Mahajan.NatGenet2018b.T2Dbmiadj.European.txt diamante.T2Dbmiadj.European.ldsc.txt.gz &

manning-summary-stats: ANALYSIS=$(DATA)/$@
manning-summary-stats:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && wget ftp://ftp.sanger.ac.uk/pub/magic/MAGIC_Manning_et_al_lnFastingInsulin_MainEffect.txt.gz && wget ftp://ftp.sanger.ac.uk/pub/magic/MAGIC_Manning_et_al_FastingGlucose_MainEffect.txt.gz
	cd $(ANALYSIS) && Rscript $(ROOT)/bin/reformat-manning-summarystats-for-ldsc.R

ukb-summary-stats: ANALYSIS = $(DATA)/$@
ukb-summary-stats:
	cd $(ANALYSIS) && wget -O ukb31063_h2_all.02Oct2019.tsv.gz https://www.dropbox.com/s/ipeqyhrpdqav5uh/ukb31063_h2_all.02Oct2019.tsv.gz?dl=1
	cd $(ANALYSIS) && Rscript $(ROOT)/bin/choose-ukb-traits-new.R ukb31063_h2_all.02Oct2019.tsv.gz manifest.tsv && echo "# drmr:job" > wget-commands.sh && cut -f4 manifest-subset.tsv | perl -pe 's/.bgz$$/.gz/' >> wget-commands.sh && drmrarray -s 20 wget-commands.sh

ldsc-data: ANALYSIS=$(DATA)/$@
ldsc-data:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baselineLD_v2.2_ldscores.tgz
	cd $(ANALYSIS) && wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_frq.tgz && tar -xvzf 1000G_Phase3_frq.tgz
	cd $(ANALYSIS) && wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz && tar -xvzf 1000G_Phase3_plinkfiles.tgz
	cd $(ANALYSIS) && wget https://data.broadinstitute.org/alkesgroup/LDSCORE/hapmap3_snps.tgz && tar -xvzf hapmap3_snps.tgv
	cd $(ANALYSIS) && wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_weights_hm3_no_MHC.tgz && tar -xvzf 1000G_Phase3_weights_hm3_no_MHC.tgz
	cd $(ANALYSIS) && wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2 && bzip2 -d w_hm3.snplist.bz2
	cd $(ANALYSIS) && wget https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz && tar -xvzf weights_hm3_no_hla.tgz
	cd $(ROOT)/bin && python make-snp-id-conversions.py

ldsc-baseline-min-maf:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume -with-trace -with-report --results $(ANALYSIS)/results --projroot $(ROOT) $(ROOT)/ldsc-update-baseline-model-min-maf.nf &

chromhmm:
	mkdir -p $(DATA)/$@
	cd $(DATA)/$@ && wget https://zenodo.org/record/3524578/files/islet-cage-zenodo.tar.gz && tar -xvzf islet-cage-zenodo.tar.gz
	cd $(DATA)/$@ && cp data/chromhmm/selected_annotated_states/files_by_state/* .
	cd $(DATA)/$@ && rename 's/cell4_11.//' *
	$(foreach t,Adipose SkeletalMuscle Liver Islets,cd $(DATA)/$@ && python $(BIN)/make-chromhmm-dense.py $(foreach s,Active_enhancer Active_TSS Bivalent_poised_TSS Flanking_TSS Genic_enhancer Quiescent_low_signal Repressed_polycomb Strong_transcription Weak_enhancer Weak_transcription Weak_TSS,$(t).$(s).bed) | sort -k1,1 -k2n,2 > $(t).dense.bed$(NL))
	cd $(DATA)/$@ && rm -r data

other-annotations:
	mkdir -p $(DATA)/$@
	ln -s /lab/work/vivekrai/2019_mohlke_adipose/work/macs2/Adipose*.broadPeak $(DATA)/$@/
	python $(ROOT)/bin/peak-sharing.py $(DATA)/$@/Adipose* | awk '$$4>=2' | cut -f1-3 | bedtools intersect -a stdin -b /lab/work/porchard/sn-muscle-project/data/mappability/hg19* -v > adipose.bed
	cat /home/vivekrai/analyses/2018_NIH_Islets.snatacseq.v2/work/2019-03-01_clustering-final/peaks/1_peaks.broadPeak.noblacklist | cut -f1-3 | sort -k1,1 -k2n,2 > $(DATA)/$@/beta_ATAC.bed
	cp /lab/work/vivekrai/2017_NIH_Islets.atacseq/work/2019-05-09_process-samples-subset/macs2/*.noblacklist $(DATA)/$@
	cd $(DATA)/$@ && rm EndoC*

roadmap-posteriors:
	cd $(DATA)/$@ && nohup bash $(ROOT)/src/download-roadmap-posteriors.sh &

orthologues:
	mkdir -p $(DATA)/$@
	cd $(DATA)/$@ && Rscript $(SRC)/make-orthologues.R

reformat-ukbb:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --projroot $(ROOT) -with-singularity $(ROOT)/singularity/general/general.simg --results $(ANALYSIS) $(ROOT)/reformat-ukbb.nf &

#### ANALYSES
snp-calling-on-bulk:
	mkdir -p $(ANALYSIS)/data
	ln -s $(PATH_TO_BULK)/prune/320-NM-1.pruned.bam $(ANALYSIS)/data/KSM1_rep1.bam
	ln -s $(PATH_TO_BULK)/prune/320-NM-2.pruned.bam $(ANALYSIS)/data/KSM1_rep2.bam
	ln -s $(PATH_TO_BULK)/prune/320-NM-3.pruned.bam $(ANALYSIS)/data/KSM2_rep1.bam
	ln -s $(PATH_TO_BULK)/prune/320-NM-4.pruned.bam $(ANALYSIS)/data/KSM2_rep2.bam
	cd $(ANALYSIS) && nohup nextflow run -resume -with-trace -with-report -with-singularity $(ROOT)/singularity/general/general.simg --results $(ANALYSIS)/results --bam_glob '$(ANALYSIS)/data/*' --fasta $(HG19_FASTA) $(ROOT)/snp-calling-on-bulk.nf &

diamante-credible-set-ppa-track:
	mkdir -p $(ANALYSIS)
	cat $(DATA)/diamante-credible-sets/genetic_credible_sets/* | cut -f2,3,4 | grep -v Pos | awk '{print($$1, $$2-1, $$2, $$3)}' | perl -pe 's/ /\t/g; s/^/chr/' | sort -k1,1 -k2n,2 | bedtools groupby -g 1,2,3 -c 4 -o max >> $(ANALYSIS)/diamante-ppa.bdg
	cat $(DATA)/diamante-credible-sets/genetic_credible_sets/g_credible_set_Eur_ARL15_5_53271420.txt | cut -f2,3,4 | grep -v Pos | awk '{print($$1, $$2-1, $$2, $$3)}' | perl -pe 's/ /\t/g; s/^/chr/' | sort -k1,1 -k2n,2 >> $(ANALYSIS)/diamante-ppa-our-loci.bdg
	cat $(DATA)/diamante-credible-sets/genetic_credible_sets/g_credible_set_Eur_ITPR2_12_26453283.txt | cut -f2,3,4 | grep -v Pos | awk '{print($$1, $$2-1, $$2, $$3)}' | perl -pe 's/ /\t/g; s/^/chr/' | sort -k1,1 -k2n,2 >> $(ANALYSIS)/diamante-ppa-our-loci.bdg
	cat $(ANALYSIS)/diamante-ppa-our-loci.bdg | sort -k1,1 -k2n,2 >> $(ANALYSIS)/diamante-ppa-our-loci.sorted
	mv $(ANALYSIS)/diamante-ppa-our-loci.sorted $(ANALYSIS)/diamante-ppa-our-loci.bdg

rnaseq:
	mkdir -p $(ANALYSIS)
	python $(CONTROL)/rnaseq/make_config.py $(ROOT) $(ROOT)/sample_info/sample_info.txt > $(CONFIG)
	cd $(ANALYSIS) && nohup nextflow run -resume --chemistry V3 -with-trace -params-file $(CONFIG) -with-singularity /lab/work/porchard/singularity/archive/snRNA/2019-11-15/snRNA.simg -with-report -qs 300 --results $(ANALYSIS)/results /home/porchard/github/snRNAseq-NextFlow/main.nf &

bulk-atacseq:
	mkdir -p $(ANALYSIS)
	python $(CONTROL)/$@/make_config.py $(ROOT) > $(ANALYSIS)/config.json
	cd $(ANALYSIS) && nohup nextflow run -resume -with-singularity /lab/work/porchard/singularity/archive/ATAC/2019-11-20/ATAC.simg -with-dag -with-timeline -with-trace -with-report -params-file $(ANALYSIS)/config.json --results $(ANALYSIS)/results /home/porchard/github/ATACseq-NextFlow/main.nf &

atacseq:
	mkdir -p $(ANALYSIS)
	python $(CONTROL)/$@/make_config.py $(ROOT) $(ROOT)/sample_info/sample_info.txt > $(CONFIG)
	cp $(ROOT)/atac.nextflow.config $(ANALYSIS)/nextflow.config
	cd $(ANALYSIS) && nohup nextflow run -resume -with-trace -params-file $(CONFIG) -with-singularity /lab/work/porchard/singularity/archive/snATAC/2019-11-07/snATAC.simg --low_read_count_threshold 1000 -with-report -qs 300 --results $(ANALYSIS)/results /home/porchard/github/snATACseq-NextFlow/main.nf &

process-as-bulk:
	mkdir -p $(ANALYSIS)/data/bams
	cp $(WORK)/atacseq/results/merge/133* $(ANALYSIS)/data/bams/
	cp $(WORK)/atacseq/results/merge/63_* $(ANALYSIS)/data/bams/
	cp $(ROOT)/atac.nextflow.config $(ANALYSIS)/nextflow.config
	cd $(ANALYSIS) && nohup nextflow run -resume --mapped_bam_glob '$(ANALYSIS)/data/bams/*.bam' -with-singularity $(ROOT)/singularity/ATAC.simg --results $(ANALYSIS)/results $(ROOT)/process-snATAC-as-bulk.nf &

counts-rna: ANALYSIS=$(WORK)/counts
counts-rna:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --projroot $(ROOT) $(ROOT)/counts.nf &

counts-atac:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --projroot $(ROOT) --bam_glob '$(ROOT)/work/atacseq/results/prune/*.bam' $(ROOT)/counts-atac.nf &

downstream-new-features:
	mkdir -p $(ANALYSIS)
	python $(ROOT)/control/$@/make_config.py $(ROOT) > $(ANALYSIS)/config.json
	cd $(ANALYSIS) && nohup nextflow run -resume -dump-hashes --results $(ANALYSIS)/results -with-dag -params-file config.json $(ROOT)/downstream-cellbender-corrected.nf &

atac-correlation-fans: LIBRARIES=$(shell seq 133151 133154)
atac-correlation-fans: BULK_LIBRARIES=320-NM-1 320-NM-2
atac-correlation-fans:
	mkdir -p $(ANALYSIS)/data/peaks
	mkdir -p $(ANALYSIS)/data/bams
	$(foreach l,$(LIBRARIES),ln -sf $(WORK)/process-as-bulk/results/prune/$(l)-hg19.pruned.bam $(ANALYSIS)/data/bams/$(NL))
	$(foreach l,$(LIBRARIES),ln -sf $(WORK)/process-as-bulk/results/macs2/$(l)-hg19_peaks.broadPeak.noblacklist $(ANALYSIS)/data/peaks/$(NL))
	$(foreach l,$(BULK_LIBRARIES),ln -sf $(PATH_TO_BULK)/macs2/$(l)_peaks.broadPeak.noblacklist $(ANALYSIS)/data/peaks/$(NL))
	$(foreach l,$(BULK_LIBRARIES),ln -sf $(PATH_TO_BULK)/prune/$(l).pruned.bam $(ANALYSIS)/data/bams/$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --bam_glob '$(ANALYSIS)/data/bams/*' --peak_glob '$(ANALYSIS)/data/peaks/*' $(ROOT)/atac-correlation.nf &

atac-correlation-loading: LIBRARIES=63_20 63_40
atac-correlation-loading: BULK_LIBRARIES=320-NM-1 320-NM-2 320-NM-3 320-NM-4
atac-correlation-loading:
	mkdir -p $(ANALYSIS)/data/peaks
	mkdir -p $(ANALYSIS)/data/bams
	$(foreach l,$(LIBRARIES),ln -sf $(WORK)/process-as-bulk/results/prune/$(l)-hg19.pruned.bam $(ANALYSIS)/data/bams/$(NL))
	$(foreach l,$(LIBRARIES),ln -sf $(WORK)/process-as-bulk/results/macs2/$(l)-hg19_peaks.broadPeak.noblacklist $(ANALYSIS)/data/peaks/$(NL))
	$(foreach l,$(BULK_LIBRARIES),ln -sf $(PATH_TO_BULK)/macs2/$(l)_peaks.broadPeak.noblacklist $(ANALYSIS)/data/peaks/$(NL))
	$(foreach l,$(BULK_LIBRARIES),ln -sf $(PATH_TO_BULK)/prune/$(l).pruned.bam $(ANALYSIS)/data/bams/$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --bam_glob '$(ANALYSIS)/data/bams/*' --peak_glob '$(ANALYSIS)/data/peaks/*' $(ROOT)/atac-correlation.nf &

logistic-regression-multiple-states-hg19:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume -with-trace -with-report --results $(ANALYSIS)/results --downstream $(DOWNSTREAM) --genome hg19 --projroot $(ROOT) $(ROOT)/logistic-regression-multiple-states.nf &

logistic-regression-multiple-states-rn6:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume -with-trace -with-report --results $(ANALYSIS)/results --downstream $(DOWNSTREAM) --genome rn6 --projroot $(ROOT) $(ROOT)/logistic-regression-multiple-states.nf &

fimo:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --meme_glob '/lab/work/porchard/fimo/work/ENCODE2013/data/*.meme' --fimo_background $(DATA)/motifs/hg19.background.txt --snp_file $(ROOT)/snps-for-manuscript.txt --results $(ANALYSIS)/results $(ROOT)/fimo.nf &

delta-svm-hg19: NUMBER_PEAKS_KEEP=40000
delta-svm-hg19: FASTA=$(HG19_FASTA)
delta-svm-hg19: CHROM_SIZES=$(HG19_CHROM_SIZES)
delta-svm-hg19: GENOME=hg19
delta-svm-hg19: BLACKLISTS=$(HG19_BLACKLISTS)
delta-svm-hg19: SNP_FILE=/lab/work/porchard/sn-muscle-project/work/1000G-SNP-vcf/1000G-snps.vcf.gz
delta-svm-hg19:
	mkdir -p $(ANALYSIS)/data
	$(foreach c,$(shell seq 0 6),cat $(DOWNSTREAM)/results/process-by-cluster-round-1/narrowpeaks/$(c)-$(GENOME)_peaks.narrowPeak.noblacklist | sort -k9n,9 | tail -n $(NUMBER_PEAKS_KEEP) > $(ANALYSIS)/data/cluster_$(c).narrow.bed$(NL))
	$(foreach c,$(shell seq 0 6),cp $(DOWNSTREAM)/results/process-by-cluster-round-1/peaks/$(c)-$(GENOME)_peaks.broadPeak.noblacklist $(ANALYSIS)/data/cluster_$(c).broad.bed$(NL))
	cd $(ANALYSIS) && cp $(ROOT)/deltasvm-config.nf nextflow.config && NXF_VER=19.10.0 nohup ~/github/snATAC_deltaSVM/main.nf -resume --repeat_bed $(DATA)/repeats/trf.bed --max_repeat_content 0.6 --results_dir $(ANALYSIS)/results --kernel 2 --peak_dir '$(ANALYSIS)/data/cluster_*.narrow.bed' --broadpeak_dir '$(ANALYSIS)/data/cluster_*.broad.bed' --exclude '$(BLACKLISTS)' --ref $(FASTA) --human_ref $(HG19_FASTA) --snp_file $(SNP_FILE) --chrom_sizes $(CHROM_SIZES) --workDir $(ANALYSIS)/work &

get-deltasvm-z-scores:
	mkdir -p $(ANALYSIS)
	cat snps-for-manuscript.txt | perl -pe 's/ /\t/g' | cut -f1,2 | grep -v BP | perl -pe 's/^/chr/' > $(ANALYSIS)/pos-file.txt
	$(foreach c,$(shell seq 0 6),python $(ROOT)/bin/get-deltaSVM-z-scores.py --deltaSVM-file $(WORK)/delta-svm-hg19/results/cluster_$(c)/snp_scores_cluster_$(c).txt --pos-file $(ANALYSIS)/pos-file.txt > $(ANALYSIS)/cluster_$(c).z-scores.txt$(NL))

gkmexplain-and-fimo:
	mkdir -p $(ANALYSIS)/data
	$(foreach c,$(shell seq 0 6),cp $(WORK)/delta-svm-hg19/results/cluster_$(c)/model.model.txt $(ANALYSIS)/data/cluster_$(c).model.txt$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --meme_glob '$(DATA)/motifs/*.meme' --fimo_background $(DATA)/motifs/hg19.background.txt --plain_motif_glob '$(DATA)/motifs/*.txt' --snp_file $(ROOT)/snps-for-manuscript.txt --results $(ANALYSIS)/results --model_glob '$(ANALYSIS)/data/*' $(ROOT)/gkmexplain-and-fimo.nf &

creatinine-snp-gb-tracks:
	mkdir -p $(ANALYSIS)
	cat $(WORK)/clump-creatine-ukb/results/ld-08/30700_irnt.ld | perl -pe 's/^\s+//; s/\s+/\t/g; s/$$/\n/' | grep -w rs227727 | cut -f1-7 > $(ANALYSIS)/ld.txt
	python $(ROOT)/bin/make-locus-bed-for-index-snp.py $(ANALYSIS)/ld.txt rs227727 > $(ANALYSIS)/locus.bed
	cat nog-snp.bed | cut -f1,3 > $(ANALYSIS)/index-snp.txt
	python $(ROOT)/bin/make-bb-from-loci-color-lead.py $(ANALYSIS)/locus.bed $(ANALYSIS)/index-snp.txt > $(ANALYSIS)/nog-locus-lead-plus-ld-buddies.bed
	Rscript $(ROOT)/bin/make-pics-bdg.R $(ANALYSIS)/ld.txt $(ROOT)/nog-pics.csv $(ANALYSIS)/nog-locus-ppa.bdg

lift-human-to-rat:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --peak_glob '$(WORK)/downstream-new-features/results/process-by-cluster-round-1/peaks/*-hg19_peaks.broadPeak.noblacklist' $(ROOT)/find-rat-peak.nf &

project-human-snps-to-rat:
	mkdir -p $(ANALYSIS)
	cat $(DATA)/diamante-credible-sets/genetic_credible_sets/* | grep -v Pos | cut -f2,3 | awk '{print($$1, $$2-20, $$2+20)}' | perl -pe 's/^/chr/; s/ /\t/g' > $(ANALYSIS)/snps-and-flanking.bed
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --bed_glob '$(ANALYSIS)/snps-and-flanking.bed' $(ROOT)/find-rat-snp.nf &

### GWAS ENRICHMENTS ###
ldsc-T2D-one-model-per-cluster-hg19: ANALYSIS=$(WORK)/ldsc/T2D/hg19/one-model-per-cluster
ldsc-T2D-one-model-per-cluster-hg19:
	mkdir -p $(ANALYSIS)/data
	$(foreach c,$(shell seq 0 6),cat $(DOWNSTREAM)/results/process-by-cluster-round-1/peaks/$(c)-hg19_peaks.broadPeak.noblacklist | sort -k1,1 -k2n,2 | cut -f1-3 > $(ANALYSIS)/data/cluster_$(c).bed$(NL))
	cp $(DATA)/other-annotations/beta_ATAC.bed $(ANALYSIS)/data/beta_cell.bed
	cp $(DATA)/other-annotations/adipose.bed $(ANALYSIS)/data/adipose.bed
	cp $(WORK)/meuleman/results/tissues/Liver.bed $(ANALYSIS)/data/liver.bed
	cd $(ANALYSIS) && nohup nextflow run -resume -with-trace -with-report --other_annotations '$(SHARED_OPEN_CHROMATIN)' --projroot $(ROOT) --results $(ANALYSIS)/results --cluster_names $(CLUSTER_NAMES) $(ROOT)/ldsc-T2D-one-model-per-cell-type.nf --peak_glob '$(ANALYSIS)/data/*.bed' &

ldsc-T2D-joint-model-hg19: ANALYSIS=$(WORK)/ldsc/T2D/hg19/joint-model
ldsc-T2D-joint-model-hg19:
	mkdir -p $(ANALYSIS)/data
	$(foreach c,$(shell seq 0 6),cat $(DOWNSTREAM)/results/process-by-cluster-round-1/peaks/$(c)-hg19_peaks.broadPeak.noblacklist | sort -k1,1 -k2n,2 | cut -f1-3 > $(ANALYSIS)/data/cluster_$(c).bed$(NL))
	cp $(DATA)/other-annotations/beta_ATAC.bed $(ANALYSIS)/data/beta_cell.bed
	cp $(DATA)/other-annotations/adipose.bed $(ANALYSIS)/data/adipose.bed
	cp $(WORK)/meuleman/results/tissues/Liver.bed $(ANALYSIS)/data/liver.bed
	cd $(ANALYSIS) && nohup nextflow run -resume -with-trace -with-report --other_annotations '$(SHARED_OPEN_CHROMATIN)' --projroot $(ROOT) --results $(ANALYSIS)/results --cluster_names $(CLUSTER_NAMES) $(ROOT)/ldsc-T2D-joint-model.nf --peak_glob '$(ANALYSIS)/data/*.bed' &

ldsc-T2D-one-model-per-cluster-rn6: ANALYSIS=$(WORK)/ldsc/T2D/rn6/one-model-per-cluster
ldsc-T2D-one-model-per-cluster-rn6:
	mkdir -p $(ANALYSIS)/data
	$(foreach c,$(shell seq 0 6),cat $(DOWNSTREAM)/results/rat-peak-liftover/$(c)-rn6_peaks_in_hg19.bed | sort -k1,1 -k2n,2 > $(ANALYSIS)/data/cluster_$(c).bed$(NL))
	cp $(DATA)/other-annotations/beta_ATAC.bed $(ANALYSIS)/data/beta_cell.bed
	cp $(DATA)/other-annotations/adipose.bed $(ANALYSIS)/data/adipose.bed
	cp $(WORK)/meuleman/results/tissues/Liver.bed $(ANALYSIS)/data/liver.bed
	cd $(ANALYSIS) && nohup nextflow run -resume -with-trace -with-report --other_annotations '$(SHARED_OPEN_CHROMATIN)' --projroot $(ROOT) --results $(ANALYSIS)/results --cluster_names $(CLUSTER_NAMES) $(ROOT)/ldsc-T2D-one-model-per-cell-type.nf --peak_glob '$(ANALYSIS)/data/*.bed' &

ldsc-T2D-joint-model-rn6: ANALYSIS=$(WORK)/ldsc/T2D/rn6/joint-model
ldsc-T2D-joint-model-rn6:
	mkdir -p $(ANALYSIS)/data
	$(foreach c,$(shell seq 0 6),cat $(DOWNSTREAM)/results/rat-peak-liftover/$(c)-rn6_peaks_in_hg19.bed | sort -k1,1 -k2n,2 > $(ANALYSIS)/data/cluster_$(c).bed$(NL))
	cp $(DATA)/other-annotations/beta_ATAC.bed $(ANALYSIS)/data/beta_cell.bed
	cp $(DATA)/other-annotations/adipose.bed $(ANALYSIS)/data/adipose.bed
	cp $(WORK)/meuleman/results/tissues/Liver.bed $(ANALYSIS)/data/liver.bed
	cd $(ANALYSIS) && nohup nextflow run -resume -with-trace -with-report --other_annotations '$(SHARED_OPEN_CHROMATIN)' --projroot $(ROOT) --results $(ANALYSIS)/results --cluster_names $(CLUSTER_NAMES) $(ROOT)/ldsc-T2D-joint-model.nf --peak_glob '$(ANALYSIS)/data/*.bed' &

# For UKB:
ldsc-ukb-hg19: OTHER_ANNOTATIONS = $(shell ls $(WORK)/meuleman/results/tissues/*.bed) $(shell ls $(DATA)/other-annotations/*.bed)
ldsc-ukb-hg19: ANALYSIS=$(WORK)/ldsc/UKB/hg19/joint-model
ldsc-ukb-hg19:
	mkdir -p $(ANALYSIS)/data
	$(foreach f,$(OTHER_ANNOTATIONS),#cp $(f) $(ANALYSIS)/data/$(NL))
	$(foreach c,$(shell seq 0 6),#cp $(DOWNSTREAM)/results/process-by-cluster-round-1/peaks/$(c)-hg19_peaks.broadPeak.noblacklist $(ANALYSIS)/data/cluster_$(c).bed$(NL))
	#cd $(ANALYSIS)/data && rm -rf Amion.bed Esophagus.bed Periodontal_Ligament.bed
	cd $(ANALYSIS) && nohup nextflow run -resume -with-trace -with-report --projroot $(ROOT) --results $(ANALYSIS)/results $(ROOT)/ldsc-UKB-new-baseline.nf --other_annotations '$(SHARED_OPEN_CHROMATIN)' --trait_list $(ROOT)/ukb-traits.txt --peak_glob '$(ANALYSIS)/data/*.bed' &

ldsc-ukb-rn6: OTHER_ANNOTATIONS = $(shell ls $(WORK)/meuleman/results/tissues/*.bed) $(shell ls $(DATA)/other-annotations/*.bed)
ldsc-ukb-rn6: ANALYSIS=$(WORK)/ldsc/UKB/rn6/joint-model
ldsc-ukb-rn6:
	mkdir -p $(ANALYSIS)/data
	$(foreach f,$(OTHER_ANNOTATIONS),cp $(f) $(ANALYSIS)/data/$(NL))
	$(foreach c,$(shell seq 0 6),cp $(DOWNSTREAM)/results/rat-peak-liftover/$(c)-rn6_peaks_in_hg19.bed $(ANALYSIS)/data/cluster_$(c).bed$(NL))
	cd $(ANALYSIS)/data && rm -rf Amion.bed Esophagus.bed Periodontal_Ligament.bed
	cd $(ANALYSIS) && nohup nextflow run -resume -with-trace -with-report --projroot $(ROOT) --results $(ANALYSIS)/results $(ROOT)/ldsc-UKB-new-baseline.nf --other_annotations '$(SHARED_OPEN_CHROMATIN)' --trait_list $(ROOT)/ukb-traits.txt --peak_glob '$(ANALYSIS)/data/*.bed' &

clump-creatine-ukb:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --trait 30700_irnt --results $(ANALYSIS)/results --projroot $(ROOT) $(ROOT)/clump-ukb.nf &

### Making figures
ldsc-annotation-table:
	zcat work/ldsc-baseline/results/baseline-model/baselineLD.10.annot.gz | head -1 | perl -pe 's/\t/\n/g' | awk 'NR>4' > ldsc-annotations.csv

manuscript-figures:
	mkdir -p $(ANALYSIS)/bin
	cp $(ROOT)/manuscript-figures.nf $(ANALYSIS)/
	cp $(ROOT)/nextflow.config $(ANALYSIS)/
	cp -uL $(ROOT)/manuscript-figure-plots/* $(ANALYSIS)/bin
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --projroot $(ROOT) manuscript-figures.nf &

diamante-locuszoom:
	mkdir -p $(ANALYSIS)
	printf "MarkerName\tP-value\n" > $(ANALYSIS)/gwas.txt
	cat /lab/work/porchard/data/gwas/diamante/Mahajan.NatGenet2018b.T2D.European.txt | cut -f1,9 | grep -v Pvalue | perl -pe 's/^5:53271420\t/rs702634\t/; s/^12:26472562\t/rs7132434\t/' >> $(ANALYSIS)/gwas.txt
	cd $(ANALYSIS) && locuszoom --metal $(ANALYSIS)/gwas.txt --refsnp rs702634 --flank 300kb --build hg19 --pop EUR --source 1000G_Nov2014 theme=publication height=6 hiStart=53269334 hiEnd=53278795
	cd $(ANALYSIS) && locuszoom --metal $(ANALYSIS)/gwas.txt --refsnp rs7132434 --flank 300kb --build hg19 --pop EUR --source 1000G_Nov2014 theme=publication height=6 hiStart=26436640 hiEnd=26491955

diamante-bmiadj-locuszoom:
	mkdir -p $(ANALYSIS)
	printf "MarkerName\tP-value\n" > $(ANALYSIS)/gwas.txt
	cat /lab/work/porchard/data/gwas/diamante/Mahajan.NatGenet2018b.T2Dbmiadj.European.txt | cut -f1,9 | grep -v Pvalue >> $(ANALYSIS)/gwas.txt
	cd $(ANALYSIS) && locuszoom --metal $(ANALYSIS)/gwas.txt --refsnp 5:53271420 --flank 300kb --build hg19 --pop EUR --source 1000G_Nov2014 theme=publication height=6 hiStart=53269334 hiEnd=53278795
	cd $(ANALYSIS) && locuszoom --metal $(ANALYSIS)/gwas.txt --refsnp 12:26472562 --flank 300kb --build hg19 --pop EUR --source 1000G_Nov2014 theme=publication height=6 hiStart=26436640 hiEnd=26491955

creatinine-locuszoom:
	mkdir -p $(ANALYSIS)
	printf "MarkerName\tP-value\n" > $(ANALYSIS)/gwas.txt
	zcat $(DATA)/ukb-summary-stats/30700_irnt.gwas.imputed_v3.both_sexes.tsv.gz | cut -f1,11 | grep -v pval | perl -pe 's/(.*):(.*):(.*):(.*)\t/$$1:$$2\t/' | perl -pe 's/^17:54776955\t/rs227727\t/' >> $(ANALYSIS)/gwas.txt
	cd $(ANALYSIS) && locuszoom --metal $(ANALYSIS)/gwas.txt --refsnp rs227727 --flank 300kb --build hg19 --pop EUR --source 1000G_Nov2014 theme=publication height=6 hiStart=54735000 hiEnd=54800000

knit: ANALYSIS=$(ROOT)/knit
knit: DOWNSTREAM=$(WORK)/downstream-new-features/results/nucleus-qc
knit: all
	# FANS vs no FANS figures
	ln -sf $(WORK)/manuscript-figures/results/fans-vs-no-fans/hqaa-vs-tss-enrichment-fans-vs-no-fans.png $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/fans-vs-no-fans/hqaa-vs-max-fraction-reads-from-single-autosome-fans-vs-no-fans.png $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/fans-vs-no-fans/fans-atac-nuclei-per-library.png $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/fans-vs-no-fans/fans-atac-nuclei-read-counts.png $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/fans-vs-no-fans/fans-aggregate-fld.pdf $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/fans-vs-no-fans/fans-aggregate-tss.pdf $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/fans-vs-no-fans/fans-aggregate-ataqv-metrics.pdf $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/fans-vs-no-fans/fans-chromhmm-overlap.pdf $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/fans-vs-no-fans/umis-vs-mitochondrial-fans-vs-no-fans.png $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/fans-vs-no-fans/fans-rna-nuclei-per-library.png $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/fans-vs-no-fans/fans-rna-correlation.png $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/fans-vs-no-fans/fans-gb.pdf $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/fans/fans-atac-correlation.png $(ANALYSIS)/
	# 20k vs 40k figures
	ln -sf $(WORK)/manuscript-figures/results/20k-vs-40k/hqaa-vs-tss-enrichment-20k-vs-40k.png $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/20k-vs-40k/hqaa-vs-max-fraction-reads-from-single-autosome-20k-vs-40k.png $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/20k-vs-40k/loading-aggregate-fld.pdf $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/20k-vs-40k/loading-aggregate-tss.pdf $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/20k-vs-40k/loading-chromhmm-overlap.pdf $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/20k-vs-40k/loading-gb.pdf $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/20k-vs-40k/umis-vs-mitochondrial-20k-vs-40k.png $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/20k-vs-40k/loading-rna-correlation.png $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/loading/loading-atac-correlation.png $(ANALYSIS)/
	# Feb. sample figures
	ln -sf $(WORK)/manuscript-figures/results/feb/hqaa-vs-tss-enrichment-feb.png $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/feb/hqaa-vs-max-fraction-reads-from-single-autosome-feb.png $(ANALYSIS)/
	ln -sf $(DOWNSTREAM)/hg19-rn6-ratio-threshold.png $(ANALYSIS)/species-assignment-feb.png
	# Libraries used downstream
	ln -sf $(WORK)/manuscript-figures/results/qc-for-downstream-libraries/*.png $(ANALYSIS)/
	# Figures for the biology part of the manuscript
	ln -sf $(WORK)/manuscript-figures/results/tables/nuclei-summary-stats.csv $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/umaps/full-umap.png $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/umaps/split-umap.png $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/marker-gene-heatmaps/marker-genes.pdf $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/per-nucleus-plots/MYH1-plus-MYH4-vs-MYH7.png $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/logistic-regression/hg19/enhancer-similarity-heatmap.pdf $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/cell-type-proportions/fraction-of-library-within-individual.jitter.pdf $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/fiber-type-rna-diffs/rubenstein-vs-our-fiber-type-lfcs.pdf $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/UKB-GWAS-enrichments/hg19/UKB-LDSC-bonferroni.pdf $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/UKB-GWAS-enrichments/hg19/UKB-LDSC-bonferroni.pdf $(ANALYSIS)/UKB-LDSC-bonferroni-hg19.pdf
	ln -sf $(WORK)/manuscript-figures/results/UKB-GWAS-enrichments/rn6/UKB-LDSC-bonferroni.pdf $(ANALYSIS)/UKB-LDSC-bonferroni-rn6.pdf
	ln -sf $(WORK)/manuscript-figures/results/UKB-GWAS-enrichments/hg19/UKB-LDSC-by.pdf $(ANALYSIS)/UKB-LDSC-by-hg19.pdf
	ln -sf $(WORK)/manuscript-figures/results/UKB-GWAS-enrichments/rn6/UKB-LDSC-by.pdf $(ANALYSIS)/UKB-LDSC-by-rn6.pdf
	ln -sf $(WORK)/manuscript-figures/results/T2D-GWAS-enrichments-no-bmiadj/hg19/T2D-FIns.pdf $(ANALYSIS)/T2D-FIns-hg19.pdf
	ln -sf $(WORK)/manuscript-figures/results/T2D-GWAS-enrichments-no-bmiadj/rn6/T2D-FIns.pdf $(ANALYSIS)/T2D-FIns-rn6.pdf
	ln -sf $(WORK)/manuscript-figures/results/UKB-GWAS-enrichments/hg19/UKB-LDSC-top.pdf $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/UKB-GWAS-enrichments/hg19/UKB-top-in-our-cell-types.pdf $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/gb-screenshots/ARL15/ARL15-other-cell-types.chr5_53269196_53277711.pdf $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/gb-screenshots/ITPR2/ITPR2-other-cell-types.chr12_26436640_26491955.pdf $(ANALYSIS)/
	ln -sf $(WORK)/manuscript-figures/results/chromatin-state-legend/chrom_state_legend.pdf $(ANALYSIS)/
	cd $(ANALYSIS) && pdflatex main.tex && pdflatex main.tex
