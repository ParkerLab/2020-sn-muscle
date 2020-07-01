// Nextflow configuration file for excecuting in Parker Lab Environment

// Universal params (independent of profile)
params {
    model {
        holdout_frac = 0.15
        l = 10
        k = 6
        seed = 57
        window_size = 150
    }
}
// profile based params
profiles {
    test {
        params.input.cluster_peaks_dir = "/home/jeremybk/repos/snATAC_deltaSVM/tests/workflow/input/*.broadPeak.noblacklist"
        params.input.bam_dir = "/home/jeremybk/repos/snATAC_deltaSVM/tests/workflow/input/*.bam"
        params.input.snp_file = "/home/jeremybk/repos/snATAC_deltaSVM/tests/workflow/input/test.vcf"
        params.input.eqtl_file = "/home/jeremybk/repos/snATAC_deltaSVM/tests/workflow/input/eqtl.tab.gz"
        params.output.results_dir  = "/home/jeremybk/repos/snATAC_deltaSVM/tests/workflow/results/"
        workDir = "/home/jeremybk/repos/snATAC_deltaSVM/tests/workflow/work"
    }
}

process {
    withLabel: cluster_submit {
        executor = "slurm"
	containerOptions='--bind "/lab:/lab" --bind "/localscratch:/tmp"'
	clusterOptions='--constraint=wolverine'
    }
    withLabel: container_conda {
        container = '/lab/work/jeremybk/snATAC_deltaSVM/images/conda_env.img'
	containerOptions='--bind "/lab:/lab" --bind "/localscratch:/tmp"'
    }
    withLabel: container_svmgkm {
        container = '/home/porchard/github/snATAC_deltaSVM/containers/gkm_svm.img'
	containerOptions='--bind "/lab:/lab" --bind "/localscratch:/tmp"'
    }
    withLabel: container_general {
        container = '/lab/work/porchard/sn-muscle-project/singularity/general/general.simg'
	containerOptions='--bind "/lab:/lab" --bind "/localscratch:/tmp"'
    }
}

singularity {
    enabled = true
    autoMounts = true
}
