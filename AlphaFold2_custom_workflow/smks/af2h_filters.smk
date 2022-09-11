rule all:
	input:
		"config_multimer_complexes.json"

rule filter_plddts:
	input:
		dirA="outdir/multimer/A",
		dirB="outdir/multimer/B",
		summary_path="outdir/multimer",
		fasta_dir="fasta_wdir/multimer/complexes",
		repo_msa="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas",
		repo_multimer="/data/dops-plant-genomics/bras4304/alphafold/alphafold_multimer"

	output:
		"config_multimer_complexes.json"

	conda:
		"envs/biopython.yaml"

	resources:
		partition="short",
		time_min=60,
		cluster="arc",
		gpus=0,
		name="filter_proteins",
		cpus=1,
		mem_mb=8000,
		email="felix.homma@bnc.ox.ac.uk"
	
	shell:
		"scripts/score_monomers.py {input.dirA} {input.dirB} {input.summary_path} {input.fasta_dir} {input.repo_msa} {input.repo_multimer}"
