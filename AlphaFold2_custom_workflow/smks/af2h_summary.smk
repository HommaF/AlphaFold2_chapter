rule all:
	input:
		"outdir/multimer/multimers_scoring.tsv"

rule summary_ptms:
	input:
		complexes = "outdir/multimer/complexes",
		summary_path = "outdir/multimer"

	output:
		"outdir/multimer/multimers_scoring.tsv"

	conda:
		"envs/biopython.yaml"

	resources:
		partition="short",
		time_min=60,
		cluster="arc",
		gpus=0,
		name="complex_scoring",
		cpus=1,
		mem_mb=8000,
		email="felix.homma@bnc.ox.ac.uk"
	
	shell:
		"scripts/score_multimers.py {input.complexes} {input.summary_path}"
