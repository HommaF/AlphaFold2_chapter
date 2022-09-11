configfile: "config_multimer_complexes.json"

rule all:
	input:
		expand("/data/dops-plant-genomics/bras4304/alphafold/alphafold_multimer/{sample}/token_compressed.txt", sample=config["samples"])


rule alphafold_multimers:
	input:
		fasta="fasta_wdir/multimer/complexes/{sample}.fasta",
		
	params:
		template_date="2021-11-02",
		db_preset="reduced_dbs",
		model_preset="multimer",
		alphafold_data_dir="/data/dops-plant-genomics/bras4304/alphafold/alphafold_data",
		outdir="outdir/multimer/complexes/",
		repo_multi="/data/dops-plant-genomics/bras4304/alphafold/alphafold_multimer/",
		repo_multi_target="/data/dops-plant-genomics/bras4304/alphafold/alphafold_multimer/{sample}",
		precomputed_msas=True

	output:
		repo_rank=protected("/data/dops-plant-genomics/bras4304/alphafold/alphafold_multimer/{sample}/ranking_debug.json")


	resources:
		partition="short",
		time_min=720,
		cluster="htc",
		gpus=1,
		name="Multi_{sample}",
		cpus=1,
		mem_mb=128000,
		email="felix.homma@bnc.ox.ac.uk"

	shell:
		"module load AlphaFold/2.1.1-fosscuda-2020b;"
		"module load TensorFlow/2.3.1-foss-2020a-Python-3.8.2;"
		"export ALPHAFOLD_DATA_DIR=/data/dops-plant-genomics/bras4304/alphafold/alphafold_data;"
		"alphafold --data_dir={params.alphafold_data_dir} --fasta_paths={input.fasta} --output_dir={params.repo_multi} --max_template_date={params.template_date} --db_preset={params.db_preset} --model_preset={params.model_preset} --use_precomputed_msas={params.precomputed_msas};"
		"ln -sf {params.repo_multi_target} {params.outdir}"

rule link_multimer:
	input:
		repo_rank="/data/dops-plant-genomics/bras4304/alphafold/alphafold_multimer/{sample}/ranking_debug.json"

	params:
		repo_multi_target="/data/dops-plant-genomics/bras4304/alphafold/alphafold_multimer/{sample}",

	output:
		outdir=directory("outdir/multimer/complexes/{sample}")

	resources:
		partition="short",
		time_min=10,
		cluster="arc",
		name="link_multi",
		cpus=1,
		mem_mb=8000,
		email="felix.homma@bnc.ox.ac.uk"

	shell:
		"ln -sf {params.repo_multi_target} {output.outdir}"


rule compress_alphafold_output:
	input:
		repo_rank="/data/dops-plant-genomics/bras4304/alphafold/alphafold_multimer/{sample}/ranking_debug.json",
		outdir="outdir/multimer/complexes/{sample}"

	output:
		repo_multi_target="/data/dops-plant-genomics/bras4304/alphafold/alphafold_multimer/{sample}/token_compressed.txt"

	params:
		repo_multi_target="/data/dops-plant-genomics/bras4304/alphafold/alphafold_multimer/{sample}",
		pkls_repo="/data/dops-plant-genomics/bras4304/alphafold/alphafold_multimer_pkls"

	conda:
		"envs/biopython.yaml"
	
	resources:
		partition="short",
		time_min=90,
		cluster="arc",
		gpus=0,
		name="compress_multi",
		cpus=1,
		mem_mb=8000,
		email="felix.homma@bnc.ox.ac.uk"

	shell:
		"scripts/compress_alphafold_output.py {params.repo_multi_target};"
		"scripts/move_pkls.py {params.repo_multi_target} {params.pkls_repo};"
		"touch {output.repo_multi_target}"
