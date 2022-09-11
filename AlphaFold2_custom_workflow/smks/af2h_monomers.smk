configfile: "config_multimer.json"

rule all:
	input:
		expand("fasta_wdir/multimer/complexes/{sample_A}_{sample_B}.fasta", sample_A=config["samples_A"], sample_B=config["samples_B"])



rule alphafold_monomer_A:
	input:
		fasta="fasta_wdir/multimer/A/{sample_A}.fasta",
		
	params:
		template_date="2021-11-02",
		db_preset="reduced_dbs",
		model_preset="monomer",
		alphafold_data_dir="/data/dops-plant-genomics/bras4304/alphafold/alphafold_data",
		out_monomer="outdir/multimer/A/{sample_A}",
		repo_monomer="/data/dops-plant-genomics/bras4304/alphafold/alphafold_monomer/{sample_A}",
		repo_msa="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/{sample_A}",
		repo_out="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas",
		msa_file="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/{sample_A}/msas",


	output:
		repo_rank=protected("/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/{sample_A}/ranking_debug.json")

	resources:
		partition="short",
		time_min=720,
		cluster="htc",
		gpus=1,
		name="Mono_{sample_A}",
		cpus=8,
		mem_mb=90000,
		email="felix.homma@bnc.ox.ac.uk"

	shell:
		"module load AlphaFold/2.1.1-fosscuda-2020b;"
		"module load TensorFlow/2.3.1-foss-2020a-Python-3.8.2;"
		"export ALPHAFOLD_DATA_DIR=/data/dops-plant-genomics/bras4304/alphafold/alphafold_data;"
		"alphafold --data_dir={params.alphafold_data_dir} --fasta_paths={input.fasta} --output_dir={params.repo_out} --max_template_date={params.template_date} --db_preset={params.db_preset} --model_preset={params.model_preset}"



rule link_monomer_A:
	input:
		"/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/{sample_A}/ranking_debug.json"

	output:
		directory("outdir/multimer/A/{sample_A}")

	params:
		repo="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/{sample_A}/"

	resources:
		partition="short",
		cluster='arc',
		time_min=10,
		name="link_A",
		cpus=1,
		mem_mb=8000,
		email="felix.homma@bnc.ox.ac.uk"
	
	shell:
		"ln -sf {params.repo} {output}"



rule jackhmmer_uniprot_A:
	input:
		repo_rank="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/{sample_A}/ranking_debug.json",
		fasta="fasta_wdir/multimer/A/{sample_A}.fasta",
		data_dir="/data/dops-plant-genomics/bras4304/alphafold/alphafold_data/uniprot/uniprot.fasta"

	output:
		uni_out="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/{sample_A}/msas/uniprot_hits.sto"

	resources:
		cpus=3,
		partition="short",
		time_min=60,
		cluster='arc',
		name="jackhmmer_uniprot",
		mem_mb=16000,
		email="felix.homma@bnc.ox.ac.uk"

	shell:
		"module load AlphaFold/2.1.1-fosscuda-2020b;"
		"scripts/jackhmmer.py {input.data_dir} {input.fasta} {output.uni_out} {resources.cpus}"

rule compress_A:
	input:
		repo_rank="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/{sample_A}/ranking_debug.json"

	output:
		gz_out="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/{sample_A}/A_compressed_token.txt"


	params:
		repo_msa="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/{sample_A}"


	conda:
		"envs/biopython.yaml"

	resources:
		partition="short",
		time_min=90,
		cluster="arc",
		gpus=0,
		name="compress_A",
		cpus=1,
		mem_mb=8000,
		email="felix.homma@bnc.ox.ac.uk"

	shell:
		"scripts/compress_alphafold_output.py {params.repo_msa};"
		"touch {output.gz_out}"


rule alphafold_monomer_B:
	input:
		fasta="fasta_wdir/multimer/B/{sample_B}.fasta",
		
	params:
		template_date="2021-11-02",
		db_preset="reduced_dbs",
		model_preset="monomer",
		alphafold_data_dir="/data/dops-plant-genomics/bras4304/alphafold/alphafold_data",
		out_monomer="outdir/multimer/B/{sample_B}",
		repo_monomer="/data/dops-plant-genomics/bras4304/alphafold/alphafold_monomer/{sample_B}",
		repo_msa="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/{sample_B}",
		repo_out="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas",
		msa_file="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/{sample_B}/msas",


	output:
		repo_rank=protected("/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/{sample_B}/ranking_debug.json")

	resources:
		partition="short",
		time_min=600,
		cluster="htc",
		gpus=1,
		name="Mono_{sample_B}",
		cpus=8,
		mem_mb=64000,
		email="felix.homma@bnc.ox.ac.uk"

	shell:
		"module load AlphaFold/2.1.1-fosscuda-2020b;"
		"module load TensorFlow/2.3.1-foss-2020a-Python-3.8.2;"
		"export ALPHAFOLD_DATA_DIR=/data/dops-plant-genomics/bras4304/alphafold/alphafold_data;"
		"alphafold --data_dir={params.alphafold_data_dir} --fasta_paths={input.fasta} --output_dir={params.repo_out} --max_template_date={params.template_date} --db_preset={params.db_preset} --model_preset={params.model_preset}"


rule link_monomer_B:
	input:
		"/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/{sample_B}/ranking_debug.json"

	output:
		directory("outdir/multimer/B/{sample_B}")
	
	params:
		repo="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/{sample_B}/"

	resources:
		partition="short",
		cluster='arc',
		time_min=10,
		name="link_B",
		cpus=1,
		mem_mb=8000,
		email="felix.homma@bnc.ox.ac.uk"
	
	shell:
		"ln -sf {params.repo} {output}"



rule jackhmmer_uniprot_B:
	input:
		repo_rank="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/{sample_B}/ranking_debug.json",
		fasta="fasta_wdir/multimer/B/{sample_B}.fasta",
		data_dir="/data/dops-plant-genomics/bras4304/alphafold/alphafold_data/uniprot/uniprot.fasta"

	output:
		uni_out="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/{sample_B}/msas/uniprot_hits.sto"

	resources:
		cpus=3,
		partition="short",
		time_min=60,
		cluster='arc',
		name="jackhmmer_uniprot",
		mem_mb=16000,
		email="felix.homma@bnc.ox.ac.uk"

	shell:
		"module load AlphaFold/2.1.1-fosscuda-2020b;"
		"scripts/jackhmmer.py {input.data_dir} {input.fasta} {output.uni_out} {resources.cpus}"

rule compress_B:
	input:
		repo_rank="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/{sample_B}/ranking_debug.json",


	output:
		gz_out="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/{sample_B}/B_compressed_token.txt"


	params:
		repo_msa="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/{sample_B}"


	conda:
		"envs/biopython.yaml"

	resources:
		partition="short",
		time_min=90,
		cluster="arc",
		gpus=0,
		name="compress_B",
		cpus=1,
		mem_mb=8000,
		email="felix.homma@bnc.ox.ac.uk"

	shell:
		"scripts/compress_alphafold_output.py {params.repo_msa};"
		"touch {output.gz_out}"


rule make_multi_fasta:
	input:
		sample_A="fasta_wdir/multimer/A/{sample_A}.fasta",
		sample_B="fasta_wdir/multimer/B/{sample_B}.fasta",
		gz_out_A="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/{sample_A}/A_compressed_token.txt",
		gz_out_B="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/{sample_B}/B_compressed_token.txt",
		uni_out_A="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/{sample_A}/msas/uniprot_hits.sto",
		uni_out_B="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/{sample_B}/msas/uniprot_hits.sto",
		link_A="outdir/multimer/A/{sample_A}",
		link_B="outdir/multimer/B/{sample_B}"


	output:
		fasta_files="fasta_wdir/multimer/complexes/{sample_A}_{sample_B}.fasta"

	params:
		fasta_path="fasta_wdir/multimer/complexes",
		outdir_path="outdir/multimer",
		repo_msas="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas"

	conda:
		"envs/biopython.yaml"

	resources:
		partition="short",
		time_min=10,
		cluster="arc",
		gpus=0,
		name="prep_mff",
		cpus=1,
		mem_mb=8000,
		email="felix.homma@bnc.ox.ac.uk"

	shell:
		"scripts/make_multifasta.py {input.sample_A} {input.sample_B} {params.fasta_path} {params.outdir_path} {params.repo_msas}"
