#!/bin/bash

snakemake -s smks/af2h_monomers.smk --profile slurm
snakemake -s smks/af2h_filters.smk --profile slurm
./prep_files.sh
snakemake -s smks/af2h_multimers.smk --profile slurm
snakemake -s smks/af2h_summary.smk --profile slurm
