#!/bin/bash

af2h_dir=$1
multifasta_dir='/data/dops-plant-genomics/bras4304/alphafold/alphafold_multifastas'

mkdir $af2h_dir

mv config_multimer_complexes.json $af2h_dir
cp config_multimer.json $af2h_dir

mv fasta_wdir/multimer/complexes/*fasta $multifasta_dir
mv fasta_wdir/ $af2h_dir
mkdir -p fasta_wdir/multimer/A
mkdir -p fasta_wdir/multimer/B
mkdir -p fasta_wdir/multimer/complexes

mv logs_slurm $af2h_dir
mkdir logs_slurm

mv outdir $af2h_dir
mkdir -p outdir/multimer/A
mkdir -p outdir/multimer/B
mkdir -p outdir/multimer/complexes

cp *sh $af2h_dir
cp -r scripts $af2h_dir
cp -r smks $af2h_dir

