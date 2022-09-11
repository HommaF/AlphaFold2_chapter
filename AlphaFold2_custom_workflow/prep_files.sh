#!/bin/bash

## For monomer pipeline

base_mono_repo="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/"
ranking_mono_file="/ranking_debug.json"
compress_file="compressed_token.txt"
uniprot_repo="/data/dops-plant-genomics/bras4304/alphafold/alphafold_msas/"
uniprot_file="/msas/uniprot_hits.sto"
link_mono_A="/data/dops-plant-genomics/bras4304/smk_pipelines/af2h_beta/outdir/multimer/A/*"
link_mono_B="/data/dops-plant-genomics/bras4304/smk_pipelines/af2h_beta/outdir/multimer/B/*"

echo "touch monomers"

ls fasta_wdir/multimer/A/*fasta | cut -d "/" -f4 | sed "s/\.fasta//g" | xargs -I@ -n 1 -P 1 bash -c "scripts/touch_ranking.sh $base_mono_repo @ $ranking_mono_file $compress_file"
ls fasta_wdir/multimer/B/*fasta | cut -d "/" -f4 | sed "s/\.fasta//g" | xargs -I@ -n 1 -P 1 bash -c "scripts/touch_ranking.sh $ranking_mono_repo @ $ranking_mono_file $compress_file"

#echo "touch compression files"

#ls fasta_wdir/multimer/A/*fasta | cut -d "/" -f4 | sed "s/\.fasta//g" | xargs -I@ -n 1 -P 1 bash -c "scripts/touch_existing.sh $compress_repo @ /A_$compress_file"
#ls fasta_wdir/multimer/B/*fasta | cut -d "/" -f4 | sed "s/\.fasta//g" | xargs -I@ -n 1 -P 1 bash -c "scripts/touch_existing.sh $compress_repo @ /B_$compress_file"

echo "touch uniprot files"

ls fasta_wdir/multimer/A/*fasta | cut -d "/" -f4 | sed "s/\.fasta//g" | xargs -I@ -n 1 -P 1 bash -c "scripts/touch_existing.sh $uniprot_repo @ $uniprot_file"
ls fasta_wdir/multimer/B/*fasta | cut -d "/" -f4 | sed "s/\.fasta//g" | xargs -I@ -n 1 -P 1 bash -c "scripts/touch_existing.sh $uniprot_repo @ $uniprot_file"

touch -h $link_mono_A
touch -h $link_mono_B

touch fasta_wdir/multimer/complexes/*fasta

## For multimer pipeline

base_multi_repo="/data/dops-plant-genomics/bras4304/alphafold/alphafold_multimer/"
ranking_file="/ranking_debug.json"
compress_file="token_compressed.txt"

link_multi="/data/dops-plant-genomics/bras4304/smk_pipelines/af2h_beta/outdir/multimer/complexes/*"



if test -f config_multimer_complexes.json; then
	echo "touching multimers"
	tail -n +3 config_multimer_complexes.json | head -n -2 | cut -d "\"" -f2 | xargs -I@ -n 1 -P 1 bash -c "scripts/touch_ranking.sh $base_multi_repo @ $ranking_file $compress_file"
	touch -h $link_multi
fi
