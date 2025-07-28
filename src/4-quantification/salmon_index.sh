#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH -c 6
#SBATCH --mem-per-cpu=12G

# activate the env
module load Miniforge3
conda activate /scratch/mcarazo/envs/salmon

mkdir h_salmon_EMN_45S_rRNA_k15
mkdir h_salmon_primary_tRNA_k15

salmon index -k 15 --gencode -t EMN_45S_rRNA.fa.gz -i h_salmon_EMN_45S_rRNA_k15
salmon index -k 15 --gencode -t Gencode.primary_tRNA_transcripts.fa -i h_salmon_primary_tRNA_k15 