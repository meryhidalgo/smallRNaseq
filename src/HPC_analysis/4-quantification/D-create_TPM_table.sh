#!/usr/bin/env bash

#SBATCH -c 2
#SBATCH --partition=biogipuzkoa-exclusive
#SBATCH --account=biogipuzkoa-exclusive
#SBATCH --mem=5G

module load Miniforge3
conda activate /scratch/mcarazo/envs/DGE_analysis

cd /scratch/mcarazo/ongoing/GB/sRNAseq_may25/

#change outputDir, gtf and outname to reuse the script for all biotypes
#EMN_45S_rRNA.gtf can be found in references directory of this repo
#the rest of gtf are generated from fasta files

outputDir="salmon_out_EMN_45S_rRNA"
metric=lengthScaledTPM
designFile="/scratch/mcarazo/ongoing/GB/sRNAseq_may25/design/design.tab"
gtf="/scratch/mcarazo/indexes/rRNA/EMN_45S_rRNA.gtf" 
outname="sRNA_rRNA_EMN_45S"

Rscript /scratch/mcarazo/ongoing/GB/sRNAseq_may25/scripts/GetSalmonTPMs_mch.R \
    $metric \
    $outputDir \
    $designFile \
    $gtf \
    $outname
