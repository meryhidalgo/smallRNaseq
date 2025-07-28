#!/usr/bin/env bash

#SBATCH -c 2
#SBATCH --partition=biogipuzkoa-exclusive
#SBATCH --account=biogipuzkoa-exclusive
#SBATCH --mem=5G

module load Miniforge3
conda activate /scratch/mcarazo/envs/DGE_analysis

cd /scratch/mcarazo/ongoing/GB/sRNAseq_may25/

outputDir="salmon_out_gencode_tRNA"
metric=lengthScaledTPM
designFile="/scratch/mcarazo/ongoing/GB/sRNAseq_may25/design/design.tab"
gtf="/scratch/mcarazo/indexes/tRNA/gencode.v48.primary_tRNA.gtf" # -> sRNA_gencode_tRNA
outname="sRNA_gencode_tRNA"

Rscript /scratch/mcarazo/ongoing/GB/sRNAseq_may25/scripts/GetSalmonTPMs.R \
    $metric \
    $outputDir \
    $designFile \
    $gtf \
    $outname

#mv sRNA_rRNA_GB_Control_COUNTS_AGGREGATED-lengthScaledTPM.tab salmon_out_rRNA/

# DiffExpressionAnalysis.R in locals
