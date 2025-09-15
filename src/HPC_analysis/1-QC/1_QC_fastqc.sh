#!/bin/bash

#SBATCH -t 8:00:00
#SBATCH -c 6
#SBATCH --mem-per-cpu=12G
#SBATCH --job-name=Qcontrol
#SBATCH --output=Qcontrol_%j.out

module load Miniforge3
conda activate /scratch/mcarazo/envs/QControl/

cd /scratch/mcarazo/ongoing/GB/sRNAseq_may25/fastqs

seqkit stat *.fastq.gz
 
mkdir fastqc_results
fastqc *.fastq.gz -o fastqc_results/
cd fastqc_results
multiqc .