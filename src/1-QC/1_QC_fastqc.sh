#!/bin/bash

#SBATCH -t 12:00:00
#SBATCH -c 6 		
#SBATCH --mem-per-cpu=12G
#SBATCH --job-name=Qcontrol
#SBATCH --output=Qcontrol_%j.out

module load Miniforge3
conda activate /scratch/mcarazo/envs/QControl/

cd /scratch/mcarazo/ongoing/GB/sRNAseq_may25/trim_galore_out

seqkit stat *.fastq.gz
 
#cd  FASTQs/fastp_out/
mkdir fastqc_results
fastqc *.fastq.gz -o fastqc_results/
cd fastqc_results
multiqc .