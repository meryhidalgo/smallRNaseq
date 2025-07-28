#!/bin/bash

#SBATCH -c 6
#SBATCH --partition=biogipuzkoa-exclusive
#SBATCH --account=biogipuzkoa-exclusive
#SBATCH --mem-per-cpu=12G
#SBATCH -o fastp_%A_%a.out
#SBATCH --array=0-85%10

module load Miniforge3
conda activate /scratch/mcarazo/envs/biotools/

cd /scratch/mcarazo/ongoing/GB/sRNAseq_may25/
mkdir -p fastp_out
cd fastp_out

# Select the fastq files
fastqs=($(realpath ../fastqs/*_R1.fastq.gz))

fq1=${fastqs[$SLURM_ARRAY_TASK_ID]}
fq2=${fq1/_R1/_R2}
base=$(basename ${fq1} _R1.fastq.gz)

fastp -x -y -i $fq1 -I $fq2 -o ${base}_trim_R1.fastq.gz -O ${base}_trim_R2.fastq.gz \
	--adapter_fasta /scratch/mcarazo/envs/trimming/share/trimmomatic-0.39-2/adapters/sRNA-adapters.fa

