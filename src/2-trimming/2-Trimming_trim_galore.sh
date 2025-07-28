#!/bin/bash
#SBATCH --time=15:00:00
#SBATCH --partition=biogipuzkoa-exclusive
#SBATCH --account=biogipuzkoa-exclusive
#SBATCH --mem-per-cpu=5G 
#SBATCH -o trim_galore_%A_%a.out
#SBATCH --array=0-85%10

module load Python
conda activate /scratch/mcarazo/envs/trimming/

cd /scratch/mcarazo/ongoing/GB/sRNAseq_may25/
mkdir -p trim_galore_out
cd trim_galore_out

# Select the fastq files
fastqs=($(realpath ../fastqs/*_R1.fastq.gz))


fq1=${fastqs[$SLURM_ARRAY_TASK_ID]}
fq2=${fq1/_R1.fastq.gz/_R2.fastq.gz}
base=$(basename ${fq1} _R1.fastq.gz)

trim_galore -q 20 --paired --retain_unpaired --small_rna --fastqc --fastqc_args "-o fastqc_results" --gzip $fq1 $fq2 

