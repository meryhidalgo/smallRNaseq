#!/bin/bash

#SBATCH -c 6
#SBATCH --partition=biogipuzkoa-exclusive
#SBATCH --account=biogipuzkoa-exclusive
#SBATCH --mem-per-cpu=12G
#SBATCH -o trimmomatic_%A_%a.out
#SBATCH --array=0-85%10

module load Miniforge3
conda activate /scratch/mcarazo/envs/biotools/

cd /scratch/mcarazo/ongoing/GB/sRNAseq_may25/
mkdir -p trimmomatic_out
cd trimmomatic_out

# Select the fastq files
fastqs=($(realpath ../fastqs/*_R1.fastq.gz))

fq1=${fastqs[$SLURM_ARRAY_TASK_ID]}
fq2=${fq1/_R1/_R2}
base=$(basename ${fq1} _R1.fastq.gz)

java -jar /scratch/mcarazo/envs/trimming/share/trimmomatic-0.39-2/trimmomatic.jar PE -phred33 ${fq1} ${fq2} \
    ""${base}_R1.fastq.gz"" ""${base}_R1_untrimmed.fastq.gz"" \
    ""${base}_R2.fastq.gz"" ""${base}_R2_untrimmed.fastq.gz"" \
    ILLUMINACLIP:/scratch/mcarazo/envs/trimming/share/trimmomatic-0.39-2/adapters/sRNA-adapters.fa:2:12:6 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15

rm ""${base}_R1_untrimmed.fastq.gz"" ""${base}_R2_untrimmed.fastq.gz""

