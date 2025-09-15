#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH -c 6
#SBATCH --mem=16G
#SBATCH --partition=biogipuzkoa-exclusive
#SBATCH --account=biogipuzkoa-exclusive
#SBATCH -o salmon_%A_%a.out
#SBATCH --array=0-85%10

# activate the env
module load Miniforge3
conda activate /scratch/mcarazo/envs/salmon

cd /scratch/mcarazo/ongoing/GB/sRNAseq_may25/
#comment or decomment as necessary to reuse the script. remember to change dir value and select appropriate index
dir="salmon_out_EMN_45S_rRNA"
mkdir -p $dir

#gencode
index=/scratch/mcarazo/indexes/h_salmon_k15 #-> salmon_out_gencode
#refseq
index=/scratch/mcarazo/indexes/h_salmon_refseq_k15 #-> salmon_out_RefSeq
#manual annotation rRNA
index=/scratch/mcarazo/indexes/rRNA/h_salmon_EMN_45S_rRNA_k15 # for rRNA quantification -> salmon_out_EMN_45S_rRNA

fastqs=($(realpath fastp_out/*_R1.fastq.gz))


fq1=${fastqs[$SLURM_ARRAY_TASK_ID]}
fq2=${fq1/_R1.fastq.gz/_R2.fastq.gz}
base=$(basename ${fq1} _trim_R1.fastq.gz)

# quantify transcripts in fastq files
salmon quant -i $index \
    -l A \
	-1 $fq1 \
	-2 $fq2 \
	-p 6 \
	--validateMappings \
    -o $dir/${base}
