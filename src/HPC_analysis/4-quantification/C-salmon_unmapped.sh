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

# SCRIPT TO BE USED AFTER rRNA STAR ALIGNMENT TO QUANTIFY UNMAPPED READS AGAINST OTHER RNA TYPES (miRNA, tRNA, snoRNA, lncRNA...)
cd /scratch/mcarazo/ongoing/GB/sRNAseq_may25/star_out_EMN_45S_rRNA
#comment or decomment as necessary to reuse the script. remember to change dir value and select appropriate index
dir="salmon_out_gencode_tRNA"
mkdir -p $dir

index=/scratch/mcarazo/indexes/miRBase/miRNA_hairpin_salmon_k15 #-> salmon_out_miRNA_hairpin
index=/scratch/mcarazo/indexes/miRBase/miRNA_mature_salmon_k15 #-> salmon_out_miRNA_mature
index=/scratch/mcarazo/indexes/h_salmon_k15 #-> salmon_out_gencode
index=/scratch/mcarazo/indexes/tRNA/GtRNAdb_salmon_k15 #-> salmon_out_GtRNAdb
index=/scratch/mcarazo/indexes/snoRNA/snoRNA_salmon_k15 #-> salmon_out_snoRNA
index=/scratch/mcarazo/indexes/lncRNA/lncRNABook_salmon_k15 #-> salmon_out_lncRNABook


# Select the fastq files
fastqs=($(realpath unmapped_fastqs/*.mate1.fastq))


fq1=${fastqs[$SLURM_ARRAY_TASK_ID]}
fq2=${fq1/.mate1.fastq/.mate2.fastq}
base=$(basename ${fq1} _trim_Unmapped.out.mate1.fastq)


# quantify transcripts in fastq files
salmon quant -i $index \
    -l A \
	-1 $fq1 \
	-2 $fq2 \
	-p 6 \
	--validateMappings \
    -o $dir/${base}
