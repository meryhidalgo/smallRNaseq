#!/bin/bash

#SBATCH -t 08:00:00
#SBATCH -c 6
#SBATCH --partition=biogipuzkoa-exclusive
#SBATCH --account=biogipuzkoa-exclusive
#SBATCH --mem=40G
#SBATCH -o STAR_R_%A_%a.out
#SBATCH --job-name=STAR
#SBATCH --array=0-85%10

module load Miniforge3
conda activate /scratch/mcarazo/envs/STAR/

# create a directory for STAR analysis
cd /scratch/mcarazo/ongoing/GB/sRNAseq_may25
mkdir -p star_out_GRCh38

#created with GRCh38.primary_assembly.genome.fa and gencode.v41.primary_assembly.annotation.gtf
index=/data/mcarazo/indexes/h_index/index_h99/ 

# Select the fastq files
fastqs=($(realpath fastp_out/*_R1.fastq.gz))

fq1=${fastqs[$SLURM_ARRAY_TASK_ID]}
fq2=${fq1/_R1/_R2}
base=$(basename ${fq1} _R1.fastq.gz)


STAR --genomeDir $index \
  --runThreadN 6 \
  --readFilesIn  $fq1 $fq2 \
  --readFilesCommand zcat \
  --outFileNamePrefix "/scratch/mcarazo/ongoing/GB/sRNAseq_may25/star_out_GRCh38/${base}_" \
  --outSAMtype BAM SortedByCoordinate \
  --outReadsUnmapped Fastx \
  --outSAMattributes Standard \
  --outFilterMultimapNmax 100 \
  --outFilterMismatchNoverLmax 0.04 \
  --alignIntronMax 1
