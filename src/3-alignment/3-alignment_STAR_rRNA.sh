#!/bin/bash

#SBATCH -t 01:00:00
#SBATCH -c 6
#SBATCH --partition=biogipuzkoa-exclusive
#SBATCH --account=biogipuzkoa-exclusive
#SBATCH --mem=64G
#SBATCH -o STAR_R_%A_%a.out
#SBATCH --job-name=STAR
#SBATCH --array=0-5%6

module load Miniforge3
conda activate /scratch/mcarazo/envs/STAR/

# create a directory for STAR analysis
cd /scratch/mcarazo/ongoing/GB/sRNAseq_may25
mkdir -p star_out_EMN_45S_rRNA

index=/scratch/mcarazo/indexes/index_rRNA
index=/scratch/mcarazo/indexes/rRNA/index_mch_main_rRNA #star_out_main_rRNA
#index=/scratch/mcarazo/indexes/rRNA/index_mch_45s_rRNA #star_out_45s_rRNA
index=/scratch/mcarazo/indexes/rRNA/index_EMN_45S_rRNA #star_out_EMN_45S_rRNA

# Select the fastq files
fastqs=($(realpath fastp_out/other/*_R1.fastq.gz))

fq1=${fastqs[$SLURM_ARRAY_TASK_ID]}
fq2=${fq1/_R1/_R2}
base=$(basename ${fq1} _R1.fastq.gz)

#parameters used for STAR MCH alignment:
STAR --genomeDir "$index" \
  --runThreadN 6 \
  --readFilesIn "$fq1" "$fq2" \
  --readFilesCommand zcat \
  --outFileNamePrefix "/scratch/mcarazo/ongoing/GB/sRNAseq_may25/star_out_EMN_45S_rRNA/${base}_" \
  --outSAMtype BAM SortedByCoordinate \
  --outReadsUnmapped Fastx \
  --outFilterMultimapNmax 100 \
  --outFilterMismatchNoverLmax 0.04 \
  --alignIntronMax 1 \
  --limitBAMsortRAM 3307885652



#  --outReadsUnmapped Fastx guarda reads no mapeados como fastq aparte
#  --outSAMunmapped Within  guarda reads no mapeados en el mismo archivo bam
#  --alignIntronMax is maximum intron size. recommended for sRNA to use 1 so NO splicing would be considered
#  --outFilterMultimapNmax is the maximum number of multiple mapping allowed. Default is 20.
#  --outFilterMismatchNmax / outFilterMismatchNoverLmax is the maximum number / proportion of mismatches allowed
#  --outFilterScoreMinOverLread 0.33 exige que la puntuación mínima sea ≥ 33% del largo del read


#https://groups.google.com/g/rna-star/c/RBWvAGFooMU/m/HAOWHxcOCgAJ
#--outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterScoreMinOverLread 0  --outFilterMatchNminOverLread 0 --alignIntronMax 1