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

# create a directory for salmon analysis
#cd /scratch/mcarazo/ongoing/GB/sRNAseq_may25/
cd /scratch/mcarazo/ongoing/GB/sRNAseq_may25/star_out_EMN_45S_rRNA
mkdir -p salmon_out_gencode_tRNA

#index=/scratch/mcarazo/indexes/h_salmon_refseq_k15
#index=/scratch/mcarazo/indexes/h_salmon_k15
index=/scratch/mcarazo/indexes/tRNA/h_salmon_primary_tRNA_k15 # after rRNA alignment with EMN + 45S -> star_out_EMN_45S_rRNA/salmon_out_gencode_tRNA
#index=/scratch/mcarazo/indexes/rRNA/h_salmon_rRNA_k15 -> salmon_out_rRNA
#index=/scratch/mcarazo/indexes/rRNA/h_salmon_mch_main_rRNA_k15 #-> salmon_out_main_mch_rRNA
#index=/scratch/mcarazo/indexes/rRNA/h_salmon_mch_45rRNA_k15 #-> salmon_out_45_mch_rRNA
#index=/scratch/mcarazo/indexes/rRNA/h_salmon_MCH_all_main_45S_rRNA_k15 #-> salmon_out_MCH_main_45S_rRNA
#index=/scratch/mcarazo/indexes/T2T/h_salmon_T2T_rRNA_k15 #-> salmon_out_T2T_rRNA

# Select the fastq files
fastqs=($(realpath *.mate1.fastq))


fq1=${fastqs[$SLURM_ARRAY_TASK_ID]}
fq2=${fq1/.mate1.fastq/.mate2.fastq}
base=$(basename ${fq1} _trim_Unmapped.out.mate1.fastq)

#mv ${fq1} ${fq1}.fastq
#mv ${fq2} ${fq2}.fastq

#fq1=${fq1}.fastq
#fq2=${fq2}.fastq

# quantify transcripts in fastq files
salmon quant -i $index \
    -l A \
	-1 $fq1 \
	-2 $fq2 \
	-p 6 \
	--validateMappings \
    -o salmon_out_gencode_tRNA/${base}
