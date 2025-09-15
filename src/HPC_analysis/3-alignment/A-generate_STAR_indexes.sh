#!/bin/bash

#SBATCH -c 6
#SBATCH --mem-per-cpu=12G

#### FOR HUMAN ####
module load Miniforge3
conda activate /scratch/mcarazo/envs/STAR/

cd /scratch/mcarazo/indexes/
mkdir index_h99

STAR --runThreadN 20 \
--runMode genomeGenerate \
--genomeDir /scratch/mcarazo/indexes/index_h99 \
--genomeFastaFiles /data/mcarazo/indexes/GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile /data/mcarazo/indexes/gencode.v41.primary_assembly.annotation.gtf \
--sjdbOverhang 99

cd /scratch/mcarazo/indexes/rRNA/
mkdir index_EMN_45S_rRNA

STAR --runThreadN 20 \
--runMode genomeGenerate \
--genomeDir /scratch/mcarazo/indexes/rRNA/index_EMN_45S_rRNA \
--genomeFastaFiles /scratch/mcarazo/indexes/rRNA/EMN_45S_rRNA.fa \
--genomeSAindexNbases 6

#EMN_45S_rRNA.fa can be found in references directory of this repo