#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH -c 6
#SBATCH --mem-per-cpu=12G

# activate the env
module load Miniforge3
conda activate /scratch/mcarazo/envs/salmon

mkdir h_salmon_EMN_45S_rRNA_k15
mkdir h_salmon_primary_tRNA_k15

salmon index -k 15 --gencode -t EMN_45S_rRNA.fa.gz -i h_salmon_EMN_45S_rRNA_k15
#EMN_45S_rRNA.fa can be found in references directory of this repo
salmon index -k 15 --gencode -t gencode.v48.transcripts.fa.gz -i h_salmon_k15 
#Gencode.primary_tRNA_transcripts.fa is taken from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.transcripts.fa.gz
salmon index -k 15 -t hg38-tRNAs.fa -i GtRNAdb_salmon_k15
#hg38-tRNAs.fa can be found in https://gtrnadb.org/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz
salmon index -k 15 -t snoDB_FASTA.fa -i snoRNA_salmon_k15
#snoDB_FASTA.fa can be found in https://bioinfo-scottgroup.med.usherbrooke.ca/snoDB/download_all/
salmon index -k 15 -t lncRNA_LncBookv2.1_GRCh38.fa -i lncRNABook_salmon_k15
#lncRNA_LncBookv2.1_GRCh38.fa can be found in https://ngdc.cncb.ac.cn/lncbook/files/lncRNA_LncBookv2.1_GRCh38.fa.gz
salmon index -k 15 -t hairpin_hsa_dna.fa -i miRNA_hairpin_salmon_k15
#hairpin_hsa_dna.fa can be found in https://www.mirbase.org/download/hairpin.fa -> only hsa sequences
salmon index -k 15 -t mature_hsa_dna.fa -i miRNA_mature_salmon_k15
#mature_hsa_dna.fa can be found in https://www.mirbase.org/download/mature.fa -> only hsa sequences
