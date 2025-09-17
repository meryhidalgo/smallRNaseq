# 🧬 smallRNaseq Pipeline

Este repositorio contiene un conjunto de scripts para el análisis de datos de **small RNA-seq**, con especial énfasis en el preprocesamiento de datos brutos, alineamiento, cuantificación, análisis diferencial y exploraciones posteriores como PCA, deconvolución, enriquecimiento funcional y análisis de supervivencia.  

El pipeline está diseñado para ejecutarse parcialmente en clústeres HPC (High Performance Computing), especialmente los pasos escritos en **Bash**, y se complementa con análisis downstream en **R** y **Python**.  

---

## 📁 Estructura del repositorio
```text
.
├── envs/                        # Entornos conda necesarios
│   └── *.yml
├── references/
│   ├── EMN_45S_rRNA.gtf
│   ├── EMN_45S_rRNA.fa / .fai
│   ├── sRNA-adapters.fa
│   └── Ruiz-Moreno_et-al_2022.xlsx
├── reports/                     
│   ├── multiqc_report_*.html
│   ├── multiqc_general_stats.txt
│   ├── fastqc_length_distribution.txt
│   └── trimmingTools/*.tsv
├── src/
│   ├── HPC_analysis/            # Scripts principales a ejecutar en HPC
│   │   ├── 1-QC/                # Control de calidad inicial (FastQC)
│   │   │   └── 1_QC_fastqc.sh
│   │   ├── 2-trimming/          # Recorte de adaptadores 
│   │   │   ├── 2-Trimming_trim_galore.sh
│   │   │   ├── 2-Trimming_trimmomatic.sh
│   │   │   └── 2-Trimming_fastp.sh
│   │   ├── 3-alignment/         # Alineamiento con STAR
│   │   │   ├── A-generate_STAR_indexes.sh
│   │   │   ├── B-alignment_STAR_rRNA.sh
│   │   │   └── B-alignment_STAR_GRCh38.sh
│   │   ├── 4-quantification/    # Cuantificación con Salmon
│   │   │   ├── A-salmon_index.sh
│   │   │   ├── B-salmon.sh
│   │   │   ├── C-salmon_unmapped.sh
│   │   │   ├── D-create_TPM_table.sh
│   │   │   ├── E-mapping_rates.sh
│   │   │   └── GetSalmonTPMs_mch.R
│   │   └── 5-reads_distribution/
│   │       └── process_depth.sh
│   └── Downstream_analysis/     # (R/Python, notebooks, visualización)
│       ├── DGE/
│       │   ├── Salmon_DEG_analysis.Rmd
│       │   └── Gene_biotype_DEG_gencode.Rmd
│       ├── GO/
│       │   └── GO_clusterprofiler.Rmd
│       ├── Figure_3/adapter_plot.ipynb
│       ├── Figure_4/General_stats.ipynb
│       ├── Figure_5/
│       │   ├── Gene_biotype_refseq.Rmd
│       │   ├── Gene_biotype_gencode.Rmd
│       │   └── plot_mappingRate_diff.ipynb
│       ├── Figure_7/plot_depth.ipynb
│       ├── Figure_12/miRNA_interaction_network.Rmd 
│       ├── Figure_16/Deconvolution.Rmd
│       └── Figure_17/
│           ├── Forest_plot.Rmd
│           └── survival_data.Rmd
```

---

## 🔁 Flujo de trabajo

1. ### Control de calidad inicial
   - Script: `src/HPC_analysis/1-QC/1_QC_fastqc.sh`  
   - Outputs: `reports/fastqc_length_distribution.txt` `reports/multiqc_general_stats.txt`
   - Visualización: `src/Downstream_analysis/Figure_4/General_stats.ipynb`  

2. ### Recorte de adaptadores
   - Scripts:
     - `src/HPC_analysis/2-trimming/2-Trimming_trim_galore.sh`
     - `src/HPC_analysis/2-trimming/2-Trimming_trimmomatic.sh`
     - `src/HPC_analysis/2-trimming/2-Trimming_fastp.sh`
   - Adaptadores: `references/sRNA-adapters.fa`
   - Output: `reports/trimmingTools/*.tsv`
   - Visualización: `src/Downstream_analysis/Figure_3/adapter_plot.ipynb` 

3. ### Alineamiento con STAR
   - Index script: `src/HPC_analysis/3-alignment/A-generate_STAR_indexes.sh`
   - Alineamiento contra rRNA: `src/HPC_analysis/3-alignment/B-alignment_STAR_rRNA.sh`
   - Alineamiento contra genoma humano (GRCh38): `src/HPC_analysis/3-alignment/B-alignment_STAR_GRCh38.sh`

4. ### Cuantificación con Salmon
   - Index script: `src/HPC_analysis/4-quantification/A-salmon_index.sh`
   - Cuantificación: `src/HPC_analysis/4-quantification/B-salmon.sh`
   - Cuantificación reads no mapeados: `src/HPC_analysis/4-quantification/C-salmon_unmapped.sh`
   - Generación de tablas TPM: `src/HPC_analysis/4-quantification/D-create_TPM_table.sh` + `GetSalmonTPMs_mch.R`
   - Porcentajes de mapeo: `src/HPC_analysis/4-quantification/E-mapping_rates.sh`
   - Visualización: `src/Downstream_analysis/Figure_5/plot_mappingRate_diff.ipynb`


5. ### Distribución de lecturas
   - Script: `src/HPC_analysis/5-reads_distribution/process_depth.sh`
   - Visualización: `src/Downstream_analysis/Figure_7/plot_depth.ipynb`

6. ### Análisis expresión diferencial
   - Script: `Downstream_analysis/DGE/Salmon_DEG_analysis.Rmd`, `Gene_biotype_DEG_gencode.Rmd`
   - Enriquecimiento funcional (GO): `Downstream_analysis/GO/GO_clusterprofiler.Rmd`
   - Red de interacción miRNA-target: `Downstream_analysis/Figure_12/miRNA_interaction_network.Rmd`. Se apoya en la base de datos miRDB. 
   - Estudio de deconvolución: `Downstream_analysis/Figure_16/Deconvolution.Rmd`. Necesita del archivo de anotaciones `reports/Ruiz-Moreno_et-al_2022.xlsx`
   - Análisis de supervivencia: `Downstream_analysis/Figure_17/survival_data.Rmd + Forest_plot.Rmd`

---

## 📊 Informes

Los informes de calidad se encuentran en la carpeta `reports/` e incluyen salidas de MultiQC y métricas derivadas de diferentes herramientas de trimming y alineamiento:

- `multiqc_report_fastp.html`
- `multiqc_report_trim_galore.html`
- `multiqc_report_trimmomatic.html`
- `multiqc_report_fastp_merged.html`
- `multiqc_general_stats.txt`
- `fastqc_length_distribution.txt`

---

## ⚙️ Dependencias

Este pipeline utiliza scripts en **Bash**, **R** y **Python**. Algunas herramientas clave utilizadas incluyen:

- `FastQC`, `Trim Galore`, `Trimmomatic`, `Fastp`
- `STAR`, `Salmon`, `MultiQC`
- `R` con paquetes para análisis de expresión, enriquecimiento y supervivencia (`DESeq2`, `ggplot2`, `clusterProfiler`, etc.)

Las dependencias están organizadas en archivos de entorno (`.yml`) dentro de la carpeta `./envs/`. Para lanzar los scripts del clúster HPC se necesitará la instalación y activación del correspondiebte. Esto garantiza la reproducibilidad en clústeres HPC.

### 📂 Entornos disponibles
- `./envs/QControl.yml` → Control de calidad (FastQC, MultiQC, etc.)
- `./envs/STAR.yml` → Alineamiento con STAR
- `./envs/salmon.yml` → Cuantificación con Salmon
- `./envs/DGE_analysis.yml` → Análisis de expresión diferencial (R con DESeq2, ggplot2, etc.)
- `./envs/biotools.yml` → Herramientas generales de bioinformática. Incluye herramientas de trimming, samtools, etc. 

### ▶️ Creación de entornos

Para crear un entorno:

```bash
conda env create -f ./envs/NOMBRE.yml
conda activate NOMBRE
```

---

## 🧑‍💻 Autora

**Maria Carazo Hidalgo**  
Bioinformática  
Correo: mariacarazohidalgo@gmail.com
