# 🧬 smallRNaseq Pipeline

Este repositorio contiene un conjunto de scripts para el análisis de datos de **small RNA-seq**, con especial énfasis en el preprocesamiento de datos brutos, alineamiento, cuantificación, análisis diferencial y exploraciones posteriores como PCA, deconvolución, enriquecimiento funcional y análisis de supervivencia.  

El pipeline está diseñado para ejecutarse parcialmente en clústeres HPC (High Performance Computing), especialmente los pasos escritos en **Bash**, y se complementa con análisis downstream en **R** y **Python**.  

---

## 📁 Estructura del repositorio
```text
.
├── design.tab                   # Tabla de diseño experimental
├── references/                  # Archivos de referencia (rRNA, índices, adaptadores, artículos previos)
│   ├── EMN_45S_rRNA.gtf
│   ├── EMN_45S_rRNA.fa / .fai
│   ├── sRNA-adapters.fa
│   └── Ruiz-Moreno_et-al_2022.xlsx
├── reports/                     # Informes de calidad y métricas de trimming/alineamiento
│   ├── multiqc_report_*.html
│   ├── multiqc_general_stats.txt
│   ├── fastqc_length_distribution.txt
│   └── trimmingTools/*.tsv
├── src/
│   ├── HPC_analysis/            # Scripts principales a ejecutar en HPC
│   │   ├── 1-QC/                # Control de calidad inicial (FastQC)
│   │   │   └── 1_QC_fastqc.sh
│   │   ├── 2-trimming/          # Recorte de adaptadores (Trim Galore, Trimmomatic, Fastp)
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
│   └── Downstream_analysis/     # Análisis posteriores (R/Python, notebooks, visualización)
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
   - `src/HPC_analysis/1-QC/1_QC_fastqc.sh`  
   - Visualización: `reports/fastqc_length_distribution.txt`

2. ### Recorte de adaptadores
   - `src/HPC_analysis/2-trimming/2-Trimming_trim_galore.sh`
   - `src/HPC_analysis/2-trimming/2-Trimming_trimmomatic.sh`
   - `src/HPC_analysis/2-trimming/2-Trimming_fastp.sh`
   - Adaptadores: `references/sRNA-adapters.fa`

3. ### Alineamiento con STAR
   - Indexado: `src/HPC_analysis/3-alignment/A-generate_STAR_indexes.sh`
   - Alineamiento contra rRNA: `src/HPC_analysis/3-alignment/B-alignment_STAR_rRNA.sh`
   - Alineamiento contra genoma humano (GRCh38): `src/HPC_analysis/3-alignment/B-alignment_STAR_GRCh38.sh`

4. ### Cuantificación con Salmon
   - Indexado: `src/HPC_analysis/4-quantification/A-salmon_index.sh`
   - Cuantificación: `src/HPC_analysis/4-quantification/B-salmon.sh`
   - Cuantificación reads no mapeados: `src/HPC_analysis/4-quantification/C-salmon_unmapped.sh`
   - Generación de tablas TPM: `src/HPC_analysis/4-quantification/D-create_TPM_table.sh`
   - Tasas de mapeo: `src/HPC_analysis/4-quantification/E-mapping_rates.sh`
   - Scripts R complementarios: `GetSalmonTPMs_mch.R`

5. ### Distribución de lecturas
   - `src/HPC_analysis/5-reads_distribution/process_depth.sh`

6. ### Análisis downstream
   - Expresión diferencial: `Downstream_analysis/DGE/Salmon_DEG_analysis.Rmd`, `Gene_biotype_DEG_gencode.Rmd`
   - Enriquecimiento funcional (GO): `Downstream_analysis/GO/GO_clusterprofiler.Rmd`
   - Visualizaciones: `Downstream_analysis/Figure_*`
     - PCA, profundidad, tasas de mapeo, interacción miRNA, deconvolución, forest plot, etc.
   - Análisis de supervivencia: `Downstream_analysis/Figure_17/survival_data.Rmd`

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

> Las dependencias detalladas y cómo instalarlas se incluirán en una futura actualización.

---

## 📌 Notas

- Este repositorio está diseñado para ejecutarse principalmente en clústeres HPC con sistema de colas.
- La tabla de diseño experimental se encuentra en `design.tab`.
- Los archivos de referencia están en la carpeta `references/`.

---

## 🧑‍💻 Autora

**Maria Carazo Hidalgo**  
Bioinformática  
Correo: MARIA.CARAZOHIDALGO@bio-gipuzkoa.eus
