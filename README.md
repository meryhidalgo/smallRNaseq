# 🧬 smallRNaseq Pipeline

Este repositorio contiene un conjunto de scripts para el análisis de datos de **small RNA-seq**, con especial énfasis en el preprocesamiento de datos brutos, alineamiento, cuantificación, análisis diferencial y exploraciones posteriores como PCA y análisis de supervivencia.

El pipeline está diseñado para ejecutarse parcialmente en clústeres HPC (High Performance Computing), especialmente los pasos escritos en **Bash**, y se complementa con análisis downstream en **R** y **Python**.

---

## 📁 Estructura del repositorio
.
├── design.tab # Tabla de diseño experimental
├── references/ # Archivos de referencia (rRNA, índices, etc.)
├── reports/ # Informes de calidad (MultiQC)
├── src/ # Scripts del pipeline organizados por pasos
│ ├── 1-QC/ # Control de calidad inicial (FastQC)
│ ├── 2-trimming/ # Recorte de adaptadores (Trim Galore, Trimmomatic, Fastp)
│ ├── 3-alignment/ # Alineamiento (STAR)
│ ├── 4-quantification/ # Cuantificación (Salmon)
│ ├── 5-DGE/ # Análisis de expresión diferencial, PCA
│ └── 6-survival/ # Análisis de supervivencia


---

## 🔁 Flujo de trabajo

1. ### Control de calidad inicial
   - `src/1-QC/1_QC_fastqc.sh`  
   - Visualización de métricas generales: `1-general_stats.ipynb`

2. ### Recorte de adaptadores
   - `src/2-trimming/2-Trimming_trim_galore.sh`
   - `src/2-trimming/2-Trimming_trimmomatic.sh`
   - `src/2-trimming/2-Trimming_fastp.sh`
   - Archivos de adaptadores: `sRNA-adapters.fa`

3. ### Alineamiento a referencias
   - Alineamiento contra rRNA: `src/3-alignment/3-alignment_STAR_rRNA.sh`
   - Alineamiento contra el genoma: `src/3-alignment/3-alignment_STAR_GRCh38.sh`

4. ### Cuantificación con Salmon
   - Indexado: `src/4-quantification/salmon_index.sh`
   - Cuantificación: `src/4-quantification/salmon.sh`, `salmon_after_rRNA.sh`
   - Generación de tablas: `create_TPM_table.sh`, `create_TPM_table_rRNA.sh`

5. ### Análisis de expresión diferencial (DGE)
   - PCA: `GetSalmon_PCA.Rmd`, `GetSalmon_PCA_rRNA.Rmd`
   - Detección por biotipo: `Gene_biotype_DEG_gencode_tRNA.Rmd`

6. ### Análisis de supervivencia
   - Script R Markdown: `src/6-survival/survival_data.Rmd`

---

## 📊 Informes

Los informes de calidad se encuentran en la carpeta `reports/` e incluyen salidas de MultiQC para diferentes herramientas de trimming:

- `multiqc_report_fastp.html`
- `multiqc_report_trim_galore.html`
- `multiqc_report_trimmomatic.html`
- `multiqc_report_fastp_merged.html`

---

## ⚙️ Dependencias

Este pipeline utiliza scripts en **Bash**, **R** y **Python**. Algunas herramientas clave utilizadas incluyen:

- `FastQC`, `Trim Galore`, `Trimmomatic`, `Fastp`
- `STAR`, `Salmon`, `MultiQC`
- `R` con paquetes para análisis de expresión y supervivencia (`DESeq2`, `ggplot2`, etc.)

> Las dependencias detalladas y cómo instalarlas se incluirán en una futura actualización.

---

## 📌 Notas

- Este repositorio está diseñado para ejecutarse principalmente en clústeres HPC con sistema de colas.
- La tabla de diseño experimental se encuentra en `design.tab`.
- Los archivos de referencia están en la carpeta `references/`.

---

## 🧑‍💻 Autora

**Mery Hidalgo**  
Bioinformática 
Correo: MARIA.CARAZOHIDALGO@bio-gipuzkoa.eus


