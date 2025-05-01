# BMP4_BMP9_adipocyte_signaling
# Distinct Signaling Dynamics of BMP4 and BMP9 in Brown versus White Adipocytes

**Authors**: Ioannis Kamzolas et al.
**Corresponding Author**: [ik352@cam.ac.uk]  
**BioStudies Accession**: [E-MTAB-14609](https://www.ebi.ac.uk/biostudies/studies/E-MTAB-14609)

---

## 🧬 Overview

This repository contains the R scripts and Dataset used by or produced by the study:
> **Title**: Distinct Signaling Dynamics of BMP4 and BMP9 in Brown versus White Adipocytes  


---

## 🗃️ Dataset

RNA-seq data (fastq files and processed files) for this study have been deposited in the BioStudies database:

- **Accession**: [E-MTAB-14609](https://www.ebi.ac.uk/biostudies/studies/E-MTAB-14609)

---

## 🔧 Methodology

### Sample Preparation
- 4 batches of **brown** and **white** pre-adipocytes were differentiated.
- Treated with **3 ng/ml BMP4 or BMP9** for 8 hours + vehicle controls.

### RNA-seq
- RNA extracted with RNeasy Plus Mini Kit (Qiagen)
- Libraries prepared and sequenced by **Novogene** on an **Illumina NovaSeq 6000**.

### Analysis Pipeline
1. **Quality Control**: `FASTQC v0.11.9`
2. **Alignment**: `Hisat2 v2.1.0` to `GRCm38` genome
3. **Counting**: `HTSeq v0.11.1`
4. **Gene Annotation**: `biomaRt`
5. **Normalization & Batch Correction**: `quantile normalization + ComBat (sva)`
6. **Differential Expression**: `DESeq2 v1.26.0` + `Benjamini-Hochberg FDR correction`
7. **Pathway Enrichment**: `FGSEA`, `EnrichR`
8. **Transcription Factor Activity**: `VIPER`

---

## 📁 Repository Contents

| Path | Description |
|------|-------------|
| `Data/` | Metadata and count matrix |
| `Scripts/` | Full analysis pipeline scripts in R |
| `Results/` | Output DEGs, enrichment results, and relevant results producing the manuscript figures |
| `Source_Files/` | Files used to run the scripts |

---

## 🖥️ How to Run

Copy the directory to a local folder and run the R scripts
