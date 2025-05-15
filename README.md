# SCG Cardiac Disease Reanalysis

This repository contains all code, Google Colab notebooks, and R scripts for the reanalysis of publicly available human RNA-seq data from ENA Project PRJNA967653.

## 🔍 Project Summary

This project presents an independent, exploratory reanalysis of human SCG bulk RNA-seq data originally published by Ziegler et al. (2023). By reprocessing and analyzing the data with updated bioinformatics tools, we aim to uncover novel transcriptomic signatures relevant to cardiac disease.

Key findings include:

- Upregulation of inflammatory and immune response genes (e.g., CXCL10, TREM2, MHCs)
- Evidence of metabolic reprogramming (e.g., oxidative phosphorylation, purine metabolism)
- Enrichment of neuroimmune and antigen presentation pathways

## 🧪 Methods Overview

- Data download and QC via FastQC (Galaxy)
- Transcript quantification using **kallisto**
- Differential expression via **DESeq2**
- Functional enrichment using **clusterProfiler**, **fgsea**
- Visualization with **ggplot2**, **EnhancedVolcano**, **pheatmap** & more

All code is written in **R (v4.4.2)** and executed in **Google Colab** and local R environments.

## 📁 Folder Structure

/FastQCs/ # Contains generated FASTQC files
/kallisto/ # Contains raw kallisto output
/Plots/ # Contans generated plots (volcano, heatmaps, etc.)
/Supplementary/ # Contains Extra files like DEGs tables, and functional enrichment analysis results

## 📬 Contact

For questions or collaborations:

- GitHub: [@yourusername]
- Email: [djosergenomics@gmail.com](mailto:djosergenomics@gmail.com)
- Website: [Djoser Genomics](https://djosergenomics.github.io)

---

**Note**: This is a reanalysis project independent of the original authors. See their publication for the primary findings and dataset generation.
