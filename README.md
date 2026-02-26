# Salmon2edgeR-Toolkit

**From Salmon Quantification to Biologically Interpretable Differential Expression**

Salmon2edgeR-Toolkit is a graphical and command-line compatible toolkit designed to streamline RNA-seq post-processing workflows based on **Salmon transcript quantification** and **edgeR differential expression analysis**.

The toolkit integrates gene-level aggregation, quality control, and biologically contextualized visualization in a reproducible and accessible framework.

---

## 🎯 Objective

Salmon2edgeR-Toolkit was developed to automate and standardize the downstream processing of RNA-seq data quantified with **Salmon**, enabling rapid preparation of gene-level matrices and biologically informed interpretation of **edgeR** differential expression results.

The toolkit bridges raw transcript-level quantification and publication-ready biological visualization.

---

## 🧬 Methodology

The toolkit supports a streamlined post-Salmon RNA-seq workflow:

1. **Transcript-to-Gene Aggregation**
   - Uses FASTA headers containing `[GeneID=XXXX]` to map transcripts to genes.
   - Aggregates transcript-level TPM and counts into gene-level matrices.

2. **Generation of edgeR-Ready Matrices**
   - Automatically prepares count matrices formatted for edgeR.
   - Supports 3–5 biological replicates per condition.

3. **Quality Control Module**
   - Log2(TPM+1) transformation  
   - Standardized PCA (PC1–PC2)  
   - Pearson correlation heatmaps  
   - Automatic outlier detection based on PCA mean distance  
   - Identification of highly discrepant genes  

4. **edgeR Result Analysis**
   - Direction-aware interpretation of log2FC (log2(A/B))  
   - Strict filtering by FDR and |logFC|  
   - Functional categorization via keyword-based biological classes  
   - Divergent logFC plots by category  
   - Separate Top-N A-up and B-up visualizations  
   - Publication-ready heatmaps (logCPM and row-wise z-score)  

> The toolkit assumes that Salmon quantification and edgeR differential expression analysis have already been performed externally.

---

## 🧰 Features

The toolkit includes three main modules:

### 1. Salmon Gene Mapper

- Converts Salmon transcript-level output to gene-level matrices  
- Supports optional gene annotation tables (NCBI GeneID / Symbol / Description)  
- Generates:
  - Wide gene-level tables  
  - edgeR-ready count matrices  
  - TPM tables  
  - TPM normalized by actin genes  

---

### 2. Salmon QC and Outlier Detection

- Performs PCA on log2(TPM+1) gene-level data  
- Generates correlation matrices  
- Detects statistical outliers automatically  
- Identifies genes contributing most to sample deviation  
- Exports:
  - PCA plots  
  - Correlation heatmaps  
  - Outlier reports  
  - Annotated gene lists  

---

### 3. edgeR Analyzer

- Imports edgeR results and corresponding count matrices  
- Automatically detects 3–5 replicates per condition  
- Applies strict significance filtering (FDR + |logFC|)  
- Assigns genes to biological categories using regex-based rules  
- Generates:
  - Divergent logFC plots by functional category  
  - Top-N A-up and B-up barplots  
  - Heatmaps (logCPM and z-score)  
  - Summary tables by category and direction  

---

## 🖥️ Requirements

- Python 3.9+

### Dependencies

- pandas  
- numpy  
- matplotlib  
- scikit-learn  
- tkinter  

Install via:

```bash
pip install pandas numpy matplotlib scikit-learn
```

---

## ▶️ How to Run

Each module can be executed independently.

### Run the Gene Mapper

```bash
python salmon_gene_mapper.py
```

### Run the QC & Outlier Module

```bash
python salmon_qc_gui_cli.py
```

CLI mode example:

```bash
python salmon_qc_gui_cli.py run \
  --inputs A.tabular B.tabular C.tabular D.tabular \
  --names A B C D \
  --fasta reference.fa \
  --gene_table genes.tsv \
  --outdir qc_out \
  --top_outlier_genes 50 \
  --min_tpm 1
```

### Run the edgeR Analyzer

```bash
python edger_analyzer.py
```

---

## 🔁 Workflow Overview

Salmon → Gene-level aggregation → QC & PCA → edgeR → Functional categorization → Biological visualization outputs

---

## 📊 Outputs

Depending on the module used, the toolkit generates:

- Gene-level count matrices
- edgeR-ready formatted tables
- TPM tables
- TPM normalized by actin genes
- PCA plots (PNG/PDF)
- Correlation heatmaps
- Outlier reports
- Annotated top discrepant genes
- Category-based divergent logFC plots
- Top-N A-up and B-up plots
- Publication-ready heatmaps
- Summary tables by category and regulation direction
- Reproducibility files (`run_manifest.json`, `reproduce_command.txt`)

---

## 🔬 Reproducibility

The QC module automatically records:

- Execution parameters  
- Input file paths  
- Python and library versions  
- Reproducible CLI command  

This ensures transparency and experiment traceability.

---

## 📜 Citation

If you use Salmon2edgeR-Toolkit in your research, please cite:

> Rocha, D. (2026). *Salmon2edgeR-Toolkit: From Salmon Quantification to Biologically Interpretable Differential Expression* (v1.0) [Computer software]. Zenodo. DOI: XXXXX

---

## 📚 References

- Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. *Nature Methods*, 14, 417–419.  
- Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*, 26(1), 139–140.  

---

## 👨‍💻 Author

**Prof. Danilo Massuia Rocha**  
ORCID: https://orcid.org/0000-0003-0059-7962  
Email: danilo.rocha@unesp.br  

Department of Biology  
School of Agricultural and Veterinary Sciences – FCAV  
São Paulo State University (UNESP)  
Jaboticabal, São Paulo, Brazil  
