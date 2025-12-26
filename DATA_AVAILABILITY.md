# Data Availability

This document describes the availability and accessibility of datasets used in the immunity manuscript. All raw and processed data are publicly available through established scientific repositories to ensure reproducibility and open science.

## Public Data Sources

### 1. Reynolds et al., Science 2021

**Study**: Single-cell and spatial immune profiling of human tissues

**Reference**: Reynolds et al. (2021). "Immunological and transcriptomic profiling reveals shared molecular and physical hubs of inflammation." *Science*, 372(6545).

**Datasets Available**:
- Raw scRNA-seq counts (10x Genomics)
- Processed expression matrices
- Cell type annotations
- Spatial transcriptomics data

**Access**:
- **GEO Accession**: GSE167107
- **URL**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167107
- **Download**: All files available through GEO FTP

**Data Types**:
- Human skin tissue samples
- Human blood immune cells
- Multiple tissue compartments
- scRNA-seq (10x v2, v3)

---

### 2. Ma et al., Nature Communications 2020

**Study**: Genome-wide association analysis and meta-analysis of human skin transcriptome

**Reference**: Ma et al. (2020). "Skin transcriptomics reveals novel fibrin-based cell state transitions in human fibroblasts." *Nature Communications*, 11, 5302.

**Datasets Available**:
- Processed scRNA-seq matrices (AnnData format)
- Cell type annotations
- Cell lineage classifications
- Quality control metrics

**Access**:
- **GEO Accession**: GSE147424
- **URL**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147424
- **Download**: H5AD files available through repository

**Data Types**:
- Human dermal tissue
- Fibroblasts, immune cells, keratinocytes
- Single-cell transcriptomics
- Cell-type-specific signatures

---

## Data Processing and Integration

This repository provides all scripts necessary to:

1. **Download** raw data from GEO/public repositories
2. **Preprocess** count matrices (quality control, normalization)
3. **Integrate** multiple datasets using batch correction algorithms
4. **Annotate** cells using reference-based and unsupervised approaches
5. **Analyze** at tissue, lineage, and system levels

**Scripts Location**: `scripts/data_processing/`

---

## Data Privacy and Ethics

### No Patient-Identifiable Information

- ✅ All data are **fully anonymized**
- ✅ No personal health information (PHI) is present
- ✅ No genetic identifiers or ancestry information
- ✅ No temporal data that could enable re-identification
- ✅ All sequences are counts/expression, not genomic coordinates for individual variants

### Ethical Oversight

- All source data collection was performed under institutional review board (IRB) approval
- Original publications contain informed consent documentation
- Data usage complies with original study agreements

---

## Accessing the Data

### Option 1: Automated Download (Recommended)

```bash
# Clone this repository
git clone https://github.com/elkinnavarro-glitch/immunity.git
cd immunity

# Run data download script
python scripts/data_processing/download_public_data.py

# This script will:
# - Download GSE167107 (Reynolds et al.)
# - Download GSE147424 (Ma et al.)
# - Verify file integrity
# - Place data in data/raw/
```

### Option 2: Manual Download

1. Visit GEO using the accession numbers above
2. Download all supplementary files
3. Place in `data/raw/` directory
4. Run preprocessing scripts

### Option 3: Use Pre-Processed Data

If you only need the processed, integrated datasets:

```bash
# Available at data/processed/
ls data/processed/

# Contains:
# - integrated_matrix.h5ad (full dataset)
# - lineage_subset_*.h5ad (lineage-specific subsets)
# - metadata_table.csv (cell annotations)
```

---

## Data Storage Specifications

### Raw Data
- **Location**: `data/raw/`
- **Size**: ~50 GB (uncompressed)
- **Format**: H5AD, CSV, TSV
- **Compression**: gzip (`.gz`)

### Processed Data
- **Location**: `data/processed/`
- **Size**: ~5 GB (uncompressed)
- **Format**: H5AD (HDF5-based)
- **Requirements**: 16GB RAM for full dataset analysis

### Intermediate Files
- **Location**: `data/intermediate/`
- **Size**: Variable (depends on analysis)
- **Format**: H5AD, pickle, numpy arrays

---

## Using the Data

### Required Software
- Python 3.9+
- scanpy >= 1.8
- anndata >= 0.8
- numpy, scipy, pandas
- (see `requirements.txt`)

### Loading Data in Python

```python
import scanpy as sc

# Load integrated dataset
adata = sc.read_h5ad('data/processed/integrated_matrix.h5ad')

# Load lineage-specific subset
adata_fib = sc.read_h5ad('data/processed/lineage_subset_fibroblasts.h5ad')

# Access cell metadata
cell_types = adata.obs['cell_type'].value_counts()
lineages = adata.obs['cell_lineage'].value_counts()
```

---

## Data Citation

When using these datasets, please cite both the original publications:

**Reynolds et al. (2021)**:
```
Reynolds, G., et al. (2021). Immunological and transcriptomic profiling 
reveals shared molecular and physical hubs of inflammation. 
Science, 372(6545), eabf7570.
```

**Ma et al. (2020)**:
```
Ma, S., et al. (2020). Skin transcriptomics reveals novel fibrin-based 
cell state transitions in human fibroblasts. 
Nature Communications, 11, 5302.
```

**This Repository**:
```
Navarro, E., et al. (2025). Critical dynamics and immune system organization. 
[Journal and DOI to be added upon publication]
```

---

## Troubleshooting

### Large File Sizes
- Data files are large; ensure adequate storage (>60GB)
- Consider working with lineage-specific subsets for faster processing
- Use cloud computing if local resources are limited

### Download Failures
- Check internet connection stability
- Verify GEO is accessible (sometimes temporary downtime)
- Re-run download script (includes retry logic)

### Data Format Issues
- Ensure scanpy is properly installed: `pip install scanpy`
- Check Python version >= 3.9
- Verify HDF5 libraries are available

---

## Questions or Issues?

For questions about data access, contact:
- **Lead**: Elkin Navarro Quiroz
- **Affiliation**: Centro de Investigaciones en Ciencias de la Vida, Universidad Simón Bolívar
- **Email**: Available in repository contact information

---

**Last Updated**: December 2025
**Data Frozen**: Yes - No further modifications to source datasets
**License**: CC BY 4.0 (for this processing pipeline)
