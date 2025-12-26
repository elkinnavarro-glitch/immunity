# Figure Reproducibility Map

This document provides a comprehensive guide to reproducing all figures in the immunity manuscript. Each figure is generated from dedicated scripts with clearly specified inputs and outputs.

## Figure 1: Cell Atlas and Lineage Composition

**Description**: Comprehensive single-cell atlas visualization and cell lineage proportions across tissue types.

- **Script**: `scripts/generate_figure1_atlas.py`
- **Input**: Processed scRNA-seq matrices (see `data/README.md` for details)
- **Output**: `figures/Figure1.pdf`
- **Dependencies**: scanpy, matplotlib, numpy, seaborn
- **Runtime**: ~5-10 minutes
- **Key Components**:
  - Panel A: UMAP projection of full dataset
  - Panel B: UMAP colored by cell lineage
  - Panel C: Barplot of lineage proportions
  - Panel D: Quality control metrics (ASW, silhouette)

---

## Figure 2: NMF-Based Tissue Programs

**Description**: Non-negative matrix factorization analysis revealing conserved tissue programs across lineages.

- **Script**: `scripts/generate_figure2_nmf.py`
- **Input**: Lineage-specific processed matrices, NMF model fits
- **Output**: `figures/Figure2.pdf`
- **Dependencies**: scikit-learn, matplotlib, numpy, pandas
- **Runtime**: ~15-20 minutes
- **Key Components**:
  - Panel A: Top genes - Fibroblast programs
  - Panel B: Top genes - Immune programs
  - Panel C: Conservation scores across lineages
  - Panel D: Pathway enrichment for selected programs

---

## Figure 3: Critical Dynamics and System Coupling

**Description**: Analysis of critical dynamics in immune cell interactions and tissue-scale coupling.

- **Script**: `scripts/generate_figure3_coupling.py`
- **Input**: TCR repertoire data, cell interaction matrices, temporal dynamics
- **Output**: `figures/Figure3.pdf`
- **Dependencies**: igraph, networkx, matplotlib, numpy, scipy
- **Runtime**: ~10-15 minutes
- **Key Components**:
  - Panel A: High-dimensional UMAP of immune cells
  - Panel B: Network topology of TCR-cognate pairs
  - Panel C: Scaling analysis (power-law relationships)
  - Panel D: Phase transition signatures in immune activation

---

## Running the Complete Pipeline

To reproduce all figures in sequence:

```bash
cd scripts/
python generate_figure1_atlas.py
python generate_figure2_nmf.py
python generate_figure3_coupling.py
```

All PDFs will be generated in the `figures/` directory.

---

## System Requirements

- Python 3.9+
- 16GB RAM (minimum for full dataset processing)
- ~5GB disk space for intermediate outputs
- See `requirements.txt` for complete dependency list

---

## Troubleshooting

- **Memory errors**: Ensure 16GB+ RAM available; consider processing by lineage separately
- **Missing dependencies**: Run `pip install -r requirements.txt`
- **Path errors**: Ensure scripts are run from the `scripts/` directory

---

## Citation

If you use these figures, please cite:

> Navarro et al. (2025). "Critical dynamics and immune system organization." [Journal]. doi:10.xxxx/xxxxx
