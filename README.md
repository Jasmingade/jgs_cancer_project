# Cancer Multi-Omics Survival Pipeline

Reproducible bioinformatics workflows for integrating transcriptomics, proteomics, mutation-derived features, and clinical metadata in cancer survival analysis.

This repository contains the computational workflows used to preprocess and harmonise large-scale cancer omics datasets, link them to clinical survival endpoints, and run feature-wise survival modelling and downstream visualisation. The project is centered on integrating transcriptomic and proteomic layers in a cancer context, with a focus on gene-level and isoform-level analyses.

## Project focus

This repository demonstrates how heterogeneous molecular datasets can be transformed into analysis-ready resources for survival modelling and biological interpretation.

Core themes include:

- multi-omics data integration
- transcriptomics and proteomics preprocessing
- clinical metadata harmonisation
- gene- and isoform-level feature modelling
- survival analysis using Cox proportional hazards models
- reproducible workflows in R, shell, Git, and Linux/HPC environments

## What this repository demonstrates

This project reflects my profile in:

- **Transcriptomics:** preprocessing, feature handling, and modelling of RNA-derived features
- **Proteomics:** preprocessing and integration of CPTAC/PDC-derived proteomics data
- **Clinical data integration:** harmonising omics data with survival endpoints and covariates
- **Survival analysis:** univariate Cox proportional hazards modelling across feature types
- **Reproducible workflows:** structured analysis pipelines, environment management, and workflow organisation
- **Cancer data science:** integrating molecular features with clinically relevant outcomes

## Repository overview

The repository is organised into two main analysis domains:

### `01_transcriptomics/`
Workflows related to transcriptomic preprocessing, harmonisation, modelling, and downstream analysis.

Typical tasks include:
- preprocessing transcriptomic input data
- feature preparation at gene and/or isoform level
- alignment with clinical covariates and survival endpoints
- running univariate survival models
- preparing outputs for downstream interpretation and visualisation

### `02_proteomics/`
Workflows related to proteomics preprocessing, integration, modelling, and downstream analysis.

Typical tasks include:
- preprocessing CPTAC/PDC proteomics data
- harmonisation with clinical metadata
- preparation of protein- or isoform-related feature matrices
- survival modelling of proteomic features
- downstream visualisation and exploratory summaries

### Root-level files
- `environment.yml` – conda environment specification
- `r_installed_packages.csv` – record of installed R packages
- `__init__.py` – project-level Python package marker
- `.gitignore` – ignored files and directories
- `LICENSE` – repository license

## Workflow summary

At a high level, the project follows this logic:

1. **Data acquisition / input preparation**  
   Omics datasets and associated metadata are prepared for analysis.

2. **Clinical harmonisation**  
   Molecular data are linked to clinical annotations and survival endpoints.

3. **Feature engineering / preprocessing**  
   Gene-, isoform-, or protein-related features are filtered, transformed, and organised into analysis-ready matrices.

4. **Survival modelling**  
   Univariate Cox proportional hazards models are run across features.

5. **Result processing and visualisation**  
   Significant associations, effect sizes, and selected feature-level results are summarised and visualised.

## Data types covered

This repository is designed around cancer datasets spanning multiple molecular layers, including:

- transcriptomics
- proteomics
- mutation-derived features
- clinical metadata
- survival endpoints

The overall goal is to make these heterogeneous layers comparable and usable in a common analytical framework.

## Methods used

The repository includes workflows centered around:

- data harmonisation
- feature-level modelling
- survival analysis
- downstream survival visualisation
- reproducible environment setup

Key methodological components include:
- **Cox proportional hazards models**
- **gene- and isoform-level feature analysis**
- **clinical covariate integration**
- **multi-omics data structuring**

## Technical environment

This repository is primarily written in **R**, with supporting **shell** components for workflow execution and environment handling. GitHub currently identifies the language split as predominantly R with a smaller shell component. :contentReference[oaicite:2]{index=2}

The project is intended to be run in a reproducible computational environment using:

- R
- shell / bash
- Git
- Linux / HPC-compatible workflows
- conda environment management via `environment.yml`
