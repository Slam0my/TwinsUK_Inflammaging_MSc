# Identification of Early Instigators of Inflammaging Within the TwinsUK Cohort
- [Introduction](#introduction) 
- [Repository Structure](#repostructure)

- 
## <a id="introduction"></a> Introduction
This repository is part of an extended research project for the MSc Applied Bioinformatics course at King's College London. 



## <a id="repostructure"></a> Repository Structure 
All Bioinformatics pipeline was conducted on the CREATE HPC. Data, results, logs and software directories have been omitted from this GitHub Repository.
- `data/`
- Directory containing all of the data used in the workflow. This includes CSVs and a copy of the raw data. 
  - `raw_olink_ip/`
    - Contains the raw, original data. Untouched, READ-ONLY.
  - `associations_test/`
    - Contains generated CSVs of the linear mixed-effects model (LMM) results for each IP-Protein pair.
  - `protein_protein_solar/`
    - Contains list of protein-protein pairs, and respective collated SOLAR results outputs.
  - `solar_analysis/`
    - Contains the SOLAR pedigree and phenotype CSVs, list of IP-Protein pairs, and respective collated SOLAR results outputs.
  - olink.test.list
    - TXT file containing a list of all the olink proteins.
- `scripts/`
- Directory for all R and shell scripts.
  - `data_exploration/`
    - The initial exploration of the Olink proteomic results containing IP-Protein Associations.
  - `lmer_testing/`
    - Scripts for performing the linear mixed-effects model (LMM).
    - `copy_helper_functions/`
      - Normalisation and lmer return regression coefficients scripts to call.
  - `protein_protein_solar/`
    - All R scripts needed to run SOLAR analysis and QUARTO (qmd) containing the analysis. 
    - `solar_chunks/`
      - TXT files containing SOLAR commands for submitting pairs.
  - `solar_analysis/`
    - From running SOLAR analysis to interpreting the results.
    - All scripts needed to run SOLAR including the initial creation of the Pedigree and Phenotype CSVs, and subsequent analysis in QUARTO (qmd).
      - `solar_chunks/`
          - TXT files containing SOLAR commands for submitting pairs.
- `results/`
- Directory designated for SOLAR analysis outputs.
  - `protein_protein_solar/`
    - job_chunk_1-200 directories with SOLAR outputs. Missing pairs from batching were performed here separately.
  - `solar_analysis/ `
    - job_chunk_1-200 directories with SOLAR outputs. Missing pairs from batching were performed here separately.
- `logs/`
- Directory designated for SOLAR analysis logs.
  - `protein_protein_solar/`
  - `solar_analysis/ `
- `software/`
- Directory for the SOLAR package installation (solar-eclipse-9.0.1-static-Linux).
- `visualisations/`
- All graphical and visual outputs from analysis.


