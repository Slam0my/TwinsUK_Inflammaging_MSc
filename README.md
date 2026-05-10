# Identification of Early Instigators of Inflammaging Within the TwinsUK Cohort

- [Background](#introduction) 
- [Repository Structure](#repostructure)

- 
## <a id="introduction"></a> Background

This repository is part of an extended research project for the MSc Applied Bioinformatics course at King's College London. This research project aims to identify the early instigators of inflammaging within the TwinsUK cohort. Inflammaging is the chronic, low grade inflammation that occurs during ageing. While current studies focus on symptomatic cohorts, investigations with a general healthy ageing cohort can enable insights into development of this intricate network where preventative interventions can be directed upon. With immunophenotypes (IPs) and immune related proteins obtained from patients, such pairs are investigated for heritability, genetic and environmental correlations utilising the Twins design. In turn, this enables filtering and research into the complex interplay, enabling such early instigators to be identified. 


## <a id="repostructure"></a> Repository Structure 

All of the bioinformatics, computational analysis pipeline was conducted on the Computational Research, Engineering and Technology Environment (CREATE) high performance computing (HPC). Data, results, logs and software directories have been omitted from this GitHub repository. 

- `data/`
- Directory containing all of the data used in the workflow. This includes a copy of the raw data, CSVs from LMM and SOLAR analysis outputs
    
- `scripts/`
- R scripts for performing linear mixed-effects model (LMM), SOLAR-Eclipse from the FDR P-value > 0.05 IP-Protein pairs, which also includes the initial creation of the Pedigree and Phenotype CSVs, and for Protein-Protein pairs.

- `notebook/`
- All qmd from initial data exploration, solar filtering analysis and within trait analysis.
  - `IP-Protein_Networks/`
    - RShiny App for an interactive network graph of IP-Protein pairs from those genetically and environmentally significant.
 
- `results/`
- Contains SOLAR analysis outputs for IP-Protein pairs and Protein-Protein pairs.

- `logs/`
- For SOLAR analysis logs.

- `software/`
- Directory for the SOLAR package installation (solar-eclipse-9.0.1-static-Linux).

- `visualisations/`
- All graphical and visual outputs from analysis. Directories are sorted based on filtering.


