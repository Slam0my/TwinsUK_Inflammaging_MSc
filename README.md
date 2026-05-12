# Identification of Early Instigators of Inflammaging Within the TwinsUK Cohort

- [Background](#introduction) 
- [Repository Structure](#repostructure)
- [References](#references)
 
## <a id="introduction"></a> Background

This repository is part of an extended research project for the MSc Applied Bioinformatics course at King's College London. This research project aims to identify the early instigators of inflammaging within the TwinsUK cohort [^1]. Inflammaging is the chronic, low grade inflammation that occurs during ageing. While current studies focus on symptomatic cohorts, investigations with a general healthy ageing cohort can enable insights into development of this intricate network where preventative interventions can be directed upon. 

With immunophenotypes (IPs) and immune related proteins obtained from patients, such pairs are investigated utilising the Twins design. Through the genetic variance components analysis tool [SOLAR-Eclipse](https://solar-eclipse-genetics.org/), heritability, genetic and environmental correlations were estimated for IP-Protein pairs. In turn, this filtering enables research into the complex interplay for identification of early instigators/


## <a id="repostructure"></a> Repository Structure 

All of the bioinformatics, computational analysis pipeline was conducted on the Computational Research, Engineering and Technology Environment (CREATE) high performance computing (HPC) [^3]. Data, results, logs and software directories have been omitted from this GitHub repository. 

- `data/`
- Directory containing all of the data used in the workflow. This includes a copy of the raw data originals (untouched), CSVs from linear mixed-effects model (LMM) and SOLAR-Eclipse analysis outputs.
    
- `scripts/`
- R scripts for performing LMM, SOLAR-Eclipse from the FDR P-value > 0.05 IP-Protein pairs, which also includes the initial creation of the Pedigree and Phenotype CSVs, and for Protein-Protein pairs.

- `notebook/`
- All qmd from initial data exploration, solar filtering analysis and within trait analysis.

- `Rshiny_app/`
- Rshiny dashboards curated for interactive visualisations of network graphs. These include IP-Protein and Protein-Protein SOLAR-Eclipse analysis pairs and filtered from those genetically and environmentally significant.

- `results/`
- Contains SOLAR analysis outputs for IP-Protein pairs and Protein-Protein pairs.

- `logs/`
- For SOLAR analysis logs.

- `software/`
- Directory for the SOLAR package installation (solar-eclipse-9.0.1-static-Linux).

- `visualisations/`
- All graphical and visual outputs from analysis. Directories are sorted based on post-SOLAR filtering and within trait analyses.




## <a id="references"></a> References

[^1]: Serena Verdi et al. TwinsUK: The UK Adult Twin Registry Update. 2019 [Read more](https://www.cambridge.org/core/ journals/twin-research-and-human-genetics/article/twinsuk-the-ukadult-twin-registry-update/87C69A6FFED4481B9E4AA8BF14F79A8C)

[^3]	King’s College London e-Research team. King’s Computational Research, Engineering and Technology Environment (CREATE). 2022. [Read more](https://docs.er.kcl.ac.uk/)
