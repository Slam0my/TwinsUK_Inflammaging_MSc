# Identification of Early Instigators of Inflammaging Within the TwinsUK Cohort

- [Background](#introduction) 
- [Repository Structure](#repostructure)
- [References](#references)
 
## <a id="introduction"></a> Background

This repository is part of an extended research project for the MSc Applied Bioinformatics course at King's College London. This research project aims to identify the early instigators of inflammaging within the TwinsUK cohort [[1]](#ref1). Inflammaging is the chronic, low grade inflammation that occurs during ageing [[2]](#ref2). While current studies focus on symptomatic cohorts, investigations with a general healthy ageing cohort can enable insights into development of this intricate network where preventative interventions can be directed upon. This study aims to expand from the previous studies by Roederer et al. (2015) [[3]](#ref3) and Mangino et al. (2017) [[4]](#ref4) where heritability and genetic architectures of the immune system were observed .


With immunophenotypes (IPs) and immune related proteins obtained from patients, such pairs are investigated utilising the Twins design. Through the genetic variance components analysis tool [SOLAR-Eclipse](https://solar-eclipse-genetics.org/), heritability, genetic and environmental correlations were estimated for IP-Protein pairs. In turn, this filtering enables research into the complex interplay for identification of early instigators.

Due to the sensitive nature of the dataset used throughout this pipeline, the raw datasets are excluded from this repository for data privacy and adherence. Access to the data may be requested via application to TwinsUK [[1]](#ref1).

## <a id="repostructure"></a> Repository Structure 

All of the bioinformatics, computational analysis pipeline was conducted on the Computational Research, Engineering and Technology Environment (CREATE) high performance computing (HPC) [[5]](#ref5). Data, results, logs and software directories have been omitted from this GitHub repository. 


- `data/` ->Directory containing all of the data used in the workflow. This includes a copy of the raw data originals (untouched), CSVs from linear mixed-effects model (LMM) and SOLAR-Eclipse analysis outputs.
    
- `scripts/` ->  R scripts for performing LMM, SOLAR-Eclipse from the FDR P-value > 0.05 IP-Protein pairs, which also includes the initial creation of the Pedigree and Phenotype CSVs, and for Protein-Protein pairs.

- `notebook/` -> All qmd from initial data exploration, solar filtering analysis and within trait analysis.

- `Rshiny_app/` -> Rshiny dashboards curated for interactive visualisations of network graphs. These include IP-Protein and Protein-Protein SOLAR-Eclipse analysis pairs and filtered from those genetically and environmentally significant. To run the RShiny app, access to the SOLAR analysis results RData is required, which may be requested via TwinsUK [[1]](#ref1).

- `results/` -> Contains SOLAR analysis outputs for IP-Protein pairs and Protein-Protein pairs.

- `logs/`-> For SOLAR analysis output logs.

- `software/`-> Directory for the SOLAR package installation (solar-eclipse-9.0.1-static-Linux).

- `visualisations/` -> All graphical and visual outputs from analysis. Directories are sorted based on post-SOLAR filtering and within trait analyses.



## <a id="references"></a> References

<a id="ref1"></a> [1] Verdi S, Abbasian G, Bowyer RCE, Lachance G, Yarand D, Christofidou P, Mangino M, Menni C, Bell JT, Falchi M, Small KS, Williams FMK, Hammond CJ, Hart DJ, Spector TD, Steves CJ. TwinsUK: The UK Adult Twin Registry Update. Twin Res Hum Genet. 2019 Dec;22(6):523-529. doi: 10.1017/thg.2019.65. Epub 2019 Sep 17. PMID: 31526404. 

<a id="ref2"></a> [2] López-Otín C, Blasco MA, Partridge L, Serrano M, Kroemer G. Hallmarks of aging: An expanding universe. Cell. 2023 Jan 19;186(2):243-278. doi: 10.1016/j.cell.2022.11.001. Epub 2023 Jan 3. PMID: 36599349. 


<a id="ref2"></a> [3]  Roederer M, Quaye L, Mangino M, Beddall MH, Mahnke Y, Chattopadhyay P, Tosi I, Napolitano L, Terranova Barberio M, Menni C, Villanova F, Di Meglio P, Spector TD, Nestle FO. The genetic architecture of the human immune system: a bioresource for autoimmunity and disease pathogenesis. Cell. 2015 Apr 9;161(2):387-403. doi: 10.1016/j.cell.2015.02.046. Epub 2015 Mar 12. PMID: 25772697; PMCID: PMC4393780.

<a id="ref2"></a> [4] Mangino M, Roederer M, Beddall MH, Nestle FO, Spector TD. Innate and adaptive immune traits are differentially affected by genetic and environmental factors. Nat Commun. 2017 Jan 5;8:13850. doi: 10.1038/ncomms13850. PMID: 28054551; PMCID: PMC5227062.

<a id="ref2"></a> [5] King’s College London e-Research team. King’s Computational Research, Engineering and Technology Environment (CREATE). 2022. [Read more](https://docs.er.kcl.ac.uk/)
