# Simulation-Based Sensitivity Analysis in Optimal Treatment Regimes and Causal Decomposition with Individualized Interventions
R codes for "Simulation-Based Sensitivity Analysis in Optimal Treatment Regimes and Causal Decomposition with Individualized Interventions"

Soojin Park<sup>1</sup>, Suyeon Kang<sup>2</sup>, and Chioun Lee<sup>3</sup>

<sup>1</sup> School of Education, University of California, Riverside  
<sup>2</sup> Department of Statistics, University of Central Florida
<sup>3</sup> Department of Sociology, University of California, Riverside


## Overview

Causal decomposition analysis aims to assess the effect of modifying risk factors on reducing social disparities in outcomes. Recently, this analysis has incorporated individual characteristics when modifying risk factors by utilizing optimal treatment regimes (OTRs). Since the newly defined individualized effects rely on the no omitted confounding assumption, developing sensitivity analyses to account for potential omitted confounding is essential. Moreover, OTRs and individualized effects are primarily based on binary risk factors, and no formal approach currently exists to benchmark the strength of omitted confounding using observed covariates for binary risk factors. To address this gap, we extend a simulation-based sensitivity analysis that simulates unmeasured confounders, addressing two sources of bias emerging from deriving OTRs and estimating individualized effects. Additionally, we propose a formal bounding strategy that benchmarks the strength of omitted confounding for binary risk factors. Using the High School Longitudinal Study 2009 (HSLS:09), we demonstrate this sensitivity analysis and benchmarking method..

For more details of our proposed methods, see [our paper](https://www.degruyter.com/document/doi/10.1515/jci-2022-0031/html). 
Here, we provide `R` codes to reproduce our simulation study and replicate our data analysis. 

## Case Study

* `data_imputed.RData` 
  
  For our case study, we used data from the Midlife Development in the U.S. (MIDUS) study. However, as the MIDUS data is restricted from circulation, and the original data can be downloaded from the MIDUS portal by clicking [here](https://www.midus.wisc.edu/data/index.php). 

* `CaseStudy_alg9_main.R` 
 
   This `R` file replicates Tables 2 and Table 3 of our study.

* `CaseStudy_alg9_Functions.R` 
 
   This `R` file contains source files used in `CaseStudy_al9_main.R`.

## Simulation Study

* `Simulation_het.R`  

   This `R` file contains the simulation codes for our propposed simulation-based sensitivity analysis for heterogenous effects. This code replicates Figures 2 and 3 of our paper.

* `Simulation_hom.R`  

   This `R` file contains the simulation codes for our propposed simulation-based sensitivity analysis for homogeneous effects. This code replicates Figures 2 and 3 of our paper.

* `functions_het.R` 
 
   This `R` file includes source functions required to run our heterogenous simulation codes.

* `functions_hom.R` 
 
   This `R` file includes source functions required to run our homogenous simulation codes. 

These supplementary materials are provided solely for the purpose of reproducibility and must be used in compliance with academic ethical guidelines. If you reference these materials in your own work, please ensure proper citation of the original sources.
