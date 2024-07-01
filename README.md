# multispecies_connectivity_risk
Analyses associated with "Multispecies migratory connectivity reveals hemispheric risk to birds"

# Code for Multispecies migratory connectivity reveals hemispheric risk to birds

## Authors

Sarah P. Saunders, William DeLuca, Brooke L. Bateman, Jill L. Deppe, Joanna Grand, Erika J. Knight, Timothy D. Meehan, Nicole L. Michel,
Nathaniel E. Seavy, Melanie A. Smith, Lotem Taylor, Chad J. Witko, Chad B. Wilsey.

## Information

Repository Type: Audubon's GitHub repository

Year of Origin:  2022

Year of Version: 2024

Version: 1.0.0

Digital Object Identifier (DOI): https://doi.org/xxxx

## Suggested Citation
Saunders, S., DeLuca, W., Bateman, B., Depe, J., Grand, J., Knight, E., Meehan, T., Michel, N., Seavy, N., Smith, M.,
Taylor, L., Witko, C., and Wilsey, C. Code for multispecies migratory connectivity reveals hemispheric risk to birds. 
Version 1.0.0; National Audubon Society. https://doi.org/xxxx

## Abstract 

This repository contains 2 R scripts in the workflow: 1_MC_ex.R; 2_MSMCrisk_ex.R. The overall objective of the project was to analyze movement data for 112 migratory bird species to estimate species-specific migratory connectivity between breeding and nonbreeding regions (hereafter MCR pairs) using an integrated model adapted from Korner-Nievergelt et al. 2017.
Then we combined species-specific estimates and associated uncertainties to quantify multispecies migratory connectivity per MCR pair by calculating weighted means and variances to derive average connectivity proportions, along with bootstrapped estimates of uncertainty, which we used to represent exposure.
Next, we combined exposure estimates with vulnerability and hazard metrics (conservation assessment scores and temperature and land-cover changes by mid-century, respectively), which represented the relative risk of future migratory bird population declines. The first script (1_MC_ex) estimates species-specific migratory connectivity. The example data provided is for Blackpoll warbler (bkpwar_data). The second
script (2_MSMCrisk_ex) estimates multispecies migratory connectivity and all relevant components of risk for an example set of 5 species for which connectivity results serve as input data (Connectivity_data_ExSp): Blackpoll Warbler, Grasshopper Sparrow, Ovenbird, Prothonotary Warbler, and Tree Swallow, along with hazard and vulnerability data (VulHaz_data). 
Raw and processed data inputs are used with data holder permissions.

### File 1: "1_MC_ex.R" 

Data inputs located in the following folder (movement data and species range shapefile): 
  bkpwar_data 

Outputs: 
  1. Model summary from JAGS (_ModelSummary.csv)
  2. Data frame of summarized results (_Results.csv) formatted for input to estimating multispecies connectivity
  in the second script of the workflow (2_MSMCrisk_ex.R script).
  3. Figure (_MigConnect.png): violin plot of estimated mean connectivity proportions and uncertainty per connection 
  with breeding MCRs as facets and nonbreeding MCR connections shown within each facet
  
### File 2: "2_MSMCrisk_ex.R" 

Data inputs located in the following folders (species connectivity results, vulnerability data, hazard data): 
  Connectivity_data_ExSp,
  VulHaz_data

Outputs: 
  1. Data frame of exposure (i.e. multispecies migratory connectivity) estimates and 
  uncertainty estimates (if applicable; only for MCR pairs with >1 species) per MCR pair (MSMC_est.csv).
  2. Example figure (MSMC_WeightMeansCI.png) of exposure and uncertainty estimates per MCR pair for those
  pairs with > 1 species (i.e. weighted means and 80% and 95% bootstrapped CIs)
  3. Data frame of risk calculations (Risk_components.csv): exposure, hazard, vulnerability, and 
  risk, both unscaled and scaled (normalized) values, for each MCR pair
  4. Model summary (Risk_regression_modsum.csv) from example post-hoc linear regression estimating relative importance of the
  three componenets of risk (exposure, vulnerability, hazard) in explaining variation in risk

Please refer to descriptions at the top of each R script for more details, including directions to run each script.

### Full_outputs: contains 3 outputs from the full analysis of 112 species:
1. MSMC_est_fullanalysis.csv: multispecies connectivity estimates and associated uncertainties (80 and 95% bootstrapped CIs for weighted means,
   i.e. connections with > 1 species; 95% CI for unweighted means, i.e. connections with 1 species)
3. Risk_components_fullanalysis.csv: estimates of risk and its components for all 921 connections
4. Studysp_groups_fullanalysis.csv: species associated with each connection (group) for the full analysis
