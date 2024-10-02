## Biochange_Paper1
This repository contains data and code relating to the study:

### **Recovery or reorganisation? Long-term increases in riverine taxonomic and functional diversity are confounded by compositional dynamics** 
### Published in Hydrobiologia (2024). DOI: https://doi.org/10.1007/s10750-024-05665-5

#### Authors:
Nathan J. Baker, Francesca Pilotto, Ellen A. R. Welti, Diana Osadčaja, Vaidas Palinauskas

## Metadata file:

#### ***AquaticMacroInvert_metadata.xls*** :
* Contains information on data providers, site characteristics, and the number of sites from each country
* Two tabs within this metadata spreedsheet contain metadata for the two csv files (in the outputs folder):
	* *All_indices_benthicMacroInverts_AllYears_alienzeros.csv* (metadata tab: metadata_siteyear)
	* *All_siteLevel.csv* (metadata tab: metadata_sites)

## Data folder:
The folder from which data is pulled for analyses

### Raw data files for calculating taxonomic and functional indices:

* ***LT_taxalist_2010-2020_long.csv*** : 
Contains raw taxonomic composition (abundance) data in long format for all included sites. Used for calculations of taxonomic diversity, XX, and XX.

* ***LT_taxalist_2010-2020_short_transposed.csv*** : 
Contains raw taxonomic composition (abundance) data in wide format for all included sites. Data are transported into a simplified site-by-species matrix. Used for calculations of functional diversity.

* ***LT_traits_2010-2020_short.csv*** : 
Contains raw trait data (extracted from freshwaterecology.info) in wide format for all included taxa. Used for calculations of functional diversity.
* ***LT_CWMs_long.csv*** : Contains site-year level calculated community-weighted means (CWMs) in long format. Used for calculations of functional turnover.

### Raw data files of environmental data and site metadata:

* ***LT_siteYr_EnvVariables.csv*** : 
Contains raw environmental data in wide format for 28 (out of 41) included sites. Used for all driver analyses, but included in main data file "LT_siteYr_AllData_wNAs_modified.csv".

* ***LT_siteYr_Factors.csv*** : 
Contains annual site-year level metadata, including sampling dates, coordinates, ecological quality ratios/classes, etc. Required for the calculation of functional turnover. 

* ***LT_site_metadata.csv*** : 
Contains condensed site-level metadata, including River sampled, coordinates, River type,	Extent of modification, etc. Required for checking metaanalytic models.

### Main data file for statistical analyses (see dedicated metadata file):

* ***LT_siteYr_AllData_wNAs_modified.csv*** : 
This file is used in all trend and statistical analyses. It contains data for each site and year on:
	* taxonomic diversity metrics
	* functional diversity metrics
	* taxonomic group specific indices
	* calculated environmental means from 12 months prior to biological sampling
	* site-level factors (coordinates, river type, extent of modification)

### TaskIDs:

* ***LT_ResponseTrends_TaskIDs.csv*** :
Task IDs for running models of taxonomic and functional indices and their model estimates

* ***LT_DriverTrends_TaskIDs.csv*** :
Task IDs for running models of environmental drivers and thier model estimates

## Outputs folder:
The folder to which data is saved following analyses.

* ***LT_siteYr_AllData.csv*** : 
Contains data for each site and year on: 
	* taxonomic diversity metrics
	* functional diversity metrics
	* taxonomic group specific indices
	* calculated environmental means from 12 months prior to biological sampling
	* site-level factors (coordinates, river type, extent of modification)

* ***glsTrends_site_level.rds*** :
Contains data for each site on slopes of functional and taxonomic diversity indices from the gls models used in "LT_macroinverts_site_trends.R". File created using "LT_combining_repsonse_trends.R". This file is necessary for illustrating the overall trends of biodiversity indices through time. Figures are generated using "LT_dentity_plots_trend_averages_and_distributions.R" and "LT_dentity_plots_trend_averages_and_distributions_wAllGroups.R". File can be viewed in R.

* ***glsTrends_site_level_drivers.rds*** :
Contains data for each site on slopes of environmental drivers from the gls models used in "LT_drivers_site_trends.R". File created using "LT_combining_driver_trends.R". This file is necessary for illustrating the overall trends of environmental drivers through time. Figures are generated using "LT_dentity_plots_trend_averages_and_distribution_drivers.R" and "LT_dentity_plots_trend_averages_and_distribution_nutrient_drivers.R". File can be viewed in R.

* ***glsTrends_drivers_site_level.rds*** :
Contains data on the responses of each taxonomic and functional index to the suite of environmental drivers from the linear mixed models used in "LT_driver_analysis_all_sites_combined.R". File created using "LT_combining_driver_estimates.R". This file is necessary for illustrating the overall responses of biodiversity indices to the the suite of environmental predictors used. Figures are generated using "LT_overall_driver_analyses_figures.R". File can be viewed in R.

* ***null_model_outputs_wNAs.csv*** :
Contains the null model outputs of the functional diversity metrics. This file is manually uploaded into the "Calculation_of_functional_indices.R" instead of having to rerun the null models every time I would like to rerun the code. PLEASE NOTE: Null model analyses are computationally exhausting. When running the null model analyses, I recommend first changing the number of iterations to 3, checking that the code runs smoothly, and thereafter increasing the number of iterations to 999. One should be aware that null model analyses can take a surpringly long time to run, for example on my machine (Dell XPS 15, 12th Gen Intel(R) Core(TM) i9-12900HK 2.50 GHz, 32.0 GB of RAM, with two graphics cards) it took ~8 hours per funcitonal index (e.g., FRic). Thus, patience is key.

* **Site_trends** folder:
Site-level GLS model outputs for biodiversity indices used in meta analytic models

* **Metaanalysis_trends** folder:
Meta analytic Bayesian mixed model outputs for biodiversity indices

* **Driver_trends** folder:
Site-level GLS model outputs for environmental drivers used in meta analytic models

* **Driver_metaanalysis_trends** folder:
Meta analytic Bayesian mixed model outputs for environmental drivers

* **Drivers_all_sites** folder:
Contains linear mixed model outputs for the overall responses of biodiversity indices to the suite of environmental drivers used.

* **outputs_movingWindow** folder:
Moving window analysis outputs

* **outputs_driver** folder:
Driver analysis outputs 

* **outputs_sensitivity** folder: 
Sensitivity analyses outputs including: 
	* high threshold moving window analysis (HTMW2)
	* one-country removal effects (Jackknife)
	* meta-analysis model comparisions (metaanalysisModelComparison)
	* site level high threshold moving window analysis (siteLevel_HTMW)
	* driver and moving window analyses only for sites with taxonomic resolution to species level (SppLevel_outputs)
	* taxonomic resolution and seasonality sensitivity analyses outputs (TaxonomicSeason)

* **TaskIDs** folder:
Task IDs for running models and model estimates

* **ClimateModel_stan** folder:
Contains outputs of models calculating temperature and precipitation trends


## R scripts folder

* **Initial_Biodiversity_FuncTrait_and_climate_calcs** folder:
Contains scripts to calculate taxonomic and functional trait metrics for each site year and calculate temperature and precipitation trends

* **stan_models** folder:
Contains scripts of all stan models used in analyses

* **trend_metaAnalysis** folder:
Contains the following scripts:
	* *HPC_macroinverts_stanmodels* : Models to calculate trends for each biodiversity metric and each of the 1816 sites
	* *HPC_Meta_analysis* : Meta-analysis models to synthesize the site-level data using the output of the previous step and get overall trend for each metric (Trend ~ 1 + (1|study_ID) + (1|Country)
	* *HPC_modelchecking* : Examining meta-analysis model fit parameters and calculating probability increasing/decreasing for each biodiversity metric

* **MovingWindow** folder: 
Contains the following scripts:
	* *HPC_MovingWindow_sites* : Models to calculate trends for each site within each window and biodiversity metric in moving window analysis
	* *HPC_Meta_analysis_movingwindow* : Models to calculate overall estimates of trends within each window and biodiverisity metric
	* *MetaMeta_movingwindowEsts_Yr* : Models to calculate overall linear year effects on biodiversity trajectories in moving window analyses

* **Driver** folder:
Contains the following scripts:
	* *HPC_Meta_analysis_drivers* : Models to caluculate driver effects on biodiversity metrics using site-level data
	* *Modelchecking_Drivers_horseshoePriors* : Examining driver model fit parameters

* **Sensitivity** folder:
Contains scripts for sensitivity checking models including: 
	* effects of trends within countries (Country_effects)
	* high threshold moving window analysis (HTMW)
	* using a one stage model rather than the two stage meta-analysis model for trend estimates (OneStage_models)
	* effects of taxonomic resolution for the moving window analysis and claculating the proportion of pos/neg trends/ window (checking_MovingWindow_parameters)
	* effects of taxonomic resolution and seasonality on trend estimates (HPC_Sensitivity_analysis)
	* examing driver model outputs for sites with taxa IDed to species level only (Modelchecking_Drivers_sppLevel_horseshoe)

* *HPC_combine_site_trends* script:
Concatenates outputs for many analyses including:  
	* biodiverity trends
	* meta-analysis
	* moving window analysis
	* driver analysis
	* sensitivity analyses

## Addtional functions folder

## Sensitivity folder

## plots folder

* ***Online Figures.docx*** : 
Includes all additional online figures not included in main text or Extended Data

* **Supplementary Table 2** : 
Contains information on time series locations, durations, and site characteristics

* **Fig1** folder: 
Contains scripts, data, and icons used to make Figure 1

* **Fig2_DensityPlots** folder : 
Contains scripts and plots related to meta-analysis trend results

* **Fig3_movingWindow** folder: 
Contains scripts and plots related to the moving window analyses

* **Fig4_drivers folder**: 
Contains scripts and plots related to driver analyses

* **descriptive_plots** folder: 
Contains scripts and plots descriping data distributions and correlations

* **Sensitivity** folder: 
Contains scripts and plots for sensitivity analyses regarding: 
	* one-country removal effects (Jackknife)
	* meta-analysis model comparisions (one stage, two stage weighted and unweighed) 
	* effects of sampling years & start year
	* effects of seasonality
	* effects of taxonomic resolution
