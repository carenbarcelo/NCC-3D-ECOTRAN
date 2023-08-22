# NCC-3D-ECOTRAN

This GitHub repository includes the code and the re-gridded ROMS driver files used for analysis in 

*Barceló et al. in review, “Non-linear and alternating spatial effects of climate change on the Northern California Current Ecosystem: Insights from an end-to-end ecosystem model”*

### Table of contents:
1.	A directory containing the full Matlab and C++ code set used in this analysis: **ECOTRAN_CodeSet** (C++ files include a version compiled for Mac OS and a raw, uncompiled version to allow compiling on other systems)
2.	A directory containing the excel VisualBasic file (.xlsm) with the mass-balanced Ecopath NCC food web and a .csv file with the full food web parameter set: **NCC_FoodWeb**.
3.	A directory containing three ROMS-ESM driver files for: ROMS-GFDL, ROMS-HAD, and ROMS-IPSL, on the Northern California Current ECOTRAN spatial grid: **ROMS_driverFiles_NCCgrid**. 
4.	A directory containing one example year of raw ROMS-ESM output and original ROMS grid information: **Raw_ROMS_driverFiles**. See Pozo Buil et al. (2021) for a full description of the original ROMS product.

### Defining food web parameters

ECOPATH-style mass-balance models are set up within Excel. Excel VisualBasic code provides the essential set of ECOPATH algorithms needed define a food web model directly without needing to use the complete ECOPATH software package. This was made available to us from Kerim Aydin (NOAA AFSC). This foodweb has been previously published in Ruzicka et al. (2016).

Please see the included Northern California Current food web parameter files, NCC2_09032022.xlsm and NCC2_09032022.csv.

### Main script m-file: 
*ECOTRANdynamic_NCC_08152023*	

* The user will need to change file directory calls to the directory structure of their local computer. 
* This calls all required sub-functions to run a complete simulation.
* The main output result of an ECOTRAN simulation is the variable **re_Y** describing the production rate of each functional group within each of 15 geographic sub-regions at each time-step (structured as a 3D matrix: time X group X sub-region).

### Main script sub-functions:
| File | Description |
| -----| ------------|
| *f_readEwEcsv_10pp_07072021* |	read the VisualBasic (ECOPATH) mass-balanced food web parameters from a .csv file as variable dat. NOTE: use for models with up to 10 (or less) primary producers| 
| *f_AggregateBiologicalModel_02052021*|	prepare ECOPATH parameter set for work in the ECOTRAN code. NOTE: you can aggregate functional groups here, if wanted| 
| *f_calcEE_12292020*	| calculate Ecotrophic Efficiency of each functional group to evaluate mass balance| 
| *f_VarianceDivision_12132018*	| calculate the variance of one variance term divided by another variance term
| *f_VarianceMultiplication_12132018*	| calculate the product of two variance terms| 
| *ECOTRANheart_09032021*	| generate an ECOTRAN model| 
| *f_ECOfunction_09032021*	| return a single ECOTRAN model from 1 "type" ECOPATH parameter set or from 1 randomly generated Monte Carlo parameter set| 
| *f_RedistributeCannibalism_11202019*	| remove cannibalism terms on diagonal of the diet matrix Dpc. NOTE: Cannibalism is directed to additional metabolism and feces production (equivalent to the cannibalism fraction of diet). This is the mechanism for reduced group Transfer Efficiency that J Steele pointed out is required if cannibalism is to be removed from diet matrix.| 
| *f_calcEE_12292020*	| calculate Ecotrophic Efficiency of each functional group to evaluate mass balance| 
| *f_calcPredationBudget_12102019*	| for each for each producer p, calculate the fraction of total predation going to each consumer c| 
| *f_E2Epedigree_08042020*	| calculate the uncertainty for every cell within the EnergyBudget (Acp) from pre-defined uncertainty values for all parameters| 
| *f_VarianceDivision_12132018*| 	calculate the variance of one variance term divided by another variance term| 
| *f_VarianceMultiplication_12132018*	| calculate the product of two variance terms| 
| *f_E2E_MonteCarlo_08042020*	| generate a stack of random EnergyBudget matrices (Acp) by drawing from a normal or a uniform distribution about each element of the EnergyBudget. The first model in the stack is the “type” model generated from the VisualBasic mass-balanced food web parameter set.| 
| *f_FunctionalResponse_MonteCarlo_05132021* | prepare array of vulnerability terms and allows for generation of random functional response terms within a predefined uncertainty level| 
| *f_InitialProductionRates_02012022* | calculate initial consumption rate q conditions| 
| *f_WebProductivity_03272019*| calculate consumption rates q of all groups under a given driver (e.g., NO3 or primary production); also accounts for defined rates of group production export when running static scenarios. NOTE: despite the function name, consumption rates are calculated, not production.| 

### Physical model functions:
| File | Description |
| -----| ------------|
| *f_OrdinalDate*	| Calculate the ordinal date from date format '01-Jan-1998'. The ordinal date is day of year with January 1 of any year = 1
| *f_ECOTRANphysics_NCC_ROMS_08152023*	| Use ROMS output to drive the ECOTRAN model through time. Sub-functions define the ECOTRAN spatial grid and aggregate ROMS volume flux output and biogeochemical model output to the ECOTRAN grid.
| *f_ROMS_GridPrep_NCC_08152023*	| Map ROMS grid to ECOTRAN grid
| *f_ROMS_FluxPrep_NCC_08152023*| 	express ROMS fluxes in DESTINY<--SOURCE format on ECOTRAN grid| 
| *f_CompactFluxTimeSeries_11182019*| 	compact physical flux time-series, arranged as 3D matrix (time X source box X destiny box). Also, eliminate source & destiny information for subregions that never communicate| 
| *f_UnCompactFluxTimeSeries_12112019*| 	UnCompact a previously compacted physical flux time-series to provide IMPORT & EXPORT fluxes for each box and the domain as a whole| 
| *f_calcNetFlux_12112019*| 	calculate net flux into and net flux out of each model box and across outer domain boundaries| 
| *f_EvaluateFluxBalance_11262021*| 	examine for flux time-series for imbalances IN & OUT of individual boxes and IN & OUT of the overall domain| 
| *f_LightIntensity_12112020*| 	Light intensity parameters: instantaneous (W/m2), daily mean averaged across 24 h (W m^-2 h^-1), & daily | integrated (W m^-2 d^-1) solar raditation at ocean surface; (vertical vector: num_t X 1). (Not actually used when driven by ROMS output)| 
| *f_ECOTRANmigration_NCC_02222022* | area of overlap of neighboring model sub-domains (for ROMS). (NOTE: this only provides the area of contact between spatial domain boxes. This information is not used by this version of the code package, but null variables are still needed for code to run)| 
| *f_DVMsinusoid_08122021*	| calculate diel vertical migration (DVM) flux rates between all model domain boxes for each functional group at each time point. Uses a daily sinusoidal migration pattern. (NOTE: Not defined for the NCC model, but null variables are still required to run the code)| 

### Physiological temperature response  functions:
| File | Description |
| -----| ------------|
| *f_physiology_Q10_12082022* | 	Calculate temperature-dependent Q10 metabolic rate scaling factors for living groups| 
| *f_TLparameterization_08182023*	| Prepare the eight required Thornton-Lessem parameters for each living group| 
| *f_ThorntonLessem_12012022*	| Calculate temperature-dependent Thornton-Lessem temperature scaling factors for ingestion rate| 

### Model solver:
The model is run by solving the system of Ordinary Differential Equations (ODE) for each functional group at each time step. 
The solver is coded in C++.

| File | Description |
| -----| ------------|
| *f_PrepMexODE_09192022*| prepare ECOTRAN variables for using C++ ODE solver mex function (mex = Matlab EXecutable function). Pack parameters and drivers along proper dimensions| 
| *f_unspoolMATRIX_04282020*| linearize (vectorize) multidimensional matrices up to 4D for use in C++.
| *mex_ECOTRANode_09182022.mexmaci64* |	use MATLAB-executable C++ function to solve the ODE for solution to functional group consumption rates at each time-point; default is for reflective boundary conditions. NOTE: this is compiled for Mac. NOTE: uncompiled code file is mex_ECOTRANode_09182022.cpp | 
