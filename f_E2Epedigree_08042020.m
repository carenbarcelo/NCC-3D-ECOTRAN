function ECOTRAN_PEDIGREE = f_E2Epedigree_08042020(ECOTRAN, PEDIGREE)
% by Jim Ruzicka
% calculate CV uncertainty terms to apply to ECOTRAN variables
% for use in Monte Carlo analyses with ECOTRAN variables (f_E2E_MonteCarlo_08032018)
%       taking advantage of fact that ECOTRAN models are by definition in balance
%       and avoiding the need to run Monte Carlo analyses at the ECOPATH (EwE) stage
%
% calls:
%       f_VarianceDivision_12132018                calculate the variance of one term divided by another term
%       f_VarianceMultiplication_12132018          calculate the variance of two products
%
% takes:
%       ECOTRAN
%       PEDIGREE
%
% returns:
%       ECOTRAN_PEDIGREE
%           biomass_CV           	CV values for biomass           (horizontal vector: 1 X num_grps)
%       	production_CV         	CV values for production rate	(horizontal vector: 1 X num_grps)
%       	pb_CV                 	CV values for pb                (horizontal vector: 1 X num_grps)
%       	qb_CV                 	CV values for qb                (horizontal vector: 1 X num_grps)
%       	pq_CV                 	CV values for pq                (horizontal vector: 1 X num_grps)
%       	ae_CV                	CV values for ae                (horizontal vector: 1 X num_grps)
%       	BioenergeticBudget_CV	CV values for each element of the BioenergeticBudget      (3D matrix: 3 X num_grps)
%                                                   1) feces
%                                                   2) metabolism
%                                                   3) production
%       	ProductionBudget_CV     CV values for each element of the ProductionBudget        (3D matrix: 5 X num_grps)
%                                                   1) eggs (reproduction)
%                                                   2) predation
%                                                   3) senescence
%                                                   4) ba (biomass accumulation)
%                                                   5) em (emigration); NOTE: negative for immigration
%       	ConsumptionBudget_CV	CV values for each element of the ConsumptionBudget       (3D matrix: 7 X num_grps)
%                                                   1) feces
%                                                   2) metabolism
%                                                   3) eggs (reproduction)
%                                                   4) predation
%                                                   5) senescence
%                                                   6) ba (biomass accumulation)
%                                                   7) em (emigration); NOTE: negative for immigration
%       	EnergyBudget_CV         CV values for each element of the EnergyBudget            (3D matrix: num_grps X num_grps)
%       	fname_E2Epedigree       name of this f_E2Epedigree function
%
% revision date: 8-4-2020

% FFF In future, could add uncertainty to detritus, metabolism fates, and retention scalers
% FFF in future, could add uncertainty of non-biological oxidation (reduction?) of NH4 to NO3


% *************************************************************************
% STEP 1: read in PEDIGREE values------------------------------------------
fname_E2Epedigree           = mfilename; % save name of this m-file to keep in saved model results
display(['Running: ' fname_E2Epedigree])


% step 1a: read in parameter PEDIGREE (CV) --------------------------------
biomass_CV                  = PEDIGREE.biomass_CV;	% (vertical vector: num_grps X 1)
pb_CV                       = PEDIGREE.pb_CV;    	% (vertical vector: num_grps X 1)
% qb_CV                       = PEDIGREE.qb_CV;     	% (vertical vector: num_grps X 1); NOTE: recalculated below as a dependent variable
pq_CV                       = PEDIGREE.pq_CV;     	% (vertical vector: num_grps X 1)
ae_CV                       = PEDIGREE.ae_CV;   	% (vertical vector: num_grps X 1)
ba_CV                       = PEDIGREE.ba_CV;    	% CV of ba absolute production rate; (vertical vector: num_grps X 1)
em_CV                       = PEDIGREE.em_CV;     	% CV of em absolute production rate; (vertical vector: num_grps X 1)


% step 1b: read in consumption & diet PEDIGREE (CV) -----------------------
CONSUMPTION_CV              = PEDIGREE.CONSUMPTION_CV;          % (2D matrix: num_grps X num_grps)
consumption_total_CV        = PEDIGREE.consumption_total_CV;	% (vertical vector: num_grps X 1)
consumption_domestic_CV     = PEDIGREE.consumption_domestic_CV;	% (vertical vector: num_grps X 1)
consumption_import_CV       = PEDIGREE.consumption_import_CV;	% (vertical vector: num_grps X 1) 
DIET_CV                     = PEDIGREE.DIET_CV;                 % (2D matrix: num_grps X num_grps)
diet_import_CV              = PEDIGREE.diet_import_CV;          % (horizontal vector: 1 X num_grps)


% step 1c: read in fleet PEDIGREE (CV) ------------------------------------
LANDINGS_CV                 = PEDIGREE.LANDINGS_CV;                 % (2D matrix: num_grps X num_fleets)
DISCARDS_CV                 = PEDIGREE.DISCARDS_CV;                 % (2D matrix: num_grps X num_fleets)
CATCH_CV                    = PEDIGREE.CATCH_CV;                    % (2D matrix: num_grps X num_fleets)
landings_TotalByFleet_CV	= PEDIGREE.landings_TotalByFleet_CV;	% (horizontal vector: 1 X num_fleets)
discards_TotalByFleet_CV    = PEDIGREE.discards_TotalByFleet_CV;	% (horizontal vector: 1 X num_fleets)
catch_TotalByFleet_CV       = PEDIGREE.catch_TotalByFleet_CV;       % (horizontal vector: 1 X num_fleets)
% landings_TotalByGroup_CV	  = PEDIGREE.landings_TotalByGroup_CV;	  % (vertical vector: num_grps X 1); NOTE: not used here
% discards_TotalByGroup_CV	  = PEDIGREE.discards_TotalByGroup_CV;	  % (vertical vector: num_grps X 1); NOTE: not used here
% catch_TotalByGroup_CV       = PEDIGREE.catch_TotalByGroup_CV;       % (vertical vector: num_grps X 1); NOTE: not used here


% step 1d: read in forced egg production & metabolism PEDIGREE (CV) -------
% EggProduction_forced_CV     = PEDIGREE.EggProduction_forced_CV;     % (vertical vector: num_grps X 1); NOTE: FFF for future
% Metabolism_forced_CV        = PEDIGREE.Metabolism_forced_CV;        % (vertical vector: num_grps X 1); NOTE: FFF for future


% step 1e: read in fates PEDIGREE (CV) ------------------------------------
%          NOTE: there is no fate_eggs_CV as egg production is 100% certain to go to a single assigned functional group
fate_feces_CV               = PEDIGREE.fate_feces_CV;       % (horizontal vector: 1 X num_AggGrps)
fate_metabolism_CV          = PEDIGREE.fate_metabolism_CV;	% (horizontal vector: 1 X num_AggGrps)
fate_senescence_CV          = PEDIGREE.fate_senescence_CV;	% (horizontal vector: 1 X num_AggGrps)


% step 1f: defined pedigrees for all model groups -------------------------
ee_eggs_CV                  = PEDIGREE.ee_eggs_CV;          % egg CV W.R.T. production budget for all groups; (CV); (scaler); FFF in futre allow for different values for each group
BacterialMTBLSM_CV          = PEDIGREE.BacterialMTBLSM_CV;	% pedigree for implicit bacterial metabolism of terminal detritus; (CV); (scaler)
Oxidation_NH4_CV         	= PEDIGREE.Oxidation_NH4_CV;	% fraction of NH4 produced oxidized directly back to NO3 abiologically; (vertical vector: num_NH4 X 1)
NutrientUptake_CV           = PEDIGREE.NutrientUptake_CV;   % pedigree for nutrient uptake by primary producers; (CV); (scaler)
% *************************************************************************





% *************************************************************************
% STEP 2: read in ECOTRAN terms--------------------------------------------

% step 2a: GroupType ------------------------------------------------------
GroupType                           = ECOTRAN.GroupType;	% includes nutrients & fleets; (vertical vector: num_grps X 1)


% step 2b: group type definitions -----------------------------------------
GroupTypeDef_ANYNitroNutr           = ECOTRAN.GroupTypeDef_ANYNitroNutr;
GroupTypeDef_NO3                    = ECOTRAN.GroupTypeDef_NO3;
GroupTypeDef_plgcNH4                = ECOTRAN.GroupTypeDef_plgcNH4;
GroupTypeDef_bnthNH4                = ECOTRAN.GroupTypeDef_bnthNH4;
GroupTypeDef_ANYPrimaryProd         = ECOTRAN.GroupTypeDef_ANYPrimaryProd;
GroupTypeDef_ANYConsumer            = ECOTRAN.GroupTypeDef_ANYConsumer;
GroupTypeDef_fleet                  = ECOTRAN.GroupTypeDef_fleet;
GroupTypeDef_eggs                   = ECOTRAN.GroupTypeDef_eggs;
GroupTypeDef_terminalPlgcDetr       = ECOTRAN.GroupTypeDef_terminalPlgcDetr;
GroupTypeDef_terminalBnthDetr       = ECOTRAN.GroupTypeDef_terminalBnthDetr;
GroupTypeDef_ANYDetritus            = ECOTRAN.GroupTypeDef_ANYDetritus;
GroupTypeDef_offal                  = ECOTRAN.GroupTypeDef_offal;
GroupTypeDef_micrograzers           = ECOTRAN.GroupTypeDef_micrograzers;
GroupTypeDef_bacteria               = ECOTRAN.GroupTypeDef_bacteria;


% step 2c: find group addresses -------------------------------------------
looky_nutrients                     = find(floor(GroupType) == GroupTypeDef_ANYNitroNutr); % row addresses of nutrients
looky_NO3                       	= find(GroupType == GroupTypeDef_NO3);
looky_plgcNH4                       = find(GroupType == GroupTypeDef_plgcNH4);
looky_bnthNH4                       = find(GroupType == GroupTypeDef_bnthNH4);
looky_NH4                           = find(GroupType == GroupTypeDef_plgcNH4 | GroupType == GroupTypeDef_bnthNH4);
looky_ANYPrimaryProd                = find(floor(GroupType) == GroupTypeDef_ANYPrimaryProd);
looky_ANYconsumer                   = find(floor(GroupType) == GroupTypeDef_ANYConsumer);
looky_bacteria                      = find(floor(GroupType) == GroupTypeDef_bacteria);
looky_eggs                          = find(GroupType == GroupTypeDef_eggs);
looky_ANYdetritus                   = find(floor(GroupType) == GroupTypeDef_ANYDetritus);
looky_terminalPLGCdetritus          = find(GroupType == GroupTypeDef_terminalPlgcDetr);
looky_terminalBNTHdetritus          = find(GroupType == GroupTypeDef_terminalBnthDetr);
looky_eggsANDdetritus               = sort([looky_ANYdetritus; looky_eggs]);
looky_fleets                        = find(floor(GroupType) == GroupTypeDef_fleet);
looky_livingANDfleets               = [looky_ANYPrimaryProd; looky_ANYconsumer; looky_bacteria; looky_fleets]; % includes primary producers & bacteria
looky_offal                         = find(GroupType == GroupTypeDef_offal);        % offal group(s)
looky_OtherDetritus                 = find((floor(GroupType) == GroupTypeDef_ANYDetritus) & (GroupType ~= GroupTypeDef_terminalPlgcDetr) & (GroupType ~= GroupTypeDef_terminalBnthDetr)); % non-terminal, non-egg detritus groups
looky_DetritusInDetritus            = find((ismember(looky_eggsANDdetritus, (looky_eggsANDdetritus(~ismember(looky_eggsANDdetritus, looky_eggs))))));    % non-egg addresses in EwE_DetritusFate matrix
looky_EggsInDetritus                = find((ismember(looky_eggsANDdetritus, (looky_eggsANDdetritus(ismember(looky_eggsANDdetritus, looky_eggs))))));     % egg addresses in EwE_DetritusFate matrix


% step 2d: number of group types ------------------------------------------
num_grps                            = ECOTRAN.num_grps;
num_nutrients                       = ECOTRAN.num_nutrients;
num_NO3                             = ECOTRAN.num_NO3;
num_NH4                             = ECOTRAN.num_NH4;
num_plgcNH4                         = ECOTRAN.num_plgcNH4;
num_bnthNH4                         = ECOTRAN.num_bnthNH4;
num_PrimaryProducers                = ECOTRAN.num_PrimaryProducers;
num_fleets                          = ECOTRAN.num_fleets;
num_nonFleetGrps                    = num_grps - num_fleets;
num_eggs                            = ECOTRAN.num_eggs;
num_ANYdetritus                     = ECOTRAN.num_ANYdetritus;
num_terminalBNTHdetritus            = ECOTRAN.num_terminalBNTHdetritus;
num_eggsANDdetritus                 = ECOTRAN.num_eggsANDdetritus;


% step 2e: read in parameter values ---------------------------------------
biomass                             = ECOTRAN.biomass;      % (t/km2); (vertical vector: num_grps X 1)
production                          = ECOTRAN.production;   % (t/km2/y); (vertical vector: num_grps X 1)
pb                                  = ECOTRAN.pb;           % (t/km2/y); (vertical vector: num_grps X 1)
qb                                  = ECOTRAN.qb;           % (1/y); (vertical vector: num_grps X 1)
pq                                  = ECOTRAN.pq;           % (1/y); (vertical vector: num_grps X 1)
ae                                  = ECOTRAN.ae;           % (dimensionless); (vertical vector: num_grps X 1)
ba                                  = ECOTRAN.ba;           % (t/km2/y); (vertical vector: num_grps X 1); NOTE: ba in absolute production rate terms
em                                  = ECOTRAN.em;           % (t/km2/y); (vertical vector: num_grps X 1); NOTE: em in absolute production rate terms


% step 2f: ecotrophic efficiency terms ------------------------------------
%          proportions of total production used within model domain
ee_eggs                             = ECOTRAN.ee_eggs;      % eggs (&/or reproduction) as proportion of production; (vertical vector: num_grps X 1)
ee_predation                        = ECOTRAN.ee_predation;	% predation as proportion of production; (vertical vector: num_grps X 1)
ee_ba                               = ECOTRAN.ee_ba;        % biomass accumulation as proportion of production; (vertical vector: num_grps X 1)
ee_em                               = ECOTRAN.ee_em;        % emigration as proportion of production; (vertical vector: num_grps X 1)
ee                                  = ECOTRAN.ee;           % proportion of production used within model domain; (vertical vector: num_grps X 1)


% step 2g: read in consumption & diet -------------------------------------
CONSUMPTION                         = ECOTRAN.CONSUMPTION;       	% (2D matrix: num_grps X num_grps)
consumption_total                   = ECOTRAN.consumption_total;  	% (vertical vector: num_grps X 1)
consumption_domestic                = ECOTRAN.consumption_domestic;	% (vertical vector: num_grps X 1)
consumption_import                  = ECOTRAN.consumption_import;	% (vertical vector: num_grps X 1)
DIET                                = ECOTRAN.DIET;                 % (2D matrix: num_grps X num_grps)
import_diet                         = ECOTRAN.import_diet;          % (horizontal vector: 1 X num_grps)


% step 2h: read in fleet values -----------------------------------
LANDINGS                            = ECOTRAN.LANDINGS;                 % (2D matrix: num_grps X num_fleets)
DISCARDS                            = ECOTRAN.DISCARDS;                 % (2D matrix: num_grps X num_fleets)
CATCH                               = ECOTRAN.CATCH;                    % (2D matrix: num_grps X num_fleets)
landings_TotalByFleet               = ECOTRAN.landings_TotalByFleet;	% (horizontal vector: 1 X num_fleets)
discards_TotalByFleet               = ECOTRAN.discards_TotalByFleet;	% (horizontal vector: 1 X num_fleets)
catch_TotalByFleet                  = ECOTRAN.catch_TotalByFleet;       % (horizontal vector: 1 X num_fleets)
DiscardFraction                     = ECOTRAN.DiscardFraction;          % fraction of catch discarded for each group; (2D matrix: num_grps X num_fleets)


% step 2i: read in fates --------------------------------------------------
fate_feces                          = ECOTRAN.fate_feces;     	% (2D matrix: num_ANYdetritus X num_grps)
fate_metabolism                     = ECOTRAN.fate_metabolism;	% (2D matrix: num_nutrients X num_grps)
fate_eggs                           = ECOTRAN.fate_eggs;      	% (2D matrix: num_eggs X num_grps)
fate_senescence                     = ECOTRAN.fate_senescence;	% (2D matrix: num_ANYdetritus X num_grps)


% step 2j: forced metabolism & egg rates ----------------------------------
% EggProduction_forced                = ECOTRAN.EggProduction_forced;     % (t/km2/y); (vertical vector: num_AggGrps X 1); NOTE: not used here; FFF for future
% Metabolism_forced                   = ECOTRAN.Metabolism_forced;     	  % (t/km2/y); (vertical vector: num_AggGrps X 1); NOTE: not used here; FFF for future


% % step 2k: retention & production loss scalers ----------------------------
% ProductionLossScaler                = ECOTRAN.ProductionLossScaler;     % (scaler: 0 to 1); (vertical vector: num_AggGrps X 1); NOTE: not used here
% RetentionScaler                     = ECOTRAN.RetentionScaler;          % (scaler: 0 to 1); (vertical vector: num_AggGrps X 1); NOTE: not used here


% % step 2l: functional response terms -------------------------------------
% FunctionalResponseParams            = ECOTRAN.FunctionalResponseParams;     % (vertical vector: num_AggGrps X 4); NOTE: not used here
% FunctionalResponse_matrix           = ECOTRAN.FunctionalResponse_matrix;	% (2D matrix: num_AggGrps X num_AggGrps); NOTE: not used here; FFF for future


% step 2m: nutrient & detritus cycling terms ------------------------------
PelagicBacterialReduction           = ECOTRAN.PelagicBacterialReduction;	% fraction of terminal pelagic detritus oxidized to pelagic NH4
BenthicBacterialReduction           = ECOTRAN.BenthicBacterialReduction;	% fraction of terminal benthic detritus oxidized to benthic NH4
Oxidation_NH4                       = ECOTRAN.Oxidation_NH4;                % fraction of NH4 produced oxidized directly back to NO3 abiologically; this should take precedence over phytoplankton uptake of NH4 in code below; (vertical vector: num_NH4 X 1)
PhytoUptake_NO3                     = ECOTRAN.PhytoUptake_NO3;              % (fraction); sums to 1; (vertical vector: num_PrimaryProducers X 1)
PhytoUptake_NH4                     = ECOTRAN.PhytoUptake_NH4;              % (fraction); sums to 1; (vertical vector: num_PrimaryProducers X 1)


% step 2n: four budget matrices -------------------------------------------
BioenergeticBudget                  = ECOTRAN.BioenergeticBudget;   % (2D matrix: 3 X num_grps)
%                                       1) feces
%                                       2) metabolism
%                                       3) production
ProductionBudget                    = ECOTRAN.ProductionBudget; % (2D matrix: 5 X num_grps)
%                                       1) eggs (reproduction)
%                                       2) predation
%                                       3) senescence
%                                       4) ba (biomass accumulation)
%                                       5) em (emigration); NOTE: negative for immigration
ConsumptionBudget                   = ECOTRAN.ConsumptionBudget; % (2D matrix: 7 X num_grps)
%                                       1) feces
%                                       2) metabolism
%                                       3) eggs (reproduction)
%                                       4) predation
%                                       5) senescence
%                                       6) ba (biomass accumulation)
%                                       7) em (emigration); NOTE: negative for immigration
EnergyBudget                        = ECOTRAN.EnergyBudget; % (2D matrix: num_grps X num_grps)
% *************************************************************************





% *************************************************************************
% STEP 3: calculate pedigree terms-----------------------------------------
%
%           RULE: VAR     = STD^2 = (mean * CV)^2
%
%           RULE: CV      = STD / mean
%
%           RULE: adding & subtracting uncertainties: 
%                       mean(A+B) = mean(A) + mean(B) OR mean(A-B) = mean(A) - mean(B) 
%                       var(A+B)  = var(A) + var(B) (variances are additive whether adding or subtracting terms)
%
%           RULE: multiplying uncertainties of INDEPENDENT variables (i.e., we can ignore covariance term): 
%                 (from: https://chem.libretexts.org/Core/Analytical_Chemistry/Quantifying_Nature/Significant_Digits/Propagation_of_Error)
%                       mean(A*B) = mean(A) * mean(B)
%                       var(A*B)  = [mean(A*B)^2] * [(var(A)/mean(A)^2) + (var(B)/mean(B)^2)]    NOTE: Goodman (1960) derives the same equation but also adds the term [... + var(A)*var(B)]
%                                                                                                      (Goodman, L.A. 1960. On the exact variance of products. Journal of the American Statistical Association. December 1960, 708?713. DOI: 10.2307/2281592)
%           RULE: dividing uncertainties of INDEPENDENT variables (i.e., we ignore covariance term):
%                 (from: Seltman, H. Approximations for Mean and Variance of a Ratio. Carnegie Mellon University)
%                       mean(A/B) = mean(A) / mean(B) 
%                       var(A/B)  = [mean(A)^2/mean(B)^2] * [(var(A)/mean(A)^2) + (var(B)/mean(B)^2)]

% step 3a: parameter variance values --------------------------------------
biomass_VAR                         = (biomass .* biomass_CV).^2;     	% variance; (vertical vector: num_grps X 1)
pb_VAR                              = (pb      .* pb_CV).^2;          	% variance; (vertical vector: num_grps X 1)
% qb_VAR                              = (qb      .* qb_CV).^2;            % variance; (vertical vector: num_grps X 1); NOTE: recalculated below as a dependent variable
pq_VAR                              = (pq      .* pq_CV).^2;         	% variance; (vertical vector: num_grps X 1)
ae_VAR                              = (ae      .* ae_CV).^2;          	% variance; (vertical vector: num_grps X 1)
ba_VAR                              = (ba      .* ba_CV).^2;          	% ba variance W.R.T. absolute production rate; (vertical vector: num_grps X 1)
em_VAR                              = (em      .* em_CV).^2;          	% em variance W.R.T. absolute production rate; (vertical vector: num_grps X 1)


% step 3b: consumption & diet variance values -----------------------------
CONSUMPTION_VAR                     = (CONSUMPTION          .* CONSUMPTION_CV).^2;        	% (2D matrix: num_grps X num_grps)
consumption_total_VAR               = (consumption_total    .* consumption_total_CV).^2;   	% (vertical vector: num_grps X 1)
consumption_domestic_VAR            = (consumption_domestic .* consumption_domestic_CV).^2;	% (vertical vector: num_grps X 1)
consumption_import_VAR              = (consumption_import   .* consumption_import_CV).^2;  	% (vertical vector: num_grps X 1)
DIET_VAR                            = (DIET                 .* DIET_CV).^2;                 % variance; (2D matrix: num_grps X num_grps)
diet_import_VAR                     = (import_diet          .* diet_import_CV).^2;          % variance; (horizontal vector: 1 X num_grps)


% step 3c: fleet variance values ------------------------------------------
LANDINGS_VAR                        = (LANDINGS              .* LANDINGS_CV).^2;	            % variance; (2D matrix: num_grps X num_fleets)
DISCARDS_VAR                        = (DISCARDS              .* DISCARDS_CV).^2;	            % variance; (2D matrix: num_grps X num_fleets)
CATCH_VAR                           = (CATCH                 .* CATCH_CV).^2;                   % variance; (2D matrix: num_grps X num_fleets)
landings_TotalByFleet_VAR           = (landings_TotalByFleet .* landings_TotalByFleet_CV).^2;	% (horizontal vector: 1 X num_fleets)
discards_TotalByFleet_VAR         	= (discards_TotalByFleet .* discards_TotalByFleet_CV).^2; 	% (horizontal vector: 1 X num_fleets)
catch_TotalByFleet_VAR              = (catch_TotalByFleet    .* catch_TotalByFleet_CV).^2;    	% (horizontal vector: 1 X num_fleets)
% landings_TotalByGroup_VAR           = (landings_TotalByGroup .* landings_TotalByGroup_CV).^2;	  % (vertical vector: num_AggGrps X 1); NOTE: not used here
% discards_TotalByGroup_VAR           = (discards_TotalByGroup .* discards_TotalByGroup_CV).^2;	  % (vertical vector: num_AggGrps X 1); NOTE: not used here
% catch_TotalByGroup_VAR              = (catch_TotalByGroup    .* catch_TotalByGroup_CV).^2;      % (vertical vector: num_AggGrps X 1); NOTE: not used here


% step 3d: forced egg production & metabolism variance values -------------
% EggProduction_forced_VAR            = (EggProduction_forced  .* EggProduction_forced_CV).^2;	  % egg production variance W.R.T. absolute egg production rates; (vertical vector: num_AggGrps X 1); NOTE: FFF for future
% Metabolism_forced_VAR               = (Metabolism_forced     .* Metabolism_forced_CV).^2;       % (vertical vector: num_AggGrps X 1); NOTE: FFF for future


% step 3e: variance for EwE eggs-as-detritus ------------------------------
%          NOTE: as of now, this value is defined as a single value for all groups
ee_eggs_CV                          = repmat(ee_eggs_CV, [num_grps, 1]);	% egg CV W.R.T. production budget for all groups; (vertical vector: num_grps X 1)
ee_eggs_VAR                         = (ee_eggs .* ee_eggs_CV).^2;           % egg variance W.R.T. production budget for all groups; (vertical vector: num_grps X 1)


% step 3f: variance for implicit bacterial metabolism of detritus ---------
plgc_BacterialMTBLSM_VAR            = (PelagicBacterialReduction .* BacterialMTBLSM_CV).^2; % (scaler); fraction of terminal pelagic detritus oxidized to pelagic NH4
bnth_BacterialMTBLSM_VAR            = (BenthicBacterialReduction .* BacterialMTBLSM_CV).^2; % (scaler); fraction of terminal pelagic detritus oxidized to benthic NH4


% step 3g: variance for abiotic oxidation of NH4 to NO3 -------------------
Oxidation_NH4_CV                    = repmat(Oxidation_NH4_CV, [num_NH4, 1]);	% CV of fraction of NH4 produced oxidized directly back to NO3 abiologically; (vertical vector: num_NH4 X 1)
Oxidation_NH4_VAR                   = (Oxidation_NH4 .* Oxidation_NH4_CV).^2;	% variance of fraction of NH4 produced oxidized directly back to NO3 abiologically; (vertical vector: num_NH4 X 1)


% step 3h: variance for E2E nutrient uptake by phytoplankton --------------
%          FFF: in future, can define different uncertainty for NO3 & NH4 uptake
% % % NO3Uptake_CV                        = repmat(NutrientUptake_CV, [num_PrimaryProducers, num_NO3]);	% NO3 uptake CV W.R.T. energy budget for all groups; (2D matrix: num_PrimaryProducers X num_NO3)
% % % NO3Uptake_VAR                       = (PhytoUptake_NO3 .* NO3Uptake_CV).^2;                         % NO3 uptake variance W.R.T. energy budget for all groups; (2D matrix: num_PrimaryProducers X num_NO3)
% % % 
% % % PhytoUptake_NH4_repmat              = repmat(PhytoUptake_NH4, [1, num_NH4]); % (2D matrix: num_PrimaryProducers X num_NH4)
% % % NH4Uptake_CV                        = repmat(NutrientUptake_CV, [num_PrimaryProducers, num_NH4]);	% NH4 uptake CV W.R.T. energy budget for all groups; (2D matrix: num_PrimaryProducers X num_NH4)
% % % NH4Uptake_VAR                       = (PhytoUptake_NH4_repmat .* NH4Uptake_CV).^2;                  % NH4 uptake variance W.R.T. energy budget for all groups; (2D matrix: num_PrimaryProducers X num_NH4)

PhytoUptake_NO3_repmat              = repmat(PhytoUptake_NO3, [1, num_NO3]);                        % (2D matrix: num_PrimaryProducers X num_NO3)
PhytoUptake_plgcNH4_repmat          = repmat(PhytoUptake_NH4, [1, num_plgcNH4]);                    % (2D matrix: num_PrimaryProducers X num_plgcNH4)
PhytoUptake_bnthNH4_repmat          = repmat(PhytoUptake_NH4, [1, num_bnthNH4]);                    % (2D matrix: num_PrimaryProducers X num_bnthNH4)

PhytoUptake_NO3_CV               	= repmat(NutrientUptake_CV, [num_PrimaryProducers, num_NO3]);	% NO3 uptake CV W.R.T. energy budget for all groups; (2D matrix: num_PrimaryProducers X num_NO3)
PhytoUptake_NO3_VAR               	= (PhytoUptake_NO3_repmat .* PhytoUptake_NO3_CV).^2;            % NO3 uptake variance W.R.T. energy budget for all groups; (2D matrix: num_PrimaryProducers X num_NO3)

PhytoUptake_plgcNH4_CV           	= repmat(NutrientUptake_CV, [num_PrimaryProducers, num_plgcNH4]);	% plgcNH4 uptake CV W.R.T. energy budget for all groups; (2D matrix: num_PrimaryProducers X num_plgcNH4)
PhytoUptake_plgcNH4_VAR          	= (PhytoUptake_plgcNH4_repmat .* PhytoUptake_plgcNH4_CV).^2;        % plgcNH4 uptake variance W.R.T. energy budget for all groups; (2D matrix: num_PrimaryProducers X num_plgcNH4)

PhytoUptake_bnthNH4_CV           	= repmat(NutrientUptake_CV, [num_PrimaryProducers, num_bnthNH4]);	% bnthNH4 uptake CV W.R.T. energy budget for all groups; (2D matrix: num_PrimaryProducers X num_bnthNH4)
PhytoUptake_bnthNH4_VAR          	= (PhytoUptake_bnthNH4_repmat .* PhytoUptake_bnthNH4_CV).^2;        % bnthNH4 uptake variance W.R.T. energy budget for all groups; (2D matrix: num_PrimaryProducers X num_bnthNH4)
% *************************************************************************





% *************************************************************************
% STEP 4: calculate derived variance terms---------------------------------

% step 4a: derived qb variance --------------------------------------------
%          NOTE: qb = pb/pq
qb_VAR                      = f_VarianceDivision_12132018(pb, pb_VAR, pq, pq_VAR); % variance; (vertical vector: num_EwEgroups X 1)
qb_VAR(looky_fleets)        = 0; % qb for fleets always = 1 beacause q = b, so qb_VAR for fleets = 0; (vertical vector: num_grps X 1)


% step 4b: derived predation variance -------------------------------------
%          NOTE: predation_total (M2 = sigma_c(D_pc * q_c))
predation               	= sum(CONSUMPTION(:, looky_livingANDfleets), 2);        % total predation on each p; (vertical vector: num_grps X 1)
repmat_predation        	= repmat(predation, [1, num_grps]);                     % (2D matrix: num_EwEgroups X num_EwEgroups)
predation_VAR           	= sum(CONSUMPTION_VAR(:, looky_livingANDfleets), 2);	% (vertical vector: num_grps X 1)
repmat_predation_VAR        = repmat(predation_VAR, [1, num_grps]);                 % (2D matrix: num_EwEgroups X num_EwEgroups)


% step 4c: derived production variance ------------------------------------
%          NOTE: production = pb * biomass
production_VAR           	= f_VarianceMultiplication_12132018(pb, pb_VAR, biomass, biomass_VAR); % variance; (vertical vector: num_EwEgroups X 1)
production_VAR(looky_nutrients)	= 0; % NO3 is a driver & NH4 pools are highly derived from the rest of the model so no production_VAR is assigned


% step 4d: derived ee variance terms --------------------------------------
%          express predation, ba, em variance W.R.T. production budget

% ee_predation = predation / production
ee_predation_VAR        	= f_VarianceDivision_12132018(predation, predation_VAR, production, production_VAR); % predation variance W.R.T. production budget; (vertical vector: num_EwEgroups X 1)

% ee_ba = ba / production
ee_ba_VAR               	= f_VarianceDivision_12132018(ba, ba_VAR, production, production_VAR); % ba variance W.R.T. production budget; (vertical vector: num_EwEgroups X 1)
                  
% ee_em = em / production
ee_em_VAR                	= f_VarianceDivision_12132018(em, em_VAR, production, production_VAR); % em variance W.R.T. production budget; (vertical vector: num_EwEgroups X 1)

% ee = e_eggs + ee_predation + ee_ba + ee_em
ee_VAR                    	= ee_eggs_VAR + ee_predation_VAR + ee_ba_VAR + ee_em_VAR;	% ee variance W.R.T. production budget; (vertical vector: num_EwEgroups X 1)
% *************************************************************************





% *************************************************************************
% STEP 5: derived BioenergeticBudget variance------------------------------
%         fraction of total consumption (Q) going to metabolism, feces (non-assimilated consumption), or production
%               1) feces
%               2) metabolism
%               3) production

% step 5a: fraction of consumption going to metabolism ---------------------
metabolism                      = ae - pq;        	% metabolism W.R.T. consumption budget (metabolism = ae - pq); (vertical vector: num_grps X 1)
metabolism_VAR                  = ae_VAR + pq_VAR;	% metabolism variance W.R.T. consumption budget; (vertical vector: num_grps X 1)
metabolism_VAR(looky_fleets)	= 0;                % fleet metabolism always = 0
metabolism_VAR(looky_terminalPLGCdetritus) = plgc_BacterialMTBLSM_VAR; % fraction of terminal pelagic detritus oxidized to pelagic NH4; (vertical vector: num_grps X 1)
metabolism_VAR(looky_terminalBNTHdetritus) = bnth_BacterialMTBLSM_VAR; % fraction of terminal pelagic detritus oxidized to benthic NH4; (vertical vector: num_grps X 1)
metabolism_VAR(looky_NH4)       = Oxidation_NH4_VAR; % variance of fraction of NH4 produced oxidized directly back to NO3 abiologically; (vertical vector: num_grps X 1)


% step 5b: fraction of consumption going to feces -------------------------
feces                           = 1 - ae;        	% feces W.R.T. consumption budget; (vertical vector: num_grps X 1)
feces_VAR                       = ae_VAR;         	% feces variance W.R.T. consumption budget; (vertical vector: num_grps X 1)

% step 5c: some error checking (should not be necessary) ------------------
BioenergeticBudget_check(1, :)	= feces';       % fraction of consumption going to feces; (horizontal vector: 1 X num_grps); NOTE transpose
BioenergeticBudget_check(2, :)	= metabolism';	% fraction of consumption going to metabolism (NH4 production); metabolism = 1 - pq - feces = (1 - (pq + (1-ae))); (horizontal vector: 1 X num_grps); NOTE transpose
BioenergeticBudget_check(3, :)  = pq';          % fraction of consumption going to production; (horizontal vector: 1 X num_grps); NOTE transpose
BioenergeticBudget_check        = BioenergeticBudget - BioenergeticBudget_check;
BioenergeticBudget_check        = sum(abs([min(min(BioenergeticBudget_check)) max(max(BioenergeticBudget_check))]));
if BioenergeticBudget_check ~= 0
    display(['WARNING: BioenergeticBudget error in ' fname_E2Epedigree]) % error checking that BioenergeticBudget is not corrupt
end


% step 5d: BioenergeticBudget variance ------------------------------------
BioenergeticBudget_VAR(1, :)    = feces_VAR';       % (horizontal vector: 1 X num_grps); NOTE transpose
BioenergeticBudget_VAR(2, :)    = metabolism_VAR';	% (horizontal vector: 1 X num_grps); NOTE transpose
BioenergeticBudget_VAR(3, :)    = pq_VAR';          % (horizontal vector: 1 X num_grps); NOTE transpose
% *************************************************************************





% *************************************************************************
% STEP 6: derived ProductionBudget variance--------------------------------
%         fraction of production going to eggs & reproduction, predation, senescence, biomass accumulation, or emigration (incl. immigration, physical transport)
%               1) eggs (reproduction)
%            	2) predation
%             	3) senescence
%               4) ba (biomass accumulation)
%              	5) em (emigration); NOTE: negative for immigration
ProductionBudget_VAR(1, :)    = ee_eggs_VAR';       % eggs (&/or reproduction) variance W.R.T. production budget; (horizontal vector: 1 X num_grps); NOTE transpose
ProductionBudget_VAR(2, :)    = ee_predation_VAR';	% predation variance W.R.T. production budget; (horizontal vector: 1 X num_grps); NOTE transpose
ProductionBudget_VAR(3, :)    = ee_VAR';            % senescence variance W.R.T. production budget; (horizontal vector: 1 X num_grps); senescence = 1 - ee, therefore variance of senescence = variance of ee; NOTE transpose
ProductionBudget_VAR(4, :)    = ee_ba_VAR';         % ba variance W.R.T. production budget; (horizontal vector: 1 X num_grps); NOTE transpose
ProductionBudget_VAR(5, :)    = ee_em_VAR';         % em variance W.R.T. production budget; (horizontal vector: 1 X num_grps); NOTE transpose
% *************************************************************************





% *************************************************************************
% STEP 7: derived ConsumptionBudget variance-------------------------------
%         fraction of consumption going to feces, metabolism, eggs, predation, senescence, biomass accumulation, or emigration (incl. immigration, physical transport)
%               1) feces
%            	2) metabolism
%            	3) eggs (reproduction)
%           	4) predation
%           	5) senescence
%            	6) ba (biomass accumulation)
%            	7) em (emigration); NOTE: negative for immigration

% step 7a: derived variance terms -----------------------------------------
%          convert variance to W.R.T. production budget to variance W.R.T. consumption budget
%          NOTE: var(A*B)  = [mean(A*B)^2] * [(var(A)/mean(A)^2) + (var(B)/mean(B)^2)]
%          NOTE: "cb" prefix stands for "consumption budget"

% cb_eggs = ee_eggs * pq
cb_eggs_VAR                 = f_VarianceMultiplication_12132018(ee_eggs',      ee_eggs_VAR',      pq', pq_VAR'); % egg (&/or reproduction) variance W.R.T. consumption budget; (horizontal vector: 1 X num_EwEgroups); NOTE transposes

% cb_predation = ee_predation * pq
cb_predation_VAR            = f_VarianceMultiplication_12132018(ee_predation', ee_predation_VAR', pq', pq_VAR'); % predation variance W.R.T. consumption budget; (horizontal vector: 1 X num_EwEgroups); NOTE transposes

% cb_senescence = (1-ee) * pq
cb_senescence_VAR           = f_VarianceMultiplication_12132018((1-ee)',       ee_VAR',           pq', pq_VAR'); % senescence variance W.R.T. consumption; (horizontal vector: 1 X num_EwEgroups); NOTE transposes

% cb_ba = ee_ba * pq
cb_ba_VAR                   = f_VarianceMultiplication_12132018(ee_ba',        ee_ba_VAR',        pq', pq_VAR'); % ba variance W.R.T. consumption; (horizontal vector: 1 X num_EwEgroups); NOTE transposes

% cb_em = ee_em * pq
cb_em_VAR                   = f_VarianceMultiplication_12132018(ee_em',        ee_em_VAR',        pq', pq_VAR'); % em variance W.R.T. consumption; (horizontal vector: 1 X num_EwEgroups); NOTE transposes


% step 7b: build ConsumptionBudget_VAR matrix -----------------------------
ConsumptionBudget_VAR(1, :)	= feces_VAR';         	% feces variance W.R.T. consumption budget; (horizontal vector: 1 X num_grps); NOTE transpose
ConsumptionBudget_VAR(2, :)	= metabolism_VAR';     	% metabolism variance W.R.T. consumption budget; (horizontal vector: 1 X num_grps); NOTE transpose
ConsumptionBudget_VAR(3, :)	= cb_eggs_VAR;        	% egg (&/or reproduction) variance W.R.T. consumption; (horizontal vector: 1 X num_grps)
ConsumptionBudget_VAR(4, :)	= cb_predation_VAR;   	% predation variance W.R.T. consumption; (horizontal vector: 1 X num_grps)
ConsumptionBudget_VAR(5, :)	= cb_senescence_VAR;	% senescence variance W.R.T. consumption; (horizontal vector: 1 X num_grps)
ConsumptionBudget_VAR(6, :)	= cb_ba_VAR;          	% ba variance W.R.T. consumption; (horizontal vector: 1 X num_grps)
ConsumptionBudget_VAR(7, :)	= cb_em_VAR;           	% em variance W.R.T. consumption; (horizontal vector: 1 X num_grps)
% *************************************************************************





% *************************************************************************
% STEP 8: derived EnergyBudget variance------------------------------------

% step 8a: initialize EnergyBudget_VAR ------------------------------------
%          NOTE: A_cp = Q_cp / M2_p
%          NOTE: PredationBudget: for each producer, p, and consumer, c: ((b_pc * Q_c) / M2_p)) = the fraction of total predation on each producer, p, going to each consumer, c.
%          NOTE: PredationBudget = CONSUMPTION ./ repmat_predation
PredationBudget_VAR     = f_VarianceDivision_12132018(CONSUMPTION, CONSUMPTION_VAR, repmat_predation, repmat_predation_VAR); % variance; (2D matrix: num_grps X num_grps)
EnergyBudget_VAR        = PredationBudget_VAR'; % initialize & transpose to A_cp orientation; (2D matrix: num_grps X num_grps)


% step 8b: scale predator rows in EnergyBudget_VAR by predation term in ConsumptionBudget
cb_predation            = ConsumptionBudget(4, :);                  % predation fraction in ConsumptionBudget; (horizontal matrix: 1 X num_grps)
repmat_cb_predation     = repmat(cb_predation, [num_grps, 1]);      % (2D matrix: num_grps X num_grps)
repmat_cb_predation_VAR = repmat(cb_predation_VAR, [num_grps, 1]);	% (2D matrix: num_grps X num_grps)
EnergyBudget_VAR        = f_VarianceMultiplication_12132018(EnergyBudget, EnergyBudget_VAR, repmat_cb_predation, repmat_cb_predation_VAR); % variance; (2D matrix: num_grps X num_grps)


% step 8c: paste nutrient flow to primary producers into EnergyBudget_VAR -
EnergyBudget_VAR(looky_ANYPrimaryProd, looky_NO3)       = PhytoUptake_NO3_VAR;
EnergyBudget_VAR(looky_ANYPrimaryProd, looky_plgcNH4)	= PhytoUptake_plgcNH4_VAR;
EnergyBudget_VAR(looky_ANYPrimaryProd, looky_bnthNH4)	= PhytoUptake_bnthNH4_VAR;


% step 8d: paste eggs (& reproduction) term from ConsumptionBudget_VAR into EnergyBudget_VAR
repmat_cb_eggs_VAR                      	= repmat(cb_eggs_VAR, [num_eggs, 1]);	% egg (&/or reproduction) variance W.R.T. consumption; (2D matrix: num_eggs X num_grps)
repmat_cb_eggs_VAR                          = repmat_cb_eggs_VAR .* fate_eggs.^2;	% redistribute among egg rows as defined by fate_eggs; NOTE: scale fate_eggs by ^2; (2D matrix: num_eggs X num_grps); NOTE: fate_eggs are usually either 0 or 1
EnergyBudget_VAR(looky_eggs, :)          	= repmat_cb_eggs_VAR;                   % paste egg terms into EnergyBudget_VAR


% step 8e: paste feces & senescence terms from ConsumptionBudget into appropriate BigMatrix detritus rows
repmat_feces_VAR                            = repmat(feces_VAR', [num_ANYdetritus, 1]);         % feces variance W.R.T. consumption budget; (2D matrix: num_ANYdetritus X num_grps); NOTE transpose
repmat_feces_VAR                            = repmat_feces_VAR .* fate_feces.^2;              	% redistribute feces among detritus rows as defined by fate_feces; NOTE: scale fate_feces by ^2; (2D matrix: num_ANYdetritus X num_grps)
repmat_cb_senescence_VAR                    = repmat(cb_senescence_VAR, [num_ANYdetritus, 1]);	% senescence variance W.R.T. consumption; (2D matrix: num_ANYdetritus X num_grps)
repmat_cb_senescence_VAR                	= repmat_cb_senescence_VAR .* fate_senescence.^2;   % redistribute senescence among detritus rows as defined by fate_senescence; NOTE: scale fate_senescence by ^2; (2D matrix: num_ANYdetritus X num_grps)
EnergyBudget_VAR(looky_ANYdetritus, :)      = repmat_feces_VAR + repmat_cb_senescence_VAR;      % paste detritus rows back into EnergyBudget_VAR; (2D matrix: num_grps X num_grps)


% step 8f: paste metabolism term from ConsumptionBudget into BigMatrix ----
%          NOTE: includes abiotic oxidation of NH4 to NO3
repmat_metabolism_VAR                       = repmat((metabolism_VAR'), [num_nutrients, 1]);	% metabolism variance W.R.T. consumption budget; (2D matrix: num_nutrients X num_grps); NOTE transpose
repmat_metabolism_VAR                       = repmat_metabolism_VAR .* fate_metabolism.^2;    	% redistribute metabolism among nutrient rows as defined by fate_metabolism; NOTE: scale fate_metabolism by ^2; (2D matrix: num_nutrients X num_grps)
EnergyBudget_VAR(looky_nutrients, :)     	= repmat_metabolism_VAR;                            % paste metabolism terms into EnergyBudget_VAR

% QQQ--->>???    BUT, still need to account for the special case of terminal
%                benthic detritus group(s) where predator rows were scaled by 
%                predation + senescence so that their column sum(s) add to 1
% *************************************************************************





% *************************************************************************
% STEP 9: derived discard fraction variance for each group and fleet-------
%         QQQ --> think I fixed this problem?-->NOTE: for later calculations with Monte Carlo models, realize
%               that the fleet pb, ae, ConsumptionBudget, and EnergyBudget
%               columns will NOT reflect the discard rates derived from the sum
%               of the discard rates of individual groups. This is because we do
%               not at this point know the production rates, and therefore we do 
%               not know the catch rates of individual groups. The
DiscardFraction_VAR     = f_VarianceDivision_12132018(DISCARDS, DISCARDS_VAR, CATCH, CATCH_VAR); % variance; (2D matrix: num_grps X num_fleets)
% *************************************************************************







% *************************************************************************
% STEP 9: express uncertainties as CV terms -------------------------------
biomass_CV                                          = (sqrt(biomass_VAR)            ./ biomass)';           % CV; (horizontal vector: 1 X num_grps); NOTE transpose
production_CV                                       = (sqrt(production_VAR)         ./ production)';        % CV; (horizontal vector: 1 X num_grps); NOTE transpose
pb_CV                                               = (sqrt(pb_VAR)                 ./ pb)';                % CV; (horizontal vector: 1 X num_grps); NOTE transpose
qb_CV                                               = (sqrt(qb_VAR)                 ./ qb)';                % CV; (horizontal vector: 1 X num_grps); NOTE transpose
pq_CV                                               = (sqrt(pq_VAR)                 ./ pq)';                % CV; (horizontal vector: 1 X num_grps); NOTE transpose
ae_CV                                               = (sqrt(ae_VAR)                 ./ ae)';                % CV; (horizontal vector: 1 X num_grps); NOTE transpose

BioenergeticBudget_CV                               = (sqrt(BioenergeticBudget_VAR) ./ BioenergeticBudget);     % CV; (3D matrix: 3 X num_grps)
ProductionBudget_CV                                 = (sqrt(ProductionBudget_VAR)   ./ abs(ProductionBudget));	% CV; (3D matrix: 5 X num_grps); NOTE use of abs() because ba & em terms may be negative
ConsumptionBudget_CV                                = (sqrt(ConsumptionBudget_VAR)  ./ abs(ConsumptionBudget));	% CV; (3D matrix: 7 X num_grps); NOTE use of abs() because ba & em terms may be negative
EnergyBudget_CV                                     = (sqrt(EnergyBudget_VAR)       ./ EnergyBudget);           % CV; (3D matrix: num_grps X num_grps)
DiscardFraction_CV                                  = (sqrt(DiscardFraction_VAR)    ./ DiscardFraction);        % CV; (2D matrix: num_grps X num_fleets)

biomass_CV(isnan(biomass_CV))                       = 0; % replace div/0 NaN values with 0; (horizontal vector: 1 X num_grps)
production_CV(isnan(production_CV))                 = 0; % replace div/0 NaN values with 0; (horizontal vector: 1 X num_grps)
pb_CV(isnan(pb_CV))                                 = 0; % replace div/0 NaN values with 0; (horizontal vector: 1 X num_grps)
qb_CV(isnan(qb_CV))                                 = 0; % replace div/0 NaN values with 0; (horizontal vector: 1 X num_grps)
pq_CV(isnan(pq_CV))                                 = 0; % replace div/0 NaN values with 0; (horizontal vector: 1 X num_grps)
ae_CV(isnan(ae_CV))                                 = 0; % replace div/0 NaN values with 0; (horizontal vector: 1 X num_grps)

BioenergeticBudget_CV(isnan(BioenergeticBudget_CV))	= 0; % replace div/0 NaN values with 0; (3D matrix: 3 X num_grps)
ProductionBudget_CV(isnan(ProductionBudget_CV))     = 0; % replace div/0 NaN values with 0; (3D matrix: 5 X num_grps)
ConsumptionBudget_CV(isnan(ConsumptionBudget_CV))	= 0; % replace div/0 NaN values with 0; (3D matrix: 7 X num_grps)
EnergyBudget_CV(isnan(EnergyBudget_CV))             = 0; % replace div/0 NaN values with 0; (3D matrix: num_grps X num_grps)
DiscardFraction_CV(isnan(DiscardFraction_CV))       = 0; % replace div/0 NaN values with 0; (3D matrix: num_grps X num_grps)

BioenergeticBudget_CV(isinf(BioenergeticBudget_CV))	= 0; % replace div/0 inf values with 0; (3D matrix: 3 X num_grps)
ProductionBudget_CV(isinf(ProductionBudget_CV))     = 0; % replace div/0 inf values with 0; (3D matrix: 5 X num_grps)
ConsumptionBudget_CV(isinf(ConsumptionBudget_CV))	= 0; % replace div/0 inf values with 0; (3D matrix: 7 X num_grps)
EnergyBudget_CV(isinf(EnergyBudget_CV))             = 0; % replace div/0 inf values with 0; (3D matrix: num_grps X num_grps)
DiscardFraction_CV(isinf(DiscardFraction_CV))       = 0; % replace div/0 inf values with 0; (3D matrix: num_grps X num_grps)
% *************************************************************************





% *************************************************************************
% STEP 10: package for export----------------------------------------------
ECOTRAN_PEDIGREE.biomass_CV             = biomass_CV;               % CV; (horizontal vector: 1 X num_grps)
ECOTRAN_PEDIGREE.production_CV          = production_CV;            % CV; (horizontal vector: 1 X num_grps)
ECOTRAN_PEDIGREE.pb_CV                  = pb_CV;                    % CV; (horizontal vector: 1 X num_grps)
ECOTRAN_PEDIGREE.qb_CV                  = qb_CV;                    % CV; (horizontal vector: 1 X num_grps)
ECOTRAN_PEDIGREE.pq_CV                  = pq_CV;                    % CV; (horizontal vector: 1 X num_grps)
ECOTRAN_PEDIGREE.ae_CV                  = ae_CV;                    % CV; (horizontal vector: 1 X num_grps)

ECOTRAN_PEDIGREE.BioenergeticBudget_CV	= BioenergeticBudget_CV;	% CV; (3D matrix: 3 X num_grps)
ECOTRAN_PEDIGREE.ProductionBudget_CV	= ProductionBudget_CV;      % CV; (3D matrix: 5 X num_grps)
ECOTRAN_PEDIGREE.ConsumptionBudget_CV	= ConsumptionBudget_CV;     % CV; (3D matrix: 7 X num_grps)
ECOTRAN_PEDIGREE.EnergyBudget_CV        = EnergyBudget_CV;          % CV; (3D matrix: num_grps X num_grps)

ECOTRAN_PEDIGREE.LANDINGS_CV                = LANDINGS_CV;                  % CV; (2D matrix: num_grps X num_fleets)
ECOTRAN_PEDIGREE.DISCARDS_CV                = DISCARDS_CV;                  % CV; (2D matrix: num_grps X num_fleets)
ECOTRAN_PEDIGREE.CATCH_CV                   = CATCH_CV;                     % CV; (2D matrix: num_grps X num_fleets)
ECOTRAN_PEDIGREE.landings_TotalByFleet_CV	= landings_TotalByFleet_CV;     % CV; (horizontal vector: 1 X num_fleets)
ECOTRAN_PEDIGREE.discards_TotalByFleet_CV	= discards_TotalByFleet_CV; 	% CV; (horizontal vector: 1 X num_fleets)
ECOTRAN_PEDIGREE.catch_TotalByFleet_CV      = catch_TotalByFleet_CV;        % CV; (horizontal vector: 1 X num_fleets)
% ECOTRAN_PEDIGREE.landings_TotalByGroup_CV	  = landings_TotalByGroup_CV;     % CV; (vertical vector: num_grps X 1); NOTE: not used here
% ECOTRAN_PEDIGREE.discards_TotalByGroup_CV	  = discards_TotalByGroup_CV;     % CV; (vertical vector: num_grps X 1); NOTE: not used here
% ECOTRAN_PEDIGREE.catch_TotalByGroup_CV      = catch_TotalByGroup_CV;        % CV; (vertical vector: num_grps X 1); NOTE: not used here
ECOTRAN_PEDIGREE.DiscardFraction_CV         = DiscardFraction_CV;           % CV; (2D matrix: num_grps X num_fleets)

ECOTRAN_PEDIGREE.fname_E2Epedigree      = fname_E2Epedigree;        % name of this f_E2Epedigree function
% *************************************************************************


% end m-file***************************************************************