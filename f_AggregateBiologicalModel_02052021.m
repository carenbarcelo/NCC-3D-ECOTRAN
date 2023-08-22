function [EwEResult, PEDIGREE] = f_AggregateBiologicalModel_02052021(dat)
% prepare EwE model and aggregate functional groups
% by Jim Ruzicka
% calculates production-weighted averages of physiological constants
% corrects for technique of defining egg production as a detritus term:
%     (since EwE does not distinguish between feces and "other mortality" detritus pools, AE must be adjusted upwards to account for eggs generated from feces detritus)
% detritus production = detritus consumption rate
% PB is recalculated based on production value that includes the detritus import rate (i.e., egg production)
%
% calls:
%       f_calcEE_12292020                   calculate Ecotrophic Efficiency
%       f_VarianceDivision_12132018         calculate the variance of one term divided by another term
%       f_VarianceMultiplication_12132018   calculate the variance of two products
%
% takes:
%		dat
%
% returns (model):
%    EwEResult.AggLabel                         group labels (vertical cell vector)
%    EwEResult.GroupType                        group type code in aggregated model (vertical vector)
%    EwEResult.FullResNums                      group numbers in fully-resolved model (vector)
%    EwEResult.EwE_code
%    EwEResult.AggNumber                        aggregate group number (vertical vector)
%    EwEResult.num_PrimaryProducers
%    EwEResult.num_eggs
%    EwEResult.num_detritus
%    EwEResult.num_EggsAndDetritus
%    EwEResult.num_NitrogenNutrients
%    EwEResult.num_fisheries
%    EwEResult.num_micrograzers
%    EwEResult.num_bacteria
%    EwEResult.Biomass                          (vertical vector)
%    EwEResult.BiomassWOegg                     (vertical vector)
%    EwEResult.Production                       summed production rate (t/km2/yr) (vertical vector)
%    EwEResult.ProductionWOegg                  (vertical vector)
%    EwEResult.PB                               (vertical vector)
%    EwEResult.QB                               (vertical vector)
%    EwEResult.PQ                               (vertical vector)
%    EwEResult.AE                               weighted mean assimilation efficiency (vertical vector)
%    EwEResult.TL                               weighted mean TL (vertical vector)
%    EwEResult.EwE_EE                           weighted mean EE (vertical vector)
%    EwEResult.BA                               (vertical vector)
%    EwEResult.EM                               (vertical vector)
%    EwEResult.Diet                             (row 1 = predator code number, column 1 = prey code number)
%    EwEResult.Consumption                      last row is import consumption (QQQ re-check this) (array) 
%    EwEResult.ImportConsumption
%    EwEResult.TotalConsumption                 total consumption by each aggregated predator group (includes import consumption) (horizontal vector)
%    EwEResult.EwE_DetritusFate
%    EwEResult.EwE_EggFractions
%    EwEResult.DetritusFate_feces
%    EwEResult.DetritusFate_senescence
%    EwEResult.ExcretionFate
%    EwEResult.ProductionLossScaler
%    EwEResult.RetentionScaler
%    EwEResult.EggProduction_forced
%    EwEResult.Metabolism_forced
%    EwEResult.PelagicBacterialReduction
%    EwEResult.BenthicBacterialReduction
%    EwEResult.Oxidation_NH4
%    EwEResult.PhytoUptake_NH4
%    EwEResult.PhytoUptake_NO3
%    EwEResult.AggLabel_nutrients
%    EwEResult.GroupType_nutrients
%    EwEResult.AggLabel_fisheries
%    EwEResult.Landings
%    EwEResult.Discards
%    EwEResult.Biomass_fisheries
%    EwEResult.BiomassWOegg_fisheries
%    EwEResult.Production_fisheries
%    EwEResult.ProductionWOegg_fisheries
%    EwEResult.PB_fisheries
%    EwEResult.QB_fisheries
%    EwEResult.PQ_fisheries
%    EwEResult.AE_fisheries
%    EwEResult.EwE_EE_fisheries
%    EwEResult.BA_fisheries
%    EwEResult.EM_fisheries
%    EwEResult.TL_fisheries
%    EwEResult.Landings_fisheries
%    EwEResult.Discards_fisheries
%    EwEResult.DetritusFate_feces_fisheries
%    EwEResult.DetritusFate_senescence_fisheries
%    EwEResult.ExcretionFate_fisheries
%    EwEResult.ProductionLossScaler_fisheries
%    EwEResult.RetentionScaler_fisheries
%    EwEResult.EggProduction_forced_fisheries
%    EwEResult.Metabolism_forced_fisheries
%    EwEResult.GroupTypeDef_ANYNitroNutr                  GroupType definitions
%    EwEResult.GroupTypeDef_NO3                           GroupType definitions
%    EwEResult.GroupTypeDef_plgcNH4                       GroupType definitions
%    EwEResult.GroupTypeDef_bnthNH4                       GroupType definitions
%    EwEResult.GroupTypeDef_ANYPrimaryProd                GroupType definitions
%    EwEResult.GroupTypeDef_LrgPhyto
%    EwEResult.GroupTypeDef_SmlPhyto
%    EwEResult.GroupTypeDef_Macrophytes
%    EwEResult.GroupTypeDef_ANYConsumer
%    EwEResult.GroupTypeDef_ConsumPlgcPlankton
%    EwEResult.GroupTypeDef_ConsumPlgcNekton
%    EwEResult.GroupTypeDef_ConsumPlgcWrmBlood
%    EwEResult.GroupTypeDef_ConsumBntcInvert
%    EwEResult.GroupTypeDef_ConsumBntcVert
%    EwEResult.GroupTypeDef_ConsumBnthWrmBlood
%    EwEResult.GroupTypeDef_eggs
%    EwEResult.GroupTypeDef_ANYDetritus
%    EwEResult.GroupTypeDef_terminalPlgcDetr
%    EwEResult.GroupTypeDef_offal
%    EwEResult.GroupTypeDef_terminalBnthDetr
%    EwEResult.GroupTypeDef_em
%    EwEResult.GroupTypeDef_ba
%    EwEResult.GroupTypeDef_fleet
%    EwEResult.GroupTypeDef_import
%    EwEResult.GroupTypeDef_micrograzers
%    EwEResult.GroupTypeDef_bacteria
% returns (Pedigree)
%    PEDIGREE.PhysiologyScaler
%    PEDIGREE.BiomassScaler
%    PEDIGREE.FisheriesScaler
%    PEDIGREE.DietScaler
%    PEDIGREE.BAScaler
%    PEDIGREE.EMScaler
%    PEDIGREE.DiscardScaler
%    PEDIGREE.minBiomass
%    PEDIGREE.minPhysiology
%    PEDIGREE.PedigreeParameters
%    PEDIGREE.PedigreeLandings
%    PEDIGREE.PedigreeDiscards
%    PEDIGREE.PedigreeDiet
%    PEDIGREE.PedigreeDietPrefRules
%    PEDIGREE.PreyGuild_name
%    PEDIGREE.PreyGuild_code
%
% revision date: 2-5-2021


% *************************************************************************
% STEP 1: get info from .dat-----------------------------------------------
fname_AggregateBiologicalModel	= mfilename; % save name of this m-file to keep in saved model results
display(['Running: ' fname_AggregateBiologicalModel])


% step 1a: get EwE biomass & physiology parameters ------------------------
%          NOTE: all in EwE model order (dimensions = num_EwEgrps = num_PreaggGrps_livingANDdetritus + num_fleets)
%          NOTE: includes fleets at end
%                       fleets have "0" place-holders (in preagg_biomass, preagg_production, preagg_pb, preagg_qb, preagg_pq, preagg_ae)
%                       fleets have values for preagg_TL
preagg_biomass                  = dat.Biomass;      % pre-aggregate EwE biomass; (t/km2);                       (vertical vector: num_EwEgrps X 1)
preagg_production               = dat.Production;	% pre-aggregate EwE production rate; (t/km2/y);             (vertical vector: num_EwEgrps X 1)
preagg_pb                       = dat.PB;           % pre-aggregate EwE growth rate; (1/y);                     (vertical vector: num_EwEgrps X 1)
preagg_qb                       = dat.QB;           % pre-aggregate EwE consumption rate; (1/y);                (vertical vector: num_EwEgrps X 1)
preagg_pq                       = dat.PQ;           % pre-aggregate EwE production efficiency; (proportion);    (vertical vector: num_EwEgrps X 1)
preagg_TL                       = dat.TL;           % pre-aggregate EwE trophic level; (unitless);              (vertical vector: num_EwEgrps X 1)
preagg_ae                       = dat.AE;           % pre-aggregate EwE assimilation efficiency; (proportion);	(vertical vector: num_EwEgrps X 1)
% -------------------------------------------------------------------------


% step 1b: get EwE diet matrix (D_pc) and import diet ---------------------
%          NOTE: in EwE order, but does not include fleets; (dimensions = num_PreaggGrps_livingANDdetritus)
%          NOTE: Diet matrix does not include import, so clms may not sum to 1; later addition of fleets and renormalization will bring sum to 1
preagg_DIET                     = dat.diet(2:end, 2:end); % pre-aggregate EwE Diet matrix (D_pc); (proportion); trim off headers;	(2D matrix: num_PreaggGrps_livingANDdetritus X num_PreaggGrps_livingANDdetritus); (NOTE: does not include import diet)
preagg_diet_import              = dat.import_diet(2:end); % pre-aggregate EwE import diet; (proportion); trim off header column;	(horizontal vector: 1 X num_PreaggGrps_livingANDdetritus)


% % step 1bb: remove cannibalism from diet matrix ---------------------------
% [preagg_DIET, toss]   	        = f_RedistributeCannibalism_11142018(preagg_DIET); % remove cannibalism terms on diagonal of DIET matrix
% display('   -->WARNING in f_AggregateResults: cannibalism is REMOVED')
% -------------------------------------------------------------------------


% step 1c: get EwE Consumption matrix (Q_pc) and import consumption -------
%          NOTE: CONSUMPTION matrix is recalculated based upon pb * biomass * DIET matrix
%                The EwE CONSUMPTION matrix is used here ONLY to get the total flow to each egg & detritus group
%          NOTE: in EwE order, but does not include fleets (dimensions = num_PreaggGrps_livingANDdetritus)
EwE_CONSUMPTION              = dat.consumption(2:end, 2:end);	% pre-aggregate EwE Consumption matrix (Q_pc); (t/km2/y); trim off headers; (2D matrix: num_PreaggGrps_livingANDdetritus X num_PreaggGrps_livingANDdetritus); (NOTE: does not include fleets nor import diet consumption)
% -------------------------------------------------------------------------


% step 1d: get EwE fleet data ---------------------------------------------
%          NOTE: all in EwE model order (dimensions = num_PreaggGrps_livingANDdetritus & num_fleets)
%          NOTE: does NOT include fleets at end
preagg_LANDINGS                 = dat.landings;                         % (t/km2/y); (2D matrix: num_PreaggGrps_livingANDdetritus X num_fleets)
preagg_DISCARDS                 = dat.discards;                         % (t/km2/y); (2D matrix: num_PreaggGrps_livingANDdetritus X num_fleets)
preagg_CATCH                    = preagg_LANDINGS + preagg_DISCARDS;    % (t/km2/y); (2D matrix: num_PreaggGrps_livingANDdetritus X num_fleets)
preagg_catch_TotalByFleet   	= sum(preagg_CATCH, 1);                 % (t/km2/y); (horizontal vector: 1 X num_fleets); (this vector will get written over after aggregation)
preagg_landings_TotalByFleet  	= sum(preagg_LANDINGS, 1);              % (t/km2/y); (horizontal vector: 1 X num_fleets); 
% -------------------------------------------------------------------------


% step 1e: get EwE detritus fates -----------------------------------------
%          NOTE: all in EwE model order (dimensions = num_EwEgrps = num_PreaggGrps_livingANDdetritus + num_fleets & num_EwE_DetritusGrps)
%          NOTE: includes fleets at end
%          NOTE: may include eggs as detritus groups
preagg_fate_EwEdetritus         = dat.EwE_DetritusFate;     % pre-aggregate EwE detritus fate definitions (not the ECOTRAN definitions); (2D matrix: num_EwEgrps X num_EwE_DetritusGrps)
% -------------------------------------------------------------------------


% step 1f: get info from "EcotranType" worksheet --------------------------
%          NOTE: all in ECOTRAN model order
preagg_E2E_NumberCode         	= dat.EcotranNumberCode;    % pre-aggregate E2E code numbers of fully resolved trophic groups;	(vertical vector: num_E2Egrps X 1)
preagg_GroupTypes             	= dat.EcotranGroupType;     % pre-aggregate E2E group type;                                     (vertical vector: num_E2Egrps X 1); (NOTE: BA & EM rows not read in from this sheet)
preagg_AggName                  = dat.AggEwEName;           % pre-aggregate E2E agg group names;                                (vertical vector: num_E2Egrps X 1)
PreaggCodes                     = dat.AggCode;              % E2E aggregate group code numbers;                                 (vertical vector: num_E2Egrps X 1); (NOTE: agg code numbers for all grps in resolved model prior to aggregation)
AggCodes                        = unique(PreaggCodes);      % unique aggregate group code numbers;                              (vertical vector: num_AggGrps X 1)
num_AggGrps                     = length(AggCodes);         % number of aggregated groups; E2E model

GroupTypeDef_ANYNitroNutr       = dat.GroupTypeDef_ANYNitroNutr;
GroupTypeDef_NO3                = dat.GroupTypeDef_NO3;
GroupTypeDef_plgcNH4            = dat.GroupTypeDef_plgcNH4;
GroupTypeDef_bnthNH4            = dat.GroupTypeDef_bnthNH4;
GroupTypeDef_ANYPrimaryProd     = dat.GroupTypeDef_ANYPrimaryProd;
GroupTypeDef_LrgPhyto           = dat.GroupTypeDef_LrgPhyto;
GroupTypeDef_SmlPhyto           = dat.GroupTypeDef_SmlPhyto;
GroupTypeDef_Macrophytes        = dat.GroupTypeDef_Macrophytes;
GroupTypeDef_ANYConsumer        = dat.GroupTypeDef_ANYConsumer;
GroupTypeDef_ConsumPlgcPlankton = dat.GroupTypeDef_ConsumPlgcPlankton;
GroupTypeDef_ConsumPlgcNekton   = dat.GroupTypeDef_ConsumPlgcNekton;
GroupTypeDef_ConsumPlgcWrmBlood = dat.GroupTypeDef_ConsumPlgcWrmBlood;
GroupTypeDef_ConsumBntcInvert   = dat.GroupTypeDef_ConsumBntcInvert;
GroupTypeDef_ConsumBntcVert     = dat.GroupTypeDef_ConsumBntcVert;
GroupTypeDef_ConsumBnthWrmBlood = dat.GroupTypeDef_ConsumBnthWrmBlood;
GroupTypeDef_eggs               = dat.GroupTypeDef_eggs;
GroupTypeDef_ANYDetritus        = dat.GroupTypeDef_ANYDetritus;
GroupTypeDef_terminalPlgcDetr   = dat.GroupTypeDef_terminalPlgcDetr;
GroupTypeDef_offal              = dat.GroupTypeDef_offal;
GroupTypeDef_terminalBnthDetr   = dat.GroupTypeDef_terminalBnthDetr;
GroupTypeDef_ba                 = dat.GroupTypeDef_BA;
GroupTypeDef_em                 = dat.GroupTypeDef_EM;
GroupTypeDef_fleet              = dat.GroupTypeDef_fishery;
GroupTypeDef_import             = dat.GroupTypeDef_import;
GroupTypeDef_micrograzers       = dat.GroupTypeDef_micrograzers;
GroupTypeDef_bacteria           = dat.GroupTypeDef_bacteria;
GroupSet_consumers              = [GroupTypeDef_ANYConsumer; GroupTypeDef_fleet];
GroupSet_detritus               = [GroupTypeDef_ANYDetritus; GroupTypeDef_eggs];
% -------------------------------------------------------------------------


% step 1g: get info from "EcotranRecycling" worksheet ---------------------
%          NOTE: in E2E ECOTRAN order with nutrients & fleets (dimensions = num_E2Egrps)
%          NOTE: detritus & metabolism fate terms are all oriented horizontally
%          NOTE: ProductionLossScaler, RetentionScaler, FunctionalResponseParams, & FunctionalResponse_matrix are all oriented vertically
preagg_ba                           = dat.BA;                           % pre-aggregate E2E biomass accumulation rate; (t/km2/y); (vertical vector: 1 X num_E2Egrps); (NOTE: absolute rate units)
preagg_em                           = dat.EM;                           % pre-aggregate E2E emigration rate; (t/km2/y); (vertical vector: 1 X num_E2Egrps); (NOTE: absolute rate units)
preagg_EggProduction_forced         = dat.EggProduction;            	% pre-aggregate E2E defined egg, gamete, or live-birth production (not defined as EwE detritus); (t/km2/y); (vertical vector: num_E2Egrps X 1); (FFF for future)
preagg_Metabolism_forced            = dat.Metabolism;               	% pre-aggregate E2E defined metabolism (not defined by EwE parameters); (t/km2/y); (vertical vector: num_E2Egrps X 1); (FFF for future)
preagg_fate_feces                   = dat.DetritusFate_feces';      	% pre-aggregate E2E feces fate; (units = 0 to 1, fractions of all detritus);            (2D matrix: num_TerminalDetritus(=2) X num_E2Egrps); NOTE transpose
preagg_fate_senescence              = dat.DetritusFate_senescence';  	% pre-aggregate E2E senescence fate; (units = 0 to 1, fractions of all detritus);       (2D matrix: num_TerminalDetritus(=2) X num_E2Egrps); NOTE transpose
preagg_fate_metabolism              = dat.ExcretionFate';           	% pre-aggregate E2E metabolic excretion fate; (units = 0 to 1, fractions of all metabolic excretion); (2D matrix: num_NH4 pools(=2) X num_E2Egrps); NOTE transpose

preagg_ProductionLossScaler         = dat.ProductionLossScaler;         % pre-aggregate E2E production loss term for steady state analyses of physical export; (scaler: 0-1); (vertical vector: num_E2Egrps X 1)
preagg_RetentionScaler              = dat.RetentionScaler;              % pre-aggregate E2E plankton retention scaler; (units = 0 for purely planktonic, 1 for purely sessile or nektonic); (vertical vector: num_E2Egrps X 1)
preagg_FunctionalResponseParams     = dat.FunctionalResponseParams;     % pre-aggregate E2E functional response parameters, 4 parameters;                       (2D matrix: num_E2Egrps X 4); (NOTE: each column could have variable definitions)
preagg_FunctionalResponse_matrix	= dat.FunctionalResponse_matrix;	% pre-aggregate E2E functional response parameters, matrix of grp-to-grp relationships; (2D matrix: num_E2Egrps X num_E2Egrps);
% -------------------------------------------------------------------------


% step 1h: get PEDIGREE info ----------------------------------------------
%          NOTE: PedigreeParameters Coefficients of Variation (CV); (2D matrix: num_PreaggGrps_livingANDdetritus X 37)
%                max EE; B (CV); PB (CV); QB (CV); PQ (CV); AE (CV); 
%                B (abs min); PB (abs min); QB (abs min); PQ (abs min); 
%                AE (abs min); B (abs max); PB (abs max); QB (abs max); 
%                PQ (abs max); AE (abs max); BA (CV); EM (CV); Egg rate (CV); 
%                Metabolism rate (CV); BA (abs min); EM (abs min); 
%                Egg rate (abs min); Metabolism rate (abs min); BA (abs max); 
%                EM (abs max); Egg rate (abs max); Metabolism rate (abs max); 
%                DetritusFate_feces (CV); DetritusFate_senescence (CV); 
%                ExcretionFate (CV); DetritusFate_feces (min); 
%                DetritusFate_senescence (min); ExcretionFate (min); 
%                DetritusFate_feces (max); DetritusFate_senescence (max); 
%                ExcretionFate (max)
%           NOTE: these terms no longer used:
%                PedigreeDietPrefRules, import_PedigreeDietPrefRules, 
%                Pedigree_PhysiologyScaler, Pedigree_BiomassScaler
%                Pedigree_FisheriesScaler, Pedigree_DietScaler
%                Pedigree_BAScaler, Pedigree_EMScaler, Pedigree_DiscardScaler,
%                Pedigree_minBiomass, Pedigree_minPhysiology
preagg_PedigreeParameters           = dat.PedigreeParameters;           % Coefficients of Variation (CV); (2D matrix: num_PreaggGrps_livingANDdetritus X 37)
preagg_PEDIGREE_LANDINGS            = dat.PedigreeLandings;             % Coefficients of Variation (CV); (2D matrix: num_PreaggGrps_livingANDdetritus X num_fleets)
preagg_PEDIGREE_DISCARDS            = dat.PedigreeDiscards;             % Coefficients of Variation (CV); (2D matrix: num_PreaggGrps_livingANDdetritus X num_fleets)
preagg_PEDIGREE_DIET                = dat.PedigreeDiet(2:end, 2:end);	% pre-aggregate diet pedigree matrix; (CV); trim off headers; (2D matrix: num_PreaggGrps_livingANDdetritus X num_PreaggGrps_livingANDdetritus); NOTE: does not include import diet
preagg_PEDIGREE_diet_import         = dat.import_PedigreeDiet(2:end);  	% pre-aggregate EwE import dietpedigree matrix; (CV); trim off header column; (horizontal vector: 1 X num_PreaggGrps_livingANDdetritus)
% *************************************************************************





% *************************************************************************
% STEP 2: put all into ECOTRAN E2E order-----------------------------------

% step 2a: addresses in pre-aggregate model -------------------------------
%          NOTE: these addresses are for E2E model with nutrients & fleets
looky_NO3                       	= find(preagg_GroupTypes == GroupTypeDef_NO3);
looky_plgcNH4                       = find(preagg_GroupTypes == GroupTypeDef_plgcNH4);
looky_bnthNH4                       = find(preagg_GroupTypes == GroupTypeDef_bnthNH4);
looky_NH4                           = find(preagg_GroupTypes == GroupTypeDef_plgcNH4 | preagg_GroupTypes == GroupTypeDef_bnthNH4);
looky_nutrients                     = find(floor(preagg_GroupTypes) == GroupTypeDef_ANYNitroNutr);	% row addresses of E2E nutrients
looky_ANYPrimaryProd                = find(floor(preagg_GroupTypes) == GroupTypeDef_ANYPrimaryProd);
looky_ANYconsumer                   = find(floor(preagg_GroupTypes) == GroupTypeDef_ANYConsumer);
looky_micrograzers                  = find(preagg_GroupTypes == GroupTypeDef_micrograzers);
looky_bacteria                      = find(floor(preagg_GroupTypes) == GroupTypeDef_bacteria);
looky_eggs                          = find(preagg_GroupTypes == GroupTypeDef_eggs);
looky_ANYdetritus                   = find(floor(preagg_GroupTypes) == GroupTypeDef_ANYDetritus);
looky_terminalPLGCdetritus          = find(preagg_GroupTypes == GroupTypeDef_terminalPlgcDetr);
looky_terminalBNTHdetritus          = find(preagg_GroupTypes == GroupTypeDef_terminalBnthDetr);
looky_terminalANYdetritus           = find(preagg_GroupTypes == GroupTypeDef_terminalPlgcDetr | preagg_GroupTypes == GroupTypeDef_terminalBnthDetr);
looky_eggsANDdetritus               = sort([looky_ANYdetritus; looky_eggs]);
looky_livingANDdetritus             = sort([looky_ANYPrimaryProd; looky_ANYconsumer; looky_bacteria; looky_eggsANDdetritus]);
looky_fleets                        = find(floor(preagg_GroupTypes) == GroupTypeDef_fleet);
looky_livingANDfleets               = [looky_ANYPrimaryProd; looky_ANYconsumer; looky_bacteria; looky_fleets]; % includes primary producers & bacteria
looky_NONnutrients                  = sort([looky_livingANDdetritus; looky_fleets]);                % addresses of all groups EXCEPT nutrients (needed to append nutrients)
% -------------------------------------------------------------------------


% step 2b: pre-aggregate counts of groups ---------------------------------
num_preagg_NO3                    	= length(looky_NO3);
num_preagg_plgcNH4                  = length(looky_plgcNH4);
num_preagg_bnthNH4                  = length(looky_bnthNH4);
num_preagg_NH4                     	= length(looky_NH4);
num_preagg_nutrients              	= length(looky_nutrients);
num_preagg_ANYdetritus            	= length(looky_ANYdetritus);                        % w/o eggs
num_preagg_eggs                    	= length(looky_eggs);
num_preagg_eggsANDdetritus         	= length(looky_eggsANDdetritus);
num_preagg_fleets                  	= length(looky_fleets);
num_preagg_livingANDdetritus    	= length(looky_livingANDdetritus);
num_EwEgrps                         = num_preagg_livingANDdetritus + num_preagg_fleets;	% PRE-aggregate count; EwE model; (count = living + detritus + fleet grps); (NO nutrients)
num_E2Egrps                         = num_EwEgrps + num_preagg_nutrients;               % PRE-aggregate count; E2E model; (count = living + detritus + fleet grps); (WITH nutrients)
num_preagg_terminalANYdetritus      = length(looky_terminalANYdetritus);                % FFF in future can specify between feces and senescence pools
% -------------------------------------------------------------------------


% step 2c: aggregation codes & counts -------------------------------
PreaggCodes_NO3                     = PreaggCodes(looky_NO3);                   % AggCodes of PRE-aggregated NO3 pools
PreaggCodes_plgcNH4                	= PreaggCodes(looky_plgcNH4);               % AggCodes of PRE-aggregated pelagic NH4 pools
PreaggCodes_bnthNH4                	= PreaggCodes(looky_bnthNH4);               % AggCodes of PRE-aggregated benthic NH4 pools
PreaggCodes_NH4                     = PreaggCodes(looky_NH4);                   % AggCodes of PRE-aggregated NH4 pools
PreaggCodes_nutrients               = PreaggCodes(looky_nutrients);             % AggCodes of PRE-aggregated nutrients pools
PreaggCodes_PrimaryProducers        = PreaggCodes(looky_ANYPrimaryProd);        % AggCodes of PRE-aggregated primary producers
PreaggCodes_consumers               = PreaggCodes(looky_ANYconsumer);           % AggCodes of PRE-aggregated consumers
PreaggCodes_micrograzers            = PreaggCodes(looky_micrograzers);          % AggCodes of PRE-aggregated micrograzers
PreaggCodes_bacteria                = PreaggCodes(looky_bacteria);              % AggCodes of PRE-aggregated bacteria
PreaggCodes_eggs                    = PreaggCodes(looky_eggs);                  % AggCodes of PRE-aggregated eggs
PreaggCodes_ANYdetritus             = PreaggCodes(looky_ANYdetritus);           % AggCodes of PRE-aggregated detritus pools
PreaggCodes_terminalPLGCdetritus	= PreaggCodes(looky_terminalPLGCdetritus);	% AggCodes of PRE-aggregated detritus pools
PreaggCodes_terminalBNTHdetritus    = PreaggCodes(looky_terminalBNTHdetritus);	% AggCodes of PRE-aggregated detritus pools
PreaggCodes_terminalANYdetritus 	= PreaggCodes(looky_terminalANYdetritus);	% AggCodes of PRE-aggregated detritus pools
PreaggCodes_eggsANDdetritus         = PreaggCodes(looky_eggsANDdetritus);       % AggCodes of PRE-aggregated detritus pools
PreaggCodes_livingANDdetritus   	= PreaggCodes(looky_livingANDdetritus);     % AggCodes of PRE-aggregated fleets
PreaggCodes_fleets                  = PreaggCodes(looky_fleets);                % AggCodes of PRE-aggregated fleets
PreaggCodes_livingANDfleets         = PreaggCodes(looky_livingANDfleets);       % AggCodes of PRE-aggregated groups excluding eggs & detritus, includes primary producers & bacteria
PreaggCodes_NONnutrients        	= PreaggCodes(looky_NONnutrients);          % AggCodes of PRE-aggregated fleets

AggCodes_NO3                        = unique(PreaggCodes_NO3);                  % AggCodes of aggregated NO3 pools
AggCodes_plgcNH4                    = unique(PreaggCodes_plgcNH4);             	% AggCodes of aggregated pelagic NH4 pools
AggCodes_bnthNH4                    = unique(PreaggCodes_bnthNH4);           	% AggCodes of aggregated benthic NH4 pools
AggCodes_NH4                        = unique(PreaggCodes_NH4);                  % AggCodes of aggregated NH4 pools
AggCodes_nutrients               	= unique(PreaggCodes_nutrients);            % AggCodes of aggregated nutrients pools
AggCodes_PrimaryProducers        	= unique(PreaggCodes_PrimaryProducers);     % AggCodes of aggregated primary producers
AggCodes_consumers                  = unique(PreaggCodes_consumers);            % AggCodes of aggregated consumers
AggCodes_micrograzers             	= unique(PreaggCodes_micrograzers);         % AggCodes of aggregated micrograzers
AggCodes_bacteria             	    = unique(PreaggCodes_bacteria);             % AggCodes of aggregated bacteria
AggCodes_eggs                       = unique(PreaggCodes_eggs);                 % AggCodes of aggregated eggs
AggCodes_ANYdetritus              	= unique(PreaggCodes_ANYdetritus);          % AggCodes of aggregated detritus pools
AggCodes_terminalPLGCdetritus    	= unique(PreaggCodes_terminalPLGCdetritus); % AggCodes of aggregated pealgic terminal detritus
AggCodes_terminalBNTHdetritus    	= unique(PreaggCodes_terminalBNTHdetritus);	% AggCodes of aggregated benthic terminal detritus
AggCodes_terminalANYdetritus     	= unique(PreaggCodes_terminalANYdetritus);  % AggCodes of aggregated ANY terminal detritus
AggCodes_eggsANDdetritus           	= unique(PreaggCodes_eggsANDdetritus);      % AggCodes of aggregated all eggs & detritus
AggCodes_livingANDdetritus       	= unique(PreaggCodes_livingANDdetritus);    % AggCodes of aggregated all primary producers + consumers + bacteria + eggs + detritus (no fleets, no nutrients)
AggCodes_fleets                     = unique(PreaggCodes_fleets);               % AggCodes of aggregated fleets
AggCodes_livingANDfleets            = unique(PreaggCodes_livingANDfleets);      % AggCodes of aggregated groups excluding eggs & detritus, includes primary producers & bacteria
AggCodes_NONnutrients             	= unique(PreaggCodes_NONnutrients);         % AggCodes of aggregated all primary producers + consumers + bacteria + eggs + detritus + fleets (no nutrients)

num_agg_NO3                         = length(AggCodes_NO3);
num_agg_plgcNH4                     = length(AggCodes_plgcNH4);
num_agg_bnthNH4                     = length(AggCodes_bnthNH4);
num_agg_NH4                         = length(AggCodes_NH4);
num_agg_nutrients                   = length(AggCodes_nutrients);
num_agg_PrimaryProducers            = length(AggCodes_PrimaryProducers);
num_agg_consumers                   = length(AggCodes_consumers);
num_agg_micrograzers                = length(AggCodes_micrograzers);
num_agg_bacteria                    = length(AggCodes_bacteria);
num_agg_eggs                        = length(AggCodes_eggs);
num_agg_ANYdetritus                 = length(AggCodes_ANYdetritus);
num_agg_terminalPLGCdetritus      	= length(AggCodes_terminalPLGCdetritus);
num_agg_terminalBNTHdetritus      	= length(AggCodes_terminalBNTHdetritus);
num_agg_terminalANYdetritus      	= length(AggCodes_terminalANYdetritus);
num_agg_eggsANDdetritus             = length(AggCodes_eggsANDdetritus);
num_agg_livingANDdetritus        	= length(AggCodes_livingANDdetritus);
num_agg_fleets                      = length(AggCodes_fleets);
num_agg_livingANDfleets             = length(AggCodes_livingANDfleets);
num_agg_NONnutrients            	= length(AggCodes_NONnutrients);

% error-checking: there must be 2 terminal detritus and at least 2 NH4 pools after aggregation
if num_agg_terminalANYdetritus ~= 2
    display('   -->WARNING: terminal detritus error:')
    display('      - ECOTRAN requires 2 terminal detritus groups in the post-aggregation model')
end

if num_agg_NH4 < 2
    display('   -->WARNING: NH4 pool error:')
    display('      - ECOTRAN requires at least 2 NH4 pools in the post-aggregation model')
end
% -------------------------------------------------------------------------


% step 2d: append nutrient addresses --------------------------------------
%          NOTE: in E2E ECOTRAN order with nutrients & fleets (dimensions = num_E2Egrps)
preagg_biomass(looky_NONnutrients)                      = preagg_biomass;                   % pre-aggregate E2E biomass; (t/km2); (vertical vector: num_E2Egrps X 1)
preagg_biomass(looky_nutrients)                         = zeros([num_preagg_nutrients 1]);
preagg_production(looky_NONnutrients)                   = preagg_production;                % pre-aggregate E2E production rate; (t/km2/y); (vertical vector: num_E2Egrps X 1)
preagg_production(looky_nutrients)                      = ones([num_preagg_nutrients 1]);   % NOTE: set to 1 for aggregation, nutrients are set to 0 after aggregation (step 10b)
preagg_pb(looky_NONnutrients)                           = preagg_pb;                        % pre-aggregate E2E biomass; (1/y); (vertical vector: num_E2Egrps X 1)
preagg_pb(looky_nutrients)                              = ones([num_preagg_nutrients 1]);
preagg_qb(looky_NONnutrients)                           = preagg_qb;                        % pre-aggregate E2E consumption rate; (1/y); (vertical vector: num_E2Egrps X 1)
preagg_qb(looky_nutrients)                              = ones([num_preagg_nutrients 1]);
preagg_pq(looky_NONnutrients)                           = preagg_pq;                        % pre-aggregate E2E production efficiency; (proportion); (vertical vector: num_E2Egrps X 1)
preagg_pq(looky_nutrients)                              = ones([num_preagg_nutrients 1]);
preagg_TL(looky_NONnutrients)                           = preagg_TL;                        % pre-aggregate E2E trophic level; (unitless); (vertical vector: num_E2Egrps X 1)
preagg_TL(looky_nutrients)                              = zeros([num_preagg_nutrients 1]);
preagg_ae(looky_NONnutrients)                           = preagg_ae;                        % pre-aggregate E2E assimilation efficiency; (proportion); (vertical vector: num_E2Egrps X 1)
preagg_ae(looky_nutrients)                              = zeros([num_preagg_nutrients 1]);

EwE_CONSUMPTION(looky_livingANDdetritus, :)          = EwE_CONSUMPTION;               % pre-aggregate EwE Consumption matrix (Q_pc); (t/km2/y); (2D matrix: num_E2Egrps X num_preagg_livingANDdetritus)
EwE_CONSUMPTION(looky_nutrients, :)                  = zeros([num_preagg_nutrients num_preagg_livingANDdetritus]);
EwE_CONSUMPTION(:, looky_livingANDdetritus)          = EwE_CONSUMPTION;
EwE_CONSUMPTION(:, looky_nutrients)                  = zeros([(num_preagg_livingANDdetritus + num_preagg_nutrients) num_preagg_nutrients]); % pre-aggregate EwE Consumption matrix (Q_pc); (t/km2/y); (2D matrix: num_E2Egrps X num_E2Egrps)

preagg_DIET(looky_livingANDdetritus, :)                 = preagg_DIET;                      % pre-aggregate EwE Diet matrix (D_pc); (proportion); (2D matrix: num_E2Egrps X num_preagg_livingANDdetritus)
preagg_DIET(looky_nutrients, :)                         = zeros([num_preagg_nutrients num_preagg_livingANDdetritus]);
preagg_DIET(:, looky_livingANDdetritus)                 = preagg_DIET;
preagg_DIET(:, looky_nutrients)                         = zeros([(num_preagg_livingANDdetritus + num_preagg_nutrients) num_preagg_nutrients]);   % pre-aggregate EwE Diet matrix (D_pc); (t/km2/y); (2D matrix: num_E2Egrps X num_E2Egrps)
preagg_diet_import(looky_livingANDdetritus)             = preagg_diet_import;               % pre-aggregate EwE import diet; (proportion); (horizontal vector: 1 X num_E2Egrps);
preagg_diet_import(looky_nutrients)                     = zeros([1 num_preagg_nutrients]);

preagg_LANDINGS(looky_livingANDdetritus, :)             = preagg_LANDINGS;                  % pre-aggregate Landings matrix; (t/km2/y); (2D matrix: num_E2Egrps X num_preagg_fleets)
preagg_LANDINGS(looky_nutrients, :)                     = zeros([num_preagg_nutrients num_preagg_fleets]);
preagg_DISCARDS(looky_livingANDdetritus, :)             = preagg_DISCARDS;                  % pre-aggregate Discard matrix; (t/km2/y); (2D matrix: num_E2Egrps X num_preagg_fleets)
preagg_DISCARDS(looky_nutrients, :)                     = zeros([num_preagg_nutrients num_preagg_fleets]);
preagg_CATCH(looky_livingANDdetritus, :)                = preagg_CATCH;                     % pre-aggregate Catch matrix; (t/km2/y); (2D matrix: num_E2Egrps X num_preagg_fleets)
preagg_CATCH(looky_nutrients, :)                        = zeros([num_preagg_nutrients num_preagg_fleets]);

% fate terms; NOTE fate_feces, fate_senescence, fate_metabolism already have nutrient addresses
preagg_fate_EwEdetritus(looky_NONnutrients, :)          = preagg_fate_EwEdetritus;          % pre-aggregate EwE detritus fates; (proportion); (2D matrix: num_E2Egrps X num_EwE_DetritusGrps)
preagg_fate_EwEdetritus(looky_nutrients, :)             = zeros([num_preagg_nutrients num_preagg_eggsANDdetritus]);

% pedigree terms; NOTE preagg_PedigreeCatch is calculated later
preagg_PedigreeParameters(looky_livingANDdetritus, :)  	= preagg_PedigreeParameters;            % pre-aggregate EwE PedigreeParameters; (CV); (2D matrix: (num_preagg_livingANDdetritus+num_preagg_nutrients) X 37)
preagg_PedigreeParameters(looky_nutrients, :)           = zeros([num_preagg_nutrients 37]);

preagg_PEDIGREE_LANDINGS(looky_livingANDdetritus, :)   	= preagg_PEDIGREE_LANDINGS;              % pre-aggregate PedigreeLandings matrix; (CV); (2D matrix: (num_preagg_livingANDdetritus+num_preagg_nutrients) X num_preagg_fleets)
preagg_PEDIGREE_LANDINGS(looky_nutrients, :)          	= zeros([num_preagg_nutrients num_preagg_fleets]);
preagg_PEDIGREE_DISCARDS(looky_livingANDdetritus, :)   	= preagg_PEDIGREE_DISCARDS;              % pre-aggregate PedigreeDiscards matrix; (CV); (2D matrix: (num_preagg_livingANDdetritus+num_preagg_nutrients) X num_preagg_fleets)
preagg_PEDIGREE_DISCARDS(looky_nutrients, :)          	= zeros([num_preagg_nutrients num_preagg_fleets]);

preagg_PEDIGREE_DIET(looky_livingANDdetritus, :)     	= preagg_PEDIGREE_DIET;                  % pre-aggregate EwE Diet pedigree (CV); (2D matrix: num_E2Egrps X num_preagg_livingANDdetritus)
preagg_PEDIGREE_DIET(looky_nutrients, :)              	= zeros([num_preagg_nutrients num_preagg_livingANDdetritus]);
preagg_PEDIGREE_DIET(:, looky_livingANDdetritus)       	= preagg_PEDIGREE_DIET;
preagg_PEDIGREE_DIET(:, looky_nutrients)              	= zeros([(num_preagg_livingANDdetritus + num_preagg_nutrients) num_preagg_nutrients]);	% pre-aggregate EwE Diet pedigree (CV); (2D matrix: num_E2Egrps X num_E2Egrps)
preagg_PEDIGREE_diet_import(looky_livingANDdetritus)  	= preagg_PEDIGREE_diet_import;           % pre-aggregate EwE import diet pedigree; (CV); (horizontal vector: 1 X num_E2Egrps);
preagg_PEDIGREE_diet_import(looky_nutrients)           	= zeros([1 num_preagg_nutrients]);
% -------------------------------------------------------------------------


% step 2e: append or insert fleet information -----------------------------
%          NOTE: fleet info for pb, qb, pq, & ae are calculated AFTER aggregation
%          NOTE: TL already contains fleet info, ee already all 0
%          NOTE: fate terms already have fleet information (both EwE & E2E)
%          NOTE: dimensions = num_E2Egrps
preagg_biomass(looky_fleets)                        = preagg_catch_TotalByFleet';                   % insert fleet "biomasses" as total_catch; (t/km2 "per year"); (vertical vector: num_E2Egrps X 1)
preagg_production(looky_fleets)                     = preagg_catch_TotalByFleet';                   % insert fleet "production" as total_catch for aggregation, fleet production will be set to agg_production = total_landings after aggregation in step 9b; (t/km2/y); (vertical vector: num_EwEgrps X 1)

preagg_LANDINGS(looky_fleets, :)                    = zeros([num_preagg_fleets num_preagg_fleets]);	% (t/km2/y); (2D matrix: num_E2Egrps X num_preagg_fleets)
preagg_DISCARDS(looky_fleets, :)                    = zeros([num_preagg_fleets num_preagg_fleets]);	% (t/km2/y); (2D matrix: num_E2Egrps X num_preagg_fleets)
preagg_CATCH(looky_fleets, :)                       = zeros([num_preagg_fleets num_preagg_fleets]);	% (t/km2/y); (2D matrix: num_E2Egrps X num_preagg_fleets)

EwE_CONSUMPTION(1:num_E2Egrps, looky_fleets)        = preagg_CATCH;                                 % (t/km2/y); (2D matrix: num_E2Egrps X num_E2Egrps); NOTE: CATCH columns appended just to keep matrix size consistent, CATCH values in EwE_CONSUMPTION are never used
EwE_CONSUMPTION_eggsANDdetritus                     = EwE_CONSUMPTION(:, looky_eggsANDdetritus);    % (t/km2/y); (2D matrix: num_E2Egrps X num_preagg_eggsANDdetritus)

preagg_DIET(1:num_E2Egrps, looky_fleets)            = preagg_CATCH;                                 % paste in fleet columns; (t/km2/y); (2D matrix: num_E2Egrps X num_E2Egrps); (NOTE: fleet columns won't sum to 1 until renormalization step)
preagg_diet_import(1, looky_fleets)                 = 0;                                            % (t/km2/y); (horizontal vector: 1 X num_E2Egrps)
preagg_DIET(1:num_E2Egrps, looky_eggsANDdetritus)	= EwE_CONSUMPTION_eggsANDdetritus;              % paste in eggs & detritus columns; (t/km2/y); (2D matrix: num_E2Egrps X num_E2Egrps); (NOTE: egg & detritus columns won't sum to 1 until renormalization step)
sum_preagg_DIET                                     = sum(preagg_DIET, 1);
preagg_DIET                                         = preagg_DIET ./ repmat(sum_preagg_DIET, [num_E2Egrps, 1]); % renormalize each column to 100% of domestic diet
preagg_DIET(isnan(preagg_DIET))                     = 0;                                            % correct div/0 errors

preagg_em(looky_fleets)                             = preagg_landings_TotalByFleet';                % fleet landings are an emigration term; (t/km2/y); (vertical vector: num_E2Egrps X 1); NOTE: absolute production rate units; NOTE transpose

preagg_PedigreeParameters(looky_fleets, :)          = zeros([num_preagg_fleets 37]);                % append dummy rows for fleets (so rows are arranged identically to physiological parameters)
preagg_PEDIGREE_LANDINGS(looky_fleets, :)         	= zeros([num_preagg_fleets num_preagg_fleets]);	% append dummy rows for fleets (so rows are arranged identically to physiological parameters)
preagg_PEDIGREE_DISCARDS(looky_fleets, :)         	= zeros([num_preagg_fleets num_preagg_fleets]);	% append dummy rows for fleets (so rows are arranged identically to physiological parameters)

preagg_PEDIGREE_DIET(1:num_E2Egrps, looky_fleets)	= zeros([num_E2Egrps num_preagg_fleets]);       % (CV); (2D matrix: num_E2Egrps X num_E2Egrps); (NOTE: fleet columns will be substituted with Catch CV below)
preagg_PEDIGREE_diet_import(1, looky_fleets)       	= 0;                                            % (CV); (horizontal vector: 1 X num_E2Egrps)
% -------------------------------------------------------------------------


% step 2f: recalculate consumption terms & CONSUMPTION matrix -------------
preagg_consumption_total                            = preagg_biomass' .* preagg_qb';                          	% pre-aggregate EwE total_consumption rate (q_c); (t/km2/y); (horizontal vector: 1 X num_E2Egrps); these values DO include import consumption
preagg_consumption_total(looky_fleets)              = sum(preagg_CATCH);                                      	% paste in fleet catches; (t/km2/y); (horizontal vector: 1 X num_E2Egrps)
preagg_consumption_total(looky_eggsANDdetritus)     = sum(EwE_CONSUMPTION_eggsANDdetritus);                   	% paste in flow to eggs & detritus; take values from EwE consumption matrix; (t/km2/y); (horizontal vector: 1 X num_E2Egrps);
preagg_consumption_import                           = preagg_consumption_total .* preagg_diet_import;          	% pre-aggregate import consumption rate (q_c); (t/km2/y); (horizontal vector: 1 X num_E2Egrps)
preagg_consumption_domestic                         = preagg_consumption_total - preagg_consumption_import;     % pre-aggregate import consumption rate (q_c); (t/km2/y); (horizontal vector: 1 X num_E2Egrps)
preagg_CONSUMPTION                                  = preagg_DIET .* repmat(preagg_consumption_domestic, [num_E2Egrps, 1]);	% pre-aggregate consumption matrix (q_c); (t/km2/y); (2D matrix: num_E2Egrps X num_E2Egrps)
% *************************************************************************





% *************************************************************************
% STEP 3: prep egg & detritus fates & import consumption as em-------------

% step 3a: emigration (em) terms ------------------------------------------
%          NOTE: absolute consumption rate units are here converted into absolute production rate units
%          NOTE: later converted to proportion of production & consumption budgets
preagg_em                                   = preagg_em - (preagg_consumption_import' .* preagg_pq);	% import diet is an import term (negative em); NOTE: conversion of import consumption to production rate terms; NOTE: no need to convert fleets because landings are already a production term; (t/km2/y); (vertical vector: num_E2Egrps X 1); NOTE transpose
% -------------------------------------------------------------------------


% step 3b: egg & detritus consumption rates into preagg_production --------
preagg_production(looky_eggsANDdetritus)	= preagg_consumption_total(looky_eggsANDdetritus)';	% insert detritus & egg "consumption" as "production"; (t/km2/y); (vertical vector: num_E2Egrps X 1)
% -------------------------------------------------------------------------


% step 3c: EwE egg fates --------------------------------------------------
preagg_GroupTypes_eggsANDdetritus	= preagg_GroupTypes(looky_eggsANDdetritus);
looky_eggsASdetritus              	= find(preagg_GroupTypes_eggsANDdetritus  == GroupTypeDef_eggs);	% egg addresses WITHIN preagg_fate_EwEdetritus 
preagg_fate_EwEeggs                 = preagg_fate_EwEdetritus(:, looky_eggsASdetritus);                 % select egg columns within EwE_DetritusFate; (2D matrix: num_E2Egrps X num_preagg_eggs)
preagg_fate_EwEdetritus(:, looky_eggsASdetritus)	= [];                                               % select detritus columns within EwE_DetritusFate; (2D matrix: num_E2Egrps X num_preagg_ANYdetritus)
% *************************************************************************





% *************************************************************************
% STEP 4: aggregate biomass, physiology, & fleets--------------------------

% step 4a: initialize aggregated model variables --------------------------
AggNumber(1:num_AggGrps, 1)    	= NaN;                                      % (vertical vector: num_AggGrps X 1)
AggLabel                       	= cell(num_AggGrps, 1);                     % (vertical vector: num_AggGrps X 1)
AggGroupType(1:num_AggGrps, 1)	= NaN;                                      % (vertical vector: num_AggGrps X 1)

agg_biomass                    	= zeros(num_AggGrps, 1);                    % (vertical vector: num_AggGrps X 1)
agg_production                	= zeros(num_AggGrps, 1);                    % (vertical vector: num_AggGrps X 1)
agg_pb                         	= zeros(num_AggGrps, 1);                    % (vertical vector: num_AggGrps X 1)
agg_qb                         	= zeros(num_AggGrps, 1);                    % (vertical vector: num_AggGrps X 1)
agg_pq                         	= zeros(num_AggGrps, 1);                    % (vertical vector: num_AggGrps X 1)
agg_TL                         	= zeros(num_AggGrps, 1);                    % (vertical vector: num_AggGrps X 1)
agg_ae                        	= zeros(num_AggGrps, 1);                    % (vertical vector: num_AggGrps X 1)
agg_ba                         	= zeros(num_AggGrps, 1);                    % (vertical vector: num_AggGrps X 1)
agg_em                         	= zeros(num_AggGrps, 1);                    % (vertical vector: num_AggGrps X 1)

agg_CONSUMPTION                	= zeros(num_AggGrps, num_AggGrps);          % (2D matrix: num_AggGrps X num_AggGrps)
agg_consumption_domestic       	= zeros(1, num_AggGrps);                    % (horizontal vector: 1 X num_AggGrps)
agg_consumption_import        	= zeros(1, num_AggGrps);                    % (horizontal vector: 1 X num_AggGrps)
agg_consumption_total         	= zeros(1, num_AggGrps);                    % (horizontal vector: 1 X num_AggGrps)

agg_DIET                     	= zeros(num_AggGrps, num_AggGrps);          % (2D matrix: num_AggGrps X num_AggGrps)
agg_diet_import               	= zeros(1, num_AggGrps);                    % (horizontal vector: 1 X num_AggGrps)
agg_DIET_method2              	= zeros(num_AggGrps, num_AggGrps);          % (2D matrix: num_AggGrps X num_AggGrps)
agg_diet_import_method2        	= zeros(1, num_AggGrps);                    % (horizontal vector: 1 X num_AggGrps)

horizontal_LANDINGS            	= zeros(num_AggGrps, num_preagg_fleets);    % for fleet aggregation step 1; (2D matrix: num_AggGrps X num_preagg_fleets)
horizontal_DISCARDS           	= zeros(num_AggGrps, num_preagg_fleets);	% for fleet aggregation step 1; (2D matrix: num_AggGrps X num_preagg_fleets)
horizontal_CATCH            	= zeros(num_AggGrps, num_preagg_fleets);	% for fleet aggregation step 1; (2D matrix: num_AggGrps X num_preagg_fleets)
agg_LANDINGS                  	= zeros(num_AggGrps, num_agg_fleets);       % (t/km2/y); (2D matrix: num_AggGrps X num_agg_fleets)
agg_DISCARDS                   	= zeros(num_AggGrps, num_agg_fleets);       % (t/km2/y); (2D matrix: num_AggGrps X num_agg_fleets)
agg_CATCH                     	= zeros(num_AggGrps, num_agg_fleets);       % (t/km2/y); (2D matrix: num_AggGrps X num_agg_fleets)

agg_EggProduction_forced       	= zeros(num_AggGrps, 1);                    % (vertical vector: num_AggGrps X 1)
agg_Metabolism_forced          	= zeros(num_AggGrps, 1);                    % (vertical vector: num_AggGrps X 1)

horizontal_fate_EwEdetritus    	= zeros(num_AggGrps, num_preagg_ANYdetritus);        	% (2D matrix: num_AggGrps X num_preagg_ANYdetritus)
horizontal_fate_EwEeggs       	= zeros(num_AggGrps, num_preagg_eggs);                  % (2D matrix: num_AggGrps X num_preagg_eggs)
agg_fate_EwEdetritus          	= zeros(num_AggGrps, num_agg_ANYdetritus);              % (2D matrix: num_AggGrps X num_agg_ANYdetritus)
agg_fate_EwEeggs              	= zeros(num_AggGrps, num_agg_eggs);                     % (2D matrix: num_AggGrps X num_agg_eggs)

agg_fate_feces                	= zeros(num_preagg_terminalANYdetritus, num_AggGrps);	% (2D matrix: num_preagg_terminalANYdetritus X num_AggGrps)
agg_fate_senescence           	= zeros(num_preagg_terminalANYdetritus, num_AggGrps);   % (2D matrix: num_preagg_terminalANYdetritus X num_AggGrps)
agg_fate_metabolism            	= zeros(num_preagg_NH4, num_AggGrps);                   % (2D matrix: num_preagg_NH4 X num_AggGrps)

agg_ProductionLossScaler      	= zeros(num_AggGrps, 1);                    % (vertical vector: num_AggGrps X 1)
agg_RetentionScaler            	= zeros(num_AggGrps, 1);                    % (vertical vector: num_AggGrps X 1)
agg_FunctionalResponseParams   	= zeros(num_AggGrps, 4);                    % (vertical vector: num_AggGrps X 4)
agg_FunctionalResponse_matrix	= zeros(num_AggGrps, num_AggGrps);          % (2D matrix: num_AggGrps X num_AggGrps)
% -------------------------------------------------------------------------


for TrophicGrp_loop = 1:num_AggGrps
    
    % step 4b: pick current groups for aggregation ------------------------
	current_AggCode           	= AggCodes(TrophicGrp_loop);               	% current aggregation group; AggCodes is the "unique" code list
    looky_currentAggCode      	= find(PreaggCodes == current_AggCode);     % row numbers of all groups included in current aggregated group
    num_currentAgg          	= length(looky_currentAggCode);         	% number of groups in the current aggregation
    % ---------------------------------------------------------------------

    
    % step 4c: aggregate group code, group label, & GroupType codes -------
    AggNumber(TrophicGrp_loop)	= current_AggCode;                        	% aggregate group number; (vertical vector: num_AggGrps X 1)
    AggLabel(TrophicGrp_loop)  	= preagg_AggName(looky_currentAggCode(1));	% aggregate group labels; (vertical vector: num_AggGrps X 1)
    current_GroupType         	= preagg_GroupTypes(looky_currentAggCode);
    % ---------------------------------------------------------------------

    
    % step 4d: error check ------------------------------------------------
    %          prevent cross-aggregation of eggs, detritus, plankton, nutrients, & consumers
    unique_GroupType        = floor(unique(current_GroupType));
    if intersect(unique_GroupType, GroupSet_consumers) & intersect(unique_GroupType, GroupTypeDef_ANYPrimaryProd)
        error('Error--> You are trying to aggregate a consumer and a primary producer')
    elseif intersect(unique_GroupType, GroupSet_consumers) & intersect(unique_GroupType, GroupTypeDef_ANYNitroNutr)
        error('Error--> You are trying to aggregate a consumer and a nutrient')
    elseif intersect(unique_GroupType, GroupSet_consumers) & intersect(unique_GroupType, GroupTypeDef_ANYDetritus)
        error('Error--> You are trying to aggregate a consumer and detritus')
    elseif intersect(unique_GroupType, GroupSet_consumers) & intersect(unique_GroupType, GroupTypeDef_eggs)
        error('Error--> You are trying to aggregate a consumer and an egg')
    elseif intersect(unique_GroupType, GroupTypeDef_ANYPrimaryProd) & intersect(unique_GroupType, GroupTypeDef_ANYNitroNutr)
        error('Error--> You are trying to aggregate a primary producer and a nutrient')
    elseif intersect(unique_GroupType, GroupTypeDef_ANYPrimaryProd) & intersect(unique_GroupType, GroupTypeDef_ANYDetritus)
        error('Error--> You are trying to aggregate a primary producer and detritus')
    elseif intersect(unique_GroupType, GroupTypeDef_ANYPrimaryProd) & intersect(unique_GroupType, GroupTypeDef_eggs)
        error('Error--> You are trying to aggregate a primary producer and an egg')
    elseif intersect(unique_GroupType, GroupTypeDef_ANYNitroNutr) & intersect(unique_GroupType, GroupSet_detritus)
        error('Error--> You are trying to aggregate a nutrient and detritus')
    elseif intersect(unique_GroupType, GroupTypeDef_ANYDetritus) & intersect(unique_GroupType, GroupTypeDef_eggs)
        error('Error--> You are trying to aggregate detritus and an egg')
    end
    % ---------------------------------------------------------------------
    
    
    % step 4e: aggregate group types --------------------------------------
    if sum(ismember(current_GroupType, 4.1)) > 0
        AggGroupType(TrophicGrp_loop, 1) = 4.1; % favor terminal pelagic detritus GroupType
    elseif sum(ismember(current_GroupType, 4.3)) > 0
        AggGroupType(TrophicGrp_loop, 1) = 4.3; % favor terminal benthic detritus GroupType
    else
        AggGroupType(TrophicGrp_loop, 1) = mode(current_GroupType); % use the most common group type in the aggregation group (in case groups of different types get aggregated together)
    end
    % ---------------------------------------------------------------------
    
    
    % step 4f: aggregation by simple sums ---------------------------------
    current_biomass                                	= preagg_biomass(looky_currentAggCode);                 % (vertical vector: num_currentAgg X 1)
    current_production                            	= preagg_production(looky_currentAggCode);              % (vertical vector: num_currentAgg X 1)
    current_ba                                     	= preagg_ba(looky_currentAggCode);                      % (t/km2/y); (vertical vector: num_currentAgg X 1)
    current_em                                    	= preagg_em(looky_currentAggCode);                      % (t/km2/y); (vertical vector: num_currentAgg X 1)
    current_consumption_domestic                  	= preagg_consumption_domestic(looky_currentAggCode);    % domestic_consumption rate (q_c); (t/km2/y); (horizontal vector: 1 X num_currentAgg);
    current_consumption_import                    	= preagg_consumption_import(looky_currentAggCode);      % import_consumption rate (q_c); (t/km2/y); (horizontal vector: 1 X num_currentAgg); 
    current_consumption_total                     	= preagg_consumption_total(looky_currentAggCode);       % total_consumption rate (q_c); (t/km2/y); (horizontal vector: 1 X num_currentAgg); NOTE: these values the same as q = biomass * q/b, but also include non-zero entries for eggs, detritus, & fleets
    current_EggProduction_forced                  	= preagg_EggProduction_forced(looky_currentAggCode);	% (t/km2/y); (vertical vector: num_currentAgg X 1)
    current_Metabolism_forced                       = preagg_Metabolism_forced(looky_currentAggCode);       % (t/km2/y); (vertical vector: num_currentAgg X 1)
    
    agg_biomass(TrophicGrp_loop, 1)                 = sum(current_biomass);                                 % aggregate biomass; (t/km2);    	(vertical vector: num_AggGrps X 1)
	agg_production(TrophicGrp_loop, 1)              = sum(current_production);                              % aggregate production; (t/km2/y);	(vertical vector: num_AggGrps X 1)
    agg_ba(TrophicGrp_loop, 1)                      = sum(current_ba);                                      % aggregate ba; (t/km2/y);        	(vertical vector: num_AggGrps X 1)
    agg_em(TrophicGrp_loop, 1)                      = sum(current_em);                                      % aggregate em; (t/km2/y);        	(vertical vector: num_AggGrps X 1)
    agg_consumption_domestic(1, TrophicGrp_loop)	= sum(current_consumption_domestic);                    % aggregate import_consumption rate (q_c); (t/km2/y); (horizontal vector: 1 X num_AggGrps); NOTE: does NOT include import consumption
    agg_consumption_import(1, TrophicGrp_loop)      = sum(current_consumption_import);                      % aggregate import_consumption rate (q_c); (t/km2/y); (horizontal vector: 1 X num_AggGrps)
    agg_consumption_total(1, TrophicGrp_loop)       = sum(current_consumption_total);                       % aggregate total_consumption rate (q_c); (t/km2/y); (horizontal vector: 1 X num_AggGrps); NOTE: DOES include import consumption
    agg_EggProduction_forced(TrophicGrp_loop)       = sum(current_EggProduction_forced);                    % aggregate EggProduction_forced; (t/km2/y); (vertical vector: num_currentAgg X 1)
    agg_Metabolism_forced(TrophicGrp_loop)          = sum(current_Metabolism_forced);                       % aggregate Metabolism_forced; (t/km2/y); (vertical vector: num_currentAgg X 1)
    % ---------------------------------------------------------------------

    
    % step 4g: aggregate fleets -------------------------------------------
    %       NOTE: sum along caught groups
    current_LANDINGS                                            = preagg_LANDINGS(looky_currentAggCode, :);	% (2D matrix: num_AggGrps X num_preagg_fleets)
    current_DISCARDS                                            = preagg_DISCARDS(looky_currentAggCode, :); % (2D matrix: num_AggGrps X num_preagg_fleets)
    current_CATCH                                               = preagg_CATCH(looky_currentAggCode, :);    % (2D matrix: num_AggGrps X num_preagg_fleets)
    
    horizontal_LANDINGS(TrophicGrp_loop, 1:num_preagg_fleets)	= sum(current_LANDINGS, 1);                 % (t/km2/y); (2D matrix: num_AggGrps X num_preagg_fleets)
    horizontal_DISCARDS(TrophicGrp_loop, 1:num_preagg_fleets)	= sum(current_DISCARDS, 1);                 % (t/km2/y); (2D matrix: num_AggGrps X num_preagg_fleets)
    horizontal_CATCH(TrophicGrp_loop, 1:num_preagg_fleets)      = sum(current_CATCH, 1);                    % (t/km2/y); (2D matrix: num_AggGrps X num_preagg_fleets)
    % ---------------------------------------------------------------------

    
    % step 4h: aggregate physiology terms ---------------------------------
    %       NOTE: (production-weighted means)
    current_pb           	= preagg_pb(looky_currentAggCode);     % (1/y);         (vertical vector: num_currentAgg X 1)
    current_qb            	= preagg_qb(looky_currentAggCode);     % (1/y);         (vertical vector: num_currentAgg X 1)
    current_pq          	= preagg_pq(looky_currentAggCode);     % (unitless);    (vertical vector: num_currentAgg X 1)
    current_TL           	= preagg_TL(looky_currentAggCode);     % (unitless);    (vertical vector: num_currentAgg X 1)
    current_ae          	= preagg_ae(looky_currentAggCode);     % (proportion);	(vertical vector: num_currentAgg X 1)

	agg_pb(TrophicGrp_loop)	= sum(current_pb .* current_production) / agg_production(TrophicGrp_loop); % weighted mean pb; (vertical vector: num_AggGrps X 1)
	agg_qb(TrophicGrp_loop)	= sum(current_qb .* current_production) / agg_production(TrophicGrp_loop); % weighted mean qb; (vertical vector: num_AggGrps X 1)
	agg_pq(TrophicGrp_loop)	= sum(current_pq .* current_production) / agg_production(TrophicGrp_loop); % weighted mean pq; (vertical vector: num_AggGrps X 1)
	agg_TL(TrophicGrp_loop)	= sum(current_TL .* current_production) / agg_production(TrophicGrp_loop); % weighted mean TL; (vertical vector: num_AggGrps X 1)
	agg_ae(TrophicGrp_loop)	= sum(current_ae .* current_production) / agg_production(TrophicGrp_loop); % weighted mean ae; (vertical vector: num_AggGrps X 1)
    % ---------------------------------------------------------------------

    
    % step 4i: aggregate fate_EwEeggs & fate_EwEdetritus ------------------
    %       NOTE: (production-weighted means)
    current_fate_EwEdetritus                      	= preagg_fate_EwEdetritus(looky_currentAggCode, :);	% (2D matrix: num_currentAgg X num_preagg_ANYdetritus)
    current_fate_EwEeggs                            = preagg_fate_EwEeggs(looky_currentAggCode, :);     % (2D matrix: num_currentAgg X num_preagg_eggs)
    
    horizontal_fate_EwEdetritus(TrophicGrp_loop, :)	= sum((current_fate_EwEdetritus .* repmat(current_production, [1 num_preagg_ANYdetritus])), 1) ./ repmat(agg_production(TrophicGrp_loop), [1 num_preagg_ANYdetritus]);	% (2D matrix: num_AggGrps X num_preagg_ANYdetritus)
    horizontal_fate_EwEeggs(TrophicGrp_loop, :)     = sum((current_fate_EwEeggs     .* repmat(current_production, [1 num_preagg_eggs])), 1)        ./ repmat(agg_production(TrophicGrp_loop), [1 num_preagg_eggs]);         % (2D matrix: num_AggGrps X num_preagg_eggs)
    % ---------------------------------------------------------------------

    
	% step 4j: aggregate fate_feces, fate_senescence, & fate_metabolism  --
    %       NOTE: (production-weighted means)
    current_fate_feces                    	= preagg_fate_feces(:, looky_currentAggCode);     	% feces fate; (units = 0 to 1, fractions of all detritus); (2D matrix: num_agg_terminalANYdetritus(=2) X num_currentAgg)
    current_fate_senescence               	= preagg_fate_senescence(:, looky_currentAggCode);	% senescence fate; (units = 0 to 1, fractions of all detritus);	(2D matrix: num_agg_terminalANYdetritus(=2) X num_currentAgg)
    current_fate_metabolism              	= preagg_fate_metabolism(:, looky_currentAggCode);	% metabolic excretion fate; (units = 0 to 1, fractions of all metabolic excretion);	(2D matrix: num_agg_NH4(=2) X num_currentAgg)
    
    agg_fate_feces(:, TrophicGrp_loop)     	= sum((current_fate_feces      .* repmat(current_production', [num_agg_terminalANYdetritus 1])), 2) ./ repmat(agg_production(TrophicGrp_loop), [num_agg_terminalANYdetritus 1]); % (2D matrix: num_agg_terminalANYdetritus(=2) X num_AggGrps)
    agg_fate_senescence(:, TrophicGrp_loop)	= sum((current_fate_senescence .* repmat(current_production', [num_agg_terminalANYdetritus 1])), 2) ./ repmat(agg_production(TrophicGrp_loop), [num_agg_terminalANYdetritus 1]); % (2D matrix: num_agg_terminalANYdetritus(=2) X num_AggGrps)
    agg_fate_metabolism(:, TrophicGrp_loop)	= sum((current_fate_metabolism .* repmat(current_production', [num_agg_NH4 1])), 2)                 ./ repmat(agg_production(TrophicGrp_loop), [num_agg_NH4 1]);                 % (2D matrix: num_agg_NH4(=2) X num_AggGrps)
    % ---------------------------------------------------------------------

    
    % step 4k: aggregate ProductionLossScaler, RetentionScaler, FunctionalResponseParams, & FunctionalResponse_matrix
    %       NOTE: (production-weighted means)
    current_ProductionLossScaler                     	= preagg_ProductionLossScaler(looky_currentAggCode);            % production loss term for steady state analyses of physical export; (units = ???);	(vertical vector: num_currentAgg X 1)
    current_RetentionScaler                          	= preagg_RetentionScaler(looky_currentAggCode);             	% plankton retention scaler; (units = 0 for purely planktonic, 1 for purely sessile or nektonic); (vertical vector: num_currentAgg X 1)
    current_FunctionalResponseParams                 	= preagg_FunctionalResponseParams(looky_currentAggCode, 1:4);	% functional response parameters, 4 parameters; (2D matrix: num_currentAgg X 4)
    current_FunctionalResponse_matrix                	= preagg_FunctionalResponse_matrix(looky_currentAggCode, :); 	% functional response parameters, matrix of grp-to-grp relationships; (2D matrix: num_currentAgg X num_currentAgg);
    
    agg_ProductionLossScaler(TrophicGrp_loop)          	= sum(current_ProductionLossScaler       .* current_production)                               / agg_production(TrophicGrp_loop);                           	% weighted mean ProductionLossScaler (INCLUDING eggs & detritus lacking physiology)
    agg_RetentionScaler(TrophicGrp_loop)             	= sum(current_RetentionScaler            .* current_production)                               / agg_production(TrophicGrp_loop);                        	% weighted mean RetentionScaler (INCLUDING eggs & detritus lacking physiology)
    agg_FunctionalResponseParams(TrophicGrp_loop, :)  	= sum((current_FunctionalResponseParams  .* repmat(current_production, [1 4])), 1)           ./ repmat(agg_production(TrophicGrp_loop), [1 4]);            	% (2D matrix: num_currentAgg X 4 (preagg))
%     agg_FunctionalResponse_matrix(TrophicGrp_loop, :)	= sum((current_FunctionalResponse_matrix .* repmat(current_production, [1 num_E2Egrps])), 1) ./ repmat(agg_production(TrophicGrp_loop), [1 num_E2Egrps]);	% (2D matrix: num_currentAgg X num_E2Egrps (preagg)) (FFF this matrix not yet used but will need a lateral aggregation); QQQ deactivated 2/5/2021
    % ---------------------------------------------------------------------

end % end TrophicGrp_loop


% step 4l: renormalize fate terms to sum to 1 -----------------------------
agg_fate_feces                                  = agg_fate_feces ./ repmat(sum(agg_fate_feces, 1), [num_preagg_terminalANYdetritus 1]);             % (2D matrix: num_agg_terminalANYdetritus(=2) X num_AggGrps)
agg_fate_feces(isnan(agg_fate_feces))           = 0;            % fix div/0 errors
agg_fate_senescence                             = agg_fate_senescence ./ repmat(sum(agg_fate_senescence, 1), [num_preagg_terminalANYdetritus 1]);	% (2D matrix: num_agg_terminalANYdetritus(=2) X num_AggGrps)
agg_fate_senescence(isnan(agg_fate_senescence)) = 0;            % fix div/0 errors
agg_fate_metabolism                             = agg_fate_metabolism ./ repmat(sum(agg_fate_metabolism, 1), [num_preagg_NH4 1]);                   % (2D matrix: num_agg_NH4(=2) X num_AggGrps)
agg_fate_metabolism(isnan(agg_fate_metabolism)) = 0;            % fix div/0 errors
% -------------------------------------------------------------------------


% step 4m: aggregate EwE egg & detritus fates horizontally ----------------
%          NOTE: no need to aggregate E2E NH4 & detritus fates (these are already defined in the ECOTRAN parameters)
% horizontal EwE detritus fate aggregation
for detritus_loop = 1:num_agg_ANYdetritus
    
    % pick current groups for aggregation -----
	current_AggCode                     = AggCodes_ANYdetritus(detritus_loop);              % current aggregation detritus group
    looky_currentAggCode                = find(PreaggCodes_ANYdetritus == current_AggCode);	% row numbers of all groups included in current aggregated detritus group
    num_currentAgg                      = length(looky_currentAggCode);                     % number of groups in the current aggregation

    % aggregation by simple sums -----
    agg_fate_EwEdetritus(:, detritus_loop)  = sum(horizontal_fate_EwEdetritus(:, looky_currentAggCode), 2);	% (t/km2/y); (2D matrix: num_AggGrps X num_agg_ANYdetritus)
    
end % detritus_loop

% horizontal EwE egg fate aggregation
for egg_loop = 1:num_agg_eggs
    
    % pick current groups for aggregation -----
	current_AggCode              	= AggCodes_eggs(egg_loop);                      % current aggregation egg group
    looky_currentAggCode          	= find(PreaggCodes_eggs == current_AggCode);	% row numbers of all groups included in current aggregated egg group
    num_currentAgg                  = length(looky_currentAggCode);                 % number of groups in the current aggregation

    % aggregation by simple sums -----
    agg_fate_EwEeggs(:, egg_loop)	= sum(horizontal_fate_EwEeggs(:, looky_currentAggCode), 2);	% (t/km2/y); (2D matrix: num_AggGrps X num_agg_eggs)
    
end % egg_loop

% renormalize agg_fate_EwEdetritus & agg_fate_EwEeggs
temp_sum                                            = [agg_fate_EwEeggs agg_fate_EwEdetritus];
temp_sum                                            = sum(temp_sum, 2);
agg_fate_EwEdetritus                                = agg_fate_EwEdetritus ./ repmat(temp_sum, [1 num_agg_ANYdetritus]);	% (2D matrix: num_AggGrps X num_agg_ANYdetritus)
agg_fate_EwEdetritus(isnan(agg_fate_EwEdetritus))   = 0;            % fix div/0 errors
agg_fate_EwEeggs                                    = agg_fate_EwEeggs ./ repmat(temp_sum, [1 num_agg_eggs]);             	% (2D matrix: num_AggGrps X num_agg_eggs)
agg_fate_EwEeggs(isnan(agg_fate_EwEeggs))           = 0;            % fix div/0 errors
% *************************************************************************





% *************************************************************************
% STEP 5: aggregate fleets horizontally------------------------------------
for fleet_loop = 1:num_agg_fleets
    
    current_AggCode_fleets    	= AggCodes_fleets(fleet_loop);
    looky_current_fleet     	= find(PreaggCodes_fleets == current_AggCode_fleets); 	% clm numbers of all fleets included in current aggregated fleet

    agg_LANDINGS(:, fleet_loop)	= sum(horizontal_LANDINGS(:, looky_current_fleet), 2);	% (t/km2/y); (2D matrix: num_AggGrps X num_agg_fleets)
    agg_DISCARDS(:, fleet_loop)	= sum(horizontal_DISCARDS(:, looky_current_fleet), 2);	% (t/km2/y); (2D matrix: num_AggGrps X num_agg_fleets)
    agg_CATCH(:, fleet_loop)  	= sum(horizontal_CATCH(:, looky_current_fleet), 2);   	% (t/km2/y); (2D matrix: num_AggGrps X num_agg_fleets)

end % fleet_loop
agg_landings_TotalByFleet   = sum(agg_LANDINGS, 1);     % (t/km2/y); (horizontal vector: 1 X num_agg_fleets)
agg_discards_TotalByFleet   = sum(agg_DISCARDS, 1);     % (t/km2/y); (horizontal vector: 1 X num_agg_fleets)
agg_catch_TotalByFleet      = sum(agg_CATCH, 1);        % (t/km2/y); (horizontal vector: 1 X num_agg_fleets)
agg_landings_TotalByGroup	= sum(agg_LANDINGS, 2);     % (t/km2/y); (vertical vector: num_agg_fleets X 1)
agg_discards_TotalByGroup   = sum(agg_DISCARDS, 2);     % (t/km2/y); (vertical vector: num_agg_fleets X 1)
agg_catch_TotalByGroup      = sum(agg_CATCH, 2);        % (t/km2/y); (vertical vector: num_agg_fleets X 1)
% *************************************************************************





% *************************************************************************
% STEP 6: aggregate CONSUMPTION matrix-------------------------------------
%         NOTE: vertical & horizontal sums
for producer_loop = 1:num_AggGrps    
	current_ProducerCode            = AggCodes(producer_loop);                              % current aggregation producer
    looky_currentProducer           = find(PreaggCodes == current_ProducerCode);            % row numbers of all producers included in current aggregated group
    num_currentProducers            = length(looky_currentProducer);                        % number of groups in the current aggregation
    current_horizontal_consumption  = sum(preagg_CONSUMPTION(looky_currentProducer, :), 1);	% sum of prey consumed; (horizontal vector: 1 X num_currentProducers)
    
    for consumer_loop  = 1:num_AggGrps
        current_ConsumerCode            = AggCodes(consumer_loop);                          % current aggregation consumer
        looky_currentConsumer           = find(PreaggCodes == current_ConsumerCode);        % column numbers of all consumers included in current aggregated group
        num_currentConsumers            = length(looky_currentConsumer);                    % number of groups in the current aggregation
        agg_CONSUMPTION(producer_loop, consumer_loop) = sum(current_horizontal_consumption(looky_currentConsumer), 2); % sum of consumer group(s) consumption
    end % consumer_loop
    
end % producer_loop
% *************************************************************************





% *************************************************************************
% STEP 7: aggregate diet---------------------------------------------------

% step 7a: aggregate DIET matrix ------------------------------------------
%          NOTE: This step does not use earlier diet variables, consumption
%                matrix is all that is needed
%          NOTE: divide consumption matrix by domestic_consumption
%          NOTE: this diet matrix does not include import diet, 
%                it refers to diet while WITHIN model domain.
%                The total consumption term in divisor does not include
%                import consumption and so does NOT equal b * q/b
%          NOTE: detritus & fleet "diets" are now included
agg_DIET                    = agg_CONSUMPTION ./ repmat(agg_consumption_domestic, [num_AggGrps 1]);
agg_DIET(isnan(agg_DIET))   = 0;        % correct for div/0 errors (nutrients & primary producers)
% -------------------------------------------------------------------------


% step 7b: aggregate import diet vector -----------------------------------
%          NOTE: this is simply the ratio of import_consumption / total_consumption
%          NOTE: total_consumption = import_consumption + domestic_consumption
agg_diet_import                         = agg_consumption_import ./ agg_consumption_total;	% aggregate import_consumption (D_ic); (percent); (horizontal vector: 1 X num_AggGrps)
agg_diet_import(isnan(agg_diet_import)) = 0;        % correct div/0 errors caused by total_consumption = 0 for nutrients & primary producers
% *************************************************************************




% QQQ deactivated 2/5/2021
% % % % *************************************************************************
% % % % STEP 8: aggregate diet METHOD 2------------------------------------------
% % % %         NOTE: this is an alternate method for scaling diet
% % % 
% % % % step 8a: aggregate Diet matrix ------------------------------------------
% % % %         NOTE: this diet calculation has all zero entrees for eggs & detritus
% % % for producer_loop = 1:num_AggGrps
% % % 	current_ProducerCode            = AggCodes(producer_loop);                              % current aggregation producer
% % %     looky_currentProducer           = find(PreaggCodes == current_ProducerCode);            % row numbers of all producers included in current aggregated group
% % %     num_currentProducers            = length(looky_currentProducer);                        % number of groups in the current aggregation
% % %     current_horizontal_diet         = sum(preagg_DIET(looky_currentProducer, :), 1);        % sum of diet fractions; (horizontal vector: 1 X num_currentProducers)
% % %     current_horizontal_diet         = current_horizontal_diet .* preagg_consumption_total;	% scale by consumption; (horizontal vector: 1 X num_currentProducers)
% % %     
% % %     for consumer_loop  = 1:num_AggGrps
% % %         current_ConsumerCode            = AggCodes(consumer_loop);                          % current aggregation consumer
% % %         looky_currentConsumer           = find(PreaggCodes == current_ConsumerCode);        % column numbers of all consumers included in current aggregated group
% % %         num_currentConsumers            = length(looky_currentConsumer);                    % number of groups in the current aggregation
% % %         agg_DIET_method2(producer_loop, consumer_loop)   = sum(current_horizontal_diet(looky_currentConsumer)) / agg_consumption_total(looky_currentConsumer); % normalize by summed agg-model consumer consumption
% % %     end % consumer_loop
% % %     
% % % end % producer_loop
% % % agg_DIET_method2                            = agg_DIET_method2 ./ repmat(sum(agg_DIET_method2, 1), [num_E2Egrps 1]); % normalize by 1, this puts fleet catch into diet form (= percent of total catch)
% % % agg_DIET_method2(isnan(agg_DIET_method2))	= 0;        % correct div/0 errors caused by total_consumption = 0 for nutrients & primary producers
% % % % -------------------------------------------------------------------------


% % % % step 8b: aggregate import diet ------------------------------------------
% % % temp_diet_import = preagg_diet_import .* preagg_consumption_total;	% scale by consumption; (horizontal vector: 1 X num_currentProducers)
% % % for consumer_loop = 1:num_AggGrps
% % %     
% % %         current_ConsumerCode            = AggCodes(consumer_loop);                   	% current aggregation consumer
% % %         looky_currentConsumer           = find(PreaggCodes == current_ConsumerCode);	% column numbers of all consumers included in current aggregated group
% % %         num_currentConsumers            = length(looky_currentConsumer);              	% number of groups in the current aggregation
% % %         agg_diet_import_method2(producer_loop, consumer_loop) = sum(temp_diet_import(looky_currentConsumer)) / agg_consumption_total(looky_currentConsumer); % normalize by summed agg-model consumer consumption
% % %                 
% % % end % consumer_loop
% % % agg_diet_import_method2(isnan(agg_diet_import_method2))	= 0;        % correct div/0 errors caused by total_consumption = 0 for nutrients & primary producers
% *************************************************************************





% *************************************************************************
% STEP 9: addresses in AGGREGATED model -----------------------------------
looky_agg_NO3                       	= find(AggGroupType == GroupTypeDef_NO3);
looky_agg_plgcNH4                       = find(AggGroupType == GroupTypeDef_plgcNH4);
looky_agg_bnthNH4                       = find(AggGroupType == GroupTypeDef_bnthNH4);
looky_agg_NH4                           = find(AggGroupType == GroupTypeDef_plgcNH4 | AggGroupType == GroupTypeDef_bnthNH4);
looky_agg_nutrients                     = find(floor(AggGroupType) == GroupTypeDef_ANYNitroNutr); % row addresses of E2E nutrients
looky_agg_ANYPrimaryProd                = find(floor(AggGroupType) == GroupTypeDef_ANYPrimaryProd);
looky_agg_ANYconsumer                   = find(floor(AggGroupType) == GroupTypeDef_ANYConsumer);
looky_agg_micrograzers                  = find(AggGroupType == GroupTypeDef_micrograzers);
looky_agg_bacteria                      = find(floor(AggGroupType) == GroupTypeDef_bacteria);
looky_agg_eggs                          = find(AggGroupType == GroupTypeDef_eggs);
looky_agg_ANYdetritus                   = find(floor(AggGroupType) == GroupTypeDef_ANYDetritus);
looky_agg_terminalPLGCdetritus          = find(AggGroupType == GroupTypeDef_terminalPlgcDetr);
looky_agg_terminalBNTHdetritus          = find(AggGroupType == GroupTypeDef_terminalBnthDetr);
looky_agg_terminalANYdetritus           = find(AggGroupType == GroupTypeDef_terminalPlgcDetr | AggGroupType == GroupTypeDef_terminalBnthDetr);
looky_agg_eggsANDdetritus               = sort([looky_agg_ANYdetritus; looky_agg_eggs]);
looky_agg_livingANDdetritus             = sort([looky_agg_ANYPrimaryProd; looky_agg_ANYconsumer; looky_agg_bacteria; looky_agg_eggsANDdetritus]);
looky_agg_fleets                        = find(floor(AggGroupType) == GroupTypeDef_fleet);
looky_agg_livingANDfleets               = [looky_agg_ANYPrimaryProd; looky_agg_ANYconsumer; looky_agg_bacteria; looky_agg_fleets]; % includes primary producers & bacteria
looky_agg_NONnutrients                  = sort([looky_agg_livingANDdetritus; looky_agg_fleets]); % addresses of all groups EXCEPT nutrients (needed to append nutrients)
% *************************************************************************





% *************************************************************************
% STEP 10: group-specific adjustments--------------------------------------

% step 10a: adjust parameters for individual groups -----------------------
%          change parameter values for phytoplankton, detritus, & eggs ----
%          NOTE: use looky_ addresses for aggregated model
agg_ae(looky_agg_nutrients)         = 1; % make ae for nutriente = 1 (no feces)
agg_ae(looky_agg_ANYPrimaryProd)   	= 1; % make ae for phytoplankton = 1 (no feces); NOTE: redundant, here just in case source model has not defined as such
agg_ae(looky_agg_ANYdetritus)      	= 1; % make ae for PelagicDetritus = 1 (no feces); NOTE: redundant, here just in case source model has not defined as such
agg_ae(looky_agg_eggs)            	= 1; % make ae for eggs = 1 (no feces); NOTE: redundant, here just in case source model has not defined as such
agg_pq(looky_agg_ANYPrimaryProd)  	= 1; % make pq for phytoplankton = 1 (all nutrient uptake goes to production)
agg_pq(looky_agg_ANYdetritus)      	= 1; % make pq for PelagicDetritus = 1 (all "consumed" goes to production)
agg_pq(looky_agg_eggs)           	= 1; % make pq for eggs = 1 (all consumed goes to production)
agg_TL(looky_agg_eggs)              = 1; % make egg TL = 1; NOTE: redundant, here just in case source model has not defined as such
% -------------------------------------------------------------------------


% step 10b: assign production values to fleets & nutrients ----------------
agg_production(looky_agg_fleets)   	= agg_landings_TotalByFleet;    % NOTE: here we use landings, while above when we aggregated groups we wanted to use catch (= landings + discards)
agg_production(looky_agg_nutrients)	= 0;                            % set nutrient production to 0, nutrients will be an external driver
% -------------------------------------------------------------------------


% step 10c: recalculate pb from definition --------------------------------
%           to correct for production rates of fleets and for groups with eggs
%           NOTE: VERY IMPORTANT LINE
%           NOTE: this also gives pb values for eggs & detritus groups
%           NOTE: egg corrections to production & parameters made in ECOTRAN code
agg_pb                              = agg_production ./ agg_biomass; % (1/y); (vertical vector: num_AggGrps X 1)
agg_pb(looky_agg_nutrients)         = 1; % nutrient pb defined as 1
% -------------------------------------------------------------------------


% step 10d: pb & pq & ae for fleets = landings / catch --------------------
agg_pq(looky_agg_fleets)            = agg_landings_TotalByFleet ./ agg_catch_TotalByFleet; % NOTE: fleet pb was corrected in step 10c
agg_ae(looky_agg_fleets)            = agg_landings_TotalByFleet ./ agg_catch_TotalByFleet; % ae of fleets = landings / catch
% -------------------------------------------------------------------------


% step 10e: change parameter values for phytoplankton, detritus, egg groups, & fleets
agg_qb(looky_agg_ANYPrimaryProd)	= agg_pb(looky_agg_ANYPrimaryProd) .* (1./agg_pq(looky_agg_ANYPrimaryProd));	% calculate qb for phytoplankton
agg_qb(looky_agg_ANYdetritus)       = agg_pb(looky_agg_ANYdetritus)    .* (1./agg_pq(looky_agg_ANYdetritus));       % calculate qb for detritus
agg_qb(looky_agg_eggs)              = agg_pb(looky_agg_eggs)           .* (1./agg_pq(looky_agg_eggs));              % calculate qb for eggs
agg_qb(looky_agg_fleets)        	= agg_pb(looky_agg_fleets)         .* (1./agg_pq(looky_agg_fleets));            % calculate qb for fleets
% -------------------------------------------------------------------------


% step 10f: phytoplankton consumption = production ------------------------
%           NOTE: agg_CONSUMPTION is not adjusted here, only total consumption by each primary producer group 
agg_consumption_domestic(looky_agg_ANYPrimaryProd)	= agg_production(looky_agg_ANYPrimaryProd); % pq = 1 for primary producers, so consumption must = production; (horizontal vector: 1 X num_AggGrps)
agg_consumption_total(looky_agg_ANYPrimaryProd)     = agg_production(looky_agg_ANYPrimaryProd); % pq = 1 for primary producers, so consumption must = production; (horizontal vector: 1 X num_AggGrps)
% *************************************************************************





% *************************************************************************
% STEP 11: calculate new ee terms------------------------------------------

% step 11a: pack values for the f_calculate_ee function -------------------
%          NOTE: these ee terms should match EwE values,
%                but they will not be the same as E2E values if there are EwE eggs defined as detritus groups
randomEwE.CONSUMPTION        	= agg_CONSUMPTION;          % (t/km2/y); (2D matrix: num_AggGrps X num_AggGrps)
randomEwE.production          	= agg_production;           % (t/km2/y); (vertical vector: num_AggGrps X 1); NOTE: fleet production is defined by landings by this point;
randomEwE.ba                 	= agg_ba;                   % (t/km2/y); (vertical vector: num_AggGrps X 1)
randomEwE.em                  	= agg_em;                   % (t/km2/y); (vertical vector: num_AggGrps X 1)
randomEwE.EggProduction     	= agg_EggProduction_forced;	% (t/km2/y); (vertical vector: num_AggGrps X 1)
randomEwE.looky_eggsANDdetritus	= looky_agg_eggsANDdetritus;
randomEwE.looky_nutrients       = looky_agg_nutrients;
% -------------------------------------------------------------------------


% step 11b: re-calculate ee based on aggregated model ---------------------
EcotrophicEfficiency            = f_calcEE_12292020(randomEwE); % ee terms are fractions of production
agg_ee                          = EcotrophicEfficiency.ee;
% -------------------------------------------------------------------------


% step 11c: some error checking, IS MODEL BALANCED? -----------------------
max_ee                          = max(agg_ee);
if max_ee > 1
    display('   -->WARNING: aggregated model is not balanced!')
end
% *************************************************************************





% *************************************************************************
% STEP 12: aggregate prey guild codes--------------------------------------
[num_guilds, toss]  = size(dat.PreyGuild_code);
agg_PreyGuild_code	= cell(num_guilds, 1);
for guild_loop = 1:num_guilds
    current_guild                   = eval(dat.PreyGuild_code{guild_loop}); % convert from strings
    agg_PreyGuild_code{guild_loop}  = unique(current_guild);                % build as vectors; filter any duplications
end
% *************************************************************************





% *************************************************************************
% STEP 13: aggregate PEDIGREE----------------------------------------------

% step 13a: parse preagg_PedigreeParameters -------------------------------

%                max EE; B (CV); PB (CV); QB (CV); PQ (CV); AE (CV); 
%                B (abs min); PB (abs min); QB (abs min); PQ (abs min); 
%                AE (abs min); B (abs max); PB (abs max); QB (abs max); 
%                PQ (abs max); AE (abs max); BA (CV); EM (CV); Egg rate (CV); 
%                Metabolism rate (CV); BA (abs min); EM (abs min); 
%                Egg rate (abs min); Metabolism rate (abs min); BA (abs max); 
%                EM (abs max); Egg rate (abs max); Metabolism rate (abs max); 
%                DetritusFate_feces (CV); DetritusFate_senescence (CV); 
%                ExcretionFate (CV); DetritusFate_feces (min); 
%                DetritusFate_senescence (min); ExcretionFate (min); 
%                DetritusFate_feces (max); DetritusFate_senescence (max); 
%                ExcretionFate (max)

% core pedigree terms (all are CV) -------
preagg_PEDIGREE_biomass                 = preagg_PedigreeParameters(:, 2);  % (CV); (vertical vector: num_E2Egrps X 1)
preagg_PEDIGREE_pb                      = preagg_PedigreeParameters(:, 3);  % (CV); (vertical vector: num_E2Egrps X 1)
% preagg_PEDIGREE_qb                      = preagg_PedigreeParameters(:, 4);  % (CV); (vertical vector: num_E2Egrps X 1); NOTE: qb uncertainty is derived below from pb & pq terms
preagg_PEDIGREE_pq                      = preagg_PedigreeParameters(:, 5);  % (CV); (vertical vector: num_E2Egrps X 1)
preagg_PEDIGREE_ae                      = preagg_PedigreeParameters(:, 6);  % (CV); (vertical vector: num_E2Egrps X 1)
preagg_PEDIGREE_ba                      = preagg_PedigreeParameters(:, 17); % (CV); (vertical vector: num_E2Egrps X 1); NOTE: must be expressed in terms of absolute production
preagg_PEDIGREE_em                      = preagg_PedigreeParameters(:, 18); % (CV); (vertical vector: num_E2Egrps X 1); NOTE: must be expressed in terms of absolute production
preagg_PEDIGREE_EggProduction_forced    = preagg_PedigreeParameters(:, 19); % (CV); (vertical vector: num_E2Egrps X 1)
preagg_PEDIGREE_Metabolism_forced       = preagg_PedigreeParameters(:, 20); % (CV); (vertical vector: num_E2Egrps X 1)
preagg_PEDIGREE_fate_feces              = preagg_PedigreeParameters(:, 29); % (CV); (vertical vector: num_E2Egrps X 1)
preagg_PEDIGREE_fate_senescence         = preagg_PedigreeParameters(:, 30); % (CV); (vertical vector: num_E2Egrps X 1)
preagg_PEDIGREE_fate_metabolism         = preagg_PedigreeParameters(:, 31);	% (CV); (vertical vector: num_E2Egrps X 1)

% other pedigree terms (min, max, etc.) -------
%       NOTE: these min & max fate pedigree terms have no meaning because
%             they must be defined for each detritus or nutrient pool and not be
%             represented as simple vectors: MinDetritusFate_feces,
%             MinDetritusFate_senescence, MinExcretionFate, MaxDetritusFate_feces,
%             MaxDetritusFate_senescence, MaxExcretionFate (FFF could correct
%             this in the future, but these terms have never been used as of 2018)
preagg_PEDIGREE_MaxEE                       = preagg_PedigreeParameters(:, 1);   % (dimensionless); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via production-weighted mean
preagg_PEDIGREE_MinBiomass                  = preagg_PedigreeParameters(:, 7);   % (t/km2); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via simple sum
preagg_PEDIGREE_MinPB                       = preagg_PedigreeParameters(:, 8);   % (1/y); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via production-weighted mean
preagg_PEDIGREE_MinQB                       = preagg_PedigreeParameters(:, 9);   % (1/y); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via production-weighted mean
preagg_PEDIGREE_MinPQ                       = preagg_PedigreeParameters(:, 10);  % (fraction: 0-1); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via production-weighted mean
preagg_PEDIGREE_MinAE                       = preagg_PedigreeParameters(:, 11);  % (fraction: 0-1); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via production-weighted mean

preagg_PEDIGREE_MaxBiomass                  = preagg_PedigreeParameters(:, 12);  % (t/km2); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via simple sum
preagg_PEDIGREE_MaxPB                       = preagg_PedigreeParameters(:, 13);  % (1/y); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via production-weighted mean
preagg_PEDIGREE_MaxQB                       = preagg_PedigreeParameters(:, 14);  % (1/y); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via production-weighted mean
preagg_PEDIGREE_MaxPQ                       = preagg_PedigreeParameters(:, 15);  % (fraction: 0-1); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via production-weighted mean
preagg_PEDIGREE_MaxAE                       = preagg_PedigreeParameters(:, 16);  % (fraction: 0-1); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via production-weighted mean

preagg_PEDIGREE_MinBA                       = preagg_PedigreeParameters(:, 21);  % (t/km2/y); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via simple sum;
preagg_PEDIGREE_MinEM                       = preagg_PedigreeParameters(:, 22);  % (t/km2/y); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via simple sum;
preagg_PEDIGREE_MinEggProduction_forced     = preagg_PedigreeParameters(:, 23);  % (t/km2/y); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via simple sum;
preagg_PEDIGREE_MinMetabolism_forced        = preagg_PedigreeParameters(:, 24);  % (t/km2/y); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via simple sum;

preagg_PEDIGREE_MaxBA                       = preagg_PedigreeParameters(:, 25);  % (t/km2/y); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via simple sum;
preagg_PEDIGREE_MaxEM                       = preagg_PedigreeParameters(:, 26);  % (t/km2/y); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via simple sum;
preagg_PEDIGREE_MaxEggProduction_forced     = preagg_PedigreeParameters(:, 27);  % (t/km2/y); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via simple sum;
preagg_PEDIGREE_MaxMetabolism_forced        = preagg_PedigreeParameters(:, 28);  % (t/km2/y); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via simple sum;
% -------------------------------------------------------------------------


% step 13b: convert from CV to variance -----------------------------------
preagg_VAR_biomass               	= (preagg_PEDIGREE_biomass	            .* preagg_biomass).^2;                  % (VAR); (vertical vector: num_E2Egrps X 1)
preagg_VAR_pb                       = (preagg_PEDIGREE_pb                   .* preagg_pb).^2;                       % (VAR); (vertical vector: num_E2Egrps X 1)
% preagg_VAR_qb                       = (preagg_PEDIGREE_qb                   .* preagg_qb).^2;                       % (VAR); (vertical vector: num_E2Egrps X 1); NOTE: qb uncertainty is derived below from pb & pq terms
preagg_VAR_pq                       = (preagg_PEDIGREE_pq                   .* preagg_pq).^2;                       % (VAR); (vertical vector: num_E2Egrps X 1)
preagg_VAR_ae                       = (preagg_PEDIGREE_ae                   .* preagg_ae).^2;                       % (VAR); (vertical vector: num_E2Egrps X 1)
preagg_VAR_ba                       = (preagg_PEDIGREE_ba                   .* preagg_ba).^2;                       % (VAR); (vertical vector: num_E2Egrps X 1); NOTE: variance of ba production rate
preagg_VAR_em                       = (preagg_PEDIGREE_em                   .* preagg_em).^2;                       % (VAR); (vertical vector: num_E2Egrps X 1); NOTE: variance of em production rate

preagg_VAR_LANDINGS              	= (preagg_PEDIGREE_LANDINGS             .* preagg_LANDINGS).^2;                 % (VAR); (2D matrix: num_E2Egrps X num_preagg_fleets)
preagg_VAR_DISCARDS             	= (preagg_PEDIGREE_DISCARDS             .* preagg_DISCARDS).^2;             	% (VAR); (2D matrix: num_E2Egrps X num_preagg_fleets)
preagg_VAR_CATCH                    = preagg_VAR_LANDINGS + preagg_VAR_DISCARDS;                                    % (VAR); (2D matrix: num_E2Egrps X num_preagg_fleets)

preagg_VAR_DIET                     = (preagg_PEDIGREE_DIET                 .* preagg_DIET).^2;                     % (VAR); (2D matrix: num_E2Egrps X num_E2Egrps);
preagg_VAR_DIET(:, looky_fleets)    = preagg_VAR_CATCH;                                                             % substitute fleet columns with CATCH VAR
preagg_VAR_diet_import              = (preagg_PEDIGREE_diet_import          .* preagg_diet_import).^2;           	% (VAR); (horizontal vector: 1 X num_E2Egrps)

preagg_VAR_EggProduction_forced     = (preagg_PEDIGREE_EggProduction_forced .* preagg_EggProduction_forced).^2;   	% (VAR); (vertical vector: num_E2Egrps X 1)
preagg_VAR_Metabolism_forced        = (preagg_PEDIGREE_Metabolism_forced    .* preagg_Metabolism_forced).^2;      	% (VAR); (vertical vector: num_E2Egrps X 1)
preagg_VAR_fate_feces               = (preagg_PEDIGREE_fate_feces           .* sum(preagg_fate_feces, 1)').^2;    	% (VAR); (vertical vector: num_E2Egrps X 1); note transpose; 1 vert vector scaled to sum of fate pools
preagg_VAR_fate_senescence          = (preagg_PEDIGREE_fate_senescence      .* sum(preagg_fate_senescence, 1)').^2;	% (VAR); (vertical vector: num_E2Egrps X 1); note transpose; 1 vert vector scaled to sum of fate pools
preagg_VAR_fate_metabolism          = (preagg_PEDIGREE_fate_metabolism      .* sum(preagg_fate_metabolism, 1)').^2;	% (VAR); (vertical vector: num_E2Egrps X 1); note transpose; 1 vert vector scaled to sum of fate pools
% -------------------------------------------------------------------------


% step 13c: weight physiological & fate VAR terms by group production -----
%           NOTE: multiplication by square of production when scaling variance terms
preagg_VAR_pb                       = preagg_VAR_pb                      .* preagg_production.^2; % (weighted VAR); (vertical vector: num_E2Egrps X 1)
% preagg_VAR_qb                       = preagg_VAR_qb                      .* preagg_production.^2; % (weighted VAR); (vertical vector: num_E2Egrps X 1); NOTE: qb uncertainty is derived below from pb & pq terms
preagg_VAR_pq                       = preagg_VAR_pq                      .* preagg_production.^2; % (weighted VAR); (vertical vector: num_E2Egrps X 1)
preagg_VAR_ae                       = preagg_VAR_ae                      .* preagg_production.^2; % (weighted VAR); (vertical vector: num_E2Egrps X 1)
preagg_VAR_EggProduction_forced     = preagg_VAR_EggProduction_forced    .* preagg_production.^2; % (weighted VAR); (vertical vector: num_E2Egrps X 1)
preagg_VAR_Metabolism_forced        = preagg_VAR_Metabolism_forced       .* preagg_production.^2; % (weighted VAR); (vertical vector: num_E2Egrps X 1)
preagg_VAR_fate_feces               = preagg_VAR_fate_feces    	         .* preagg_production.^2; % (weighted VAR); (vertical vector: num_E2Egrps X 1)
preagg_VAR_fate_senescence          = preagg_VAR_fate_senescence         .* preagg_production.^2; % (weighted VAR); (vertical vector: num_E2Egrps X 1)
preagg_VAR_fate_metabolism          = preagg_VAR_fate_metabolism         .* preagg_production.^2; % (weighted VAR); (vertical vector: num_E2Egrps X 1)

preagg_PEDIGREE_MaxEE             	= preagg_PEDIGREE_MaxEE              .* preagg_production;  % weighted; (dimensionless); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via production-weighted mean
preagg_PEDIGREE_MinPB            	= preagg_PEDIGREE_MinPB              .* preagg_production;  % weighted; (1/y); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via production-weighted mean
preagg_PEDIGREE_MinQB              	= preagg_PEDIGREE_MinQB              .* preagg_production;  % weighted; (1/y); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via production-weighted mean
preagg_PEDIGREE_MinPQ              	= preagg_PEDIGREE_MinPQ              .* preagg_production;  % weighted; (fraction: 0-1); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via production-weighted mean
preagg_PEDIGREE_MinAE              	= preagg_PEDIGREE_MinAE              .* preagg_production;  % weighted; (fraction: 0-1); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via production-weighted mean
preagg_PEDIGREE_MaxPB              	= preagg_PEDIGREE_MaxPB              .* preagg_production;  % weighted; (1/y); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via production-weighted mean
preagg_PEDIGREE_MaxQB             	= preagg_PEDIGREE_MaxQB              .* preagg_production;  % weighted; (1/y); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via production-weighted mean
preagg_PEDIGREE_MaxPQ              	= preagg_PEDIGREE_MaxPQ              .* preagg_production;  % weighted; (fraction: 0-1); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via production-weighted mean
preagg_PEDIGREE_MaxAE             	= preagg_PEDIGREE_MaxAE              .* preagg_production;  % weighted; (fraction: 0-1); (vertical vector: num_E2Egrps X 1); NOTE: aggregate via production-weighted mean
% -------------------------------------------------------------------------


% step 13d: initialize aggregated pedigree terms --------------------------
agg_VAR_biomass                         = zeros(num_AggGrps, 1);                    % (VAR); (vertical vector: num_AggGrps X 1)
agg_VAR_pb                              = zeros(num_AggGrps, 1);                    % (VAR); (vertical vector: num_AggGrps X 1)
% agg_VAR_qb                              = zeros(num_AggGrps, 1);                    % (VAR); (vertical vector: num_AggGrps X 1); NOTE: qb uncertainty is derived below from pb & pq terms
agg_VAR_pq                              = zeros(num_AggGrps, 1);                    % (VAR); (vertical vector: num_AggGrps X 1)
agg_VAR_ae                              = zeros(num_AggGrps, 1);                    % (VAR); (vertical vector: num_AggGrps X 1)
agg_VAR_ba                              = zeros(num_AggGrps, 1);                    % (VAR); (vertical vector: num_AggGrps X 1)
agg_VAR_em                              = zeros(num_AggGrps, 1);                    % (VAR); (vertical vector: num_AggGrps X 1)

horizontal_VAR_LANDINGS                 = zeros(num_AggGrps, num_preagg_fleets);    % (VAR); (2D matrix: num_AggGrps X num_preagg_fleets)
horizontal_VAR_DISCARDS             	= zeros(num_AggGrps, num_preagg_fleets);	% (VAR); (2D matrix: num_AggGrps X num_preagg_fleets)
horizontal_VAR_CATCH                  	= zeros(num_AggGrps, num_preagg_fleets);	% (VAR); (2D matrix: num_AggGrps X num_preagg_fleets)
agg_VAR_LANDINGS                        = zeros(num_AggGrps, num_agg_fleets);       % (VAR); (2D matrix: num_AggGrps X num_agg_fleets)
agg_VAR_DISCARDS                        = zeros(num_AggGrps, num_agg_fleets);       % (VAR); (2D matrix: num_AggGrps X num_agg_fleets)
agg_VAR_CATCH                           = zeros(num_AggGrps, num_agg_fleets);       % (VAR); (2D matrix: num_AggGrps X num_agg_fleets)

horizontal_VAR_DIET                 	= zeros(num_AggGrps, num_E2Egrps);          % (VAR); (2D matrix: num_AggGrps X num_E2Egrps)
agg_VAR_DIET                            = zeros(num_AggGrps, num_AggGrps);          % (VAR); (2D matrix: num_AggGrps X num_AggGrps)
agg_VAR_diet_import                     = zeros(1, num_AggGrps);                    % (VAR); (horizontal vector: 1 X num_AggGrps)

agg_VAR_EggProduction_forced            = zeros(num_AggGrps, 1);                    % (VAR); (vertical vector: num_AggGrps X 1)
agg_VAR_Metabolism_forced               = zeros(num_AggGrps, 1);                    % (VAR); (vertical vector: num_AggGrps X 1)
agg_VAR_fate_feces                      = zeros(num_AggGrps, 1);                    % (VAR); (vertical vector: num_AggGrps X 1)
agg_VAR_fate_senescence                 = zeros(num_AggGrps, 1);                    % (VAR); (vertical vector: num_AggGrps X 1)
agg_VAR_fate_metabolism                 = zeros(num_AggGrps, 1);                    % (VAR); (vertical vector: num_AggGrps X 1)

agg_PEDIGREE_MinBiomass                 = zeros(num_AggGrps, 1);                    % (t/km2); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MinBA                      = zeros(num_AggGrps, 1);                    % (t/km2); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MinEM                      = zeros(num_AggGrps, 1);                    % (t/km2); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MinEggProduction_forced	= zeros(num_AggGrps, 1);                    % (t/km2); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MinMetabolism_forced       = zeros(num_AggGrps, 1);                    % (t/km2); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MaxBiomass                 = zeros(num_AggGrps, 1);                    % (t/km2); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MaxBA                      = zeros(num_AggGrps, 1);                    % (t/km2/y); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MaxEM                      = zeros(num_AggGrps, 1);                    % (t/km2/y); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MaxEggProduction_forced	= zeros(num_AggGrps, 1);                    % (t/km2); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MaxMetabolism_forced       = zeros(num_AggGrps, 1);                    % (t/km2); (vertical vector: num_AggGrps X 1)

agg_PEDIGREE_MaxEE                      = zeros(num_AggGrps, 1);                    % (dimensionless); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MinPB                      = zeros(num_AggGrps, 1);                    % (1/y); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MinQB                      = zeros(num_AggGrps, 1);                    % (1/y); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MinPQ                      = zeros(num_AggGrps, 1);                    % (fraction: 0-1); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MinAE                      = zeros(num_AggGrps, 1);                    % (fraction: 0-1); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MaxPB                      = zeros(num_AggGrps, 1);                    % (1/y); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MaxQB                      = zeros(num_AggGrps, 1);                    % (1/y); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MaxPQ                      = zeros(num_AggGrps, 1);                    % (fraction: 0-1); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MaxAE                      = zeros(num_AggGrps, 1);                    % (fraction: 0-1); (vertical vector: num_AggGrps X 1)
% -------------------------------------------------------------------------


% step 13e: aggregate first set of terms ----------------------------------
% NOTE: RULES FOR AGGREGATING UNCERTAINTIES--------------------------------
%
%           RULE: VAR     = STD^2 = (mean * CV)^2
%
%           RULE: CV      = STD / mean
%
%           RULE: multiply variance by a constant
%                       var(C*A)  = C^2 * var(A) (Multiplying a random variable by a constant increases the variance by the square of the constant)
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

for TrophicGrp_loop = 1:num_AggGrps
    
    % pick current groups for aggregation -------
	current_AggCode                                 = AggCodes(TrophicGrp_loop);               	% current aggregation group; AggCodes is the "unique" code list
    looky_currentAggCode                            = find(PreaggCodes == current_AggCode);     % row numbers of all groups included in current aggregated group
    num_currentAgg                                  = length(looky_currentAggCode);         	% number of groups in the current aggregation
    
    % summed variance terms -------
    current_VAR_biomass                             = preagg_VAR_biomass(looky_currentAggCode);	% (VAR); (vertical vector: num_currentAgg X 1)
    current_VAR_ba                                  = preagg_VAR_ba(looky_currentAggCode);     	% (VAR); (vertical vector: num_currentAgg X 1); NOTE: variance of ba production rate
    current_VAR_em                                  = preagg_VAR_em(looky_currentAggCode);     	% (VAR); (vertical vector: num_currentAgg X 1); NOTE: variance of em production rate
    
    agg_VAR_biomass(TrophicGrp_loop)                = sum(current_VAR_biomass);                 % (VAR); (vertical vector: num_AggGrps X 1)
    agg_VAR_ba(TrophicGrp_loop)                     = sum(current_VAR_ba);                      % (VAR); (vertical vector: num_AggGrps X 1); NOTE: variance of ba production rate
    agg_VAR_em(TrophicGrp_loop)                     = sum(current_VAR_em);                      % (VAR); (vertical vector: num_AggGrps X 1); NOTE: variance of em production rate
    
    % production-weighted variance terms -------
    current_VAR_pb                                  = preagg_VAR_pb(looky_currentAggCode);                    	% (weighted VAR); (vertical vector: num_currentAgg X 1)
%     current_VAR_qb                                  = preagg_VAR_qb(looky_currentAggCode);                     	% (weighted VAR); (vertical vector: num_currentAgg X 1); NOTE: qb uncertainty is derived below from pb & pq terms
    current_VAR_pq                                  = preagg_VAR_pq(looky_currentAggCode);                    	% (weighted VAR); (vertical vector: num_currentAgg X 1)
    current_VAR_ae                                  = preagg_VAR_ae(looky_currentAggCode);                    	% (weighted VAR); (vertical vector: num_currentAgg X 1)
    current_VAR_EggProduction_forced                = preagg_VAR_EggProduction_forced(looky_currentAggCode);  	% (weighted VAR); (vertical vector: num_currentAgg X 1)
    current_VAR_Metabolism_forced                   = preagg_VAR_Metabolism_forced(looky_currentAggCode);    	% (weighted VAR); (vertical vector: num_currentAgg X 1)
    current_VAR_fate_feces                          = preagg_VAR_fate_feces(looky_currentAggCode);              % (weighted VAR); (vertical vector: num_currentAgg X 1)
    current_VAR_fate_senescence                     = preagg_VAR_fate_senescence(looky_currentAggCode);         % (weighted VAR); (vertical vector: num_currentAgg X 1)
    current_VAR_fate_metabolism                  	= preagg_VAR_fate_metabolism(looky_currentAggCode);      	% (weighted VAR); (vertical vector: num_currentAgg X 1)

	agg_VAR_pb(TrophicGrp_loop)                  	= sum(current_VAR_pb)                   / agg_production(TrophicGrp_loop).^2; % weighted mean variance; (vertical vector: num_AggGrps X 1)
% 	agg_VAR_qb(TrophicGrp_loop)                   	= sum(current_VAR_qb)                   / agg_production(TrophicGrp_loop).^2; % weighted mean variance; (vertical vector: num_AggGrps X 1); NOTE: qb uncertainty is derived below from pb & pq terms
	agg_VAR_pq(TrophicGrp_loop)                   	= sum(current_VAR_pq)                   / agg_production(TrophicGrp_loop).^2; % weighted mean variance; (vertical vector: num_AggGrps X 1)
	agg_VAR_ae(TrophicGrp_loop)                    	= sum(current_VAR_ae)                   / agg_production(TrophicGrp_loop).^2; % weighted mean variance; (vertical vector: num_AggGrps X 1)
	agg_VAR_EggProduction_forced(TrophicGrp_loop)	= sum(current_VAR_EggProduction_forced) / agg_production(TrophicGrp_loop).^2; % weighted mean variance; (vertical vector: num_AggGrps X 1)
	agg_VAR_Metabolism_forced(TrophicGrp_loop)    	= sum(current_VAR_Metabolism_forced)    / agg_production(TrophicGrp_loop).^2; % weighted mean variance; (vertical vector: num_AggGrps X 1)
    agg_VAR_fate_feces(TrophicGrp_loop)             = sum(current_VAR_fate_feces)           / agg_production(TrophicGrp_loop).^2; % weighted mean variance; (vertical vector: num_AggGrps X 1)
	agg_VAR_fate_senescence(TrophicGrp_loop)        = sum(current_VAR_fate_senescence)      / agg_production(TrophicGrp_loop).^2; % weighted mean variance; (vertical vector: num_AggGrps X 1)
	agg_VAR_fate_metabolism(TrophicGrp_loop)     	= sum(current_VAR_fate_metabolism)      / agg_production(TrophicGrp_loop).^2; % weighted mean variance; (vertical vector: num_AggGrps X 1)
    
    % simple summed fleet pedigree -------
    current_VAR_LANDINGS                            = preagg_VAR_LANDINGS(looky_currentAggCode, :);	% (VAR); (2D matrix: num_currentAgg X num_preagg_fleets)
    current_VAR_DISCARDS                            = preagg_VAR_DISCARDS(looky_currentAggCode, :);	% (VAR); (2D matrix: num_currentAgg X num_preagg_fleets)
    current_VAR_CATCH                               = preagg_VAR_CATCH(looky_currentAggCode, :);	% (VAR); (2D matrix: num_currentAgg X num_preagg_fleets)
    
    horizontal_VAR_LANDINGS(TrophicGrp_loop, :)     = sum(current_VAR_LANDINGS, 1);                 % (VAR); (2D matrix: num_AggGrps X num_preagg_fleets)
    horizontal_VAR_DISCARDS(TrophicGrp_loop, :)     = sum(current_VAR_DISCARDS, 1);                 % (VAR); (2D matrix: num_AggGrps X num_preagg_fleets)
    horizontal_VAR_CATCH(TrophicGrp_loop, :)        = sum(current_VAR_CATCH, 1);                    % (VAR); (2D matrix: num_AggGrps X num_preagg_fleets)
    
    % simple summed diet pedigree -------
    current_VAR_DIET                                = preagg_VAR_DIET(looky_currentAggCode, :);     % (VAR); (2D matrix: num_currentAgg X num_E2Egrps)
    horizontal_VAR_DIET(TrophicGrp_loop, :)         = sum(current_VAR_DIET, 1);                     % (VAR); (2D matrix: num_AggGrps X num_E2Egrps)
    
    % simple summed other pedigree terms -------
    current_PEDIGREE_MinBiomass                           	= preagg_PEDIGREE_MinBiomass(looky_currentAggCode);                 % (t/km2); (vertical vector: num_currentAgg X 1); NOTE: aggregate via simple sum
    current_PEDIGREE_MinBA                               	= preagg_PEDIGREE_MinBA(looky_currentAggCode);                      % (t/km2/y); (vertical vector: num_currentAgg X 1); NOTE: aggregate via simple sum;
    current_PEDIGREE_MinEM                                 	= preagg_PEDIGREE_MinEM(looky_currentAggCode);                      % (t/km2/y); (vertical vector: num_currentAgg X 1); NOTE: aggregate via simple sum;
    current_PEDIGREE_MinEggProduction_forced              	= preagg_PEDIGREE_MinEggProduction_forced(looky_currentAggCode);	% (t/km2/y); (vertical vector: num_currentAgg X 1); NOTE: aggregate via simple sum;
    current_PEDIGREE_MinMetabolism_forced                 	= preagg_PEDIGREE_MinMetabolism_forced(looky_currentAggCode);       % (t/km2/y); (vertical vector: num_currentAgg X 1); NOTE: aggregate via simple sum;
    current_PEDIGREE_MaxBiomass                           	= preagg_PEDIGREE_MaxBiomass(looky_currentAggCode);                 % (t/km2); (vertical vector: num_currentAgg X 1); NOTE: aggregate via simple sum
    current_PEDIGREE_MaxBA                                	= preagg_PEDIGREE_MaxBA(looky_currentAggCode);                      % (t/km2/y); (vertical vector: num_currentAgg X 1); NOTE: aggregate via simple sum;
    current_PEDIGREE_MaxEM                               	= preagg_PEDIGREE_MaxEM(looky_currentAggCode);                      % (t/km2/y); (vertical vector: num_currentAgg X 1); NOTE: aggregate via simple sum;
    current_PEDIGREE_MaxEggProduction_forced               	= preagg_PEDIGREE_MaxEggProduction_forced(looky_currentAggCode);    % (t/km2/y); (vertical vector: num_currentAgg X 1); NOTE: aggregate via simple sum;
    current_PEDIGREE_MaxMetabolism_forced                   = preagg_PEDIGREE_MaxMetabolism_forced(looky_currentAggCode);       % (t/km2/y); (vertical vector: num_currentAgg X 1); NOTE: aggregate via simple sum;
    
    agg_PEDIGREE_MinBiomass(TrophicGrp_loop)              	= sum(current_PEDIGREE_MinBiomass);                                 % (t/km2); (vertical vector: num_AggGrps X 1); NOTE: aggregate via simple sum
    agg_PEDIGREE_MinBA(TrophicGrp_loop)                   	= sum(current_PEDIGREE_MinBA);                                      % (t/km2/y); (vertical vector: num_AggGrps X 1); NOTE: aggregate via simple sum
    agg_PEDIGREE_MinEM(TrophicGrp_loop)                    	= sum(current_PEDIGREE_MinEM);                                      % (t/km2/y); (vertical vector: num_AggGrps X 1); NOTE: aggregate via simple sum
    agg_PEDIGREE_MinEggProduction_forced(TrophicGrp_loop)	= sum(current_PEDIGREE_MinEggProduction_forced);                    % (t/km2); (vertical vector: num_AggGrps X 1); NOTE: aggregate via simple sum
    agg_PEDIGREE_MinMetabolism_forced(TrophicGrp_loop)    	= sum(current_PEDIGREE_MinMetabolism_forced);                       % (t/km2); (vertical vector: num_AggGrps X 1); NOTE: aggregate via simple sum
    agg_PEDIGREE_MaxBiomass(TrophicGrp_loop)             	= sum(current_PEDIGREE_MaxBiomass);                                 % (t/km2); (vertical vector: num_AggGrps X 1); NOTE: aggregate via simple sum
    agg_PEDIGREE_MaxBA(TrophicGrp_loop)                   	= sum(current_PEDIGREE_MaxBA);                                      % (t/km2/y); (vertical vector: num_AggGrps X 1); NOTE: aggregate via simple sum
    agg_PEDIGREE_MaxEM(TrophicGrp_loop)                   	= sum(current_PEDIGREE_MaxEM);                                      % (t/km2/y); (vertical vector: num_AggGrps X 1); NOTE: aggregate via simple sum
    agg_PEDIGREE_MaxEggProduction_forced(TrophicGrp_loop)	= sum(current_PEDIGREE_MaxEggProduction_forced);                    % (t/km2); (vertical vector: num_AggGrps X 1); NOTE: aggregate via simple sum
    agg_PEDIGREE_MaxMetabolism_forced(TrophicGrp_loop)    	= sum(current_PEDIGREE_MaxMetabolism_forced);                       % (t/km2); (vertical vector: num_AggGrps X 1); NOTE: aggregate via simple sum

    % production-weighted other pedigree terms -------
    current_PEDIGREE_MaxEE                      = preagg_PEDIGREE_MaxEE(looky_currentAggCode);                     	% weighted; (dimensionless); (vertical vector: num_currentAgg X 1); NOTE: aggregate via production-weighted mean
    current_PEDIGREE_MinPB                      = preagg_PEDIGREE_MinPB(looky_currentAggCode);                     	% weighted; (1/y); (vertical vector: num_currentAgg X 1); NOTE: aggregate via production-weighted mean
    current_PEDIGREE_MinQB                      = preagg_PEDIGREE_MinQB(looky_currentAggCode);                     	% weighted; (1/y); (vertical vector: num_currentAgg X 1); NOTE: aggregate via production-weighted mean
    current_PEDIGREE_MinPQ                      = preagg_PEDIGREE_MinPQ(looky_currentAggCode);                    	% weighted; (fraction: 0-1); (vertical vector: num_currentAgg X 1); NOTE: aggregate via production-weighted mean
    current_PEDIGREE_MinAE                      = preagg_PEDIGREE_MinAE(looky_currentAggCode);                    	% weighted; (fraction: 0-1); (vertical vector: num_currentAgg X 1); NOTE: aggregate via production-weighted mean
    current_PEDIGREE_MaxPB                      = preagg_PEDIGREE_MaxPB(looky_currentAggCode);                    	% weighted; (1/y); (vertical vector: num_currentAgg X 1); NOTE: aggregate via production-weighted mean
    current_PEDIGREE_MaxQB                      = preagg_PEDIGREE_MaxQB(looky_currentAggCode);                    	% weighted; (1/y); (vertical vector: num_currentAgg X 1); NOTE: aggregate via production-weighted mean
    current_PEDIGREE_MaxPQ                      = preagg_PEDIGREE_MaxPQ(looky_currentAggCode);                    	% weighted; (fraction: 0-1); (vertical vector: num_currentAgg X 1); NOTE: aggregate via production-weighted mean
    current_PEDIGREE_MaxAE                      = preagg_PEDIGREE_MaxAE(looky_currentAggCode);                   	% weighted; (fraction: 0-1); (vertical vector: num_currentAgg X 1); NOTE: aggregate via production-weighted mean
    
	agg_PEDIGREE_MaxEE(TrophicGrp_loop)         = sum(current_PEDIGREE_MaxEE) / agg_production(TrophicGrp_loop);	% weighted; (dimensionless); (vertical vector: num_AggGrps X 1); NOTE: aggregate via production-weighted mean
    agg_PEDIGREE_MinPB(TrophicGrp_loop)         = sum(current_PEDIGREE_MinPB) / agg_production(TrophicGrp_loop);	% weighted; (1/y); (vertical vector: num_AggGrps X 1); NOTE: aggregate via production-weighted mean
    agg_PEDIGREE_MinQB(TrophicGrp_loop)         = sum(current_PEDIGREE_MinQB) / agg_production(TrophicGrp_loop);	% weighted; (1/y); (vertical vector: num_AggGrps X 1); NOTE: aggregate via production-weighted mean
    agg_PEDIGREE_MinPQ(TrophicGrp_loop)         = sum(current_PEDIGREE_MinPQ) / agg_production(TrophicGrp_loop);	% weighted; (fraction: 0-1); (vertical vector: num_AggGrps X 1); NOTE: aggregate via production-weighted mean
    agg_PEDIGREE_MinAE(TrophicGrp_loop)         = sum(current_PEDIGREE_MinAE) / agg_production(TrophicGrp_loop);    % weighted; (fraction: 0-1); (vertical vector: num_AggGrps X 1); NOTE: aggregate via production-weighted mean
    agg_PEDIGREE_MaxPB(TrophicGrp_loop)         = sum(current_PEDIGREE_MaxPB) / agg_production(TrophicGrp_loop);    % weighted; (1/y); (vertical vector: num_AggGrps X 1); NOTE: aggregate via production-weighted mean
    agg_PEDIGREE_MaxQB(TrophicGrp_loop)         = sum(current_PEDIGREE_MaxQB) / agg_production(TrophicGrp_loop);    % weighted; (1/y); (vertical vector: num_AggGrps X 1); NOTE: aggregate via production-weighted mean
    agg_PEDIGREE_MaxPQ(TrophicGrp_loop)         = sum(current_PEDIGREE_MaxPQ) / agg_production(TrophicGrp_loop);    % weighted; (fraction: 0-1); (vertical vector: num_AggGrps X 1); NOTE: aggregate via production-weighted mean
    agg_PEDIGREE_MaxAE(TrophicGrp_loop)         = sum(current_PEDIGREE_MaxAE) / agg_production(TrophicGrp_loop);    % weighted; (fraction: 0-1); (vertical vector: num_AggGrps X 1); NOTE: aggregate via production-weighted mean
    
end % end TrophicGrp_loop

% fix div/0 errors
agg_PEDIGREE_MinBiomass(isnan(agg_PEDIGREE_MinBiomass))                             = 0;	% fix div/0 errors; (t/km2); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MinBA(isnan(agg_PEDIGREE_MinBA))                                       = 0;	% fix div/0 errors; (t/km2); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MinEM(isnan(agg_PEDIGREE_MinEM))                                       = 0;  	% fix div/0 errors; (t/km2); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MinEggProduction_forced(isnan(agg_PEDIGREE_MinEggProduction_forced))   = 0;	% fix div/0 errors; (t/km2); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MinMetabolism_forced(isnan(agg_PEDIGREE_MinMetabolism_forced))         = 0;  	% fix div/0 errors; (t/km2); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MaxBiomass(isnan(agg_PEDIGREE_MaxBiomass))                             = 0;  	% fix div/0 errors; (t/km2); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MaxBA(isnan(agg_PEDIGREE_MaxBA))                                       = 0;  	% fix div/0 errors; (t/km2); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MaxEM(isnan(agg_PEDIGREE_MaxEM))                                       = 0; 	% fix div/0 errors; (t/km2); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MaxEggProduction_forced(isnan(agg_PEDIGREE_MaxEggProduction_forced))	= 0;    % fix div/0 errors; (t/km2); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MaxMetabolism_forced(isnan(agg_PEDIGREE_MaxMetabolism_forced))         = 0;  	% fix div/0 errors; (t/km2); (vertical vector: num_AggGrps X 1)

agg_PEDIGREE_MaxEE(isnan(agg_PEDIGREE_MaxEE))                                       = 0;	% fix div/0 errors; (dimensionless); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MinPB(isnan(agg_PEDIGREE_MinPB))                                       = 0;    % fix div/0 errors; (1/y); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MinQB(isnan(agg_PEDIGREE_MinQB))                                       = 0;    % fix div/0 errors; (1/y); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MinPQ(isnan(agg_PEDIGREE_MinPQ))                                       = 0;    % fix div/0 errors; (fraction: 0-1); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MinAE(isnan(agg_PEDIGREE_MinAE))                                       = 0;    % fix div/0 errors; (fraction: 0-1); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MaxPB(isnan(agg_PEDIGREE_MaxPB))                                       = 0;    % fix div/0 errors; (1/y); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MaxQB(isnan(agg_PEDIGREE_MaxQB))                                       = 0;    % fix div/0 errors; (1/y); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MaxPQ(isnan(agg_PEDIGREE_MaxPQ))                                       = 0;    % fix div/0 errors; (fraction: 0-1); (vertical vector: num_AggGrps X 1)
agg_PEDIGREE_MaxAE(isnan(agg_PEDIGREE_MaxAE))                                       = 0;    % fix div/0 errors; (fraction: 0-1); (vertical vector: num_AggGrps X 1)

agg_VAR_pb(isnan(agg_VAR_pb))                                                       = 0;    % fix div/0 errors;  weighted mean variance; (vertical vector: num_AggGrps X 1)
% agg_VAR_qb(isnan(agg_VAR_qb))                                                       = 0;    % fix div/0 errors;  weighted mean variance; (vertical vector: num_AggGrps X 1); NOTE: qb uncertainty is derived below from pb & pq terms
agg_VAR_pq(isnan(agg_VAR_pq))                                                       = 0;    % fix div/0 errors;  weighted mean variance; (vertical vector: num_AggGrps X 1)
agg_VAR_ae(isnan(agg_VAR_ae))                                                       = 0;    % fix div/0 errors;  weighted mean variance; (vertical vector: num_AggGrps X 1)
agg_VAR_EggProduction_forced(isnan(agg_VAR_EggProduction_forced))                   = 0;    % fix div/0 errors;  weighted mean variance; (vertical vector: num_AggGrps X 1)
agg_VAR_Metabolism_forced(isnan(agg_VAR_Metabolism_forced))                         = 0;    % fix div/0 errors;  weighted mean variance; (vertical vector: num_AggGrps X 1)
agg_VAR_fate_feces(isnan(agg_VAR_fate_feces))                                       = 0;    % fix div/0 errors;  weighted mean variance; (vertical vector: num_AggGrps X 1)
agg_VAR_fate_senescence(isnan(agg_VAR_fate_senescence))                          	= 0;    % fix div/0 errors;  weighted mean variance; (vertical vector: num_AggGrps X 1)
agg_VAR_fate_metabolism(isnan(agg_VAR_fate_metabolism))                             = 0;    % fix div/0 errors;  weighted mean variance; (vertical vector: num_AggGrps X 1)
% -------------------------------------------------------------------------


% step 13f: aggregate fleet pedigrees horizontally ------------------------
for fleet_loop = 1:num_agg_fleets
    
    % pick current fleets for aggregation ---------
    current_AggCode_fleets        	= AggCodes_fleets(fleet_loop);
    looky_current_fleet         	= find(PreaggCodes_fleets == current_AggCode_fleets); % clm numbers of all fleets included in current aggregated fleet

    % aggregate fleet pedigree ---------
    agg_VAR_LANDINGS(:, fleet_loop)	= sum(horizontal_VAR_LANDINGS(:, looky_current_fleet), 2);	% (VAR); (2D matrix: num_AggGrps X num_agg_fleets)
    agg_VAR_DISCARDS(:, fleet_loop)	= sum(horizontal_VAR_DISCARDS(:, looky_current_fleet), 2); 	% (VAR); (2D matrix: num_AggGrps X num_agg_fleets)
    agg_VAR_CATCH(:, fleet_loop) 	= sum(horizontal_VAR_CATCH(:, looky_current_fleet), 2);  	% (VAR); (2D matrix: num_AggGrps X num_agg_fleets)
    
end % fleet_loop

agg_VAR_landings_TotalByFleet   = sum(agg_VAR_LANDINGS, 1);  % (VAR); (horizontal vector: 1 X num_agg_fleets)
agg_VAR_discards_TotalByFleet   = sum(agg_VAR_DISCARDS, 1);  % (VAR); (horizontal vector: 1 X num_agg_fleets)
agg_VAR_catch_TotalByFleet      = sum(agg_VAR_CATCH, 1);     % (VAR); (horizontal vector: 1 X num_agg_fleets)
agg_VAR_landings_TotalByGroup   = sum(agg_VAR_LANDINGS, 2);  % (VAR); (vertical vector: num_AggGrps X 1)
agg_VAR_discards_TotalByGroup   = sum(agg_VAR_DISCARDS, 2);  % (VAR); (vertical vector: num_AggGrps X 1)
agg_VAR_catch_TotalByGroup      = sum(agg_VAR_CATCH, 2);     % (VAR); (vertical vector: num_AggGrps X 1)
% -------------------------------------------------------------------------


% step 13g: adjustment of specific groups ---------------------------------
%          NOTE: use looky_ addresses for aggregated model

% calculate derived agg_VAR_qb --------
agg_VAR_qb                              = f_VarianceDivision_12132018(agg_pb, agg_VAR_pb, agg_pq, agg_VAR_pq); % qb = pb / pq; (VAR); (horizontal vector: 1 X num_agg_fleets); var(A/B) = [mean(A)^2/mean(B)^2] * [(var(A)/mean(A)^2) + (var(B)/mean(B)^2)]

% uncertainties of fleet terms --------
agg_VAR_biomass(looky_agg_fleets)       = agg_VAR_catch_TotalByFleet'; % (VAR); (vertical vector: num_AggGrps X 1); NOTE transpose
VAR_fleet_pb                            = f_VarianceDivision_12132018(agg_landings_TotalByFleet, agg_VAR_landings_TotalByFleet, agg_catch_TotalByFleet, agg_VAR_catch_TotalByFleet); % fleet pb = landings / catch; (VAR); (horizontal vector: 1 X num_agg_fleets); var(A/B) = [mean(A)^2/mean(B)^2] * [(var(A)/mean(A)^2) + (var(B)/mean(B)^2)]
VAR_fleet_pq                            = VAR_fleet_pb;         % fleet pq = pb = landings / catch
VAR_fleet_ae                            = VAR_fleet_pb;         % fleet ae = pb = landings / catch
fleet_pb                                = agg_pb(looky_agg_fleets)';
fleet_pq                                = agg_pq(looky_agg_fleets)';
agg_VAR_pb(looky_agg_fleets)            = VAR_fleet_pb';        % (vertical vector: num_EwEgroups X 1)
agg_VAR_qb(looky_agg_fleets)            = 0;                    % qb for fleets are always 1 beacause q = b, so agg_VAR_qb for fleets = 0; (vertical vector: num_EwEgroups X 1)
agg_VAR_pq(looky_agg_fleets)            = VAR_fleet_pq';        % (vertical vector: num_EwEgroups X 1)
agg_VAR_ae(looky_agg_fleets)            = VAR_fleet_ae';        % (vertical vector: num_EwEgroups X 1)

% uncertainties of detritus & egg terms --------
agg_VAR_pb(looky_agg_eggsANDdetritus)	= 0; % defined as 0; detritus biomass is arbitrary and production is is estimated as functions of all other groups; consumption variance for detritus is set to 0, so the values of VAR_pb & VAR_qb for eggs & detritus play no role later on
agg_VAR_qb(looky_agg_eggsANDdetritus)	= 0; % qb = pb * (1/pq) and pq = 1 for detritus
agg_VAR_pq(looky_agg_eggsANDdetritus)	= 0; % pq = 1 for all detritus so variance = 0
agg_VAR_ae(looky_agg_eggsANDdetritus)	= 0;

% uncertainties of primary producer terms --------
agg_VAR_qb(looky_agg_ANYPrimaryProd)    = agg_VAR_pb(looky_agg_ANYPrimaryProd); % pb = qb for primary producers
% -------------------------------------------------------------------------


% step 13h: aggregate diet & import_diet pedigrees horizontally -----------
%         NOTE: weight diet variance by consumer consumption
%         NOTE: multiplication by square of consumption when scaling variance terms
%         NOTE: use domestic_consumption for Diet
%         NOTE: use import_consumption for import_diet
horizontal_VAR_DIET     = horizontal_VAR_DIET    .* (repmat(preagg_consumption_domestic, [num_AggGrps 1])).^2;	% (weighted VAR); (2D matrix: num_AggGrps X num_E2Egrps)
preagg_VAR_diet_import	= preagg_VAR_diet_import .* preagg_consumption_import.^2;                               % (weighted VAR); (horizontal vector: 1 X num_E2Egrps)

for consumer_loop  = 1:num_AggGrps
    
	% pick current consumers for aggregation ---------
    current_ConsumerCode                    = AggCodes(consumer_loop);                          % current aggregation consumer
	looky_currentConsumer                   = find(PreaggCodes == current_ConsumerCode);        % column numbers of all consumers included in current aggregated group
	num_currentConsumers                    = length(looky_currentConsumer);                    % number of groups in the current aggregation
        
    % aggregate diet pedigree ---------
    current_VAR_DIET                        = horizontal_VAR_DIET(:, looky_currentConsumer);    % (weighted VAR); (2D matrix: num_AggGrps X num_currentConsumers)
    current_VAR_diet_import                 = preagg_VAR_diet_import(looky_currentConsumer);    % (weighted VAR); (horizontal vector: 1 X num_currentConsumers)
    
	agg_VAR_DIET(:, consumer_loop)          = sum(current_VAR_DIET, 2)    ./ repmat(agg_consumption_domestic(consumer_loop), [num_AggGrps 1]).^2;	% weighted mean variance; (2D matrix: num_AggGrps X num_AggGrps)
    agg_VAR_diet_import(consumer_loop)      = sum(current_VAR_diet_import) / agg_consumption_import(consumer_loop)^2;                             % weighted mean variance; (horizontal vector: 1 X num_AggGrps)

end % consumer_loop

agg_VAR_DIET(isnan(agg_VAR_DIET))               	= 0;                    % correct for div/0 errors
agg_VAR_diet_import(isnan(agg_VAR_diet_import))     = 0;                    % correct for div/0 errors
agg_VAR_diet_domestic                               = agg_VAR_diet_import;	% this is the same as import diet VAR and is not sum(agg_VAR_DIET, 1) since everything NOT import must be domestic; use for consumption pedigree; (horizontal vector: 1 X num_AggGrps)

% NOTE: uncertainty of egg & detritus diet and consumption matrices should all be set to zero. 
%       These egg & detritus columns do not contribute to the EnergyBudget matrix.
agg_VAR_DIET(:, looky_agg_eggsANDdetritus)          = 0;
agg_VAR_diet_import(looky_agg_eggsANDdetritus)      = 0;
agg_VAR_diet_domestic(looky_agg_eggsANDdetritus)	= 0;

% input fleet diets; 
%        scale VAR_CATCH by total fleet catch (fleet biomasses), same as
%        done to calculate fleet diets (fleet_DIET = CATCH / catch_TotalByFleet)
repmat_agg_catch_TotalByFleet       = repmat(agg_catch_TotalByFleet, [num_AggGrps, 1]);     % (t/km2/y); (2D matrix: num_AggGrps X num_agg_fleets)
repmat_agg_VAR_catch_TotalByFleet	= repmat(agg_VAR_catch_TotalByFleet, [num_AggGrps, 1]);	% (VAR); (2D matrix: num_AggGrps X num_agg_fleets)
VAR_fleet_DIET                      = f_VarianceDivision_12132018(agg_CATCH, agg_VAR_CATCH, repmat_agg_catch_TotalByFleet, repmat_agg_VAR_catch_TotalByFleet);
agg_VAR_DIET(:, looky_agg_fleets)   = VAR_fleet_DIET; % (VAR); (2D matrix: num_AggGrps X num_AggGrps)
% -------------------------------------------------------------------------


% step 13i: calculate variance of consumption matrix and vectors ----------
%           NOTE: reminder on variance math rules
%                 var(A*B)  = [mean(A*B)^2] * [(var(A)/mean(A)^2) + (var(B)/mean(B)^2)]
%                 var(A/B)  = [mean(A)^2/mean(B)^2] * [(var(A)/mean(A)^2) + (var(B)/mean(B)^2)]
%           NOTE: agg_VAR_pb & agg_VAR_qb could have non-zero entries for detritus (pb = production/biomass), but leaving these as zero because detritus production will vary with the physiologies of consumers & bacteria
%           NOTE: agg_VAR_pq = 0 for primary producers & detritus because pq = 1 for primary producers & detritus
%           NOTE: agg_VAR_ae = 0 for detritus because ae = 0 for detritus; agg_VAR_ae should = 0 for primary producers because ae defined as 1 (doesn't = 1 because entered as non-0 value in .xlsm file)
%           NOTE: use consumption_domestic instead of consumption_total to calculate VAR_CONSUMPTION (since this matrix only considers domestic consumption and excludes import consumption)
%           NOTE: detritus entries have no meaning because the biomass value, and thus the pb value, are arbitrary
agg_VAR_consumption_total       	= f_VarianceMultiplication_12132018(agg_biomass', agg_VAR_biomass', agg_qb', agg_VAR_qb');                                          % (VAR); (horizontal vector: 1 X num_AggGrps); total consumption = b * qb; var(A*B) = [mean(A*B)^2] * [(var(A)/mean(A)^2) + (var(B)/mean(B)^2)]; NOTE transposes
agg_VAR_consumption_import        	= f_VarianceMultiplication_12132018(agg_consumption_total, agg_VAR_consumption_total, agg_diet_import, agg_VAR_diet_import);        % (VAR); (horizontal vector: 1 X num_AggGrps); import consumption = total consumption * import diet; var(A*B) = [mean(A*B)^2] * [(var(A)/mean(A)^2) + (var(B)/mean(B)^2)]
agg_VAR_consumption_domestic     	= f_VarianceMultiplication_12132018(agg_consumption_total, agg_VAR_consumption_total, (1-agg_diet_import), agg_VAR_diet_domestic);	% (VAR); (horizontal vector: 1 X num_AggGrps); domestic consumption = total consumption * total domestic diet fractions; var(A*B) = [mean(A*B)^2] * [(var(A)/mean(A)^2) + (var(B)/mean(B)^2)]

repmat_agg_consumption_domestic     = repmat(agg_consumption_domestic, [num_AggGrps, 1]);       % (t/km2/y); (2D matrix: num_AggGrps X num_AggGrps)
repmat_agg_VAR_consumption_domestic	= repmat(agg_VAR_consumption_domestic, [num_AggGrps, 1]);	% (VAR); (2D matrix: num_AggGrps X num_AggGrps)
agg_VAR_CONSUMPTION                 = f_VarianceMultiplication_12132018(repmat_agg_consumption_domestic, repmat_agg_VAR_consumption_domestic, agg_DIET, agg_VAR_DIET); % (VAR); agg_CONSUMPTION = agg_consumption_domestic * agg_DIET; (2D matrix: num_AggGrps X num_AggGrps);

% NOTE: uncertainty of egg & detritus diet and consumption matrices should all be set to zero.
%       NOTE: Values of b and qb used to calculate agg_VAR_CONSUMPTION for eggs & detritus are arbitrary
%       NOTE: egg & detritus columns in agg_CONSUMPTION & agg_VAR_CONSUMPTION do not contribute to the EnergyBudget matrix.
agg_VAR_CONSUMPTION(:, looky_agg_eggsANDdetritus)	= 0;


% plug catch VAR back into agg_VAR_CONSUMPTION fleet clms 
%       Only use the total-catch scaled catch in the DIET pedigree (step 13h) 
%       and not in the CONSUMPTION because this unnecessarily 
%       increases catch uncertainty in CONSUMPTION over that which was already pre-defined
%       (change made 11/26/2019)
agg_VAR_CONSUMPTION(:, looky_agg_fleets) = agg_VAR_CATCH;
% -------------------------------------------------------------------------


% step 13j: adjustment of em VAR ------------------------------------------
%           add fleet landing (em) variance
%           add import consumption variance (NOTE: scaling of consumption by pq to convert to production term)
agg_VAR_em(looky_agg_fleets)      	= agg_VAR_landings_TotalByFleet';           % fleet em = landings; (VAR); (vertical vector: num_AggGrps X 1); NOTE transpose; NOTE: use of landings means that fleet VAR is expressed in terms of production;
temp_VAR_consumption_import         = f_VarianceMultiplication_12132018(agg_consumption_import', agg_VAR_consumption_import', agg_pq, agg_VAR_pq); % import consumption variance in terms of production; (VAR); (vertical vector: num_AggGrps X 1); NOTE transpose
agg_VAR_em                          = agg_VAR_em + temp_VAR_consumption_import;	% add import consumption VAR; (VAR); (vertical vector: num_AggGrps X 1)
% -------------------------------------------------------------------------


% step 13k: convert from VAR to CV ----------------------------------------
PEDIGREE_biomass                    = sqrt(agg_VAR_biomass)               ./ agg_biomass;                 	% (CV); (vertical vector: num_AggGrps X 1)
PEDIGREE_pb                         = sqrt(agg_VAR_pb)                    ./ agg_pb;                      	% (CV); (vertical vector: num_AggGrps X 1)
PEDIGREE_qb                         = sqrt(agg_VAR_qb)                    ./ agg_qb;                      	% (CV); (vertical vector: num_AggGrps X 1)
PEDIGREE_pq                         = sqrt(agg_VAR_pq)                    ./ agg_pq;                      	% (CV); (vertical vector: num_AggGrps X 1)
PEDIGREE_ae                         = sqrt(agg_VAR_ae)                    ./ agg_ae;                      	% (CV); (vertical vector: num_AggGrps X 1)
PEDIGREE_ba                         = sqrt(agg_VAR_ba)                    ./ agg_ba;                       	% (CV); (vertical vector: num_AggGrps X 1); NOTE: CV of ba production rate
PEDIGREE_em                         = sqrt(agg_VAR_em)                    ./ abs(agg_em);                   % (CV); (vertical vector: num_AggGrps X 1); NOTE: CV of em production rate

PEDIGREE_CONSUMPTION                = sqrt(agg_VAR_CONSUMPTION)           ./ agg_CONSUMPTION;               % (CV); (2D matrix: num_AggGrps X num_AggGrps)
PEDIGREE_consumption_total          = sqrt(agg_VAR_consumption_total)     ./ agg_consumption_total;         % (CV); (horizontal vector: 1 X num_AggGrps)
PEDIGREE_consumption_domestic       = sqrt(agg_VAR_consumption_domestic)  ./ agg_consumption_domestic;      % (CV); (horizontal vector: 1 X num_AggGrps)
PEDIGREE_consumption_import         = sqrt(agg_VAR_consumption_import)    ./ agg_consumption_import;        % (CV); (horizontal vector: 1 X num_AggGrps)

PEDIGREE_LANDINGS                   = sqrt(agg_VAR_LANDINGS)              ./ agg_LANDINGS;                 	% (CV); (2D matrix: num_AggGrps X num_agg_fleets)
PEDIGREE_DISCARDS                   = sqrt(agg_VAR_DISCARDS)              ./ agg_DISCARDS;                  	% (CV); (2D matrix: num_AggGrps X num_agg_fleets)
PEDIGREE_CATCH                      = sqrt(agg_VAR_CATCH)                 ./ agg_CATCH;                   	% (CV); (2D matrix: num_AggGrps X num_agg_fleets)
PEDIGREE_landings_TotalByFleet      = sqrt(agg_VAR_landings_TotalByFleet) ./ agg_landings_TotalByFleet;     % (CV); (horizontal vector: 1 X num_agg_fleets)
PEDIGREE_discards_TotalByFleet      = sqrt(agg_VAR_discards_TotalByFleet) ./ agg_discards_TotalByFleet;     % (CV); (horizontal vector: 1 X num_agg_fleets)
PEDIGREE_catch_TotalByFleet         = sqrt(agg_VAR_catch_TotalByFleet)    ./ agg_catch_TotalByFleet;        % (CV); (horizontal vector: 1 X num_agg_fleets)
PEDIGREE_landings_TotalByGroup      = sqrt(agg_VAR_landings_TotalByGroup) ./ agg_landings_TotalByGroup;     % (CV); (vertical vector: num_AggGrps X 1)
PEDIGREE_discards_TotalByGroup      = sqrt(agg_VAR_discards_TotalByGroup) ./ agg_discards_TotalByGroup;     % (CV); (vertical vector: num_AggGrps X 1)
PEDIGREE_catch_TotalByGroup         = sqrt(agg_VAR_catch_TotalByGroup)    ./ agg_catch_TotalByGroup;        % (CV); (vertical vector: num_AggGrps X 1)

PEDIGREE_DIET                       = sqrt(agg_VAR_DIET)                 ./ agg_DIET;                   	% (CV); (2D matrix: num_AggGrps X num_AggGrps)
PEDIGREE_diet_import                = sqrt(agg_VAR_diet_import)          ./ agg_diet_import;             	% (CV); (horizontal vector: 1 X num_AggGrps)

PEDIGREE_EggProduction_forced       = sqrt(agg_VAR_EggProduction_forced) ./ agg_EggProduction_forced;     	% (CV); (vertical vector: num_AggGrps X 1)
PEDIGREE_Metabolism_forced          = sqrt(agg_VAR_Metabolism_forced)    ./ agg_Metabolism_forced;       	% (CV); (vertical vector: num_AggGrps X 1)

PEDIGREE_fate_feces                 = sqrt(agg_VAR_fate_feces)           ./ sum(agg_fate_feces, 1)';        % (CV); (vertical vector: num_AggGrps X 1); NOTE transpose
PEDIGREE_fate_senescence            = sqrt(agg_VAR_fate_senescence)      ./ sum(agg_fate_senescence, 1)';	% (CV); (vertical vector: num_AggGrps X 1); NOTE transpose
PEDIGREE_fate_metabolism         	= sqrt(agg_VAR_fate_metabolism)      ./ sum(agg_fate_metabolism, 1)';  	% (CV); (vertical vector: num_AggGrps X 1); NOTE transpose
% -------------------------------------------------------------------------


% step 13l: filter div/0 NaN errors ---------------------------------------
PEDIGREE_biomass(isnan(PEDIGREE_biomass))                        	= 0; % correct for div/0 errors
PEDIGREE_pb(isnan(PEDIGREE_pb))                                  	= 0; % correct for div/0 errors
PEDIGREE_qb(isnan(PEDIGREE_qb))                                   	= 0; % correct for div/0 errors
PEDIGREE_pq(isnan(PEDIGREE_pq))                                  	= 0; % correct for div/0 errors
PEDIGREE_ae(isnan(PEDIGREE_ae))                                 	= 0; % correct for div/0 errors
PEDIGREE_ba(isnan(PEDIGREE_ba))                                  	= 0; % correct for div/0 errors
PEDIGREE_em(isnan(PEDIGREE_em))                                   	= 0; % correct for div/0 errors

PEDIGREE_CONSUMPTION(isnan(PEDIGREE_CONSUMPTION))                   = 0; % correct for div/0 errors; (CV); (2D matrix: num_AggGrps X num_AggGrps)
PEDIGREE_consumption_total(isnan(PEDIGREE_consumption_total))     	= 0; % correct for div/0 errors; (CV); (horizontal vector: 1 X num_AggGrps)
PEDIGREE_consumption_domestic(isnan(PEDIGREE_consumption_domestic))	= 0; % correct for div/0 errors; (CV); (horizontal vector: 1 X num_AggGrps)
PEDIGREE_consumption_import(isnan(PEDIGREE_consumption_import))   	= 0; % correct for div/0 errors; (CV); (horizontal vector: 1 X num_AggGrps)
PEDIGREE_consumption_total(isinf(PEDIGREE_consumption_total))     	= 0; % correct for div/0 errors; (CV); (horizontal vector: 1 X num_AggGrps)
PEDIGREE_consumption_domestic(isinf(PEDIGREE_consumption_domestic))	= 0; % correct for div/0 errors; (CV); (horizontal vector: 1 X num_AggGrps)
PEDIGREE_consumption_import(isinf(PEDIGREE_consumption_import))   	= 0; % correct for div/0 errors; (CV); (horizontal vector: 1 X num_AggGrps)

PEDIGREE_LANDINGS(isnan(PEDIGREE_LANDINGS))                             = 0; % correct for div/0 errors
PEDIGREE_DISCARDS(isnan(PEDIGREE_DISCARDS))                             = 0; % correct for div/0 errors
PEDIGREE_CATCH(isnan(PEDIGREE_CATCH))                                   = 0; % correct for div/0 errors
PEDIGREE_landings_TotalByFleet(isnan(PEDIGREE_landings_TotalByFleet))	= 0; % correct for div/0 errors
PEDIGREE_discards_TotalByFleet(isnan(PEDIGREE_discards_TotalByFleet))	= 0; % correct for div/0 errors
PEDIGREE_catch_TotalByFleet(isnan(PEDIGREE_catch_TotalByFleet))         = 0; % correct for div/0 errors
PEDIGREE_landings_TotalByGroup(isnan(PEDIGREE_landings_TotalByGroup))   = 0; % correct for div/0 errors
PEDIGREE_discards_TotalByGroup(isnan(PEDIGREE_discards_TotalByGroup))   = 0; % correct for div/0 errors
PEDIGREE_catch_TotalByGroup(isnan(PEDIGREE_catch_TotalByGroup))         = 0; % correct for div/0 errors

PEDIGREE_DIET(isinf(PEDIGREE_DIET))                             	= 0; % correct for div/0 errors
PEDIGREE_DIET(isnan(PEDIGREE_DIET))                             	= 0; % correct for div/0 errors
PEDIGREE_diet_import(isnan(PEDIGREE_diet_import))                 	= 0; % correct for div/0 errors

PEDIGREE_EggProduction_forced(isnan(PEDIGREE_EggProduction_forced))	= 0; % correct for div/0 errors
PEDIGREE_Metabolism_forced(isnan(PEDIGREE_Metabolism_forced))     	= 0; % correct for div/0 errors

PEDIGREE_fate_feces(isnan(PEDIGREE_fate_feces))                     = 0; % correct for div/0 errors
PEDIGREE_fate_senescence(isnan(PEDIGREE_fate_senescence))           = 0; % correct for div/0 errors
PEDIGREE_fate_metabolism(isnan(PEDIGREE_fate_metabolism))        	= 0; % correct for div/0 errors
% *************************************************************************





% *************************************************************************
% STEP 14: pack EwEResult for export---------------------------------------
EwEResult.fname_AggregateBiologicalModel	= fname_AggregateBiologicalModel;

% step 14a: group type definitions ----------------------------------------
EwEResult.GroupTypeDef_ANYNitroNutr         = dat.GroupTypeDef_ANYNitroNutr;
EwEResult.GroupTypeDef_NO3                  = dat.GroupTypeDef_NO3;
EwEResult.GroupTypeDef_plgcNH4              = dat.GroupTypeDef_plgcNH4;
EwEResult.GroupTypeDef_bnthNH4              = dat.GroupTypeDef_bnthNH4;
EwEResult.GroupTypeDef_ANYPrimaryProd       = dat.GroupTypeDef_ANYPrimaryProd;
EwEResult.GroupTypeDef_LrgPhyto             = dat.GroupTypeDef_LrgPhyto;
EwEResult.GroupTypeDef_SmlPhyto             = dat.GroupTypeDef_SmlPhyto;
EwEResult.GroupTypeDef_Macrophytes          = dat.GroupTypeDef_Macrophytes;
EwEResult.GroupTypeDef_ANYConsumer          = dat.GroupTypeDef_ANYConsumer;
EwEResult.GroupTypeDef_ConsumPlgcPlankton   = dat.GroupTypeDef_ConsumPlgcPlankton;
EwEResult.GroupTypeDef_ConsumPlgcNekton     = dat.GroupTypeDef_ConsumPlgcNekton;
EwEResult.GroupTypeDef_ConsumPlgcWrmBlood   = dat.GroupTypeDef_ConsumPlgcWrmBlood;
EwEResult.GroupTypeDef_ConsumBntcInvert     = dat.GroupTypeDef_ConsumBntcInvert;
EwEResult.GroupTypeDef_ConsumBntcVert       = dat.GroupTypeDef_ConsumBntcVert;
EwEResult.GroupTypeDef_ConsumBnthWrmBlood   = dat.GroupTypeDef_ConsumBnthWrmBlood;
EwEResult.GroupTypeDef_fleet                = dat.GroupTypeDef_fishery;
EwEResult.GroupTypeDef_eggs                 = dat.GroupTypeDef_eggs;
EwEResult.GroupTypeDef_offal                = dat.GroupTypeDef_offal;
EwEResult.GroupTypeDef_terminalPlgcDetr     = dat.GroupTypeDef_terminalPlgcDetr;
EwEResult.GroupTypeDef_terminalBnthDetr     = dat.GroupTypeDef_terminalBnthDetr;
EwEResult.GroupTypeDef_ANYDetritus          = dat.GroupTypeDef_ANYDetritus;
EwEResult.GroupTypeDef_micrograzers         = dat.GroupTypeDef_micrograzers;
EwEResult.GroupTypeDef_bacteria             = dat.GroupTypeDef_bacteria;
EwEResult.GroupTypeDef_ba                   = dat.GroupTypeDef_BA;
EwEResult.GroupTypeDef_em                   = dat.GroupTypeDef_EM;
EwEResult.GroupTypeDef_import               = dat.GroupTypeDef_import;            % FFF unused group type for the future
% -------------------------------------------------------------------------


% step 14b: numbers of group types ----------------------------------------
EwEResult.num_grps                          = num_AggGrps;                  % number of aggregated groups, including nutrients & fleets
EwEResult.num_NO3                           = num_agg_NO3;
EwEResult.num_plgcNH4                       = num_agg_plgcNH4;
EwEResult.num_bnthNH4                       = num_agg_bnthNH4;
EwEResult.num_NH4                           = num_agg_NH4;
EwEResult.num_nutrients                     = num_agg_nutrients;
EwEResult.num_PrimaryProducers              = num_agg_PrimaryProducers;
EwEResult.num_consumers                     = num_agg_consumers;
EwEResult.num_micrograzers                  = num_agg_micrograzers;
EwEResult.num_bacteria                      = num_agg_bacteria;
EwEResult.num_eggs                          = num_agg_eggs;
EwEResult.num_ANYdetritus                 	= num_agg_ANYdetritus;          % ANY detritus groups (w/o eggs)
EwEResult.num_terminalPLGCdetritus          = num_agg_terminalPLGCdetritus;
EwEResult.num_terminalBNTHdetritus          = num_agg_terminalBNTHdetritus;
EwEResult.num_terminalANYdetritus        	= num_agg_terminalANYdetritus;
EwEResult.num_eggsANDdetritus               = num_agg_eggsANDdetritus;
EwEResult.num_livingANDdetritus             = num_agg_livingANDdetritus;
EwEResult.num_fleets                        = num_agg_fleets;
EwEResult.num_livingANDfleets               = num_agg_livingANDfleets;
EwEResult.num_NONnutrients                  = num_agg_NONnutrients;
% -------------------------------------------------------------------------


% step 14c: labels, GroupType, & number codes -----------------------------
EwEResult.label                             = AggLabel;                 % text cells; (vertical vector: num_AggGrps X 1)
EwEResult.GroupType                         = AggGroupType;             % includes nutrients & fleets; (vertical vector: num_AggGrps X 1)
EwEResult.CodeNumber                      	= AggNumber;                % (vertical vector: num_AggGrps X 1)
% -------------------------------------------------------------------------


% step 14d: parameters ----------------------------------------------------
%          NOTE: include nutrients & fleets
EwEResult.biomass                           = agg_biomass;       % (vertical vector: num_AggGrps X 1)
EwEResult.production                        = agg_production;    % (vertical vector: num_AggGrps X 1)
EwEResult.pb                                = agg_pb;            % (1/y); (vertical vector: num_AggGrps X 1)
EwEResult.qb                                = agg_qb;            % (vertical vector: num_AggGrps X 1)
EwEResult.pq                                = agg_pq;            % (vertical vector: num_AggGrps X 1)
EwEResult.ae                                = agg_ae;            % (vertical vector: num_AggGrps X 1)
EwEResult.TL                                = agg_TL;            % (vertical vector: num_AggGrps X 1)
EwEResult.ee                                = agg_ee;            % (vertical vector: num_AggGrps X 1)
EwEResult.ba                                = agg_ba;            % (t/km2/y); (vertical vector: num_AggGrps X 1); NOTE: ba production rate
EwEResult.em                                = agg_em;            % (t/km2/y); (vertical vector: num_AggGrps X 1); NOTE: em production rate
% -------------------------------------------------------------------------


% step 14e: diet & consumption --------------------------------------------
%          NOTE: include nutrients & fleets
EwEResult.CONSUMPTION                       = agg_CONSUMPTION;           % (t/km2/y); (2D matrix: num_AggGrps X num_AggGrps)
EwEResult.consumption_domestic              = agg_consumption_domestic'; % (t/km2/y); (vertical vector: num_AggGrps X 1); NOTE transpose;
EwEResult.consumption_import                = agg_consumption_import';   % (t/km2/y); (vertical vector: num_AggGrps X 1); NOTE transpose
EwEResult.consumption_total                 = agg_consumption_total';    % (t/km2/y); (vertical vector: num_AggGrps X 1); NOTE transpose
EwEResult.DIET                              = agg_DIET;                  % (2D matrix: num_AggGrps X num_AggGrps)
EwEResult.import_diet                       = agg_diet_import;           % (horizontal vector: 1 X num_AggGrps)
% -------------------------------------------------------------------------


% step 14f: fleet info ----------------------------------------------------
EwEResult.LANDINGS                        	= agg_LANDINGS;                 % (t/km2/y); (2D matrix: num_AggGrps X num_agg_fleets)
EwEResult.DISCARDS                        	= agg_DISCARDS;                 % (t/km2/y); (2D matrix: num_AggGrps X num_agg_fleets)
EwEResult.CATCH                        	    = agg_CATCH;                    % (t/km2/y); (2D matrix: num_AggGrps X num_agg_fleets)
EwEResult.landings_TotalByFleet           	= agg_landings_TotalByFleet;	% (t/km2/y); (horizontal vector: 1 X num_agg_fleets)
EwEResult.discards_TotalByFleet         	= agg_discards_TotalByFleet;	% (t/km2/y); (horizontal vector: 1 X num_agg_fleets)
EwEResult.catch_TotalByFleet             	= agg_catch_TotalByFleet;       % (t/km2/y); (horizontal vector: 1 X num_agg_fleets)
EwEResult.landings_TotalByGroup             = agg_landings_TotalByGroup;	% (t/km2/y); (vertical vector: num_agg_fleets X 1)
EwEResult.discards_TotalByGroup             = agg_discards_TotalByGroup;	% (t/km2/y); (vertical vector: num_agg_fleets X 1)
EwEResult.catch_TotalByGroup                = agg_catch_TotalByGroup;       % (t/km2/y); (vertical vector: num_agg_fleets X 1)
% -------------------------------------------------------------------------


% step 14g: detritus & egg fates, forced metabolism & egg rates -----------
EwEResult.fate_EwEdetritus                  = agg_fate_EwEdetritus;         % (proportions); (2D matrix: num_AggGrps X num_agg_ANYdetritus); NOTE: does not include EwE egg fate columns
EwEResult.fate_eggs                         = agg_fate_EwEeggs;             % (proportions); (2D matrix: num_AggGrps X num_agg_eggs)
EwEResult.fate_feces                        = agg_fate_feces;               % (proportions); (2D matrix: num_TerminalDetritus X num_AggGrps) QQQ confirm rows are aggregated
EwEResult.fate_senescence                   = agg_fate_senescence;          % (proportions); (2D matrix: num_TerminalDetritus X num_AggGrps) QQQ confirm rows are aggregated
EwEResult.fate_metabolism                   = agg_fate_metabolism;          % (proportions); (2D matrix: num_NH4 X num_AggGrps) QQQ confirm rows are aggregated
% -------------------------------------------------------------------------


% step 14h: forced metabolism & egg rates ---------------------------------
EwEResult.EggProduction_forced              = agg_EggProduction_forced;	% (t/km2/y); (vertical vector: num_AggGrps X 1)
EwEResult.Metabolism_forced                 = agg_Metabolism_forced;   	% (t/km2/y); (vertical vector: num_AggGrps X 1)
% -------------------------------------------------------------------------


% step 14i: retention & production loss scalers ---------------------------
EwEResult.ProductionLossScaler              = agg_ProductionLossScaler;      % (scaler: 0 to 1); (vertical vector: num_AggGrps X 1)
EwEResult.RetentionScaler                   = agg_RetentionScaler;           % (scaler: 0 to 1); (vertical vector: num_AggGrps X 1)
% -------------------------------------------------------------------------


% step 14j: functional response terms -------------------------------------
EwEResult.FunctionalResponseParams          = agg_FunctionalResponseParams;  % (vertical vector: num_AggGrps X 4)
EwEResult.FunctionalResponse_matrix         = agg_FunctionalResponse_matrix; % (2D matrix: num_AggGrps X num_AggGrps); FFF for future
% -------------------------------------------------------------------------


% step 14k: nutrient & detritus cycling terms -----------------------------
%          NOTE: QQQ I've not yet any aggregation options for these
EwEResult.PelagicBacterialReduction         = dat.PelagicBacterialReduction;
EwEResult.BenthicBacterialReduction         = dat.BenthicBacterialReduction;
EwEResult.Oxidation_NH4                     = dat.Oxidation_NH4;
EwEResult.PhytoUptake_NH4                   = dat.PhytoUptake_NH4; % or use = dat.PhytoUptake_NH4(~isnan(dat.PhytoUptake_NH4))
EwEResult.PhytoUptake_NO3                   = dat.PhytoUptake_NO3; % or use = dat.PhytoUptake_NO3(~isnan(dat.PhytoUptake_NO3))
% *************************************************************************





% *************************************************************************
% STEP 15: pack PEDIGREE info for export-----------------------------------
PEDIGREE.biomass_CV                 = PEDIGREE_biomass;                     % (CV); (vertical vector: num_AggGrps X 1)
PEDIGREE.pb_CV                      = PEDIGREE_pb;                          % (CV); (vertical vector: num_AggGrps X 1)
PEDIGREE.qb_CV                      = PEDIGREE_qb;                          % (CV); (vertical vector: num_AggGrps X 1)
PEDIGREE.pq_CV                      = PEDIGREE_pq;                          % (CV); (vertical vector: num_AggGrps X 1)
PEDIGREE.ae_CV                      = PEDIGREE_ae;                          % (CV); (vertical vector: num_AggGrps X 1)
PEDIGREE.ba_CV                      = PEDIGREE_ba;                          % (CV); (vertical vector: num_AggGrps X 1); NOTE: CV of ba production rate
PEDIGREE.em_CV                      = PEDIGREE_em;                          % (CV); (vertical vector: num_AggGrps X 1); NOTE: CV of em production rate

PEDIGREE.CONSUMPTION_CV             = PEDIGREE_CONSUMPTION;                 % (CV); (2D matrix: num_AggGrps X num_AggGrps)
PEDIGREE.consumption_total_CV       = PEDIGREE_consumption_total';      	% (CV); (vertical vector: num_AggGrps X 1); NOTE transpose
PEDIGREE.consumption_domestic_CV	= PEDIGREE_consumption_domestic';    	% (CV); (vertical vector: num_AggGrps X 1); NOTE transpose
PEDIGREE.consumption_import_CV      = PEDIGREE_consumption_import';         % (CV); (vertical vector: num_AggGrps X 1); NOTE transpose

PEDIGREE.DIET_CV                    = PEDIGREE_DIET;                        % (CV); (2D matrix: num_AggGrps X num_AggGrps)
PEDIGREE.diet_import_CV             = PEDIGREE_diet_import;                 % (CV); (horizontal vector: 1 X num_AggGrps)

PEDIGREE.LANDINGS_CV                = PEDIGREE_LANDINGS;                    % (CV); (2D matrix: num_AggGrps X num_agg_fleets)
PEDIGREE.DISCARDS_CV                = PEDIGREE_DISCARDS;                    % (CV); (2D matrix: num_AggGrps X num_agg_fleets)
PEDIGREE.CATCH_CV                   = PEDIGREE_CATCH;                       % (CV); (2D matrix: num_AggGrps X num_agg_fleets)
PEDIGREE.landings_TotalByFleet_CV	= PEDIGREE_landings_TotalByFleet;       % (CV); (horizontal vector: 1 X num_agg_fleets)
PEDIGREE.discards_TotalByFleet_CV	= PEDIGREE_discards_TotalByFleet;       % (CV); (horizontal vector: 1 X num_agg_fleets)
PEDIGREE.catch_TotalByFleet_CV      = PEDIGREE_catch_TotalByFleet;          % (CV); (horizontal vector: 1 X num_agg_fleets)
PEDIGREE.landings_TotalByGroup_CV	= PEDIGREE_landings_TotalByGroup;       % (CV); (vertical vector: num_AggGrps X 1)
PEDIGREE.discards_TotalByGroup_CV	= PEDIGREE_discards_TotalByGroup;       % (CV); (vertical vector: num_AggGrps X 1)
PEDIGREE.catch_TotalByGroup_CV      = PEDIGREE_catch_TotalByGroup;          % (CV); (vertical vector: num_AggGrps X 1)

PEDIGREE.EggProduction_forced_CV    = PEDIGREE_EggProduction_forced;        % (CV); (vertical vector: num_AggGrps X 1)
PEDIGREE.Metabolism_forced_CV       = PEDIGREE_Metabolism_forced;           % (CV); (vertical vector: num_AggGrps X 1)

PEDIGREE.fate_feces_CV              = PEDIGREE_fate_feces';                 % (CV); (horizontal vector: 1 X num_AggGrps); NOTE transpose
PEDIGREE.fate_senescence_CV         = PEDIGREE_fate_senescence';            % (CV); (horizontal vector: 1 X num_AggGrps); NOTE transpose
PEDIGREE.fate_metabolism_CV       	= PEDIGREE_fate_metabolism';            % (CV); (horizontal vector: 1 X num_AggGrps); NOTE transpose

PEDIGREE.MinBiomass                 = agg_PEDIGREE_MinBiomass;              % (t/km2); (vertical vector: num_AggGrps X 1)
PEDIGREE.MinBA                      = agg_PEDIGREE_MinBA;                   % (t/km2/y); (vertical vector: num_AggGrps X 1)
PEDIGREE.MinEM                      = agg_PEDIGREE_MinEM;                   % (t/km2/y); (vertical vector: num_AggGrps X 1)
PEDIGREE.MinEggProduction_forced	= agg_PEDIGREE_MinEggProduction_forced;	% (t/km2); (vertical vector: num_AggGrps X 1)
PEDIGREE.MinMetabolism_forced       = agg_PEDIGREE_MinMetabolism_forced;   	% (t/km2); (vertical vector: num_AggGrps X 1)
PEDIGREE.MaxBiomass                 = agg_PEDIGREE_MaxBiomass;              % (t/km2); (vertical vector: num_AggGrps X 1)
PEDIGREE.MaxBA                      = agg_PEDIGREE_MaxBA;                   % (t/km2/y); (vertical vector: num_AggGrps X 1)
PEDIGREE.MaxEM                      = agg_PEDIGREE_MaxEM;                   % (t/km2/y); (vertical vector: num_AggGrps X 1)
PEDIGREE.MaxEggProduction_forced	= agg_PEDIGREE_MaxEggProduction_forced;	% (t/km2); (vertical vector: num_AggGrps X 1)
PEDIGREE.MaxMetabolism_forced       = agg_PEDIGREE_MaxMetabolism_forced;   	% (t/km2); (vertical vector: num_AggGrps X 1)

PEDIGREE.MaxEE                      = agg_PEDIGREE_MaxEE;                   % (dimensionless); (vertical vector: num_AggGrps X 1)
PEDIGREE.MinPB                      = agg_PEDIGREE_MinPB;                   % (1/y); (vertical vector: num_AggGrps X 1)
PEDIGREE.MinQB                      = agg_PEDIGREE_MinQB;                   % (1/y); (vertical vector: num_AggGrps X 1)
PEDIGREE.MinPQ                      = agg_PEDIGREE_MinPQ;                   % (fraction: 0-1); (vertical vector: num_AggGrps X 1)
PEDIGREE.MinAE                      = agg_PEDIGREE_MinAE;                   % (fraction: 0-1); (vertical vector: num_AggGrps X 1)
PEDIGREE.MaxPB                      = agg_PEDIGREE_MaxPB;                   % (1/y); (vertical vector: num_AggGrps X 1)
PEDIGREE.MaxQB                      = agg_PEDIGREE_MaxQB;                   % (1/y); (vertical vector: num_AggGrps X 1)
PEDIGREE.MaxPQ                      = agg_PEDIGREE_MaxPQ;                   % (fraction: 0-1); (vertical vector: num_AggGrps X 1)
PEDIGREE.MaxAE                      = agg_PEDIGREE_MaxAE;                   % (fraction: 0-1); (vertical vector: num_AggGrps X 1)

PEDIGREE.PreyGuild_name          	= dat.PreyGuild_name;
PEDIGREE.PreyGuild_code          	= agg_PreyGuild_code;
% *************************************************************************


% end m-file***************************************************************