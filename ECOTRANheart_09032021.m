function [ECOTRAN] = ECOTRANheart_09032021(EwEResult, MonteCarloStore)
% convert EwE results into a bottom-up (donor-driven), bioenergetic budget (= consumption budget)
% by Jim Ruzicka
%
% takes:
%       EwEResult           the base EwE model as prepared by f_AggregateBiologicalModel_11182019.m
%       MonteCarloStore     retained Monte Carlo models
%
% calls:
%       f_ECOfunction_09032021                      returns a single ECOTRAN model for 1 "type" EwE model or 1 MonteCarlo EwE model
%           f_RedistributeCannibalism_11202019      remove cannibalism terms on diagonal of matrix EwE_diet
%           f_calcEE_12292020                       calculate Ecotrophic Efficiency 
%           f_calcPredationBudget_12102019      	for each producer, p, and consumer, c: ((b_pc * Q_c) / M2_p)) = the fraction of total predation on each producer, p, going to each consumer, c; (2D matrix: num_grps X num_grps; consumers X producers)
%
% returns:
%   ECOTRAN.EnergyBudget_matrix
%       producer groups (p) run across columns
%       rows:
%           metabolism term (2 rows, pelagic & benthic excretion)
%           predation by individual predators (n rows)
%           fisheries (n rows)
%           terminal pelagic detritus (1 row)
%           terminal benthic detritus (1 row)
%           BA = biomass accumulation (1 row)
%           EM = emigration (1 row) (immigration is negative)
%       NOTE: each column sums to 1 except for difference due to BA & EM (fishery clms do not sum to 1; landings are removed from system as EM); Terminal Benthic Detritus column sums to 1 BUT does NOT include OtherMortality nor EM
%       NOTE: BA & EM are handled in BioenergeticBudget_ProductionDetail, NOT in EnergyBudget_matrix
%   ECOTRAN.BioenergeticBudget (sum of each clm = 1, except for fisheries)
%           row 1 = BioenergeticBudget_production
%           row 2 = BioenergeticBudget_feces
%           row 3 = BioenergeticBudget_metabolism
%   ECOTRAN.BioenergeticBudget_ProductionDetail (sum of each clm = ECOTRAN.BioenergeticBudget_production
%           row 1 = BioenergeticBudget_SumPredation (total consumption of each group p going to all its consumers)
%           row 2 = BioenergeticBudget_eggs (total consumption of each group p going to eggs, gametes, or live births)
%           row 3 = BioenergeticBudget_OtherMortality (total consumption of each group p going to 
%           row 4 = BioenergeticBudget_BA (total consumption of each group p going to biomass accumulation)
%           row 5 = BioenergeticBudget_EM (total consumption of each group p going to emigration)
%   ECOTRAN.BioenergeticBudget_OtherMortDetail
%           row 1 = flow to terminal pelagic detritus
%           row 2 = flow to terminal benthic detritus
%   ECOTRAN.BioenergeticBudget_FecesDetail
%           row 1 = flow to terminal pelagic detritus
%           row 2 = flow to terminal benthic detritus
%   ECOTRAN.ProductionBudget (each clm sums to 1; negative values are allowed)
%           row 1 = ProductionBudget_predation (total production of each group p going to the sum of all predation pressure on each group p)
%           row 2 = ProductionBudget_eggs (total production of each group p going to eggs, gametes, or live births)
%           row 3 = ProductionBudget_OtherMortality
%           row 4 = ProductionBudget_BA
%           row 5 = ProductionBudget_EM
%   ECOTRAN.ProductionBudget_OtherMortDetail
%           row 1 = flow to terminal pelagic detritus
%           row 2 = flow to terminal benthic detritus
%
% revision date: 9-3-2021
%       12/29/2020 making adjustments to ConsumptionBudget & EnergyBudget to account for cannibalism
%       2/6/2021 setting detritus & egg pb, qb, pq values = 1
%       9/3/2021 corrected errors in cannibalism correction within f_ECOfunction_ and tested output


% *************************************************************************
% STEP 1: compile ModelDefinitions-----------------------------------------
%         NOTE: this info to be used by f_ECOfunction and is identical for all Monte Carlo models
fname_ECOTRANheart           = mfilename; % save name of this m-file to keep in saved model results
display(['Running: ' fname_ECOTRANheart])
ECOTRAN.fname_ECOTRANheart	= fname_ECOTRANheart;


% step 1a: GroupType ------------------------------------------------------
ModelDefinitions.GroupType                      = EwEResult.GroupType;                  % includes nutrients & fleets; (vertical vector: num_AggGrps X 1)


% step 1b: group type definitions -----------------------------------------
ModelDefinitions.GroupTypeDef_ANYNitroNutr   	= EwEResult.GroupTypeDef_ANYNitroNutr;
ModelDefinitions.GroupTypeDef_NO3             	= EwEResult.GroupTypeDef_NO3;
ModelDefinitions.GroupTypeDef_plgcNH4         	= EwEResult.GroupTypeDef_plgcNH4;
ModelDefinitions.GroupTypeDef_bnthNH4          	= EwEResult.GroupTypeDef_bnthNH4;
ModelDefinitions.GroupTypeDef_ANYPrimaryProd   	= EwEResult.GroupTypeDef_ANYPrimaryProd;
ModelDefinitions.GroupTypeDef_ANYConsumer    	= EwEResult.GroupTypeDef_ANYConsumer;
ModelDefinitions.GroupTypeDef_eggs           	= EwEResult.GroupTypeDef_eggs;
ModelDefinitions.GroupTypeDef_ANYDetritus     	= EwEResult.GroupTypeDef_ANYDetritus;
ModelDefinitions.GroupTypeDef_terminalPlgcDetr	= EwEResult.GroupTypeDef_terminalPlgcDetr;
ModelDefinitions.GroupTypeDef_terminalBnthDetr	= EwEResult.GroupTypeDef_terminalBnthDetr;
ModelDefinitions.GroupTypeDef_fleet           	= EwEResult.GroupTypeDef_fleet;
ModelDefinitions.GroupTypeDef_micrograzers    	= EwEResult.GroupTypeDef_micrograzers;
ModelDefinitions.GroupTypeDef_bacteria    	    = EwEResult.GroupTypeDef_bacteria;


% step 1c: numbers of group types -----------------------------------------
ModelDefinitions.num_grps                       = EwEResult.num_grps;
ModelDefinitions.num_NO3                        = EwEResult.num_NO3;
ModelDefinitions.num_plgcNH4                    = EwEResult.num_plgcNH4;
ModelDefinitions.num_bnthNH4                    = EwEResult.num_bnthNH4;
ModelDefinitions.num_NH4                        = EwEResult.num_NH4;
ModelDefinitions.num_nutrients                  = EwEResult.num_nutrients;
ModelDefinitions.num_eggs                     	= EwEResult.num_eggs;
ModelDefinitions.num_ANYdetritus               	= EwEResult.num_ANYdetritus;
ModelDefinitions.num_terminalANYdetritus        = EwEResult.num_terminalANYdetritus;
ModelDefinitions.num_livingANDfleets            = EwEResult.num_livingANDfleets;


% step 1d: nutrient fates -------------------------------------------------
ModelDefinitions.Oxidation_NH4               	= EwEResult.Oxidation_NH4;           	% fraction of NH4 produced oxidized directly back to NO3 abiologically; this should take precedence over phytoplankton uptake of NH4 in code below; (vertical vector)
ModelDefinitions.PhytoUptake_NO3             	= EwEResult.PhytoUptake_NO3;         	% row 1 = fraction of NO3 to small phytoplankton; row 2 = fraction of NO3 to large phytoplankton; (vertical vector)
ModelDefinitions.PhytoUptake_NH4              	= EwEResult.PhytoUptake_NH4;            % row 1 = fraction of NH4 to small phytoplankton; row 2 = fraction of NH4 to large phytoplankton; (vertical vector)


% step 1e: implicit bacterial metabolism of detritus ----------------------
%          NOTE: applies to TERMINAL detritus pools
ModelDefinitions.PelagicBacterialReduction    	= EwEResult.PelagicBacterialReduction;	% fraction of pelagic detritus returned to DIN (scaler)
ModelDefinitions.BenthicBacterialReduction     	= EwEResult.BenthicBacterialReduction;	% fraction of benthic detritus returned to DIN (assumed no direct return to surface mixed layer) (scaler)


% step 1f: metabolism, egg, feces, & senescence fates ---------------------
ModelDefinitions.fate_eggs                      = EwEResult.fate_eggs;                	% EwE egg fate; (proportions); (2D matrix: num_grps X num_eggs); NOTE: does not sum to 1 (eggs + detritus sum to 1)
ModelDefinitions.fate_EwEdetritus           	= EwEResult.fate_EwEdetritus;           % EwE detritus fate (NOT the ECOTRAN parameters); (2D matrix: num_grps X num_ANYdetritus); NOTE: does not include egg columns
ModelDefinitions.fate_metabolism              	= EwEResult.fate_metabolism;            % ECOTRAN metabolism fate; (proportions); (2D matrix: num_NH4 X num_grps)
ModelDefinitions.fate_feces                     = EwEResult.fate_feces;                 % ECOTRAN feces detritus fate; (proportions); (2D matrix: num_terminalANYdetritus X num_grps)
ModelDefinitions.fate_senescence                = EwEResult.fate_senescence;            % ECOTRAN senescence detritus fate; (proportions); (2D matrix: num_terminalANYdetritus X num_grps)


% step 1g: size terms needed to initialize ECOTRAN variables --------------
num_grps                                        = EwEResult.num_grps;
num_nutrients                                   = EwEResult.num_nutrients;
num_eggs                                        = EwEResult.num_eggs;
num_ANYdetritus                                 = EwEResult.num_ANYdetritus;
num_livingANDfleets                             = EwEResult.num_livingANDfleets;
% *************************************************************************





% *************************************************************************
% STEP 2: prep MonteCarloStore structure variable from .mat----------------

% step 2a: either use Monte Carlo set OR the single "type" model ----------
if isempty(MonteCarloStore) % if Monte Carlo models are NOT used, use the single "type" model

    num_MC                          = 1; % use the single "type" model in EwEResult
    
    % parameters -----
    biomass                         = EwEResult.biomass;                % (vertical vector: num_grps X 1)
    production                      = EwEResult.production;             % (vertical vector: num_grps X 1)
    pb                              = EwEResult.pb;                     % (1/y); (vertical vector: num_grps X 1)
    qb                              = EwEResult.qb;                     % (vertical vector: num_grps X 1)
    pq                              = EwEResult.pq;                     % (vertical vector: num_grps X 1)
    ae                              = EwEResult.ae;                     % (vertical vector: num_grps X 1)
    TL                              = EwEResult.TL;                     % (vertical vector: num_grps X 1)
    ee                              = EwEResult.ee;                     % (vertical vector: num_grps X 1)
    ba                              = EwEResult.ba;                     % (t/km2/y); (vertical vector: num_grps X 1); NOTE: ba production rate
    em                              = EwEResult.em;                     % (t/km2/y); (vertical vector: num_grps X 1); NOTE: em production rate
    
    % diet & consumption -----
    CONSUMPTION                     = EwEResult.CONSUMPTION;            % (t/km2/y); (2D matrix: num_grps X num_grps)
    consumption_domestic            = EwEResult.consumption_domestic;	% (t/km2/y); (vertical vector: num_grps X 1)
    consumption_import              = EwEResult.consumption_import;     % (t/km2/y); (vertical vector: num_grps X 1)
    consumption_total               = EwEResult.consumption_total;      % (t/km2/y); (vertical vector: num_grps X 1)
    DIET                            = EwEResult.DIET;                   % (2D matrix: num_grps X num_grps); NOTE: at this point DIET columns are normalized to 1 and do not include import diet
    import_diet                     = EwEResult.import_diet;            % (horizontal vector: 1 X num_grps)

    % fleet info -----
    LANDINGS                        = EwEResult.LANDINGS;              	% (t/km2/y); (2D matrix: num_grps X num_fleets)
    DISCARDS                        = EwEResult.DISCARDS;            	% (t/km2/y); (2D matrix: num_grps X num_fleets)
    CATCH                           = EwEResult.CATCH;                	% (t/km2/y); (2D matrix: num_grps X num_fleets)
    landings_TotalByFleet           = EwEResult.landings_TotalByFleet;	% (t/km2/y); (horizontal vector: 1 X num_fleets)
    discards_TotalByFleet           = EwEResult.discards_TotalByFleet;	% (t/km2/y); (horizontal vector: 1 X num_fleets)
    catch_TotalByFleet              = EwEResult.catch_TotalByFleet;    	% (t/km2/y); (horizontal vector: 1 X num_fleets)
    landings_TotalByGroup           = EwEResult.landings_TotalByGroup;	% (t/km2/y); (vertical vector: num_fleets X 1)
    discards_TotalByGroup           = EwEResult.discards_TotalByGroup;	% (t/km2/y); (vertical vector: num_fleets X 1)
    catch_TotalByGroup              = EwEResult.catch_TotalByGroup;    	% (t/km2/y); (vertical vector: num_fleets X 1)
    LandingsFraction                = LANDINGS ./ CATCH;                % fraction of catch landed for each group; (2D matrix: num_grps X num_fleets)
    DiscardFraction                 = DISCARDS ./ CATCH;                % fraction of catch discarded for each group; (2D matrix: num_grps X num_fleets)
    LandingsFraction(isnan(LandingsFraction))	= 0;                    % correct NaNs from div/0 errors; (2D matrix: num_grps X num_fleets)
    DiscardFraction(isnan(DiscardFraction))     = 0;                    % correct NaNs from div/0 errors; (2D matrix: num_grps X num_fleets)
    
    % % forced metabolism & egg rates -----
    % EggProduction_forced            = EwEResult.EggProduction_forced;	% (t/km2/y); (vertical vector: num_grps X 1); FFF for future
    % Metabolism_forced               = EwEResult.Metabolism_forced;   	% (t/km2/y); (vertical vector: num_grps X 1); FFF for future
    
elseif ~isempty(MonteCarloStore) % if Monte Carlo models ARE used, use MonteCarloStore

    [toss, num_MC]                  = size(MonteCarloStore.biomass);

    % parameters -----
    biomass                         = MonteCarloStore.biomass;                  % (2D matrix: num_grps X num_MC)
    production                      = MonteCarloStore.production;               % (2D matrix: num_grps X num_MC)
    pb                              = MonteCarloStore.pb;                       % (1/y); (2D matrix: num_grps X num_MC)
    qb                              = MonteCarloStore.qb;                       % (2D matrix: num_grps X num_MC)
    pq                              = MonteCarloStore.pq;                       % (2D matrix: num_grps X num_MC)
    ae                              = MonteCarloStore.ae;                       % (2D matrix: num_grps X num_MC)
    TL                              = MonteCarloStore.TL;                       % (2D matrix: num_grps X num_MC)
    ee                              = MonteCarloStore.ee;                       % (2D matrix: num_grps X num_MC)
    ba                              = MonteCarloStore.ba;                       % (t/km2/y); (2D matrix: num_grps X num_MC); NOTE: ba production rate
    em                              = MonteCarloStore.em;                       % (t/km2/y); (2D matrix: num_grps X num_MC); NOTE: em production rate
    
    % diet & consumption -----
    CONSUMPTION                     = MonteCarloStore.CONSUMPTION;              % (t/km2/y); (3D matrix: num_grps X num_grps X num_MC)
    consumption_domestic            = MonteCarloStore.consumption_domestic;     % (t/km2/y); (2D matrix: num_grps X num_MC)
    consumption_import              = MonteCarloStore.consumption_import;       % (t/km2/y); (2D matrix: num_grps X num_MC)
    consumption_total               = MonteCarloStore.consumption_total;        % (t/km2/y); (2D matrix: num_grps X num_MC)
    DIET                            = MonteCarloStore.DIET;                     % (3D matrix: num_grps X num_grps X num_MC)
    import_diet                     = MonteCarloStore.import_diet;              % (3D matrix: 1 X num_grps X num_MC)

    % fleet info -----
    LANDINGS                        = MonteCarloStore.LANDINGS;              	% (t/km2/y); (3D matrix: num_grps X num_fleets X num_MC)
    DISCARDS                        = MonteCarloStore.DISCARDS;                 % (t/km2/y); (3D matrix: num_grps X num_fleets X num_MC)
    CATCH                           = MonteCarloStore.CATCH;                	% (t/km2/y); (3D matrix: num_grps X num_fleets X num_MC)
    landings_TotalByFleet           = MonteCarloStore.landings_TotalByFleet;	% (t/km2/y); (3D matrix: 1 X num_fleets X num_MC)
    discards_TotalByFleet           = MonteCarloStore.discards_TotalByFleet;	% (t/km2/y); (3D matrix: 1 X num_fleets X num_MC)
    catch_TotalByFleet              = MonteCarloStore.catch_TotalByFleet;    	% (t/km2/y); (3D matrix: 1 X num_fleets X num_MC)
    landings_TotalByGroup           = MonteCarloStore.landings_TotalByGroup;	% (t/km2/y); (2D matrix: num_fleets X num_MC)
    discards_TotalByGroup           = MonteCarloStore.discards_TotalByGroup;	% (t/km2/y); (2D matrix: num_fleets X num_MC)
    catch_TotalByGroup              = MonteCarloStore.catch_TotalByGroup;    	% (t/km2/y); (2D matrix: num_fleets X num_MC)
    LandingsFraction                             	= LANDINGS ./ CATCH;        % fraction of catch landed for each group; (3D matrix: num_grps X num_fleets X num_MC)
    DiscardFraction                             	= DISCARDS ./ CATCH;        % fraction of catch discarded for each group; (3D matrix: num_grps X num_fleets X num_MC)
    LandingsFraction(isnan(LandingsFraction_type))	= 0;                        % correct NaNs from div/0 errors; (3D matrix: num_grps X num_fleets X num_MC)
    DiscardFraction(isnan(DiscardFraction_type))	= 0;                        % correct NaNs from div/0 errors; (3D matrix: num_grps X num_fleets X num_MC)
    
    % % forced metabolism & egg rates -----
    % EggProduction_forced            = MonteCarloStore.EggProduction_forced;     % (t/km2/y); (2D matrix: num_grps X num_MC); FFF for future
    % Metabolism_forced               = MonteCarloStore.Metabolism_forced;        % (t/km2/y); (2D matrix: num_grps X num_MC); FFF for future
    
end
% *************************************************************************





% *************************************************************************
% STEP 3: define OR initialize ECOTRAN results structure-------------------

% step 3a: group type definitions -----------------------------------------
ECOTRAN.GroupTypeDef_ANYNitroNutr       = EwEResult.GroupTypeDef_ANYNitroNutr;
ECOTRAN.GroupTypeDef_NO3                = EwEResult.GroupTypeDef_NO3;
ECOTRAN.GroupTypeDef_plgcNH4            = EwEResult.GroupTypeDef_plgcNH4;
ECOTRAN.GroupTypeDef_bnthNH4            = EwEResult.GroupTypeDef_bnthNH4;
ECOTRAN.GroupTypeDef_ANYPrimaryProd     = EwEResult.GroupTypeDef_ANYPrimaryProd;
ECOTRAN.GroupTypeDef_LrgPhyto           = EwEResult.GroupTypeDef_LrgPhyto;
ECOTRAN.GroupTypeDef_SmlPhyto           = EwEResult.GroupTypeDef_SmlPhyto;
ECOTRAN.GroupTypeDef_Macrophytes        = EwEResult.GroupTypeDef_Macrophytes;
ECOTRAN.GroupTypeDef_ANYConsumer        = EwEResult.GroupTypeDef_ANYConsumer;
ECOTRAN.GroupTypeDef_ConsumPlgcPlankton = EwEResult.GroupTypeDef_ConsumPlgcPlankton;
ECOTRAN.GroupTypeDef_ConsumPlgcNekton   = EwEResult.GroupTypeDef_ConsumPlgcNekton;
ECOTRAN.GroupTypeDef_ConsumPlgcWrmBlood	= EwEResult.GroupTypeDef_ConsumPlgcWrmBlood;
ECOTRAN.GroupTypeDef_ConsumBntcInvert   = EwEResult.GroupTypeDef_ConsumBntcInvert;
ECOTRAN.GroupTypeDef_ConsumBntcVert     = EwEResult.GroupTypeDef_ConsumBntcVert;
ECOTRAN.GroupTypeDef_ConsumBnthWrmBlood = EwEResult.GroupTypeDef_ConsumBnthWrmBlood;
ECOTRAN.GroupTypeDef_fleet              = EwEResult.GroupTypeDef_fleet;
ECOTRAN.GroupTypeDef_eggs               = EwEResult.GroupTypeDef_eggs;
ECOTRAN.GroupTypeDef_offal              = EwEResult.GroupTypeDef_offal;
ECOTRAN.GroupTypeDef_terminalPlgcDetr   = EwEResult.GroupTypeDef_terminalPlgcDetr;
ECOTRAN.GroupTypeDef_terminalBnthDetr   = EwEResult.GroupTypeDef_terminalBnthDetr;
ECOTRAN.GroupTypeDef_ANYDetritus        = EwEResult.GroupTypeDef_ANYDetritus;
ECOTRAN.GroupTypeDef_micrograzers       = EwEResult.GroupTypeDef_micrograzers;
ECOTRAN.GroupTypeDef_bacteria           = EwEResult.GroupTypeDef_bacteria;
ECOTRAN.GroupTypeDef_ba                 = EwEResult.GroupTypeDef_ba;
ECOTRAN.GroupTypeDef_em                 = EwEResult.GroupTypeDef_em;
ECOTRAN.GroupTypeDef_import             = EwEResult.GroupTypeDef_import;


% step 3b: numbers of group types -----------------------------------------
ECOTRAN.num_grps                        = EwEResult.num_grps;                  % number of aggregated groups, including nutrients & fleets
ECOTRAN.num_MC                          = num_MC;
ECOTRAN.num_NO3                         = EwEResult.num_NO3;
ECOTRAN.num_plgcNH4                     = EwEResult.num_plgcNH4;
ECOTRAN.num_bnthNH4                     = EwEResult.num_bnthNH4;
ECOTRAN.num_NH4                         = EwEResult.num_NH4;
ECOTRAN.num_nutrients                   = EwEResult.num_nutrients;
ECOTRAN.num_PrimaryProducers            = EwEResult.num_PrimaryProducers;
ECOTRAN.num_consumers                   = EwEResult.num_consumers;
ECOTRAN.num_micrograzers                = EwEResult.num_micrograzers;
ECOTRAN.num_bacteria                    = EwEResult.num_bacteria;
ECOTRAN.num_eggs                        = EwEResult.num_eggs;
ECOTRAN.num_ANYdetritus                 = EwEResult.num_ANYdetritus;            % ANY detritus groups (w/o eggs)
ECOTRAN.num_terminalPLGCdetritus        = EwEResult.num_terminalPLGCdetritus;
ECOTRAN.num_terminalBNTHdetritus        = EwEResult.num_terminalBNTHdetritus;
ECOTRAN.num_terminalANYdetritus         = EwEResult.num_terminalANYdetritus;
ECOTRAN.num_eggsANDdetritus             = EwEResult.num_eggsANDdetritus;
ECOTRAN.num_livingANDdetritus           = EwEResult.num_livingANDdetritus;
ECOTRAN.num_fleets                      = EwEResult.num_fleets;
ECOTRAN.num_livingANDfleets             = EwEResult.num_livingANDfleets;
ECOTRAN.num_NONnutrients                = EwEResult.num_NONnutrients;


% step 3c: labels, GroupType, & number codes ------------------------------
ECOTRAN.label                           = EwEResult.label;      	% text cells; (vertical vector: num_AggGrps X 1)
ECOTRAN.GroupType                       = EwEResult.GroupType;      % includes nutrients & fleets; (vertical vector: num_AggGrps X 1)
ECOTRAN.CodeNumber                      = EwEResult.CodeNumber;     % (vertical vector: num_AggGrps X 1)


% step 3d: parameters -----------------------------------------------------
ECOTRAN.biomass                         = biomass;                  % (2D matrix: num_grps X num_MC)
ECOTRAN.production                      = production;               % (2D matrix: num_grps X num_MC)
ECOTRAN.pb                              = pb;                       % (1/y); (2D matrix: num_grps X num_MC)
ECOTRAN.qb                              = qb;                       % (2D matrix: num_grps X num_MC)
ECOTRAN.pq                              = pq;                       % (2D matrix: num_grps X num_MC); NOTE: will be redefined by f_ECOfunction
ECOTRAN.ae                              = ae;                       % (2D matrix: num_grps X num_MC); NOTE: will be redefined by f_ECOfunction
ECOTRAN.TL                              = TL;                       % (2D matrix: num_grps X num_MC)
ECOTRAN.ba                              = ba;                       % (t/km2/y); (2D matrix: num_grps X num_MC); NOTE: ba production rate
ECOTRAN.em                              = em;                       % (t/km2/y); (2D matrix: num_grps X num_MC); NOTE: em production rate


% step 3e: diet & consumption ---------------------------------------------
ECOTRAN.CONSUMPTION                     = CONSUMPTION;              % (t/km2/y); (3D matrix: num_grps X num_grps X num_MC)
ECOTRAN.consumption_domestic            = consumption_domestic;     % (t/km2/y); (2D matrix: num_grps X num_MC)
ECOTRAN.consumption_import              = consumption_import;       % (t/km2/y); (2D matrix: num_grps X num_MC)
ECOTRAN.consumption_total               = consumption_total;        % (t/km2/y); (2D matrix: num_grps X num_MC)
ECOTRAN.DIET                            = DIET;                     % (3D matrix: num_grps X num_grps X num_MC)
ECOTRAN.import_diet                     = import_diet;              % (3D matrix: 1 X num_grps X num_MC)


% step 3f: fleet info -----------------------------------------------------
ECOTRAN.LANDINGS                        = LANDINGS;              	% (t/km2/y); (3D matrix: num_grps X num_fleets X num_MC)
ECOTRAN.DISCARDS                        = DISCARDS;                 % (t/km2/y); (3D matrix: num_grps X num_fleets X num_MC)
ECOTRAN.CATCH                           = CATCH;                	% (t/km2/y); (3D matrix: num_grps X num_fleets X num_MC)
ECOTRAN.landings_TotalByFleet           = landings_TotalByFleet;	% (t/km2/y); (3D matrix: 1 X num_fleets X num_MC)
ECOTRAN.discards_TotalByFleet           = discards_TotalByFleet;	% (t/km2/y); (3D matrix: 1 X num_fleets X num_MC)
ECOTRAN.catch_TotalByFleet              = catch_TotalByFleet;    	% (t/km2/y); (3D matrix: 1 X num_fleets X num_MC)
ECOTRAN.landings_TotalByGroup           = landings_TotalByGroup;	% (t/km2/y); (2D matrix: num_fleets X num_MC)
ECOTRAN.discards_TotalByGroup           = discards_TotalByGroup;	% (t/km2/y); (2D matrix: num_fleets X num_MC)
ECOTRAN.catch_TotalByGroup              = catch_TotalByGroup;    	% (t/km2/y); (2D matrix: num_fleets X num_MC)
ECOTRAN.LandingsFraction             	= LandingsFraction;         % fraction of catch landed for each group; (3D matrix: num_grps X num_fleets X num_MC)
ECOTRAN.DiscardFraction              	= DiscardFraction;          % fraction of catch discarded for each group; (3D matrix: num_grps X num_fleets X num_MC)


% step 3g: metabolism, egg, feces, & senescence fates ---------------------
%          NOTE: each MonteCarlo layer should be identical QQQ confirm they are identical and cut to 1 layer
ECOTRAN.fate_metabolism                 = zeros(num_nutrients,       num_grps, num_MC);	% (3D matrix: num_nutrients X num_grps X num_MC)
ECOTRAN.fate_eggs                       = zeros(num_eggs,            num_grps, num_MC);	% (3D matrix: num_eggs X num_grps X num_MC)
ECOTRAN.fate_feces                      = zeros(num_ANYdetritus,     num_grps, num_MC);	% (3D matrix: num_ANYdetritus X num_grps X num_MC)
ECOTRAN.fate_senescence                 = zeros(num_ANYdetritus,     num_grps, num_MC);	% (3D matrix: num_ANYdetritus X num_grps X num_MC)
ECOTRAN.fate_predation                  = zeros(num_livingANDfleets, num_grps, num_MC); % (3D matrix: num_livingANDfleets X num_grps X num_MC)


% % step 3h: forced metabolism & egg rates ----------------------------------
% ECOTRAN.EggProduction_forced            = EggProduction_forced;     % (t/km2/y); (2D matrix: num_grps X num_MC); FFF for future
% ECOTRAN.Metabolism_forced               = Metabolism_forced;        % (t/km2/y); (2D matrix: num_grps X num_MC); FFF for future


% step 3i: retention & production loss scalers ----------------------------
ECOTRAN.ProductionLossScaler            = EwEResult.ProductionLossScaler;       % (scaler: 0 to 1); (vertical vector: num_AggGrps X 1)
ECOTRAN.RetentionScaler                 = EwEResult.RetentionScaler;            % (scaler: 0 to 1); (vertical vector: num_AggGrps X 1)


% step 3j: functional response terms -------------------------------------
ECOTRAN.FunctionalResponseParams        = EwEResult.FunctionalResponseParams;   % (vertical vector: num_AggGrps X 4)
ECOTRAN.FunctionalResponse_matrix       = EwEResult.FunctionalResponse_matrix;	% (2D matrix: num_AggGrps X num_AggGrps); FFF for future


% step 3k: nutrient fates & detritus cycling terms ------------------------
ECOTRAN.Oxidation_NH4                   = EwEResult.Oxidation_NH4;
ECOTRAN.PhytoUptake_NH4                 = EwEResult.PhytoUptake_NH4;
ECOTRAN.PhytoUptake_NO3                 = EwEResult.PhytoUptake_NO3;
ECOTRAN.PelagicBacterialReduction       = EwEResult.PelagicBacterialReduction;
ECOTRAN.BenthicBacterialReduction       = EwEResult.BenthicBacterialReduction;


% step 3l: BUDGET MATRICES ------------------------------------------------
ECOTRAN.BioenergeticBudget              = zeros(3, num_grps, num_MC);           % (3D matrix: 3 X num_grps X num_MC)
ECOTRAN.ProductionBudget                = zeros(5, num_grps, num_MC);           % (3D matrix: 5 X num_grps X num_MC)
ECOTRAN.ConsumptionBudget               = zeros(7, num_grps, num_MC);        	% (3D matrix: 7 X num_grps X num_MC)
ECOTRAN.EnergyBudget                    = zeros(num_grps, num_grps, num_MC);	% (3D matrix: num_grps X num_grps X num_MC)


% step 3m: ecotrophic efficiency terms ------------------------------------
ECOTRAN.ee_eggs                         = zeros(num_grps, num_MC); % (2D matrix: num_grps X num_MC)
ECOTRAN.ee_predation                    = zeros(num_grps, num_MC); % (2D matrix: num_grps X num_MC)
ECOTRAN.ee_ba                           = zeros(num_grps, num_MC); % (2D matrix: num_grps X num_MC)
ECOTRAN.ee_em                           = zeros(num_grps, num_MC); % (2D matrix: num_grps X num_MC)
ECOTRAN.ee                              = zeros(num_grps, num_MC); % (2D matrix: num_grps X num_MC)
ECOTRAN.TransferEfficiency              = zeros(num_grps, num_MC); % (2D matrix: num_grps X num_MC)


% step 3n: diet, consumption, & predation terms with/without cannibalism --
%          NOTE: QQQ not certain these are needed to be passed along
ECOTRAN.DIET_NoCannibalism                  = zeros(num_grps, num_grps, num_MC);	% (3D matrix: num_grps X num_grps X num_MC)
ECOTRAN.diet_cannibalism                    = zeros(1, num_grps, num_MC);       	% (3D matrix: 1 X num_grps X num_MC)
ECOTRAN.consumption_domestic_cannibalism    = zeros(num_grps, num_MC);              % cannibalism rates; (t/km2/y); (2D matrix: num_grps X num_MC)
ECOTRAN.CONSUMPTION_NoCannibalism           = zeros(num_grps, num_grps, num_MC);	% (3D matrix: num_grps X num_grps X num_MC)
ECOTRAN.predation_total                     = zeros(num_grps, num_MC);              % total predation ON each producer group, p; includes cannibalism; (M2_p); (t/km2/yr); (2D matrix: num_grps X num_MC)
ECOTRAN.predation_total_NoCannibalism       = zeros(num_grps, num_MC);              % total predation ON each producer group, p; does NOT include cannibalism; (M2_p); (t/km2/yr); (2D matrix: num_grps X num_MC)
ECOTRAN.consumption_domestic_NOcannibalism  = zeros(num_grps, num_MC);              % total consumption BY each consumer c; does NOT include cannibalism; (t/km2/yr); (2D matrix: num_grps X num_MC)
ECOTRAN.consumption_domestic_cannibalism    = zeros(num_grps, num_MC);              % total cannibalism BY each consumer c; (t/km2/yr); (2D matrix: num_grps X num_MC)
% *************************************************************************





% *************************************************************************
% STEP 4: perform ECOTRAN conversion for each MonteCarlo model-------------
for layer_loop = 1:num_MC
    
    % step 4a: select current MonteCarlo model ----------------------------
    CurrentModel.production             = production(:, layer_loop);            % (t/km2/y); (vertical vector: num_grps X 1)
    CurrentModel.pb                     = pb(:, layer_loop);                    % (1/y); (vertical vector: num_grps X 1)
    CurrentModel.qb                     = qb(:, layer_loop);                    % (1/y); (vertical vector: num_grps X 1)
    CurrentModel.pq                     = pq(:, layer_loop);                    % (dimensionless); (vertical vector: num_grps X 1)
    CurrentModel.ae                     = ae(:, layer_loop);                    % (dimensionless); (vertical vector: num_grps X 1)
    CurrentModel.ba                     = ba(:, layer_loop);                    % (t/km2/y); (vertical vector: num_grps X 1); NOTE: still in absolute production rate units at this point
    CurrentModel.em                     = em(:, layer_loop);                    % (t/km2/y); (vertical vector: num_grps X 1); NOTE: still in absolute production rate units at this point
    CurrentModel.CONSUMPTION            = CONSUMPTION(:, :, layer_loop);
    CurrentModel.consumption_total      = consumption_total(:, layer_loop);     % (t/km2/y); (vertical vector: num_grps X 1)
    CurrentModel.DIET                   = DIET(:, :, layer_loop);
    % CurrentModel.EggProduction_forced	  = EggProduction_forced(:, layer_loop);	% (t/km2/y); (vertical vector: num_grps X 1); FFF for future
    
    
    % step 4b: run it through ECOTRAN -------------------------------------
    [SINGLE_ECO]                        = f_ECOfunction_09032021(ModelDefinitions, CurrentModel); % returns a single ECOTRAN model for 1 "type" EwE model or 1 MonteCarlo EwE model
    
    
    % step 4c: build-up set of ECOTRAN models -----------------------------
    % budget terms -----
    ECOTRAN.BioenergeticBudget(1:3, 1:num_grps, layer_loop)                 = SINGLE_ECO.BioenergeticBudget;	% (3D matrix: 3 X num_grps X num_MC)
    ECOTRAN.ProductionBudget(1:5, 1:num_grps, layer_loop)                   = SINGLE_ECO.ProductionBudget;      % (3D matrix: 5 X num_grps X num_MC)
    ECOTRAN.ConsumptionBudget(1:7, 1:num_grps, layer_loop)                  = SINGLE_ECO.ConsumptionBudget;     % (3D matrix: 7 X num_grps X num_MC)
    ECOTRAN.EnergyBudget(1:num_grps, 1:num_grps, layer_loop)                = SINGLE_ECO.EnergyBudget;          % (3D matrix: num_grps X num_grps X num_MC)
    
    % fate terms -----
    ECOTRAN.fate_metabolism(1:num_nutrients, 1:num_grps, layer_loop)        = SINGLE_ECO.fate_metabolism;       % (3D matrix: num_nutrients X num_grps X num_MC)
    ECOTRAN.fate_eggs(1:num_eggs, 1:num_grps, layer_loop)                   = SINGLE_ECO.fate_eggs;             % (3D matrix: num_eggs X num_grps X num_MC)
    ECOTRAN.fate_feces(1:num_ANYdetritus, 1:num_grps, layer_loop)           = SINGLE_ECO.fate_feces;            % (3D matrix: num_ANYdetritus X num_grps X num_MC)
    ECOTRAN.fate_senescence(1:num_ANYdetritus, 1:num_grps, layer_loop)      = SINGLE_ECO.fate_senescence;       % (3D matrix: num_ANYdetritus X num_grps X num_MC)
    ECOTRAN.fate_predation(1:num_livingANDfleets, 1:num_grps, layer_loop)	= SINGLE_ECO.fate_predation;        % (2D matrix: num_livingANDfleets X num_grps)
    
    % ecotrophic efficiency terms -----
    ECOTRAN.ee_eggs(1:num_grps, layer_loop)                                 = SINGLE_ECO.ee_eggs;               % (2D matrix: num_grps X num_MC)
    ECOTRAN.ee_predation(1:num_grps, layer_loop)                            = SINGLE_ECO.ee_predation;          % (2D matrix: num_grps X num_MC)
    ECOTRAN.ee_ba(1:num_grps, layer_loop)                                   = SINGLE_ECO.ee_ba;                 % (2D matrix: num_grps X num_MC)
    ECOTRAN.ee_em(1:num_grps, layer_loop)                                   = SINGLE_ECO.ee_em;                 % (2D matrix: num_grps X num_MC)
    ECOTRAN.ee(1:num_grps, layer_loop)                                      = SINGLE_ECO.ee;                    % (2D matrix: num_grps X num_MC)
    ECOTRAN.TransferEfficiency(1:num_grps, layer_loop)                      = SINGLE_ECO.TransferEfficiency;	% (2D matrix: num_grps X num_MC)

    % physiology parameters -----
    ECOTRAN.ae(1:num_grps, layer_loop)                                      = SINGLE_ECO.ae;                    % (2D matrix: num_grps X num_MC)
    ECOTRAN.pb(1:num_grps, layer_loop)                                      = SINGLE_ECO.pb;                    % (2D matrix: num_grps X num_MC)
    ECOTRAN.qb(1:num_grps, layer_loop)                                      = SINGLE_ECO.qb;                    % (2D matrix: num_grps X num_MC)
    ECOTRAN.pq(1:num_grps, layer_loop)                                      = SINGLE_ECO.pq;                    % (2D matrix: num_grps X num_MC)

    % diet & consumption terms -----
    %      NOTE: I don't think it is necessary to pass along these terms)
    ECOTRAN.DIET_NoCannibalism(1:num_grps, 1:num_grps, layer_loop)        	= SINGLE_ECO.DIET_NoCannibalism;                    % (3D matrix: num_grps X num_grps X num_MC)
    ECOTRAN.diet_cannibalism(1, 1:num_grps, layer_loop)                     = SINGLE_ECO.diet_cannibalism;                      % (3D matrix: 1 X num_grps X num_MC)
    ECOTRAN.consumption_domestic_cannibalism(1:num_grps, layer_loop)        = SINGLE_ECO.consumption_domestic_cannibalism;      % cannibalism rates; (t/km2/y); (2D matrix: num_grps X num_MC)
    ECOTRAN.CONSUMPTION_NoCannibalism(1:num_grps, 1:num_grps, layer_loop)	= SINGLE_ECO.CONSUMPTION_NoCannibalism;             % (3D matrix: num_grps X num_grps X num_MC)
    ECOTRAN.predation_total(1:num_grps, layer_loop)                         = SINGLE_ECO.predation_total;                       % total predation ON each producer group, p; includes cannibalism; (M2_p); (t/km2/yr); (2D matrix: num_grps X num_MC)
    ECOTRAN.predation_total_NoCannibalism(1:num_grps, layer_loop)           = SINGLE_ECO.predation_total_NoCannibalism;         % total predation ON each producer group, p; does NOT include cannibalism; (M2_p); (t/km2/yr); (2D matrix: num_grps X num_MC)
    ECOTRAN.consumption_domestic_NOcannibalism(1:num_grps, layer_loop)      = SINGLE_ECO.consumption_domestic_NOcannibalism;	% total consumption BY each consumer c; does NOT include cannibalism; (t/km2/yr); (2D matrix: num_grps X num_MC)
    ECOTRAN.consumption_domestic_cannibalism(1:num_grps, layer_loop)        = SINGLE_ECO.consumption_domestic_cannibalism;      % total cannibalism BY each consumer c; (t/km2/yr); (2D matrix: num_grps X num_MC)
    
end

ECOTRAN.fname_ECOTRANheart              = fname_ECOTRANheart;
ECOTRAN.fname_ECOfunction               = SINGLE_ECO.fname_ECOfunction;     % name of this version of f_ECOfunction
% ECOTRAN.fname_RedistributeCannibalism	  = fname_RedistributeCannibalism;	  % FFF this function does not yet pass its name along; name of this version of f_RedistributeCannibalism
ECOTRAN.fname_calcEE                    = SINGLE_ECO.fname_calcEE;          % file name of this f_calcEEsub-function
% ECOTRAN.fname_CalcPredationBudget       = fname_CalcPredationBudget;        % FFF this function does not yet pass its name along; file name of the sub-function f_CalcPredationBudget
% *************************************************************************


% end m-file***************************************************************