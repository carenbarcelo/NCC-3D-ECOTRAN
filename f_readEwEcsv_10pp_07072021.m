function [dat] = f_readEwEcsv_10pp_07072021(readFile)
% load EwE model from named .csv file and prepare data
% includes conversion from mortality matrix in K-path to consumption matrix
% NOTE: this code allows for 10 primary producers
% by Jim Ruzicka
% calls:
%	none
% takes:
%	readFile                            concatenated directory & file name
% returns:
% number codes & names
%         dat.EcotranNumberCode
%         dat.EcotranGroupType
%         dat.AggCode
%         dat.EwENumberCode
%         dat.EwEName
%         dat.EcotranName
%         dat.AggEwEName
% parameters (vectors)        
%         dat.Biomass
%         dat.Production
%         dat.PB
%         dat.QB
%         dat.EE
%         dat.PQ
%         dat.AE
%         dat.TL
%         dat.landings                  (t km^-2 y^-1)
%         dat.discards                  (t km^-2 y^-1)
%         dat.BA                        (t km^-2 y^-1)
%         dat.EM                        (t km^-2 y^-1)  (-) for immigration; (+) for emigration
% diet & consumption (matrices)
%         dat.diet                      (column 1 & row 1 are trophic group number headers)
%         dat.import_diet               (column 1 is header)
%         dat.consumption               (t km^-2 y^-1) 	(array)
%         dat.import_consumption
%         dat.total_consumption
%         dat.mortality
% ECOTRAN recycling definitions & production parameters
%         dat.DetritusImport
%         dat.EwE_DetritusFate
%         dat.EggProduction             (forced egg, gamete, or live-birth production --- not derived as EwE detritus)
%         dat.Metabolism                (forced metabolism --- not derived from EwE parameters)
%         dat.DetritusFate_feces
%         dat.DetritusFate_senescence
%         dat.ExcretionFate
%         dat.ProductionLossScaler
%         dat.RetentionScaler
%         dat.PelagicBacterialReduction
%         dat.BenthicBacterialReduction
%         dat.Oxidation_NH4
%         dat.PhytoUptake_NH4
%         dat.PhytoUptake_NO3
% numbers of various group types
%         dat.num_living
%         dat.num_detritus
%         dat.num_gear
%         dat.num_EwEGroups
%         dat.num_LivingAndDetritus
%         dat.num_NitrogenNutrients
%         dat.num_NonNitrogenNutrients
%         dat.num_EcotranGroups
%         dat.num_eggs
%         dat.num_EggsAndDetritus
%         dat.num_micrograzers
%         dat.num_bacteria
% group type code definitions
%         dat.GroupTypeDef_ANYNitroNutr
%         dat.GroupTypeDef_NO3
%         dat.GroupTypeDef_plgcNH4
%         dat.GroupTypeDef_bnthNH4
%         dat.GroupTypeDef_ANYPrimaryProd
%         dat.GroupTypeDef_LrgPhyto
%         dat.GroupTypeDef_SmlPhyto
%         dat.GroupTypeDef_Macrophytes
%         dat.GroupTypeDef_ANYConsumer
%         dat.GroupTypeDef_ConsumPlgcPlankton
%         dat.GroupTypeDef_ConsumPlgcNekton
%         dat.GroupTypeDef_ConsumPlgcWrmBlood
%         dat.GroupTypeDef_ConsumBntcInvert
%         dat.GroupTypeDef_ConsumBntcVert
%         dat.GroupTypeDef_ConsumBnthWrmBlood
%         dat.GroupTypeDef_eggs
%         dat.GroupTypeDef_ANYDetritus
%         dat.GroupTypeDef_terminalPlgcDetr
%         dat.GroupTypeDef_offal
%         dat.GroupTypeDef_terminalBnthDetr
%         dat.GroupTypeDef_BA
%         dat.GroupTypeDef_EM
%         dat.GroupTypeDef_fishery
%         dat.GroupTypeDef_import
%         dat.GroupTypeDef_micrograzers
%         dat.GroupTypeDef_bacteria
% pedigree information
%         dat.PreyGuild_name
%         dat.PreyGuild_code
%         dat.PedigreeParameters
%         dat.PedigreeLandings
%         dat.PedigreeDiscards
%         dat.PedigreeDiet
%         dat.import_PedigreeDiet
%         dat.PedigreeDietPrefRules
%         dat.import_PedigreeDietPrefRules
%         dat.Pedigree_PhysiologyScaler
%         dat.Pedigree_BiomassScaler
%         dat.Pedigree_FisheriesScaler
%         dat.Pedigree_DietScaler
%         dat.Pedigree_BAScaler
%         dat.Pedigree_EMScaler
%         dat.Pedigree_DiscardScaler
%         dat.Pedigree_minBiomass
%         dat.Pedigree_minPhysiology
% functional response parameters
%         dat.FunctionalResponseParams
%         dat.FunctionalResponse_matrix
%
% revision date: 7-7-2021
%       7/7/2021 fixed omission where I did not account for detritus flows
%       between detritus pools (previous omission could allow for negative
%       values in detritus-to-detritus cells in consumption matrix)


% *************************************************************************
% STEP 1: read model time series .csv file---------------------------------

% step 1a: save function name for run log ---------------------------------
fname_ReadEwE       = mfilename; % save name of this m-file to keep in saved model results
display(['Running: ' fname_ReadEwE])
dat.fname_ReadEwE	= fname_ReadEwE;
% -------------------------------------------------------------------------

% step 1b: read in first column (which contains clm address info for parameters)
fid         = fopen(readFile);
[ModelInfo] = textscan(fid, '%f %*[^\r]', 'delimiter', ',');
fclose(fid);
a           = cell2mat(ModelInfo);
numcols     = sum(a(30:43)) + 13*1; % add 13 to account for 1 blank clm btn each parameter type (14 parameter classes)

fid         = fopen(readFile);
textformat  = ['%f %*f %s %s %s %s %s %*f ', repmat('%f', 1, numcols), '%*[^\r]']; % clm space string string string string space parameters...
indata      = textscan(fid, textformat, 'delimiter', ',', 'CollectOutput', 1);
fclose(fid);
% *************************************************************************





% *************************************************************************
% STEP 2: parse group counts-----------------------------------------------
dat.num_living               = a(1);
dat.num_detritus             = a(2);
dat.num_gear                 = a(3);
dat.num_EwEGroups            = dat.num_living + dat.num_detritus + dat.num_gear;
dat.num_LivingAndDetritus    = dat.num_living + dat.num_detritus;
dat.num_NitrogenNutrients    = a(4);
dat.num_NonNitrogenNutrients = a(5);
dat.num_EcotranGroups        = dat.num_living + dat.num_detritus + dat.num_gear + dat.num_NitrogenNutrients;
% *************************************************************************





% *************************************************************************
% STEP 3: initial parsing of .csv info-------------------------------------
% step 3a: get parameter column addresses start with scalers of parameter type clm length
clms_Parameters            = a(30);
clms_landings              = a(31);
clms_discards              = a(32);
clms_Diet                  = a(33);
clms_Mortalities           = a(34);
clms_DetritusFate          = a(35);
clms_EcotranType           = a(36);
clms_EcotranRecycling      = a(37);
clms_PedigreeParameters    = a(38);
clms_PedigreeLandings      = a(39);
clms_PedigreeDiscards      = a(40);
clms_PedigreeDiet          = a(41);
clms_PedigreeDietPrefRules = a(42);
clms_FunctionalResponse    = a(43);
% -------------------------------------------------------------------------

% step 3b: convert to clm address vectors (note sequential conversion of scalers to vectors)
clms_Parameters            = 1:clms_Parameters;
clms_landings              = (max(clms_Parameters)+2):(max(clms_Parameters)+1+clms_landings);
clms_discards              = (max(clms_landings)+2):(max(clms_landings)+1+clms_discards);
clms_Diet                  = (max(clms_discards)+2):(max(clms_discards)+1+clms_Diet);
clms_Mortalities           = (max(clms_Diet)+2):(max(clms_Diet)+1+clms_Mortalities);
clms_DetritusFate          = (max(clms_Mortalities)+2):(max(clms_Mortalities)+1+clms_DetritusFate);
clms_EcotranType           = (max(clms_DetritusFate)+2):(max(clms_DetritusFate)+1+clms_EcotranType);
clms_EcotranRecycling      = (max(clms_EcotranType)+2):(max(clms_EcotranType)+1+clms_EcotranRecycling);
clms_PedigreeParameters    = (max(clms_EcotranRecycling)+2):(max(clms_EcotranRecycling)+1+clms_PedigreeParameters);
clms_PedigreeLandings      = (max(clms_PedigreeParameters)+2):(max(clms_PedigreeParameters)+1+clms_PedigreeLandings);
clms_PedigreeDiscards      = (max(clms_PedigreeLandings)+2):(max(clms_PedigreeLandings)+1+clms_PedigreeDiscards);
clms_PedigreeDiet          = (max(clms_PedigreeDiscards)+2):(max(clms_PedigreeDiscards)+1+clms_PedigreeDiet);
clms_PedigreeDietPrefRules = (max(clms_PedigreeDiet)+2):(max(clms_PedigreeDiet)+1+clms_PedigreeDietPrefRules);
clms_FunctionalResponse    = (max(clms_PedigreeDietPrefRules)+2):(max(clms_PedigreeDietPrefRules)+1+clms_FunctionalResponse);
% -------------------------------------------------------------------------

% step 3c: parse parameter group types
data_to_parse              = cell2mat(indata(3));
looky_NaN                  = find(data_to_parse == 99999); % find blank cells via blank cell code
data_to_parse(looky_NaN)   = NaN; % replace blank cells with NaN
parameters                 = data_to_parse(1:dat.num_EwEGroups, clms_Parameters);
landings                   = data_to_parse(1:dat.num_LivingAndDetritus, clms_landings);
discards                   = data_to_parse(1:dat.num_LivingAndDetritus, clms_discards);
diet                       = data_to_parse(1:(dat.num_LivingAndDetritus+1), clms_Diet);
mortalities                = data_to_parse(1:dat.num_LivingAndDetritus, clms_Mortalities);
EwE_DetritusFate           = data_to_parse(1:dat.num_EwEGroups, clms_DetritusFate); % note that detritus export column is not read in
EcotranType                = data_to_parse(1:dat.num_EcotranGroups, clms_EcotranType); % note that EM & BA rows are not exported from .xlsm file
EcotranRecycling           = data_to_parse(1:dat.num_EcotranGroups, clms_EcotranRecycling);
PedigreeParameters         = data_to_parse(1:dat.num_LivingAndDetritus, clms_PedigreeParameters);
PedigreeLandings           = data_to_parse(1:dat.num_LivingAndDetritus, clms_PedigreeLandings);
PedigreeDiscards           = data_to_parse(1:dat.num_LivingAndDetritus, clms_PedigreeDiscards);
PedigreeDiet               = data_to_parse(1:(dat.num_LivingAndDetritus+1), clms_PedigreeDiet);
PedigreeDietPrefRules      = data_to_parse(1:(dat.num_LivingAndDetritus+1), clms_PedigreeDietPrefRules);
FunctionalResponseParams   = data_to_parse(1:(dat.num_EcotranGroups), clms_FunctionalResponse);
% *************************************************************************





% *************************************************************************
% STEP 4: parse EcotranType and Aggregation groups-------------------------
dat.EcotranNumberCode = EcotranType(:, 1);
dat.EcotranGroupType  = EcotranType(:, 2);
dat.AggCode           = EcotranType(:, 3);
% *************************************************************************





% *************************************************************************
% STEP 5: parse EcotranRecycling-------------------------------------------
dat.BA                                          = EcotranRecycling(:, 1);
dat.EM                                          = EcotranRecycling(:, 2);
dat.EggProduction                               = EcotranRecycling(:, 3); % defined egg, gamete, or live young production (as opposed to use of eggs as EwE detritus method)
dat.Metabolism                                  = EcotranRecycling(:, 4); % defined metabolism (as opposed to derivation from EwE params (QQQ for future options)
dat.DetritusFate_feces                          = EcotranRecycling(:, 5:6);
dat.DetritusFate_senescence                     = EcotranRecycling(:, 7:8);
dat.ExcretionFate                               = EcotranRecycling(:, 9:10);
dat.ProductionLossScaler                        = EcotranRecycling(:, 11);
dat.RetentionScaler                             = EcotranRecycling(:, 12);
dat.PelagicBacterialReduction                   = a(6);
dat.BenthicBacterialReduction                   = a(7);
temp_Oxidation_NH4                              = a(8:9);
temp_Oxidation_NH4(isnan(temp_Oxidation_NH4))   = 0; % convert any NaN to a 0 (3/12/15)
dat.Oxidation_NH4                               = temp_Oxidation_NH4;
temp_PhytoUptake_NH4                            = a(10:19);
temp_PhytoUptake_NO3                            = a(20:29);
dat.PhytoUptake_NH4                             = temp_PhytoUptake_NH4(~isnan(temp_PhytoUptake_NH4)); % trim out NaNs formed by blank cells (3/12/15)
dat.PhytoUptake_NO3                             = temp_PhytoUptake_NO3(~isnan(temp_PhytoUptake_NO3)); % trim out NaNs formed by blank cells (3/12/15)
% *************************************************************************





% *************************************************************************
% STEP 6: parse GroupTypes definitions-------------------------------------
dat.GroupTypeDef_ANYNitroNutr       = a(45);
dat.GroupTypeDef_NO3                = a(46);
dat.GroupTypeDef_plgcNH4            = a(47);
dat.GroupTypeDef_bnthNH4            = a(48);
dat.GroupTypeDef_ANYPrimaryProd     = a(49);
dat.GroupTypeDef_LrgPhyto           = a(50);
dat.GroupTypeDef_SmlPhyto           = a(51);
dat.GroupTypeDef_Macrophytes        = a(52);
dat.GroupTypeDef_ANYConsumer        = a(53);
dat.GroupTypeDef_ConsumPlgcPlankton = a(54);
dat.GroupTypeDef_ConsumPlgcNekton   = a(55);
dat.GroupTypeDef_ConsumPlgcWrmBlood = a(56);
dat.GroupTypeDef_ConsumBntcInvert   = a(57);
dat.GroupTypeDef_ConsumBntcVert     = a(58);
dat.GroupTypeDef_ConsumBnthWrmBlood = a(59);
dat.GroupTypeDef_eggs               = a(60);
dat.GroupTypeDef_ANYDetritus        = a(61); % (not including eggs)
dat.GroupTypeDef_terminalPlgcDetr   = a(62);
dat.GroupTypeDef_offal              = a(63);
dat.GroupTypeDef_terminalBnthDetr   = a(64);
dat.GroupTypeDef_BA                 = a(65); % (biomass accumulation)
dat.GroupTypeDef_EM                 = a(66); % (emigration)
dat.GroupTypeDef_fishery            = a(67);
dat.GroupTypeDef_import             = a(68); % assigned in matlab code later (import??)
dat.GroupTypeDef_micrograzers       = a(69);
dat.GroupTypeDef_bacteria           = a(70);
% *************************************************************************





% *************************************************************************
% STEP 7: parse additional group counts------------------------------------
dat.num_eggs                 = length(find(dat.EcotranGroupType  == dat.GroupTypeDef_eggs));
dat.num_detritus             = length(find(floor(dat.EcotranGroupType) == dat.GroupTypeDef_ANYDetritus));
dat.num_EggsAndDetritus      = dat.num_eggs + dat.num_detritus;
dat.num_micrograzers         = length(find(dat.EcotranGroupType  == dat.GroupTypeDef_micrograzers));
dat.num_bacteria             = length(find(dat.EcotranGroupType  == dat.GroupTypeDef_bacteria));
% *************************************************************************





% *************************************************************************
% STEP 8: parse parameters-------------------------------------------------
dat.EwENumberCode                             = parameters(1:dat.num_EwEGroups, 1);  % fully resolved group code number
dat.TL                                        = parameters(1:dat.num_EwEGroups, 2);
dat.Biomass                                   = parameters(1:dat.num_EwEGroups, 3);
dat.PB                                        = parameters(1:dat.num_EwEGroups, 4);
dat.QB                                        = parameters(1:dat.num_EwEGroups, 5);
dat.EE                                        = parameters(1:dat.num_EwEGroups, 6);
dat.PQ                                        = parameters(1:dat.num_EwEGroups, 7);
dat.AE                                        = 1 - parameters(1:dat.num_EwEGroups, 8); % assimilation efficiency = 1 - NonAssimilation
toss_EwE_BA                                   = parameters(1:dat.num_EwEGroups, 9); % this is BA from EwE main sheet, use definition on EcotranRecycling sheet instead
dat.Production                                = (dat.PB .* dat.Biomass);
dat.Production(isnan(dat.Production))         = 0; % convert any NaNs (detritus groups) into 0; probably not needed
dat.DetritusImport                            = parameters(1:dat.num_EwEGroups, 10); % (t km^-2 y^-1)
dat.DetritusImport(isnan(dat.DetritusImport)) = 0; % convert any NaNs into 0; probably not needed
% *************************************************************************





% *************************************************************************
% STEP 9: parse landings & discards----------------------------------------
dat.landings = landings(1:dat.num_LivingAndDetritus, :);
dat.discards = discards(1:dat.num_LivingAndDetritus, :);
% *************************************************************************





% *************************************************************************
% STEP 10: parse detritus fate----------------------------------------------
dat.EwE_DetritusFate = EwE_DetritusFate(1:dat.num_EwEGroups, :);
% *************************************************************************





% *************************************************************************
% STEP 11: parse diet------------------------------------------------------
dat.diet                            = zeros((dat.num_LivingAndDetritus+1), (dat.num_LivingAndDetritus+1)); % initialize (add 1 row and 1 clm for headers)
dat.import_diet                     = zeros(1, (dat.num_LivingAndDetritus+1)); % initialize (add 1 clm for header)
dat.diet(2:end, 1)                  = dat.EwENumberCode(1:dat.num_LivingAndDetritus); % producer group header
dat.diet(1, 2:end)                  = dat.EwENumberCode(1:dat.num_LivingAndDetritus)'; % consumer group header
dat.diet(1, 1)                      = NaN;
dat.import_diet(1)                  = NaN;
dat.diet(2:(dat.num_LivingAndDetritus+1), 2:(dat.num_LivingAndDetritus+1)) = diet(1:dat.num_LivingAndDetritus, 1:dat.num_LivingAndDetritus); % this final dat.diet does not include import diet
dat.import_diet(2:(dat.num_LivingAndDetritus+1))                           = diet((dat.num_LivingAndDetritus+1), 1:dat.num_LivingAndDetritus); % next last row is import diet
% *************************************************************************





% *************************************************************************
% STEP 12: parse mortalities-----------------------------------------------
dat.mortality                       = zeros((dat.num_LivingAndDetritus+1), (dat.num_LivingAndDetritus+1)); % initialize
dat.mortality(2:end, 1)             = dat.EwENumberCode(1:(dat.num_LivingAndDetritus)); % producer group header
dat.mortality(1, 2:end)             = dat.EwENumberCode(1:(dat.num_LivingAndDetritus))'; % consumer group header
dat.mortality(1, 1)                 = NaN;
dat.mortality(2:(dat.num_LivingAndDetritus+1), 2:(dat.num_living+1)) = mortalities(1:dat.num_LivingAndDetritus, (dat.num_gear+1):(dat.num_gear+dat.num_living));
% *************************************************************************





% *************************************************************************
% STEP 13: create consumption matrix---------------------------------------
% step 13a: calculate consumption rate (t/km2/y) for all living groups ----
%           NOTE: consumption = mortality (1/y) * biomass (t/km2)
dat.consumption         = dat.mortality; % initialize; (2D matrix: (num_LivingAndDetritus+1) X (num_LivingAndDetritus+1)); has a row & clm header
dat.consumption(2:(dat.num_living+1), 2:(dat.num_living+1)) = ...
                          dat.mortality(2:(dat.num_living+1), 2:(dat.num_living+1)) .* repmat(dat.Biomass(1:dat.num_living), [1, dat.num_living]); % (2D matrix: producers (num_LivingAndDetritus+1) X consumers (num_LivingAndDetritus+1))
% -------------------------------------------------------------------------

% step 13b: calculate production, feces, & senescence rates for all living groups ----
production              = dat.Biomass(1:dat.num_living) .* dat.PB(1:dat.num_living); % production by each living group; (vertical vector: num_living X 1)
total_consumption       = dat.Biomass(1:dat.num_living) .* dat.QB(1:dat.num_living); % total consumption by each living group; (vertical vector: num_living X 1)
feces_production        = total_consumption .* (1 - dat.AE(1:dat.num_living));       % feces production by each living group; (vertical vector: num_living X 1)
total_landings          = sum(dat.landings(1:dat.num_living, :), 2);                 % total landings from each living group by all fleets combined; (vertical vector: num_living X 1)
total_predation         = sum(dat.consumption(2:(dat.num_living+1), 2:(dat.num_living+1)), 2); % total predation upon each living group; (vertical vector: num_living X 1)
senescence_production	= production - (total_predation + total_landings);           % (vertical vector: num_living X 1)
detritus_production     = feces_production + senescence_production;                  % detritus production by each living group; (vertical vector: num_living X 1)

total_consumed_detritus	= sum(dat.consumption((dat.num_living+2):(dat.num_LivingAndDetritus+1), 2:(dat.num_living+1)), 2); % detritus (incl. eggs) that is consumed by living groups; (vertical vector: num_EggsAndDetritus X 1)
% -------------------------------------------------------------------------


% step 13c: calculate flow to detritus ------------------------------------
%           flow from living; (t/km2/y)
%           flow from detritus groups is still 0
Flow2Detritus           = zeros(dat.num_LivingAndDetritus, dat.num_EggsAndDetritus); % initialize; (2D matrix: num_LivingAndDetritus X num_EggsAndDetritus)
Flow2Detritus(1:dat.num_living, 1:dat.num_EggsAndDetritus) = ...
                        dat.EwE_DetritusFate(1:dat.num_living, 1:dat.num_EggsAndDetritus) .* repmat(detritus_production, [1, dat.num_EggsAndDetritus]); % (2D matrix: num_LivingAndDetritus X num_EggsAndDetritus)
% -------------------------------------------------------------------------


% step 13d: calculate fishery discard -------------------------------------
%           flow from each living group due to fleet discards; (t/km2/y)
discard_detritus        = zeros(dat.num_LivingAndDetritus, dat.num_EggsAndDetritus); % initialize; (2D matrix: num_LivingAndDetritus X num_EggsAndDetritus)
for fleet_loop = 1:dat.num_gear
    
    current_Discard     = dat.discards(:, fleet_loop); % discard of each live group by the current fleet; (vertical vector: num_LivingAndDetritus X 1)
    repmat_Discard      = repmat(current_Discard, [1, dat.num_EggsAndDetritus]); % (3D matrix: num_LivingAndDetritus X num_EggsAndDetritus)
    
    current_DiscardFate	= dat.EwE_DetritusFate((dat.num_LivingAndDetritus + fleet_loop), :); % (horizontal vector: 1 X num_EggsAndDetritus); NOTE: fleets in rows of model groups must be in same order as across clms in discards;
    repmat_DiscardFate	= repmat(current_DiscardFate, [dat.num_LivingAndDetritus, 1]); % (3D matrix: num_LivingAndDetritus X num_EggsAndDetritus)

    discard_detritus = discard_detritus + (repmat_Discard .* repmat_DiscardFate); % summed discard rate across all fleets; (t/km2/y); (2D matrix: num_LivingAndDetritus X num_EggsAndDetritus)
end
    
Flow2Detritus = Flow2Detritus + discard_detritus; % (2D matrix: num_LivingAndDetritus X num_EggsAndDetritus)
% -------------------------------------------------------------------------


% step 13e: calculate detritus flow to other detritus groups --------------
total_Flow2Detritus         = sum(Flow2Detritus, 1); % (horizontal vector: 1 X num_EggsAndDetritus)
surplus_detritus            = total_Flow2Detritus' - total_consumed_detritus; % (vertical vector: num_EggsAndDetritus X 1); NOTE: may have negative values at this point until detritus flow to detritus is included
DetritusToDetritusFate      = dat.EwE_DetritusFate((dat.num_living+1):dat.num_LivingAndDetritus, :); % (2D matrix: num_EggsAndDetritus (source) X num_EggsAndDetritus (destiny)); NOTE: this assumes detritus grps are located after living groups and before fleets

DetritusToDetritus_array    = repmat(surplus_detritus, [1, dat.num_EggsAndDetritus]) .* DetritusToDetritusFate; % flow between detritus groups; (t/km2/y); (2D matrix: num_EggsAndDetritus (source) X num_EggsAndDetritus (destiny))

% account for flow between detritus groups
sum_DetritusToDetritus      = sum(DetritusToDetritus_array);              % (horizontal vector: 1 X num_EggsAndDetritus (destiny))
sum_DetritusToDetritus(sum_DetritusToDetritus<0) = 0;                     % filter out negative values (negative flows will be error-checked below)
adjusted_surplus_detritus	= surplus_detritus + sum_DetritusToDetritus'; % (vertical vector: num_EggsAndDetritus (source) X 1)
DetritusToDetritus_array	= repmat(adjusted_surplus_detritus, [1, dat.num_EggsAndDetritus]) .* DetritusToDetritusFate; % flow between detritus groups; (t/km2/y); (2D matrix: num_EggsAndDetritus (source) X num_EggsAndDetritus (destiny))
% -------------------------------------------------------------------------


% step 13f: put final consumption matrix together -------------------------
dat.consumption(2:(dat.num_LivingAndDetritus+1), (dat.num_living+2):(dat.num_LivingAndDetritus+1))                  = Flow2Detritus;
dat.consumption((dat.num_living+2):(dat.num_LivingAndDetritus+1), (dat.num_living+2):(dat.num_LivingAndDetritus+1)) = DetritusToDetritus_array;
% -------------------------------------------------------------------------


% step 13g: calculate import consumption (total consumption .* import_diet)
dat.total_consumption(1, 1:(dat.num_LivingAndDetritus+1)) = NaN; % initialize, add 1 cell up front for header
dat.total_consumption(1, 2:end)                           = dat.Biomass(1:dat.num_LivingAndDetritus) .* dat.QB(1:dat.num_LivingAndDetritus); % ignore fisheries at the end
dat.import_consumption                                    = dat.total_consumption .* dat.import_diet; % (import_consumption = total consumption .* import_diet) (horizontal vector; first cell is a header)
% *************************************************************************





% *************************************************************************
% STEP 14: parse group names (& PreyGuild strings)-------------------------
names_to_parse                        = indata(2);
names_to_parse                        = names_to_parse{1};
dat.EwEName                           = names_to_parse(1:dat.num_EwEGroups, 1); % formal EwE trophic group names
dat.EcotranName                       = names_to_parse(1:dat.num_EcotranGroups, 2); % formal Ecotran trophic group names
dat.AggEwEName                        = names_to_parse(1:dat.num_EcotranGroups, 3); % formal AggEwE trophic group names
temp_PreyGuild_name                   = names_to_parse(1:20, 4); % excel EwE macro allows only 20 PreyGuilds right now
temp_PreyGuild_code                   = names_to_parse(1:20, 5); % excel EwE macro allows only 20 PreyGuilds right now
looky_emptyGuild                      = find(strcmp(temp_PreyGuild_code, '0') == 1 | strcmp(temp_PreyGuild_code, '') == 1);
temp_PreyGuild_name(looky_emptyGuild) = []; % remove empty ('' or '0') rows (3/12/15)
temp_PreyGuild_code(looky_emptyGuild) = []; % remove empty ('' or '0') rows (3/12/15)
dat.PreyGuild_name                    = temp_PreyGuild_name;
dat.PreyGuild_code                    = temp_PreyGuild_code;
% *************************************************************************





% *************************************************************************
% STEP 15: parse pedigree--------------------------------------------------
% step 15a: pedigrees for parameters and fisheries ------------------------
dat.PedigreeParameters        = PedigreeParameters;
dat.PedigreeLandings          = PedigreeLandings;
dat.PedigreeDiscards          = PedigreeDiscards;
% -------------------------------------------------------------------------


% step 15b: PedigreeDiet --------------------------------------------------
dat.PedigreeDiet                    = zeros((dat.num_LivingAndDetritus+1), (dat.num_LivingAndDetritus+1)); % initialize (add 1 row and 1 clm for headers)
dat.import_PedigreeDiet             = zeros(1, (dat.num_LivingAndDetritus+1)); % initialize (add 1 clm for header)
dat.PedigreeDiet(2:end, 1)          = dat.EwENumberCode(1:dat.num_LivingAndDetritus); % producer group header
dat.PedigreeDiet(1, 2:end)          = dat.EwENumberCode(1:dat.num_LivingAndDetritus)'; % consumer group header
dat.PedigreeDiet(1, 1)              = NaN;
dat.import_PedigreeDiet(1)          = NaN;
dat.PedigreeDiet(2:(dat.num_LivingAndDetritus+1), 2:(dat.num_LivingAndDetritus+1)) = PedigreeDiet(1:dat.num_LivingAndDetritus, 1:dat.num_LivingAndDetritus); % this final dat.PedigreeDiet does not include import diet
dat.import_PedigreeDiet(2:(dat.num_LivingAndDetritus+1))                           = PedigreeDiet((dat.num_LivingAndDetritus+1), 1:dat.num_LivingAndDetritus); % next last row is import PedigreeDiet
% -------------------------------------------------------------------------


% step 15c: PedigreeDietPrefRules -----------------------------------------
dat.PedigreeDietPrefRules                    = zeros((dat.num_LivingAndDetritus+1), (dat.num_LivingAndDetritus+1)); % initialize (add 1 row and 1 clm for headers)
dat.import_PedigreeDietPrefRules             = zeros(1, (dat.num_LivingAndDetritus+1)); % initialize (add 1 clm for header)
dat.PedigreeDietPrefRules(2:end, 1)          = dat.EwENumberCode(1:dat.num_LivingAndDetritus); % producer group header
dat.PedigreeDietPrefRules(1, 2:end)          = dat.EwENumberCode(1:dat.num_LivingAndDetritus)'; % consumer group header
dat.PedigreeDietPrefRules(1, 1)              = NaN;
dat.import_PedigreeDietPrefRules(1)          = NaN;
dat.PedigreeDietPrefRules(2:(dat.num_LivingAndDetritus+1), 2:(dat.num_LivingAndDetritus+1)) = PedigreeDietPrefRules(1:dat.num_LivingAndDetritus, 1:dat.num_LivingAndDetritus); % this final dat.PedigreeDietPrefRules does not include import diet
dat.import_PedigreeDietPrefRules(2:(dat.num_LivingAndDetritus+1))                           = PedigreeDietPrefRules((dat.num_LivingAndDetritus+1), 1:dat.num_LivingAndDetritus); % next last row is import PedigreeDietPrefRules
% -------------------------------------------------------------------------


% step 15d: pedigree scaling terms ----------------------------------------
dat.Pedigree_PhysiologyScaler = a(75);
dat.Pedigree_BiomassScaler    = a(76);
dat.Pedigree_FisheriesScaler  = a(77);
dat.Pedigree_DietScaler       = a(78);
dat.Pedigree_BAScaler         = a(79);
dat.Pedigree_EMScaler         = a(80);
dat.Pedigree_DiscardScaler    = a(81);
dat.Pedigree_minBiomass       = a(82);
dat.Pedigree_minPhysiology    = a(83);
% *************************************************************************





% *************************************************************************
% STEP 16: parse functional response parameters----------------------------
%   NOTE: leaving parsing for future. clms 1-4 apply to all prey for each predator
%         clms 5-end are specific params for each trophic link (not yet defined)
dat.FunctionalResponseParams  = FunctionalResponseParams(:, 1:4);
dat.FunctionalResponse_matrix = FunctionalResponseParams(:, 5:end);
% *************************************************************************





% *************************************************************************
% STEP 17: some error-checking---------------------------------------------
% step 17a: check number of primary producers against nutrient uptake definitions
num_PrimaryProducer = length(find(fix(dat.EcotranGroupType)  == dat.GroupTypeDef_ANYPrimaryProd));
dat.PhytoUptake_NH4 = temp_PhytoUptake_NH4(1:num_PrimaryProducer); % (2/2/16)
dat.PhytoUptake_NO3 = temp_PhytoUptake_NO3(1:num_PrimaryProducer); % (2/2/16)

num_PhytoUptake_NH4 = length(dat.PhytoUptake_NH4);
num_PhytoUptake_NO3 = length(dat.PhytoUptake_NO3);

if num_PhytoUptake_NH4 ~= num_PhytoUptake_NO3
    error('ECOTRAN setup ERROR: number of NO3 and number of NO3 uptake definitions do not match')
end

if num_PrimaryProducer ~= num_PhytoUptake_NO3
    error('ECOTRAN setup ERROR: number of primary producers and nutrient uptake definitions do not match')
end
% -------------------------------------------------------------------------


% step 17b: check for negative metabolism ---------------------------------
BioenergeticBudget_metabolism = (dat.AE - dat.PQ); % fraction of consumption going to metabolism (NH4 production); Metabolism = 1 - PQ - feces = (1 - (EwE_PQ + (1-EwE_AE)))
looky_negativeMetabolism = find(BioenergeticBudget_metabolism <= 0);
if ~isempty(looky_negativeMetabolism)
    display('   -->EwE parameter WARNING: parameters imply negative or zero')
    display('      metabolism for at least 1 group (metabolism = ae - pq)')
    display(['      - bad groups: ' num2str(looky_negativeMetabolism(:)')])
end
% -------------------------------------------------------------------------


% % step 17c: check for ECOTRAN feces fate entries for all groups with feces
% FFF will need to sync rows of dat.AE and at.DetritusFate_feces for this test to work
% looky_feces_producers = find(dat.AE < 1);
% sum_FecesFate         = sum(dat.DetritusFate_feces, 2);
% FecesProduced         = sum_FecesFate(looky_feces_producers);
% looky_noFecesFate     = find(FecesProduced <= 0);
% if ~isempty(looky_noFecesFate)
%     error(['ECOTRAN setup error: Feces fate not defined for at least one consumer group; bad groups: ' num2str(looky_noFecesFate(:)')])
% end
% % -------------------------------------------------------------------------


% step 17d: check for full accounting of EwE detritus ---------------------
sum_EwE_DetritusFate                                = sum(dat.EwE_DetritusFate, 2);
looky_terminalBenthicDetritus                       = find(dat.EcotranGroupType == dat.GroupTypeDef_terminalBnthDetr);
looky_terminalBenthicDetritus                       = looky_terminalBenthicDetritus - dat.num_NitrogenNutrients;
sum_EwE_DetritusFate(looky_terminalBenthicDetritus) = []; % remove terminal benthic detritus
looky_ExportedDetritus                              = find(sum_EwE_DetritusFate  < 1);
if ~isempty(looky_ExportedDetritus)
    display('   -->EwE parameter WARNING: EwE detritus fates do not sum to 1.')
    display('      - These terms will reflect that loss of feces and senescence from the system:')
    display('         EnergyBudget_matrix, BioenergeticBudget_ProductionDetail,')
    display('         BioenergeticBudget_OtherMortDetail, BioenergeticBudget_FecesDetail,')
    display('         & ProductionBudget_OtherMortDetail.') 
    display('      - The BioenergeticBudget & ProductioBudget will NOT reflect this exported feces')
    display('         and senescence detritus.')
end
% -------------------------------------------------------------------------

% step 17e: check consumption matrix & balance between INFLOW & OUTFLOW ---
test_inflow     = sum(dat.consumption(2:end, 2:end)); % consumption INTO each living and detritus group
test_outflow	= sum(dat.consumption(2:end, 2:end), 2)'; % predation OUT of each living and detritus group; NOTE transpose
test_balance	= test_inflow - test_outflow;

[~, looky_NegativeConsumer] = find(dat.consumption(2:end, 2:end) < 0);
if ~isempty(looky_NegativeConsumer)
    display('   -->consumption matrix WARNING: consumer(s) with negative consumption.')
    display(['      - bad groups: ' num2str(looky_NegativeConsumer)'])
end

% [~, looky_ConsumerImbalance] = find(test_balance < 0); % need to ignore primary producers for this test to work
% *************************************************************************


% end function*************************************************************