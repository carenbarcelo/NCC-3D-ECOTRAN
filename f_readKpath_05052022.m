function [dat] = f_readKpath_05052022(ReadFile)
% load EwE model from named Kpath .xlsm file and prepare data
% includes conversion from mortality matrix in K-path to consumption matrix
% NOTE: this code allows for 10 primary producers
% by Jim Ruzicka
% calls:
%	none
% takes:
%	ReadFile                            concatenated directory & file name
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
% revision date: 5-4-2022
%       7/7/2021 fixed omission where I did not account for detritus flows
%       between detritus pools (previous omission could allow for negative
%       values in detritus-to-detritus cells in consumption matrix)


% ReadFile = '/Users/jimsebi/Documents/11_FoodWeb_models/2_GoMexOcn/GoMexOcn_12282021_ZZH.xlsm';


% *************************************************************************
% STEP 1: read model time series .csv file---------------------------------

% step 1a: save function name for run log ---------------------------------
fname_readKpath         = mfilename; % save name of this m-file to keep in saved model results
display(['Running: ' fname_readKpath])
% -------------------------------------------------------------------------

% step 1b: read and parse info from 'Main' tab ----------------------------
ModelInfo1              = readtable(ReadFile, 'sheet', 'Main', 'range', 'B3:B5', 'ReadVariableNames', false);

num_living              = ModelInfo1{1, 1};
num_detritus_EwE      	= ModelInfo1{2, 1};
num_gear                = ModelInfo1{3, 1};

num_EwEGroups           = num_living + num_detritus_EwE + num_gear;
num_LivingAndDetritus	= num_living + num_detritus_EwE;

range_Main              = ['E2:N' num2str(num_EwEGroups + 1)];
ModelInfo2              = readtable(ReadFile, 'sheet', 'Main', 'range', range_Main, 'ReadVariableNames', false);

text_Main               = ModelInfo2{:, 1};
double_Main             = ModelInfo2{:, 2:10};
double_Main(isnan(double_Main)) = 0; % convert missing parameter values (NaN) to 0

EwEName                 = text_Main;
BA                      = double_Main(:, 7);
AE                      = 1 - double_Main(:, 8);
DetritusImport          = double_Main(:, 9);

EwENumberCode           = (1:num_EwEGroups)';

% -------------------------------------------------------------------------


% step 1c: read and parse info from 'MainOutputs' tab ---------------------
range_MainOutputs	= ['C2:H' num2str(num_EwEGroups + 1)];
ModelInfo3          = readtable(ReadFile, 'sheet', 'MainOutputs', 'range', range_MainOutputs, 'ReadVariableNames', false);
double_MainOutputs 	= ModelInfo3{:, 1:end};
double_MainOutputs(isnan(double_MainOutputs)) = 0; % convert missing parameter values (NaN) to 0

TL                  = double_MainOutputs(:, 1);
Biomass         	= double_MainOutputs(:, 2);
PB                  = double_MainOutputs(:, 3);
QB                  = double_MainOutputs(:, 4);
EE                  = double_MainOutputs(:, 5);
PQ                  = double_MainOutputs(:, 6);
% -------------------------------------------------------------------------


% step 1d: read and parse info from 'Diets' tab ---------------------------
leader_clms     = 1;
second_letter	= char(mod((leader_clms + num_LivingAndDetritus - 1), ('Z' - 'A' + 1)) + 'A');

if num_LivingAndDetritus > (26-leader_clms)
    first_letter	= char('A' + mod(floor((leader_clms + num_LivingAndDetritus - 1 - 26)/26), 'A'));
    column_letters	= [first_letter second_letter];
else
    column_letters	= second_letter;
end

leader_rows             = 1;
range_Diets             = ['B2:' column_letters num2str(num_LivingAndDetritus + leader_rows + 1)];
ModelInfo3              = readtable(ReadFile, 'sheet', 'Diets', 'range', range_Diets, 'ReadVariableNames', false);

diet                    = zeros((num_LivingAndDetritus + 1), (num_LivingAndDetritus + 1)); % initialize with blank row and column (FFF eventually need to get rid of these null headers)
diet(2:end, 1)          = EwENumberCode(1:num_LivingAndDetritus); % producer group header
diet(1, 2:end)          = EwENumberCode(1:num_LivingAndDetritus)'; % consumer group header
import_diet             = zeros(1, (num_LivingAndDetritus + 1)); % initialize with blank column (FFF eventually need to get rid of these null headers)

diet(2:end, 2:end)      = ModelInfo3{1:num_LivingAndDetritus, 1:num_LivingAndDetritus};
import_diet(1, 2:end)	= ModelInfo3{end, :};
% -------------------------------------------------------------------------


% step 1e: read and parse info from 'Detritus' tab ------------------------
leader_clms             = 2;
second_letter           = char(mod((leader_clms + num_detritus_EwE - 1), ('Z' - 'A' + 1)) + 'A');

if num_detritus_EwE > (26-leader_clms)
    first_letter	= char('A' + mod(floor((leader_clms + num_detritus_EwE - 1 - 26)/26), 'A'));
    column_letters	= [first_letter second_letter];
else
    column_letters	= second_letter;
end

leader_rows             = 1;
range_Detritus          = ['C2:' column_letters num2str(num_EwEGroups + leader_rows)];
ModelInfo4              = readtable(ReadFile, 'sheet', 'Detritus', 'range', range_Detritus, 'ReadVariableNames', false);
EwE_DetritusFate        = ModelInfo4{:, :};
% -------------------------------------------------------------------------


% step 1f: read and parse info from 'Fishing' & "Discards" tabs -----------
leader_clms             = 2;
second_letter           = char(mod((leader_clms + num_gear - 1), ('Z' - 'A' + 1)) + 'A');

if num_gear > (26-leader_clms)
    first_letter	= char('A' + mod(floor((leader_clms + num_gear - 1 - 26)/26), 'A'));
    column_letters	= [first_letter second_letter];
else
    column_letters	= second_letter;
end

leader_rows             = 1;
range_gear              = ['C2:' column_letters num2str(leader_rows + num_LivingAndDetritus)]; % NOTE: code can only handle 26 detritus groups
ModelInfo5              = readtable(ReadFile, 'sheet', 'Fishing', 'range', range_gear, 'ReadVariableNames', false);
landings                = ModelInfo5{:, :};
ModelInfo6              = readtable(ReadFile, 'sheet', 'Discards', 'range', range_gear, 'ReadVariableNames', false);
discards                = ModelInfo6{:, :};
% -------------------------------------------------------------------------


% step 1g: read and parse info from 'Mortalities" tab ---------------------
leader_clms             = 3;
second_letter           = char(mod((leader_clms + num_gear - 1 + 1), ('Z' - 'A' + 1)) + 'A');

if num_gear > (26-leader_clms)
    first_letter	= char('A' + mod(floor((leader_clms + num_gear - 1 - 26)/26), 'A'));
    column_letters	= [first_letter second_letter];
else
    column_letters	= second_letter;
end

leader_rows             = 1;
start_corner            = [column_letters num2str(leader_rows + 1)];


leader_clms             = 3;
second_letter           = char(mod((leader_clms + (num_living + num_gear) - 1), ('Z' - 'A' + 1)) + 'A');

if (num_living + num_gear) > (26-leader_clms)
    first_letter	= char('A' + mod(floor((leader_clms + (num_living + num_gear) - 1 - 26)/26), 'A'));
    column_letters	= [first_letter second_letter];
else
    column_letters	= second_letter;
end

leader_rows             = 1;
range_Mortalities     	= [start_corner ':' column_letters num2str(leader_rows + num_LivingAndDetritus)];
ModelInfo7            	= readtable(ReadFile, 'sheet', 'Mortalities', 'range', range_Mortalities, 'ReadVariableNames', false);

mortality             	= zeros((num_LivingAndDetritus+1), (num_LivingAndDetritus+1)); % initialize

mortality(2:end, 1)   	= EwENumberCode(1:num_LivingAndDetritus); % producer group header
mortality(1, 2:end)   	= EwENumberCode(1:num_LivingAndDetritus)'; % consumer group header
mortality(2:(num_LivingAndDetritus+1), 2:(num_living+1)) = ModelInfo7{:, :};
mortality(isnan(mortality))	= 0;
% -------------------------------------------------------------------------


% step 1h: read and parse info from 'EcotranType" tab ---------------------
ModelInfo8                  = readtable(ReadFile, 'sheet', 'EcotranType', 'range', 'B2:B3', 'ReadVariableNames', false);

num_NitrogenNutrients       = ModelInfo8{1, 1};
num_NonNitrogenNutrients	= ModelInfo8{2, 1};
num_EcotranGroups           = num_EwEGroups + num_NitrogenNutrients + num_NonNitrogenNutrients;

range_EcotranType       	= ['D2:H' num2str(num_EcotranGroups + 1)];
ModelInfo9                  = readtable(ReadFile, 'sheet', 'EcotranType', 'range', range_EcotranType, 'ReadVariableNames', false);

EcotranName                 = ModelInfo9{:, 1};
EcotranNumberCode           = ModelInfo9{:, 2};
EcotranGroupType            = ModelInfo9{:, 3};
AggCode                     = ModelInfo9{:, 4};
AggEwEName                  = ModelInfo9{:, 5};

range_GrpTypeDef            = ['A7:B32'];
ModelInfo10              	= readtable(ReadFile, 'sheet', 'EcotranType', 'range', range_GrpTypeDef, 'ReadVariableNames', false);

GrpTypeDef_names                = ModelInfo10{:, 1};
GroupTypeDef_ANYNitroNutr       = ModelInfo10{1, 2};
GroupTypeDef_NO3                = ModelInfo10{2, 2};
GroupTypeDef_plgcNH4            = ModelInfo10{3, 2};
GroupTypeDef_bnthNH4            = ModelInfo10{4, 2};
GroupTypeDef_ANYPrimaryProd     = ModelInfo10{5, 2};
GroupTypeDef_LrgPhyto           = ModelInfo10{6, 2};
GroupTypeDef_SmlPhyto           = ModelInfo10{7, 2};
GroupTypeDef_Macrophytes        = ModelInfo10{8, 2};
GroupTypeDef_ANYConsumer        = ModelInfo10{9, 2};
GroupTypeDef_ConsumPlgcPlankton = ModelInfo10{10, 2};
GroupTypeDef_ConsumPlgcNekton   = ModelInfo10{11, 2};
GroupTypeDef_ConsumPlgcWrmBlood	= ModelInfo10{12, 2};
GroupTypeDef_ConsumBntcInvert   = ModelInfo10{13, 2};
GroupTypeDef_ConsumBntcVert     = ModelInfo10{14, 2};
GroupTypeDef_ConsumBnthWrmBlood = ModelInfo10{15, 2};
GroupTypeDef_eggs               = ModelInfo10{16, 2};
GroupTypeDef_ANYDetritus        = ModelInfo10{17, 2};
GroupTypeDef_terminalPlgcDetr   = ModelInfo10{18, 2};
GroupTypeDef_offal              = ModelInfo10{19, 2};
GroupTypeDef_terminalBnthDetr   = ModelInfo10{20, 2};
GroupTypeDef_BA                 = ModelInfo10{21, 2};
GroupTypeDef_EM                 = ModelInfo10{22, 2};
GroupTypeDef_fishery            = ModelInfo10{23, 2};
GroupTypeDef_import             = ModelInfo10{24, 2};
GroupTypeDef_micrograzers       = ModelInfo10{25, 2};
GroupTypeDef_bacteria           = ModelInfo10{26, 2};

range_PreyGuild             = ['J2:K21'];
ModelInfo11              	= readtable(ReadFile, 'sheet', 'EcotranType', 'range', range_PreyGuild, 'ReadVariableNames', false);

PreyGuild_name              = table2cell(ModelInfo11(:, 1));
PreyGuild_code              = table2cell(ModelInfo11(:, 2));
looky_NaN                   = isnan([PreyGuild_name{:}]);
PreyGuild_name(looky_NaN)	= [];
PreyGuild_code(looky_NaN)	= [];
% -------------------------------------------------------------------------


% step 1i: read and parse info from 'EcotranRecycling" tab ----------------
ModelInfo12                         = readtable(ReadFile, 'sheet', 'EcotranRecycling', 'range', 'B2:B28', 'ReadVariableNames', false);

PelagicBacterialReduction           = ModelInfo12{1, 1};
BenthicBacterialReduction           = ModelInfo12{2, 1};
PelagicBacterialReduction(isnan(PelagicBacterialReduction))	= 0; % convert any NaN to 0
BenthicBacterialReduction(isnan(BenthicBacterialReduction))	= 0; % convert any NaN to 0

Oxidation_NH4                       = ModelInfo12{4:5, 1};
Oxidation_NH4(isnan(Oxidation_NH4))	= 0; % convert any NaN to 0

PhytoUptake_NH4                 	= ModelInfo12{7:16, 1};
PhytoUptake_NO3                   	= ModelInfo12{18:27, 1};
PhytoUptake_NH4                    	= PhytoUptake_NH4(~isnan(PhytoUptake_NH4)); % trim out NaNs formed by blank cells
PhytoUptake_NO3                  	= PhytoUptake_NO3(~isnan(PhytoUptake_NO3)); % trim out NaNs formed by blank cells


range_EcotranRecycling              = ['G3:R' num2str(num_EcotranGroups + 2)];
ModelInfo13                         = readtable(ReadFile, 'sheet', 'EcotranRecycling', 'range', range_EcotranRecycling, 'ReadVariableNames', false);

BA                                  = ModelInfo13{:, 1};
EM                                  = ModelInfo13{:, 2};
EggProduction                       = ModelInfo13{:, 3}; % defined egg, gamete, or live young production (as opposed to use of eggs as EwE detritus method)
Metabolism                          = ModelInfo13{:, 4}; % defined metabolism (as opposed to derivation from EwE params (QQQ for future options)
DetritusFate_feces               	= ModelInfo13{:, 5:6};
DetritusFate_senescence          	= ModelInfo13{:, 7:8};
ExcretionFate                       = ModelInfo13{:, 9:10};
ProductionLossScaler              	= ModelInfo13{:, 11};
RetentionScaler                     = ModelInfo13{:, 12};
% -------------------------------------------------------------------------


% step 1j: read and parse info from 'FunctionalResponse' tab --------------
leader_clms         = 1;
second_letter       = char(mod((leader_clms + (num_EcotranGroups+4) - 1), ('Z' - 'A' + 1)) + 'A');

if (num_EcotranGroups+4) > (26-leader_clms)
    first_letter	= char('A' + mod(floor((leader_clms + (num_EcotranGroups+4) - 1 - 26)/26), 'A'));
    column_letters	= [first_letter second_letter];
else
    column_letters	= second_letter;
end

leader_rows                 = 1;
range_FunctionalResponse	= ['B2:' column_letters num2str(leader_rows + num_EcotranGroups)];
ModelInfo14             	= readtable(ReadFile, 'sheet', 'FunctionalResponse', 'range', range_FunctionalResponse, 'ReadVariableNames', false);

%   NOTE: leaving parsing for future. clms 1-4 apply to all prey for each predator
%         clms 5-end are specific params for each trophic link (not yet defined)
FunctionalResponseParams	= ModelInfo14{:, 1:4};
FunctionalResponse_matrix	= ModelInfo14{:, 5:end};

FunctionalResponseParams(isnan(FunctionalResponseParams))	= 0;
FunctionalResponse_matrix(isnan(FunctionalResponse_matrix))	= 0;
% -------------------------------------------------------------------------


% step 1k: read and parse info from 'pedigree_Parameters' tab --------------
range_pedigree_Parameters        	= ['F2:AP'  num2str(num_LivingAndDetritus + 1)];

ModelInfo15                     	= readtable(ReadFile, 'sheet', 'pedigree_Parameters', 'range', range_pedigree_Parameters, 'ReadVariableNames', false);

PedigreeParameters                  = ModelInfo15{:, :};

ModelInfo20                     	= readtable(ReadFile, 'sheet', 'pedigree_Parameters', 'range', 'B2:B10', 'ReadVariableNames', false);
Pedigree_PhysiologyScaler           = ModelInfo20{1, 1};
Pedigree_BiomassScaler              = ModelInfo20{2, 1};
Pedigree_FisheriesScaler            = ModelInfo20{3, 1};
Pedigree_DietScaler                 = ModelInfo20{4, 1};
Pedigree_BAScaler                   = ModelInfo20{5, 1};
Pedigree_EMScaler                   = ModelInfo20{6, 1};
Pedigree_DiscardScaler              = ModelInfo20{7, 1};
Pedigree_minBiomass                 = ModelInfo20{8, 1};
Pedigree_minPhysiology              = ModelInfo20{9, 1};
% -------------------------------------------------------------------------


% step 1l: read and parse info from 'pedigree_Diet' & 'pedigree_DietPrefRules' tabs
ModelInfo16                  	= readtable(ReadFile, 'sheet', 'pedigree_Diet', 'range', range_Diets, 'ReadVariableNames', false);

PedigreeDiet                    = zeros((num_LivingAndDetritus+1), (num_LivingAndDetritus+1)); % initialize (add 1 row and 1 clm for headers)
import_PedigreeDiet             = zeros(1, (num_LivingAndDetritus+1)); % initialize (add 1 clm for header)
PedigreeDiet(2:end, 1)          = EwENumberCode(1:num_LivingAndDetritus); % producer group header
PedigreeDiet(1, 2:end)          = EwENumberCode(1:num_LivingAndDetritus)'; % consumer group header
PedigreeDiet(1, 1)              = NaN;
import_PedigreeDiet(1)          = NaN;
PedigreeDiet(2:(num_LivingAndDetritus+1), 2:(num_LivingAndDetritus+1))	= ModelInfo16{1:num_LivingAndDetritus, 1:num_LivingAndDetritus}; % this final PedigreeDiet does not include import diet
import_PedigreeDiet(2:(num_LivingAndDetritus+1))                    	= ModelInfo16{(num_LivingAndDetritus+1), 1:num_LivingAndDetritus}; % next last row is import PedigreeDiet


ModelInfo17                   	= readtable(ReadFile, 'sheet', 'pedigree_DietPrefRules', 'range', range_Diets, 'ReadVariableNames', false);

PedigreeDietPrefRules       	= zeros((num_LivingAndDetritus+1), (num_LivingAndDetritus+1)); % initialize (add 1 row and 1 clm for headers)
import_PedigreeDietPrefRules  	= zeros(1, (num_LivingAndDetritus+1)); % initialize (add 1 clm for header)
PedigreeDietPrefRules(2:end, 1)	= EwENumberCode(1:num_LivingAndDetritus); % producer group header
PedigreeDietPrefRules(1, 2:end)	= EwENumberCode(1:num_LivingAndDetritus)'; % consumer group header
PedigreeDietPrefRules(1, 1)    	= NaN;
import_PedigreeDietPrefRules(1)	= NaN;
PedigreeDietPrefRules(2:(num_LivingAndDetritus+1), 2:(num_LivingAndDetritus+1))	= ModelInfo17{1:num_LivingAndDetritus, 1:num_LivingAndDetritus}; % this final PedigreeDiet does not include import diet
import_PedigreeDietPrefRules(2:(num_LivingAndDetritus+1))                    	= ModelInfo17{(num_LivingAndDetritus+1), 1:num_LivingAndDetritus}; % next last row is import PedigreeDiet
% -------------------------------------------------------------------------


% step 1m: read and parse info from 'pedigree_Landings' & 'pedigree_Discards' tabs
ModelInfo18                     = readtable(ReadFile, 'sheet', 'pedigree_Landings', 'range', range_gear, 'ReadVariableNames', false);
PedigreeLandings                = ModelInfo18{:, :};
ModelInfo19                     = readtable(ReadFile, 'sheet', 'pedigree_Discards', 'range', range_gear, 'ReadVariableNames', false);
PedigreeDiscards            	= ModelInfo19{:, :};
% -------------------------------------------------------------------------










% *************************************************************************
% STEP 2: calculate derived parameters-------------------------------------

% step 2a: calculate additional group counts ------------------------------
num_detritus            = length(find(floor(EcotranGroupType) == GroupTypeDef_ANYDetritus));
num_eggs                = length(find(EcotranGroupType  == GroupTypeDef_eggs));
num_EggsAndDetritus     = num_eggs + num_detritus;
num_micrograzers        = length(find(EcotranGroupType  == GroupTypeDef_micrograzers));
num_bacteria            = length(find(EcotranGroupType  == GroupTypeDef_bacteria));
num_PrimaryProducer     = length(find(fix(EcotranGroupType)  == GroupTypeDef_ANYPrimaryProd));

if num_EggsAndDetritus ~= num_detritus_EwE
    error('Number of ECOTRAN-defined detritus + eggs groups does not equal number of detritus groups provided on K-path Main tab')
end
% -------------------------------------------------------------------------


% step 2b: EwE Production rate (t/km2/y) for all living groups ------------
Production                  = (PB .* Biomass);
% -------------------------------------------------------------------------






% *************************************************************************
% STEP 3: create consumption matrix----------------------------------------
% step 3a: calculate consumption rate (t/km2/y) for all living groups -----
%           NOTE: consumption = mortality (1/y) * biomass (t/km2)
consumption                 = mortality; % initialize; (2D matrix: (num_LivingAndDetritus+1) X (num_LivingAndDetritus+1)); has a row & clm header
consumption(2:(num_living+1), 2:(num_living+1)) = mortality(2:(num_living+1), 2:(num_living+1)) .* repmat(Biomass(1:num_living), [1, num_living]); % consumption rate; (t/km2/y); (2D matrix: producers (num_LivingAndDetritus+1) X consumers (num_LivingAndDetritus+1))
% -------------------------------------------------------------------------


% step 3b: calculate production, feces, & senescence rates for all living groups ----
production              = Biomass(1:num_living) .* PB(1:num_living); % production by each living group; (vertical vector: num_living X 1)
total_consumption_temp	= Biomass(1:num_living) .* QB(1:num_living); % total consumption by each living group; (vertical vector: num_living X 1)
feces_production        = total_consumption_temp .* (1 - AE(1:num_living));       % feces production by each living group; (vertical vector: num_living X 1)
total_landings          = sum(landings(1:num_living, :), 2);                 % total landings from each living group by all fleets combined; (vertical vector: num_living X 1)
total_predation         = sum(consumption(2:(num_living+1), 2:(num_living+1)), 2); % total predation upon each living group; (vertical vector: num_living X 1)
senescence_production	= production - (total_predation + total_landings);           % (vertical vector: num_living X 1)
detritus_production     = feces_production + senescence_production;                  % detritus production by each living group; (vertical vector: num_living X 1)

total_consumed_detritus	= sum(consumption((num_living+2):(num_LivingAndDetritus+1), 2:(num_living+1)), 2); % detritus (incl. eggs) that is consumed by living groups; (vertical vector: num_EggsAndDetritus X 1)
% -------------------------------------------------------------------------


% step 3c: calculate flow to detritus -------------------------------------
%           flow from living; (t/km2/y)
%           flow from detritus groups is still 0
Flow2Detritus               = zeros(num_LivingAndDetritus, num_EggsAndDetritus); % initialize; (2D matrix: num_LivingAndDetritus X num_EggsAndDetritus)
Flow2Detritus(1:num_living, 1:num_EggsAndDetritus)	= EwE_DetritusFate(1:num_living, 1:num_EggsAndDetritus) .* repmat(detritus_production, [1, num_EggsAndDetritus]); % (2D matrix: num_LivingAndDetritus X num_EggsAndDetritus)
% -------------------------------------------------------------------------


% step 3d: calculate fleet discard ----------------------------------------
%           flow from each living group due to fleet discards; (t/km2/y)
discard_detritus            = zeros(num_LivingAndDetritus, num_EggsAndDetritus); % initialize; (2D matrix: num_LivingAndDetritus X num_EggsAndDetritus)

for fleet_loop = 1:num_gear
    current_Discard     = discards(:, fleet_loop);	% discard of each live group by the current fleet; (vertical vector: num_LivingAndDetritus X 1)
    repmat_Discard      = repmat(current_Discard, [1, num_EggsAndDetritus]);         % (3D matrix: num_LivingAndDetritus X num_EggsAndDetritus)
    
    current_DiscardFate	= EwE_DetritusFate((num_LivingAndDetritus + fleet_loop), :); % (horizontal vector: 1 X num_EggsAndDetritus); NOTE: fleets in rows of model groups must be in same order as across clms in discards;
    repmat_DiscardFate	= repmat(current_DiscardFate, [num_LivingAndDetritus, 1]);   % (3D matrix: num_LivingAndDetritus X num_EggsAndDetritus)

    discard_detritus	= discard_detritus + (repmat_Discard .* repmat_DiscardFate);    % summed discard rate across all fleets; (t/km2/y); (2D matrix: num_LivingAndDetritus X num_EggsAndDetritus)
end % fleet_loop
    
Flow2Detritus               = Flow2Detritus + discard_detritus;	% (2D matrix: num_LivingAndDetritus X num_EggsAndDetritus)
% -------------------------------------------------------------------------


% step 3e: calculate detritus flow to other detritus groups --------------
total_Flow2Detritus         = sum(Flow2Detritus, 1); % (horizontal vector: 1 X num_EggsAndDetritus)
surplus_detritus            = total_Flow2Detritus' - total_consumed_detritus; % (vertical vector: num_EggsAndDetritus X 1); NOTE: may have negative values at this point until detritus flow to detritus is included
DetritusToDetritusFate      = EwE_DetritusFate((num_living+1):num_LivingAndDetritus, :); % (2D matrix: num_EggsAndDetritus (source) X num_EggsAndDetritus (destiny)); NOTE: this assumes detritus grps are located after living groups and before fleets

DetritusToDetritus_array    = repmat(surplus_detritus, [1, num_EggsAndDetritus]) .* DetritusToDetritusFate; % flow between detritus groups; (t/km2/y); (2D matrix: num_EggsAndDetritus (source) X num_EggsAndDetritus (destiny))

% account for flow between detritus groups
sum_DetritusToDetritus      = sum(DetritusToDetritus_array);              % (horizontal vector: 1 X num_EggsAndDetritus (destiny))
sum_DetritusToDetritus(sum_DetritusToDetritus<0) = 0;                     % filter out negative values (negative flows will be error-checked below)
adjusted_surplus_detritus	= surplus_detritus + sum_DetritusToDetritus'; % (vertical vector: num_EggsAndDetritus (source) X 1)
DetritusToDetritus_array	= repmat(adjusted_surplus_detritus, [1, num_EggsAndDetritus]) .* DetritusToDetritusFate; % flow between detritus groups; (t/km2/y); (2D matrix: num_EggsAndDetritus (source) X num_EggsAndDetritus (destiny))
% -------------------------------------------------------------------------


% step 3f: put final consumption matrix together -------------------------
consumption(2:(num_LivingAndDetritus+1), (num_living+2):(num_LivingAndDetritus+1))              = Flow2Detritus;
consumption((num_living+2):(num_LivingAndDetritus+1), (num_living+2):(num_LivingAndDetritus+1))	= DetritusToDetritus_array;
% -------------------------------------------------------------------------


% step 3g: calculate total_consumption & import_consumption (total_consumption .* import_diet)
total_consumption                       = zeros(1, (num_LivingAndDetritus+1)); % initialize, add 1 cell up front for header
total_consumption(1)                    = NaN;
total_consumption(1, 2:end)           	= Biomass(1:num_LivingAndDetritus) .* QB(1:num_LivingAndDetritus); % ignore fisheries at the end
import_consumption                   	= total_consumption .* import_diet; % (import_consumption = total consumption .* import_diet) (horizontal vector; first cell is a header)
% *************************************************************************




% *************************************************************************
% STEP 4: some error-checking----------------------------------------------
% step 4a: check number of primary producers against nutrient uptake definitions

if length(PhytoUptake_NO3) ~= num_PrimaryProducer
    error('ECOTRAN setup ERROR: number of primary producers and NO3 uptake definitions do not match')
end

if length(PhytoUptake_NH4) ~= num_PrimaryProducer
    error('ECOTRAN setup ERROR: number of primary producers and NH4 uptake definitions do not match')
end

% -------------------------------------------------------------------------


% step 4b: check for negative metabolism ---------------------------------
BioenergeticBudget_metabolism = (AE - PQ); % fraction of consumption going to metabolism (NH4 production); Metabolism = 1 - PQ - feces = (1 - (EwE_PQ + (1-EwE_AE)))
looky_negativeMetabolism = find(BioenergeticBudget_metabolism <= 0);
if ~isempty(looky_negativeMetabolism)
    display('   -->EwE parameter WARNING: parameters imply negative or zero')
    display('      metabolism for at least 1 group (metabolism = ae - pq)')
    display(['      - bad groups: ' num2str(looky_negativeMetabolism(:)')])
end
% -------------------------------------------------------------------------


% % step 4c: check for ECOTRAN feces fate entries for all groups with feces
% FFF will need to sync rows of dat.AE and at.DetritusFate_feces for this test to work
% looky_feces_producers = find(dat.AE < 1);
% sum_FecesFate         = sum(dat.DetritusFate_feces, 2);
% FecesProduced         = sum_FecesFate(looky_feces_producers);
% looky_noFecesFate     = find(FecesProduced <= 0);
% if ~isempty(looky_noFecesFate)
%     error(['ECOTRAN setup error: Feces fate not defined for at least one consumer group; bad groups: ' num2str(looky_noFecesFate(:)')])
% end
% % -------------------------------------------------------------------------


% step 4d: check for full accounting of EwE detritus ----------------------
sum_EwE_DetritusFate                                = sum(EwE_DetritusFate, 2);
looky_terminalBenthicDetritus                       = find(EcotranGroupType == GroupTypeDef_terminalBnthDetr);
looky_terminalBenthicDetritus                       = looky_terminalBenthicDetritus - num_NitrogenNutrients;
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


% step 4e: check consumption matrix & balance between INFLOW & OUTFLOW ----
test_inflow     = sum(consumption(2:end, 2:end)); % consumption INTO each living and detritus group
test_outflow	= sum(consumption(2:end, 2:end), 2)'; % predation OUT of each living and detritus group; NOTE transpose
test_balance	= test_inflow - test_outflow;

[~, looky_NegativeConsumer] = find(consumption(2:end, 2:end) < 0);
if ~isempty(looky_NegativeConsumer)
    display('   -->consumption matrix WARNING: consumer(s) with negative consumption.')
    display(['      - bad groups: ' num2str(looky_NegativeConsumer)'])
end

% [~, looky_ConsumerImbalance] = find(test_balance < 0); % need to ignore primary producers for this test to work
% *************************************************************************





% *************************************************************************
% STEP 3: pack up model parameters-----------------------------------------
dat.fname_ReadEwE                      = fname_readKpath;

dat.num_EwEGroups                      = num_EwEGroups;
dat.num_EcotranGroups                  = num_EcotranGroups;
dat.num_NitrogenNutrients              = num_NitrogenNutrients;
dat.num_NonNitrogenNutrients           = num_NonNitrogenNutrients;
dat.num_living                         = num_living;
dat.num_PrimaryProducer                = num_PrimaryProducer;
dat.num_micrograzers                   = num_micrograzers;
dat.num_bacteria                       = num_bacteria;
dat.num_eggs                           = num_eggs;
dat.num_detritus                       = num_detritus;
dat.num_gear                           = num_gear;
dat.num_LivingAndDetritus              = num_LivingAndDetritus;
dat.num_EggsAndDetritus                = num_EggsAndDetritus;

dat.EwENumberCode                      = EwENumberCode;
dat.EcotranNumberCode                  = EcotranNumberCode;
dat.AggCode                            = AggCode;

dat.EwEName                            = EwEName;
dat.AggEwEName                         = AggEwEName;
dat.EcotranName                        = EcotranName;

dat.EcotranGroupType                   = EcotranGroupType;

dat.GroupTypeDef_ANYNitroNutr          = GroupTypeDef_ANYNitroNutr;
dat.GroupTypeDef_NO3                   = GroupTypeDef_NO3;
dat.GroupTypeDef_plgcNH4               = GroupTypeDef_plgcNH4;
dat.GroupTypeDef_bnthNH4               = GroupTypeDef_bnthNH4;
dat.GroupTypeDef_ANYPrimaryProd        = GroupTypeDef_ANYPrimaryProd;
dat.GroupTypeDef_LrgPhyto              = GroupTypeDef_LrgPhyto;
dat.GroupTypeDef_SmlPhyto              = GroupTypeDef_SmlPhyto;
dat.GroupTypeDef_Macrophytes           = GroupTypeDef_Macrophytes;
dat.GroupTypeDef_ANYConsumer           = GroupTypeDef_ANYConsumer;
dat.GroupTypeDef_ConsumPlgcPlankton    = GroupTypeDef_ConsumPlgcPlankton;
dat.GroupTypeDef_ConsumPlgcNekton      = GroupTypeDef_ConsumPlgcNekton;
dat.GroupTypeDef_ConsumPlgcWrmBlood	= GroupTypeDef_ConsumPlgcWrmBlood;
dat.GroupTypeDef_ConsumBntcInvert      = GroupTypeDef_ConsumBntcInvert;
dat.GroupTypeDef_ConsumBntcVert        = GroupTypeDef_ConsumBntcVert;
dat.GroupTypeDef_ConsumBnthWrmBlood    = GroupTypeDef_ConsumBnthWrmBlood;
dat.GroupTypeDef_eggs                  = GroupTypeDef_eggs;
dat.GroupTypeDef_ANYDetritus           = GroupTypeDef_ANYDetritus;
dat.GroupTypeDef_offal                 = GroupTypeDef_offal;
dat.GroupTypeDef_terminalPlgcDetr      = GroupTypeDef_terminalPlgcDetr;
dat.GroupTypeDef_terminalBnthDetr      = GroupTypeDef_terminalBnthDetr;
dat.GroupTypeDef_BA                    = GroupTypeDef_BA;
dat.GroupTypeDef_EM                    = GroupTypeDef_EM;
dat.GroupTypeDef_fishery               = GroupTypeDef_fishery;
dat.GroupTypeDef_import                = GroupTypeDef_import;
dat.GroupTypeDef_micrograzers          = GroupTypeDef_micrograzers;
dat.GroupTypeDef_bacteria              = GroupTypeDef_bacteria;

% EwE parameters
dat.TL                                 = TL;               % (vertical vector: num_grps X 1)
dat.Biomass                            = Biomass;        	% (vertical vector: num_grps X 1)
dat.PB                                 = PB;               % (vertical vector: num_grps X 1)
dat.QB                                 = QB;               % (vertical vector: num_grps X 1)
dat.PQ                                 = PQ;               % (vertical vector: num_grps X 1)
dat.AE                                 = AE;               % (vertical vector: num_grps X 1)
dat.EE                                 = EE;               % (vertical vector: num_grps X 1)
dat.DetritusImport                     = DetritusImport;	% (vertical vector: num_grps X 1)

dat.Production                         = Production;       % (vertical vector: num_grps X 1)

dat.diet                               = diet;             % (2D matrix: num_grps X num_grps)
dat.import_diet                        = import_diet;      % (horizontal vector: 1 X num_grps)

dat.EwE_DetritusFate                   = EwE_DetritusFate; % (2D matrix: num_grps X num_EggsAndDetritus)

dat.landings                           = landings;         % (2D matrix: num_LivingAndDetritus X num_gear)
dat.discards                           = discards;         % (2D matrix: num_LivingAndDetritus X num_gear)


% nutrient & detritus uptake & recycling parameters
dat.PelagicBacterialReduction          = PelagicBacterialReduction; % (scaler)
dat.BenthicBacterialReduction          = BenthicBacterialReduction; % (scaler)
dat.Oxidation_NH4                      = Oxidation_NH4;    % (vertical vector: 2 X 1)
dat.PhytoUptake_NH4                    = PhytoUptake_NH4;  % (vertical vector: num_PrimaryProducer X 1)
dat.PhytoUptake_NO3                    = PhytoUptake_NO3;  % (vertical vector: num_PrimaryProducer X 1)

% ECOTRAN recycling parameters
dat.BA                             	= BA;                           % (vertical vector: num_EcotranGroups X 1)
dat.EM                              	= EM;                           % (vertical vector: num_EcotranGroups X 1)
dat.EggProduction                  	= EggProduction;                % FFF for future; defined egg, gamete, or live young production (as opposed to use of eggs as EwE detritus method); (vertical vector: num_EcotranGroups X 1)
dat.Metabolism                       	= Metabolism;                   % FFF for future; defined metabolism (as opposed to derivation from EwE params; (vertical vector: num_EcotranGroups X 1)
dat.DetritusFate_feces               	= DetritusFate_feces;           % (2D matrix: num_EcotranGroups X 2)
dat.DetritusFate_senescence          	= DetritusFate_senescence;      % (2D matrix: num_EcotranGroups X 2)
dat.ExcretionFate                    	= ExcretionFate;                % (2D matrix: num_EcotranGroups X 2)
dat.ProductionLossScaler              	= ProductionLossScaler;         % (vertical vector: num_EcotranGroups X 1)
dat.RetentionScaler                  	= RetentionScaler;              % (vertical vector: num_EcotranGroups X 1)

% mortality & consumption parameters
dat.mortality                          = mortality;                    % (2D matrix: num_EwEGroups X num_EwEGroups)
dat.consumption                        = consumption;                  % (2D matrix: num_EwEGroups X num_EwEGroups)
dat.total_consumption                  = total_consumption;            % (horizontal vector: 1 X num_EwEGroups)
dat.import_consumption                 = import_consumption;           % (horizontal vector: 1 X num_EwEGroups)

% functional response terms
dat.FunctionalResponseParams           = FunctionalResponseParams;     % (2D matrix: num_EcotranGroups X 4)
dat.FunctionalResponse_matrix          = FunctionalResponse_matrix;    % (2D matrix: num_EcotranGroups X num_EcotranGroups)

% pedigree terms
dat.PreyGuild_name                     = PreyGuild_name;
dat.PreyGuild_code                     = PreyGuild_code;

dat.PedigreeParameters                 = PedigreeParameters;           % (2D matrix: num_LivingAndDetritus X 37)
dat.PedigreeLandings                 	= PedigreeLandings;             % (2D matrix: num_LivingAndDetritus X num_gear)
dat.PedigreeDiscards                 	= PedigreeDiscards;             % (2D matrix: num_LivingAndDetritus X num_gear)
dat.PedigreeDiet                     	= PedigreeDiet;                 % (2D matrix: num_EwEGroups X num_EwEGroups)
dat.import_PedigreeDiet               	= import_PedigreeDiet;          % (horizontal matrix: 1 X num_EwEGroups)
dat.PedigreeDietPrefRules            	= PedigreeDietPrefRules;        % (2D matrix: num_EwEGroups X num_EwEGroups)
dat.import_PedigreeDietPrefRules      	= import_PedigreeDietPrefRules; % (horizontal matrix: 1 X num_EwEGroups)
dat.Pedigree_PhysiologyScaler         	= Pedigree_PhysiologyScaler;    % (scaler)
dat.Pedigree_BiomassScaler           	= Pedigree_BiomassScaler;       % (scaler)
dat.Pedigree_FisheriesScaler         	= Pedigree_FisheriesScaler;     % (scaler)
dat.Pedigree_DietScaler              	= Pedigree_DietScaler;          % (scaler)
dat.Pedigree_BAScaler                 	= Pedigree_BAScaler;            % (scaler)
dat.Pedigree_EMScaler                 	= Pedigree_EMScaler;         	% (scaler)
dat.Pedigree_DiscardScaler           	= Pedigree_DiscardScaler;       % (scaler)
dat.Pedigree_minBiomass               	= Pedigree_minBiomass;          % (scaler)
dat.Pedigree_minPhysiology            	= Pedigree_minPhysiology;       % (scaler)
% *************************************************************************


% end function*************************************************************