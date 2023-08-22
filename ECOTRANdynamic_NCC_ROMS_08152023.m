% ECOTRANdynamic_NCC_ROMS_08152023
% run a dynamic model over time
% by Jim Ruzicka
%
% calls:
%       f_readEwEcsv_10pp_07072021                	read in ECOPATH (EwE) model from VisualBasic .csv file and store as variable 'dat'; (use for VisualBasic food web files with 10 primary producers)
%       f_AggregateBiologicalModel_02052021        	prepare EwE model for use by ECOTRAN; also, aggregate functional groups here if wanted
%           f_calcEE_12292020                          calculate Ecotrophic Efficiency
%           f_VarianceDivision_12132018                calculate the variance of one term divided by another term
%           f_VarianceMultiplication_12132018          calculate the variance of two products
%
%       ECOTRANheart_09032021                       generate an ECOTRAN (E2E) model; the heart of ECOTRAN
%           f_ECOfunction_09032021                     returns a single ECOTRAN model for 1 "type" EwE model or 1 MonteCarlo EwE model
%               f_RedistributeCannibalism_11202019       remove cannibalism terms on diagonal of matrix EwE_diet
%               f_calcEE_12292020                        calculate Ecotrophic Efficiency 
%               f_calcPredationBudget_12102019      	 for each producer, p, and consumer, c: ((b_pc * Q_c) / M2_p)) = the fraction of total predation on each producer, p, going to each consumer, c; (2D matrix: num_grps X num_grps; consumers X producers)
%
%       f_E2Epedigree_08042020                      calculate uncertainty terms for all elements of the ECOTRAN EnergyBudget (A_cp) as Coefficients of Variation (CV)
%           f_VarianceDivision_12132018                calculate the variance of one term divided by another term
%           f_VarianceMultiplication_12132018          calculate the variance of two products
%       f_E2E_MonteCarlo_08042020                   generate a set of randomly generated ECOTRAN models, multiple alternate versions of the EnergyBudget (A_cp); model 1 of the stack is the "type" model defined by the parameters of the VisualBasic .csv file
%
%       f_OrdinalDate                               calculate ordinal dates; day of year with January 1 of any year = 1
%
%       ALTERNATE PHYSICAL MODELS:
%       2D cross-shelf:
%           OPTIONAL METHOD                         (not included in this code package)
%
%       3D ROMS:
%             f_ECOTRANphysics_NCC_ROMS_08152023
%               f_ROMS_GridPrep_NCC_08152023            Map ROMS grid to ECOTRAN grid
%               f_ROMS_FluxPrep_NCC_08152023            express ROMS fluxes in DESTINY<--SOURCE format on ECOTRAN grid
%               f_CompactFluxTimeSeries_11182019        compact physical flux time-series, arranged as 3D matrix (time X source box X destiny box), eliminate source & destiny infomration for any subregions that never communicate
%               f_UnCompactFluxTimeSeries_12112019      UnCompact a previously compacted physical flux time-series to provide IMPORT & EXPORT fluxes for each box and the domain as a whole
%                   f_calcNetFlux_12112019                  calculate net flux into and net flux out of each model box and across outer domain boundaries
%               f_EvaluateFluxBalance_11262021          examine for flux time-series for imbalances IN & OUT of individual boxes and IN & OUT of the overall domain
%               f_LightIntensity_12112020               instantaneous (W/m2), daily mean averaged across 24 h (W m^-2 h^-1), & daily integrated (W m^-2 d^-1) solar raditation at ocean surface; (vertical vector: num_t X 1)
%
%       f_ECOTRANmigration_NCC_02222022             area of overlap of neighboring model sub-domains (for ROMS); (not used by this code package but null variables are still needed for code to run)
%       f_CompactFluxTimeSeries_11182019            compact physical flux time-series, arranged as 3D matrix (time X source box X destiny box), eliminate source & destiny infomration for any subregions that never communicate
%       f_DVMsinusoid_08122021                      Calculate DVM flux rates between all model domain boxes for each functional group at each time point; Uses a daily sinusoidal migration pattern
%
%       f_FunctionalResponse_MonteCarlo_05132021    prepare array of vulnerability terms and allows for random generation of functional response terms within a predefined uncertainty level
%
%       TEMPERATURE RESPONSE PHYSIOLOGY FUNCTIONS:
%           f_physiology_Q10_12082022               Calculate temperature-dependent Q10 metabolic rate scaling factors for living groups
%           f_TLparameterization_08182023           Prepare the eight required Thornton-Lessem parameters for each living group
%           f_ThorntonLessem_12012022               Calculate temperature-dependent Thornton-Lessem temperature scaling factors for ingestion rate
%
%       f_InitialProductionRates_02012022          	calculate initial or mean production conditions
%           f_WebProductivity_03272019                  calculate production rates of all groups under a given driver (e.g., NO3 or primary production); also accounts for defined rates of group production export when running static scenarios
%
%       ODE SOLVER:%
%           f_PrepMexODE_09192022                       prepare ECOTRAN variables for using C++ ODE solver mex function for non-shelf 2D & 3D ROMS cases. Pack parameters & drivers along proper dimensions
%               f_unspoolMATRIX_04282020                    linearize multidimenional matrices up to 4-D for use in C++
%           mex_ECOTRANode_09182022                  solve the ecosystem ODE in C++ for 2D NON-shelf & 3D ROMS cases
%
% returns:
%       store_T                     time-series of day umbers
%       store_ProductionRates       time-series of production rates for each functional group; (t/km2/d) 
%
% NOTE: code cannot currently accomodate seasonal changes in flows to non-terminal detritus pools (e.g., fisheries and detritus columns of EnergyBudget)
%
% revision date: 8-15-2023


%% *************************************************************************
% STEP 1: load & aggregate EwE results-------------------------------------
% step 1a: set operating conditions ---------------------------------------
fname_ECOTRANdynamic	= 'ECOTRANdynamic_NCC_ROMS_08152023'; % save name of this m-file to keep in saved model results

switch_MonteCarlo           = 'MonteCarlo_build';	% generate (and optionally save) a stack of MonteCarlo food webs
% switch_MonteCarlo           = 'MonteCarlo_load';	% load a saved stack of MonteCarlo food webs
% switch_MonteCarlo           = 'MonteCarlo_TypeModel';	% use NO MonteCarlo food webs

switch_FunctionalResponse	= 'NonLinear_default';	% NonLinear_default functional response
% switch_FunctionalResponse	= 'Linear';             % linear functional response; NOTE: STRICTLY DONER-DRIVEN DYNAMICS (rate of consumption by each consumer is a direct proportion of the production by each of its prey groups)
% switch_FunctionalResponse	= 'NonLinear_alt';      % alternate NonLinear functional response (user must define)

switch_Q10                  = 'Q10_ON';
% switch_Q10                  = 'Q10_OFF';

switch_ThorntonLessem       = 'ThorntonLessem_ON';
% switch_ThorntonLessem       = 'ThorntonLessem_OFF';

switch_SubModel             = 'identical_SubModel';                 % OPTION 1: use the same food web across the shelf
% switch_SubModel             = 'independent_SubModel';               % OPTION 2: use independently defined webs for each shelf zone
% -------------------------------------------------------------------------


% step 1b: define food web model to use -----------------------------------
ReadFile_directory      = '/Volumes/Fortress/Build_Package/For_GitHub/FoodWeb_models/';	% SSS change to local directory

BiologicalModel_name	= 'NCC_11242020.csv'; % Jim's pre-heatwave model

readFile            	= [ReadFile_directory BiologicalModel_name];
% -------------------------------------------------------------------------


% step 1c: load ECOPATH (EwE) model from Aydin VisualBasic file (.csv format)
dat                  	= f_readEwEcsv_10pp_07072021(readFile);	% use for models with up to 10 primary producers
% -------------------------------------------------------------------------


% step 1d: aggregate model results & prep EwEResult for analysis ----------
[EwEResult, PEDIGREE] 	= f_AggregateBiologicalModel_02052021(dat);
% -------------------------------------------------------------------------


% step 1e: define filename and directory for saving results ---------------
SaveFile_directory      = '/Volumes/Fortress/Build_Package/For_GitHub/Result_files/';	% SSS change to local directory

SaveFile_label          = 'NCC_TestRun_';
% *************************************************************************





%% ************************************************************************
% STEP 2: ECOTRAN conversion-----------------------------------------------
MonteCarloStore         = [];
[ECOTRAN]            	= ECOTRANheart_09032021(EwEResult, MonteCarloStore);
% *************************************************************************





%% *************************************************************************
% STEP 3: Generate E2E Monte Carlo models based on ECOTRAN EnergyBudget----
%         Start with the one original ECOTRAN base model and generate a set of 
%           Monte Carlo models from the ECOTRAN EnergyBudget & ConsumptionBudget
%           matrices using predefined CV values

switch switch_MonteCarlo

	case 'MonteCarlo_build' % generate (and save) a stack of Monte Carlo models
    
        num_MC              = 3;        % SSS set this value
        
        disp(['MonteCarlo: building stack of ' num2str(num_MC) ' food webs'])
        
        ECOTRAN.num_MC      = num_MC;

        PEDIGREE.ee_eggs_CV                               = 0.01; % SSS egg pedigree W.R.T. production budget for all groups; (CV); (scaler)
        PEDIGREE.BacterialMTBLSM_CV                       = 0.01; % SSS pedigree for implicit bacterial metabolism of terminal detritus (CV)
        PEDIGREE.Oxidation_NH4_CV                         = 0.01;	% fraction of NH4 produced oxidized directly back to NO3 abiologically; QQQ scaler?? (vertical vector: num_NH4 X 1)??
        PEDIGREE.NutrientUptake_CV                        = 0.01; % SSS pedigree for nutrient uptake by primary producers (CV)
        ECOTRAN_PEDIGREE                               	  = f_E2Epedigree_08042020(ECOTRAN, PEDIGREE); % NEW!!!

                % SSS use for standardized pedigree
                %     overwrite the pedigree values from the ECOPATH (EwE) model from VisualBasic file (.csv format)
                [rows, clms]                                = size(ECOTRAN_PEDIGREE.EnergyBudget_CV);
                ECOTRAN_PEDIGREE.EnergyBudget_CV            = 0.001 * ones(rows, clms); % changes how important predators are relative to eachother
                ECOTRAN_PEDIGREE.ConsumptionBudget_CV       = zeros(7, clms);
                ECOTRAN_PEDIGREE.ConsumptionBudget_CV(1, :) = 0.05; % feces
                ECOTRAN_PEDIGREE.ConsumptionBudget_CV(2, :) = 0.05; % metabolism
                ECOTRAN_PEDIGREE.ConsumptionBudget_CV(3, :) = 0.05; % eggs
                ECOTRAN_PEDIGREE.ConsumptionBudget_CV(4, :) = 0.05;	% predation
                ECOTRAN_PEDIGREE.ConsumptionBudget_CV(5, :) = 0.05;  % senescence
                ECOTRAN_PEDIGREE.ConsumptionBudget_CV(6, :) = 0.05;  % ba
                ECOTRAN_PEDIGREE.ConsumptionBudget_CV(7, :) = 0.05;  % em
                ECOTRAN_PEDIGREE.DiscardFraction_CV         = 0.05 * ECOTRAN_PEDIGREE.DiscardFraction_CV; % QQQ reduce DiscardFraction_CV

        MonteCarloConditions.num_MC               	    = num_MC;	% number of random Monte Carlo models to generate
        MonteCarloConditions.DistributionType         	= 'normal';	% SSS distribution to draw random values from ('normal' or 'uniform'); (NOTE: code not fully proofed for uniform)

        ECOTRAN_MC                                      = f_E2E_MonteCarlo_08042020(MonteCarloConditions, ECOTRAN, ECOTRAN_PEDIGREE);
        EnergyBudget_MC                                 = ECOTRAN_MC.EnergyBudget_MC;
        ConsumptionBudget_MC                            = ECOTRAN_MC.ConsumptionBudget_MC;
        DiscardFraction_MC                              = ECOTRAN_MC.DiscardFraction_MC;
        ECOTRAN_MC.num_MC                               = num_MC;

%         % activate to save this stack of Monte Carlo food webs
%         filename_MC      = [SaveFile_directory 'MonteCarlo_NCC_stack_' date '.mat'];
%         disp(['MonteCarlo: SAVING stack: ' filename_MC])
%         save(filename_MC, 'ECOTRAN_MC')

    % end (case 'build_MonteCarlo') --------------

    
    case 'MonteCarlo_load' % load a set of Monte Carlo models
        filename_MC = [SaveFile_directory 'MonteCarlo_NCC_stack_27-Jun-2022.mat']; % SSS be sure to give correct saved file name here
        disp(['MonteCarlo: LOADING stack: ' filename_MC])
        load(filename_MC, 'ECOTRAN_MC')
        
        num_MC                            	= ECOTRAN_MC.num_MC;
        EnergyBudget_MC                     = ECOTRAN_MC.EnergyBudget_MC;
        ConsumptionBudget_MC                = ECOTRAN_MC.ConsumptionBudget_MC;
        DiscardFraction_MC                  = ECOTRAN_MC.DiscardFraction_MC;
    % end (case 'load_MonteCarlo') --------------
        
    
    case 'MonteCarlo_TypeModel'
        
        disp('MonteCarlo: using the defining TypeModel')
    
        num_MC                          = 1;       % only the "type" model is used
        EnergyBudget_MC                 = ECOTRAN.EnergyBudget;
        ConsumptionBudget_MC            = ECOTRAN.ConsumptionBudget;
        DiscardFraction_MC              = ECOTRAN.DiscardFraction;

        ECOTRAN_MC.num_MC               = num_MC;
        ECOTRAN_MC.EnergyBudget_MC      = ECOTRAN.EnergyBudget; % (3D matrix: num_grps (consumers) X num_grps (producers) X num_MC (1))
        ECOTRAN_MC.ConsumptionBudget_MC	= ECOTRAN.ConsumptionBudget; % (3D matrix: 7 X num_grps (producers) X num_MC (1))
        ECOTRAN_MC.DiscardFraction_MC	= ECOTRAN.DiscardFraction;
        
        ECOTRAN_MC.fate_metabolism      = ECOTRAN.fate_metabolism;	% (3D matrix: num_nutrients X num_grps (producers) X num_MC (1))
        ECOTRAN_MC.fate_eggs            = ECOTRAN.fate_eggs;        % (3D matrix: num_eggs X num_grps (producers) X num_MC (1))
        ECOTRAN_MC.fate_feces           = ECOTRAN.fate_feces;       % (3D matrix: num_ANYdetritus X num_grps (producers) X num_MC (1))
        ECOTRAN_MC.fate_senescence      = ECOTRAN.fate_senescence;  % (3D matrix: num_ANYdetritus X num_grps (producers) X num_MC (1))
        ECOTRAN_MC.fate_predation       = ECOTRAN.fate_predation;   % (3D matrix: num_livingANDfleets X num_grps (producers) X num_MC (1))
	% end (case 'MonteCarlo_TypeModel') --------------

end % (switch_MonteCarlo) -------------------------------------------------
% *************************************************************************





%% ************************************************************************
% STEP 4: prep and pack ECOTRAN model parameters---------------------------

% step 4a: read in ECOTRAN structure variables ----------------------------
%          (so that no changes are made to original values)
GroupType                           = ECOTRAN.GroupType;
label                            	= ECOTRAN.label;
% EnergyBudget_MC                 	  = ECOTRAN.EnergyBudget;
biomass                          	= ECOTRAN.biomass;              % (vertical vector: num_grps X 1); note inclusion of separately constructed regional models
pb                               	= ECOTRAN.pb;                   % (vertical vector: num_grps X 1)
qb                               	= ECOTRAN.qb;                   % (vertical vector: num_grps X 1)
fate_feces                       	= ECOTRAN.fate_feces;
fate_metabolism                  	= ECOTRAN.fate_metabolism;
fate_eggs                        	= ECOTRAN.fate_eggs;
fate_senescence                  	= ECOTRAN.fate_senescence;
ProductionLossScaler             	= ECOTRAN.ProductionLossScaler;	% (vertical vector: num_grps X 1)
RetentionScaler                 	= ECOTRAN.RetentionScaler;      % sensitivity to advection & mixing (0 = more advection <--> less advection =1); (vertical vector: num_grps X 1)
FunctionalResponseParams         	= ECOTRAN.FunctionalResponseParams;
num_grps                            = ECOTRAN.num_grps;             % number of model groups
% num_MC                              = ECOTRAN.num_MC;               % number of Monte Carlo models
% TransferEfficiency                  = ECOTRAN.TransferEfficiency;	  % gets redefined manually below
% -------------------------------------------------------------------------


% step 4b: find detritus, nutrients, ba & em ------------------------------
%           row addresses in EnergyBudget_MC
%           NOTE: ba = biomass accumulation term, em = emigration term

looky_NO3                       	= find(GroupType        == ECOTRAN.GroupTypeDef_NO3);
looky_plgcNH4                       = find(GroupType        == ECOTRAN.GroupTypeDef_plgcNH4);
looky_bnthNH4                       = find(GroupType        == ECOTRAN.GroupTypeDef_bnthNH4);
looky_NH4                           = find(GroupType        == ECOTRAN.GroupTypeDef_plgcNH4 | GroupType == ECOTRAN.GroupTypeDef_bnthNH4);
looky_nutrients                     = find(floor(GroupType)	== ECOTRAN.GroupTypeDef_ANYNitroNutr);	% row addresses of nutrients
looky_ANYPrimaryProducer        	= find(floor(GroupType) == ECOTRAN.GroupTypeDef_ANYPrimaryProd);
looky_macroalgae                  	= find(GroupType        == ECOTRAN.GroupTypeDef_Macrophytes);
looky_phytoplankton                 = looky_ANYPrimaryProducer(~ismember(looky_ANYPrimaryProducer, looky_macroalgae));
looky_ANYconsumer                 	= find(floor(GroupType) == ECOTRAN.GroupTypeDef_ANYConsumer);
% looky_PLGCbacteria                  = find(GroupType        == ECOTRAN.GroupTypeDef_plgcBacteria);
% looky_BNTHbacteria                  = find(GroupType        == ECOTRAN.GroupTypeDef_bnthBacteria);
% looky_ANYbacteria                   = find(GroupType        == ECOTRAN.GroupTypeDef_plgcBacteria | GroupType == ECOTRAN.GroupTypeDef_bnthBacteria);
looky_micrograzers                  = find(GroupType        == ECOTRAN.GroupTypeDef_micrograzers);
looky_bacteria                      = find(floor(GroupType) == ECOTRAN.GroupTypeDef_bacteria);
looky_eggs                          = find(GroupType        == ECOTRAN.GroupTypeDef_eggs);
looky_ANYdetritus                   = find(floor(GroupType) == ECOTRAN.GroupTypeDef_ANYDetritus);
looky_terminalPLGCdetritus          = find(GroupType        == ECOTRAN.GroupTypeDef_terminalPlgcDetr);
looky_terminalBNTHdetritus          = find(GroupType        == ECOTRAN.GroupTypeDef_terminalBnthDetr);
looky_eggsANDdetritus               = sort([looky_ANYdetritus; looky_eggs]);
looky_livingANDdetritus             = sort([looky_ANYPrimaryProducer; looky_ANYconsumer; looky_bacteria; looky_eggsANDdetritus]);
looky_terminalANYdetritus           = find(GroupType == ECOTRAN.GroupTypeDef_terminalPlgcDetr | GroupType == ECOTRAN.GroupTypeDef_terminalBnthDetr);
looky_fleets                        = find(floor(GroupType) == ECOTRAN.GroupTypeDef_fleet);
looky_livingANDfleets               = [looky_ANYPrimaryProducer; looky_ANYconsumer; looky_bacteria; looky_fleets]; % includes primary producers & bacteria
looky_NONnutrients                  = sort([looky_livingANDdetritus; looky_fleets]);	% addresses of all groups EXCEPT nutrients (needed to append nutrients)
looky_nonNO3                     	= 1:num_grps;
looky_nonNO3(looky_NO3)          	= [];

num_nutrients                   	= length(looky_nutrients);
num_NO3                             = length(looky_NO3);
num_NH4                          	= length(looky_NH4);
num_plgcNH4                         = length(looky_plgcNH4); % FFF replace mat_num_plgcNH4 & mat_num_bnthNH4 with simple mat_num_NH4
num_bnthNH4                         = length(looky_bnthNH4); % FFF replace mat_num_plgcNH4 & mat_num_bnthNH4 with simple mat_num_NH4
num_ANYPrimaryProd                	= length(looky_ANYPrimaryProducer);
num_phytoplankton                	= length(looky_phytoplankton);
num_macroalgae                   	= length(looky_macroalgae);
num_ANYconsumers                	= length(looky_ANYconsumer);
num_fleets                          = length(looky_fleets);
num_predators                       = num_ANYconsumers + num_fleets;
num_livingANDfleets                 = length(looky_livingANDfleets);
num_eggs                         	= length(looky_eggs);
num_ANYdetritus                   	= length(looky_ANYdetritus);
% -------------------------------------------------------------------------


% step 4c: pack variables for ODE solver ----------------------------------
ODEinput.looky_nutrients         	= looky_nutrients;
ODEinput.looky_NO3                	= looky_NO3;
ODEinput.looky_NH4                  = looky_NH4;
ODEinput.looky_plgcNH4            	= looky_plgcNH4;
ODEinput.looky_bnthNH4            	= looky_bnthNH4;
ODEinput.looky_ANYPrimaryProducer	= looky_ANYPrimaryProducer;
ODEinput.looky_phytoplankton      	= looky_phytoplankton;
ODEinput.looky_macroalgae       	= looky_macroalgae;
ODEinput.looky_fleets             	= looky_fleets;
ODEinput.looky_ANYconsumer       	= looky_ANYconsumer;
ODEinput.looky_livingANDfleets      = looky_livingANDfleets;
ODEinput.looky_eggs             	= looky_eggs;
ODEinput.looky_terminalPLGCdetritus	= looky_terminalPLGCdetritus;
ODEinput.looky_terminalBNTHdetritus	= looky_terminalBNTHdetritus;
ODEinput.looky_ANYdetritus      	= looky_ANYdetritus;
ODEinput.looky_nonNO3             	= looky_nonNO3;

ODEinput.num_grps                   = num_grps;
ODEinput.num_nutrients              = num_nutrients;
ODEinput.num_NO3                    = num_NO3;
ODEinput.num_NH4                    = num_NH4;
ODEinput.num_plgcNH4              	= num_plgcNH4; % FFF replace mat_num_plgcNH4 & mat_num_bnthNH4 with simple mat_num_NH4
ODEinput.num_bnthNH4              	= num_bnthNH4; % FFF replace mat_num_plgcNH4 & mat_num_bnthNH4 with simple mat_num_NH4
ODEinput.num_ANYPrimaryProd     	= num_ANYPrimaryProd;
ODEinput.num_phytoplankton       	= num_phytoplankton;
ODEinput.num_macroalgae          	= num_macroalgae;
ODEinput.num_ANYconsumers           = num_ANYconsumers;
ODEinput.num_predators           	= num_predators;
ODEinput.num_livingANDfleets        = num_livingANDfleets;
ODEinput.num_eggs                 	= num_eggs;
ODEinput.num_ANYdetritus            = num_ANYdetritus;
% -------------------------------------------------------------------------


% step 4e: matrix addresses of ECOTRAN generic model groups ---------------
%           These variables are simply used as references to the row &
%              column positions of each functional group in trophic matrix
%              A_cg.

%          NCC_11242020
rc_NO3                      = 1;
rc_plgcNH4                  = 2;
rc_bnthNH4                  = 3;
rc_lrg_phyto                = 4;
rc_sml_phyto                = 5;
rc_microZoop                = 6;
rc_lrg_copepods          	= 7;
rc_sml_copepods             = 8;
rc_sml_invertLarvae         = 9;
rc_pteropods                = 10;
rc_plgc_amphipods           = 11;
rc_plgcShrimp               = 12; % (Sergestidae Panaeidae)
rc_macroZoop                = 13;
rc_jellyfish_netFeeders     = 14;
rc_jellyfish_carnivores     = 15;
rc_lrg_jellyfish            = 16;
rc_Epacifica                = 17;
rc_Tspinifera               = 18;
rc_smallSquid               = 19;
rc_humboldtSquid            = 20;
rc_smelt                    = 21;
rc_shad                     = 22;
rc_sardine                  = 23;
rc_herring                  = 24;
rc_anchovy                  = 25;
rc_saury                    = 26;
rc_cohoYrlg_wild            = 27;
rc_cohoYrlg_hatchery        = 28;
rc_ChinookSubYrlg_wild      = 29;
rc_ChinookSubYrlg_hatchery  = 30;
rc_ChinookYrlg_wild         = 31;
rc_ChinookYrlg_hatchery     = 32;
rc_otherJuvSalmon           = 33;
rc_mesopelagicFish          = 34;
rc_planktivorous_rockfish   = 35;
rc_coho                     = 36;
rc_Chinook                  = 37;
rc_otherSalmon              = 38;
rc_sharks                   = 39;
rc_jackMackerel             = 40;
rc_PacificMackerel          = 41;
rc_piscivorousRockfish      = 42;
rc_dogfish                  = 43;
rc_hake                     = 44;
rc_tuna                     = 45;
rc_sablefish                = 46;
rc_Hexagrammidae            = 47;
rc_flatfish_plgcFeeders     = 48;
rc_skates_rays              = 49;
rc_sml_bnthFish             = 50;
rc_benthivorousRockfish     = 51;
rc_Gadidae                  = 52;
rc_flatfish_bnthFeeders     = 53;
rc_sml_flatfish             = 54;
rc_grenadier                = 55;
rc_juv_rockfish             = 56;
rc_juv_fish_other           = 57;
rc_juv_fish_chondrichthys   = 58;
rc_infauna                  = 59;
rc_Pandalus                 = 60;
rc_epibenthic_shrimp        = 61; % (Caridea)
rc_mysids                   = 62;
rc_echinoderms              = 63;
rc_bnth_amphipods           = 64; % (isopods & cumaceans)
rc_bivalves                 = 65;
rc_epifauna_suspnsnFeeders  = 66;
rc_DungenessCrab            = 67;
rc_TannerCrab               = 68;
rc_epifauna_carnivores      = 69;
rc_sootyShearwaters         = 70;
rc_commonMurre              = 71;
rc_gulls_terns              = 72;
rc_alcids                   = 73;
rc_lrg_plgcSeabirds         = 74; % (albatross, jaegers, 
rc_other_plgcSeabirds       = 75;
rc_coastalSeabirds_divers   = 76; % (cormorants)
rc_stormPetrels             = 77;
rc_grayWhales               = 78;
rc_baleenWhales             = 79;
rc_sml_pinnipeds            = 80;
rc_lrg_pinnipeds            = 81;
rc_sml_odontocetes          = 82;
rc_lrg_odontocetes          = 83;
rc_orcas                    = 84;
rc_invert_eggs              = 85;
rc_fish_eggs                = 86;
rc_plgc_detritus            = 87;
rc_fishery_offal            = 88;
rc_bnth_detritus            = 89;
rc_fleet_longline           = 90;
rc_fleet_troll              = 91;
rc_fleet_bottomfish_troll   = 92;
rc_fleet_hook_line          = 93; % (except troll)
rc_fleet_offshore_hook_line = 94; % (except troll)
rc_fleet_trawls             = 95; % (except shrimp trawls)
rc_fleet_shrimp_trawls      = 96;
rc_fleet_plgc_netGear       = 97; % (except trawls)
rc_fleet_midWater_trawls    = 98;
rc_fleet_gillNets           = 99;
rc_fleet_seine              = 100;
rc_fleet_fishPot            = 101;
rc_fleet_crabPot            = 102;
rc_fleet_other_pot          = 103; % (and trap gear)
rc_fleet_diving             = 104;
rc_fleet_other_known_gear   = 105;
rc_fleet_recreational       = 106;

% -------------------------------------------------------------------------


% step 4f: define functional group life-history, physiology, DVM (diel vertical migration), and other model tuning parameters --------------------
%          1) set depth distribution parameters HERE (for vertically-integrated food web technique)
%          2) set particle sinking speeds in STEP 6f
%          3) set Diel Vertical Migration (DVM) speeds in STEP 6g (not applied in CAFA ecophysiology project, and all are set to 1)
%          4) set production loss fractions in STEP 6h; NOTE: ALL USUALLY SET TO ZERO
%          5) set physiology parameters in STEP 10 (Q10 and Thornton-Lessem

% DEPTH DISTRIBUTION TERMS HERE (to define water-column positions and assignment of Thornton-Lessem temperature response parameter sets)
depthdist = [0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 ...
             1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 2 2 2 1 3 3 2 3 3 ...
             3 3 3 3 3 1 1 2 3 3 3 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
             1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; % 0 = no depth assignment, 1 = surface, 2 = mid-water, 3 = benthic;
 % *************************************************************************





%% *************************************************************************
% STEP 5: define TransferEfficiency terms----------------------------------
%         SET MANUALLY: TransferEfficiency = 1 for all groups because we are
%                       working with the consumption budget matrix that tracks the fate of ALL consumption (not just predation)
%                       system losses (groups where TransferEfficiency < 1 include benthic detritus, fishery, (others?)
%         NOTE: terminal benthic detritus can be removed via TE < 1 or by a tuned sequestration (senescence) rate in STEP 7
TransferEfficiency      = ones(num_MC, num_grps);	% re-initialize as rows of horizontal vector of ones
% *************************************************************************




% NOTE: DEACTIVATE THIS FULL VERSION OF STEP 6 IF LOADING PRE-AGGREGATED ROMS OUTPUT USING THE "OPTIONAL STEP 6 version" (below)
% % % %% *************************************************************************
% % % % STEP 6: prepare physical parameters--------------------------------------
% % % 
% % % % step 6a: Define time-frame and time-step info
% % % disp('--->>> 3D ROMS physics')
% % % 
% % % % SSS -->define directories holding ROMS files
% % % PHYSICSinput.ROMSfile_directory     = '/Volumes/Fortress/Build_Package/For_GitHub/ROMS_driverFiles/'; % directory to ROMS time-series
% % % 
% % % % ROMS product netcdf files
% % % PHYSICSinput.filename_list          = {
% % %                                        'wc12_avg_gfdl_1980_trimmed.nc'
% % %                                       }; % list of ROMS netcdf files to read (one file for each year)
% % %                                   
% % % % files needed to re-map ROMS to ECOTRAN grid
% % % PHYSICSinput.readFile_DepthLevels	= [PHYSICSinput.ROMSfile_directory 'depth_levels_trimmed.nc'];
% % % PHYSICSinput.readFile_grid        	= [PHYSICSinput.ROMSfile_directory 'wc12_grd.nc'];
% % % PHYSICSinput.readFile_ExampleYear	= [PHYSICSinput.ROMSfile_directory PHYSICSinput.filename_list{1}]; % 2008 time-series (used wc12_avg_2008_trimmed.nc)
% % % 
% % % datestart                       = datenum('01-Jan-1980'); % SSS --> enter starting date
% % % dateend                         = datenum('31-Dec-1980'); % SSS --> enter ending date (default for dynamic runs tests ('31-Dec-2020'))
% % % dt                              = 24/24; % t-step; (days); (dt = 24/24 = 1 d; dt = 3/24 = 3 hours)
% % %                                   % NOTE: take care to select good dt values for diel vertical migration 
% % %                                   %       (other values do not scale well between 1 & -1 in sin diel cycle (probably due to rounding error of pi() function)
% % % 
% % % datestart_OrdinalDate           = f_OrdinalDate(datestart);
% % % min_t                           = datestart_OrdinalDate;
% % % max_t                           = (dateend - datestart) + datestart_OrdinalDate;
% % % t_grid                          = linspace(min_t, (max_t+1-dt), ((dateend - datestart + 1)/dt))'; % QQQ NEW VERSION!!!; (vertical vector); t_grid runs from datestart to dateend and is inclusive of dateend; intervals = dt
% % % num_t                           = length(t_grid); % length of t_grid; (scaler)
% % % PHYSICSinput.datestart          = datestart;
% % % PHYSICSinput.dateend            = dateend;
% % % PHYSICSinput.dt                 = dt;
% % % PHYSICSinput.t_grid             = t_grid;
% % % calendar_time                   = datevec(datestart + t_grid - 1); % the actual calendar time for t_grid
% % % 
% % % PHYSICSinput.t_grid_real     	= linspace(datestart, (dateend+1-dt), ((dateend - datestart + 1)/dt))'; % QQQ NEW VERSION!!!; (vertical vector); t_grid runs from datestart to dateend and is inclusive of dateend; intervals = dt
% % % % -------------------------------------------------------------------------
% % % 
% % % 
% % % % step 6b: prepare advection & mixing time-series for each model box ------
% % % % 3D ROMS driver
% % % ECOTRANphysics                  = f_ECOTRANphysics_NCC_ROMS_08152023(PHYSICSinput); % use for NCC ROMS
% % % % -------------------------------------------------------------------------
% % % 
% % % 
% % % % step 6c: migration of each group (shared box face areas) ----------------
% % % %          NOTE: code does not handle migration flux uncertainty (nor physical flux uncertainty)
% % % ECOTRANmigration                = f_ECOTRANmigration_NCC_02222022(ECOTRANphysics);  % SSS; use for NCC ROMS
% % % 
% % % % ODEinput.biomass_migrator       = ECOTRANmigration.biomass_migrator;  % SSS special definition of boundary biomasses for migrators; (mmoles N/m3); (3D matrix: num_t X num_grps X num_boxes); NOTE: de-comment migrator biomass lines within ODE code
% % % % -------------------------------------------------------------------------
% % % 
% % % 
% % % % step 6d: unpack & process physics variables -----------------------------
% % % CompactFlux_ADVECTION           = ECOTRANphysics.CompactFlux_ADVECTION; % (structure)
% % % CompactFlux_HORIZONTALMIXING	= ECOTRANphysics.CompactFlux_HORIZONTALMIXING; % (structure)
% % % CompactFlux_VERTICALMIXING      = ECOTRANphysics.CompactFlux_VERTICALMIXING; % (structure)
% % % CompactFlux_SINKING             = ECOTRANphysics.CompactFlux_SINKING; % (structure); compact SINKING as box floor areas and connectivity information; apply functional group sinking speeds in ECOTRANdynamic_ code
% % % CompactFlux_MIGRATION           = ECOTRANmigration.CompactFlux_MIGRATION; % (structure)
% % % 
% % % num_boxes                       = ECOTRANphysics.num_boxes;
% % % num_domains                     = ECOTRANphysics.num_domains;	% number of geographic domains (does not consider depth layers)
% % % num_z                           = ECOTRANphysics.num_z;         % number of depth layers
% % % 
% % % BoxVolume                    	= ECOTRANphysics.BoxVolume;     	% (m3); (2D matrix: num_t X num_boxes)
% % % % BoxLength                    	= ECOTRANphysics.BoxLength;        	% (m); (2D matrix: num_t X num_boxes)
% % % % BoxHeight                  	= ECOTRANphysics.BoxHeight;       	% (m); (2D matrix: num_t X num_boxes)
% % % % BoxWidth                     	= ECOTRANphysics.BoxWidth;       	% (m); (2D matrix: num_t X num_boxes)
% % % 
% % % % ROMS BioGeoChemical (BGC) model info
% % % ROMS_temperature_initial     	= ECOTRANphysics.ROMS_temperature_initial;          % (deg C); (horizontal vector: 1 X num_boxes)
% % % ROMS_diatom                     = ECOTRANphysics.diatom_timeseries;                 % (mmole N/m3); (2D matrix: num_t X num_boxes)
% % % ROMS_nanophytoplankton          = ECOTRANphysics.nanophytoplankton_timeseries;      % (mmole N/m3); (2D matrix: num_t X num_boxes)
% % % ROMS_diatom_initial          	= ECOTRANphysics.ROMS_diatom_initial;               % (mmoles N/m3); (horizontal vector: 1 X num_boxes)
% % % ROMS_nanophytoplankton_initial	= ECOTRANphysics.ROMS_nanophytoplankton_initial;	% (mmoles N/m3); (horizontal vector: 1 X num_boxes)
% % % 
% % % % light intensity parameters (not actually used in CAFA project analyses)
% % % Io                            	= ECOTRANphysics.Io;             	% time-series of surface PAR light intensity; daily integrated solar raditation at ocean surface averaged across 24 hours; (W m^-2 h^-1); (vertical vector: num_t X 1)
% % % current_light                   = ECOTRANphysics.current_light;     % time-series of surface solar raditation at current time (W/m2); NOTE: current time = midnight when dt = 1 day; (vertical vector: num_t X 1)
% % % sunrise                         = ECOTRANphysics.sunrise;           % time of sunrise (time in h from midnight; 12 = noon, set as a default)
% % % sunset                          = ECOTRANphysics.sunset;            % time of sunset  (time in h from midnight; 12 = noon, set as a default)
% % % Kw                            	= ECOTRANphysics.Kw;              	% Light attenuation_seawater; (scalar)
% % % Kp                             	= ECOTRANphysics.Kp;              	% Light attenuation_phytoplankton (Newberger et al., 2003); (m2/mmol N); (scalar)
% % % MLD                           	= ECOTRANphysics.MLD;            	% mixed-layer depth; (m); (vertical vector: num_t X 1)
% % % EuphoticDepth                  	= repmat(MLD, [1 num_boxes]);     	% FFF (eventually move to f_physics code); depth of euphotic zone, used when converting the vertically-integrated EwE primary producer biomass to biomass/volume; depth; (m); (2D matrix: num_t X num_boxes)
% % % 
% % % WWT_to_C                      	= ECOTRANphysics.WWT_to_C;              % (scalar)
% % % atomic_mass_C                  	= ECOTRANphysics.atomic_mass_C;         % (scalar)
% % % C_to_N_phytoplankton          	= ECOTRANphysics.C_to_N_phytoplankton;  % (scalar)
% % % 
% % % grp_row                      	= 1:num_grps;
% % % 
% % % spatial_BiomassScalers         	= ones(1, num_boxes);	% NCC scalers for estimating initial (or mean) primary producer biomasses across model domain; NOTE: these values are assumed; NOTE: x2 in Box I used to compensate for 30m depth relative to 15 m depths in Boxes II & IV; FFF apply NPZD scalers here
% % % % -------------------------------------------------------------------------
% % % 
% % % 
% % % %% step 6e: prep compact fluxes for the ODE solver-------------------------
% % % %          remove information defining non-existing box links
% % % %       CompactFlux
% % % %           compact_flux            (2D matrix: num_t X num_fluxes)
% % % %                                       NOTE: fluxes include all linked boxes +1 for external links
% % % %           looky_flux              (2D matrix: num_fluxes X 3)
% % % %                                       clm 1: (destiny box) = list of boxes importing water volume (+ boundary)
% % % %                                       clm 2: (source box) = list of boxes exporting water volume (+ boundary)
% % % %                                       clm 3: flux address in non-compacted flux 2D matrix: destiny X source
% % % %                                       NOTE: fluxes include all linked boxes +1 for external links
% % % %                                       NOTE: values constant independent of t
% % % %           looky_boundary_import	(2D matrix: num_fluxes_BoundaryImport X 3)
% % % %                                       clm 1: (destiny box) = identity of boxes importing water volume
% % % %                                       clm 2: (source box) = identity of boxes exporting water volume (always the boundary flux number)
% % % %                                       clm 3: (import flux address) = addresses of import flux clm in compact_flux)
% % % %           looky_boundary_export   (2D matrix: num_fluxes_BoundaryExport X 3)
% % % %                                       clm 1: (destiny box) = identity of boxes importing water volume (always the boundary flux number)
% % % %                                       clm 2: (source box) = identity of boxes exporting water volume
% % % %                                       clm 3: (export flux address) = addresses of export fluxes clm in compact_flux
% % % %           unique_source           (vertical vector: list of source boxes (+ boundary))
% % % %           unique_destiny          (vertical vector: list of destiny boxes (+ boundary))
% % % %           num_fluxes              number of realized fluxes between boxes (or boundary) over full time-series
% % % %           fname_CompactFlux       name of this FuncName_CompactFlux function
% % % 
% % % % ADVECTION -----
% % % ODEinput.num_fluxes_advection                   = CompactFlux_ADVECTION.num_fluxes;
% % % ODEinput.ADVECTION_compact                      = CompactFlux_ADVECTION.compact_flux;                         % (m3/d); (2D matrix: num_t X num_fluxes_advection)
% % % ODEinput.looky_AdvectionFlux                    = CompactFlux_ADVECTION.looky_flux;                           % (2D matrix: num_fluxes_advection X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
% % % ODEinput.looky_AdvectionBoundary_import         = CompactFlux_ADVECTION.looky_boundary_import;                % (2D matrix: num_fluxes_BoundaryImport X [(DESTINY box) (SOURCE box) (flux clm)]);
% % % ODEinput.looky_AdvectionBoundary_export         = CompactFlux_ADVECTION.looky_boundary_export;                % (2D matrix: num_fluxes_BoundaryExport X [(DESTINY box) (SOURCE box) (flux clm)]);
% % % 
% % % repmat_GrpRow_ADVECTION                         = repmat(grp_row', [1, CompactFlux_ADVECTION.num_fluxes]);        % addressses; (2D matrix: num_grps X num_fluxes); NOTE transpose
% % % repmat_looky_ADVECTION_source                   = repmat(CompactFlux_ADVECTION.looky_flux(:, 2)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_advection); NOTE transpose
% % % repmat_looky_ADVECTION_destiny                  = repmat(CompactFlux_ADVECTION.looky_flux(:, 1)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_advection); NOTE transpose
% % % ODEinput.repmat_looky_ADVECTION_source          = [repmat_GrpRow_ADVECTION(:) repmat_looky_ADVECTION_source(:)];  % vectorize for accumarray function with the ODE; (2D matrix: (num_grps * num_fluxes_advection) X 2)
% % % ODEinput.repmat_looky_ADVECTION_destiny         = [repmat_GrpRow_ADVECTION(:) repmat_looky_ADVECTION_destiny(:)]; % vectorize for accumarray function with the ODE; (2D matrix: (num_grps * num_fluxes_advection) X 2)
% % % % ----------------------
% % % 
% % % % HORIZONTALMIXING -----
% % % ODEinput.num_fluxes_HorizontalMixing            = CompactFlux_HORIZONTALMIXING.num_fluxes;
% % % ODEinput.HORIZONTALMIXING_compact               = CompactFlux_HORIZONTALMIXING.compact_flux;                         % (m3/d); (2D matrix: num_t X num_fluxes_HorizontalMixing)
% % % ODEinput.looky_HorizontalMixingFlux             = CompactFlux_HORIZONTALMIXING.looky_flux;                           % (2D matrix: num_fluxes_HorizontalMixing X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
% % % ODEinput.looky_HorizontalMixingBoundary_import	= CompactFlux_HORIZONTALMIXING.looky_boundary_import;                % (2D matrix: num_fluxes_BoundaryImport X [(DESTINY box) (SOURCE box) (flux clm)]);
% % % ODEinput.looky_HorizontalMixingBoundary_export	= CompactFlux_HORIZONTALMIXING.looky_boundary_export;                % (2D matrix: num_fluxes_BoundaryExport X [(DESTINY box) (SOURCE box) (flux clm)]);
% % % 
% % % repmat_GrpRow_HORIZONTALMIXING                	= repmat(grp_row', [1, CompactFlux_HORIZONTALMIXING.num_fluxes]);        % addressses; (2D matrix: num_grps X num_fluxes_HorizontalMixing); NOTE transpose
% % % repmat_looky_HORIZONTALMIXING_source          	= repmat(CompactFlux_HORIZONTALMIXING.looky_flux(:, 2)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_HorizontalMixing); NOTE transpose
% % % repmat_looky_HORIZONTALMIXING_destiny         	= repmat(CompactFlux_HORIZONTALMIXING.looky_flux(:, 1)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_HorizontalMixing); NOTE transpose
% % % ODEinput.repmat_looky_HORIZONTALMIXING_source	= [repmat_GrpRow_HORIZONTALMIXING(:) repmat_looky_HORIZONTALMIXING_source(:)];  % vectorize for accumarray function with the ODE; (2D matrix: (num_grps * num_fluxes_HorizontalMixing) X 2)
% % % ODEinput.repmat_looky_HORIZONTALMIXING_destiny	= [repmat_GrpRow_HORIZONTALMIXING(:) repmat_looky_HORIZONTALMIXING_destiny(:)]; % vectorize for accumarray function with the ODE; (2D matrix: (num_grps * num_fluxes_HorizontalMixing) X 2)
% % % % ----------------------
% % % 
% % % % VERTICALMIXING -------
% % % ODEinput.num_fluxes_VerticalMixing              = CompactFlux_VERTICALMIXING.num_fluxes;
% % % ODEinput.VERTICALMIXING_compact                 = CompactFlux_VERTICALMIXING.compact_flux;                         % (m3/d); (2D matrix: num_t X num_fluxes_VerticalMixing)
% % % ODEinput.looky_VerticalMixingFlux               = CompactFlux_VERTICALMIXING.looky_flux;                           % (2D matrix: num_fluxes_VerticalMixing X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
% % % ODEinput.looky_VerticalMixingBoundary_import	= CompactFlux_VERTICALMIXING.looky_boundary_import;                % (2D matrix: num_fluxes_BoundaryImport X [(DESTINY box) (SOURCE box) (flux clm)]);
% % % ODEinput.looky_VerticalMixingBoundary_export	= CompactFlux_VERTICALMIXING.looky_boundary_export;                % (2D matrix: num_fluxes_BoundaryExport X [(DESTINY box) (SOURCE box) (flux clm)]);
% % % 
% % % repmat_GrpRow_VERTICALMIXING                 	= repmat(grp_row', [1, CompactFlux_VERTICALMIXING.num_fluxes]);        % addressses; (2D matrix: num_grps X num_fluxes_VerticalMixing); NOTE transpose
% % % repmat_looky_VERTICALMIXING_source          	= repmat(CompactFlux_VERTICALMIXING.looky_flux(:, 2)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_VerticalMixing); NOTE transpose
% % % repmat_looky_VERTICALMIXING_destiny         	= repmat(CompactFlux_VERTICALMIXING.looky_flux(:, 1)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_VerticalMixing); NOTE transpose
% % % ODEinput.repmat_looky_VERTICALMIXING_source  	= [repmat_GrpRow_VERTICALMIXING(:) repmat_looky_VERTICALMIXING_source(:)];  % vectorize for accumarray function with the ODE; (2D matrix: (num_grps * num_fluxes_VerticalMixing) X 2)
% % % ODEinput.repmat_looky_VERTICALMIXING_destiny	= [repmat_GrpRow_VERTICALMIXING(:) repmat_looky_VERTICALMIXING_destiny(:)]; % vectorize for accumarray function with the ODE; (2D matrix: (num_grps * num_fluxes_VerticalMixing) X 2)
% % % % ----------------------
% % % 
% % % % SINKING --------------
% % % ODEinput.num_fluxes_sinking                     = CompactFlux_SINKING.num_fluxes;
% % % num_fluxes_sinking                              = CompactFlux_SINKING.num_fluxes;
% % % SinkingArea_compact                             = CompactFlux_SINKING.compact_flux;                         % box floor area between sinking SOURCE box and DESTINY box; (m2); (2D matrix: num_t X num_fluxes_sinking)
% % % SinkingArea_compact                             = reshape(SinkingArea_compact, [num_t, 1, num_fluxes_sinking]); % (m2); (3D matrix: num_t X 1 X num_fluxes_sinking)
% % % SinkingArea_compact                             = repmat(SinkingArea_compact, [1, num_grps, 1]);	% (m2); (3D matrix: num_t X num_grps X num_fluxes_sinking)
% % % % ODEinput.SINKING_compact <<--sinking as volumetric flux (m3/d) for each functional group is calculated below
% % % ODEinput.looky_SinkingFlux                      = CompactFlux_SINKING.looky_flux;                           % (2D matrix: num_fluxes_sinking X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
% % % ODEinput.looky_SinkingBoundary_import           = CompactFlux_SINKING.looky_boundary_import;                % (2D matrix: num_fluxes_BoundaryImport X [(DESTINY box) (SOURCE box) (flux clm)]);
% % % ODEinput.looky_SinkingBoundary_export           = CompactFlux_SINKING.looky_boundary_export;                % (2D matrix: num_fluxes_BoundaryExport X [(DESTINY box) (SOURCE box) (flux clm)]);
% % % 
% % % repmat_GrpRow_SINKING                           = repmat(grp_row', [1, CompactFlux_SINKING.num_fluxes]);        % addressses; (2D matrix: num_grps X num_fluxes_sinking); NOTE transpose
% % % repmat_looky_SINKING_source                     = repmat(CompactFlux_SINKING.looky_flux(:, 2)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_sinking); NOTE transpose
% % % repmat_looky_SINKING_destiny                    = repmat(CompactFlux_SINKING.looky_flux(:, 1)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_sinking); NOTE transpose
% % % ODEinput.repmat_looky_SINKING_source            = [repmat_GrpRow_SINKING(:) repmat_looky_SINKING_source(:)];  % vectorize for accumarray function with the ODE; (2D matrix: (num_grps * num_fluxes_sinking) X 2)
% % % ODEinput.repmat_looky_SINKING_destiny           = [repmat_GrpRow_SINKING(:) repmat_looky_SINKING_destiny(:)]; % vectorize for accumarray function with the ODE; (2D matrix: (num_grps * num_fluxes_sinking) X 2)
% % % % ----------------------
% % % 
% % % % MIGRATION ------------
% % % ODEinput.num_fluxes_migration                 	= CompactFlux_MIGRATION.num_fluxes;                           % migration fluxes (potential number of fluxes for each group); includes external fluxes (external flux count defaults to 2 if all external fluxes equal 0)
% % % num_fluxes_migration                            = CompactFlux_MIGRATION.num_fluxes;                           % migration fluxes (potential number of fluxes for each group); includes external fluxes (external flux count defaults to 2 if all external fluxes equal 0)
% % % MigrationArea_compact                           = CompactFlux_MIGRATION.compact_flux;                         % migration as boundary area between boxes; (m2); (3D matrix: num_t X num_fluxes_migration)
% % % % ODEinput.MIGRATION_compact <<--migration as volumetric flux (m3/d) for each functional group is calculated below
% % % looky_MigrationFlux                             = CompactFlux_MIGRATION.looky_flux;                           % (2D matrix: num_fluxes_migration X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
% % % ODEinput.looky_MigrationFlux                    = CompactFlux_MIGRATION.looky_flux;                           % (2D matrix: num_fluxes_migration X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
% % % ODEinput.looky_MigrationBoundary_import         = CompactFlux_MIGRATION.looky_boundary_import;                % (2D matrix: num_fluxes_BoundaryImport X [(DESTINY box) (SOURCE box) (flux clm)]);
% % % ODEinput.looky_MigrationBoundary_export         = CompactFlux_MIGRATION.looky_boundary_export;                % (2D matrix: num_fluxes_BoundaryExport X [(DESTINY box) (SOURCE box) (flux clm)]);
% % % 
% % % repmat_GrpRow_MIGRATION                         = repmat(grp_row', [1, CompactFlux_MIGRATION.num_fluxes]);        % addressses; (2D matrix: num_grps X num_fluxes_migration); NOTE transpose
% % % repmat_looky_MIGRATION_source                   = repmat(CompactFlux_MIGRATION.looky_flux(:, 2)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_migration); NOTE transpose
% % % repmat_looky_MIGRATION_destiny                  = repmat(CompactFlux_MIGRATION.looky_flux(:, 1)', [num_grps, 1]);	% (2D matrix: num_grps X num_fluxes_migration); NOTE transpose
% % % ODEinput.repmat_looky_MIGRATION_source          = [repmat_GrpRow_MIGRATION(:) repmat_looky_MIGRATION_source(:)];  % vectorize for accumarray function with the ODE; (2D matrix: (num_grps * num_fluxes_migration) X 2)
% % % ODEinput.repmat_looky_MIGRATION_destiny         = [repmat_GrpRow_MIGRATION(:) repmat_looky_MIGRATION_destiny(:)]; % vectorize for accumarray function with the ODE; (2D matrix: (num_grps * num_fluxes_migration) X 2)
% % % % -------------------------------------------------------------------------
% % % 
% % % 
% % % %% step 6f: sinking rate of each group (m/d) -------------------------------
% % % %          NOTE: sinking is treated like a physical term (i.e., not incorporated into EnergyBudget)
% % % %          NOTE: apply this factor whether or not using Michaelis-Menten for primary producers
% % % SinkingSpeed                                = zeros(num_t, num_grps, num_fluxes_sinking);	% initialze sinking speed time-series; (m/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)
% % % 
% % % SinkingSpeed(:, rc_plgc_detritus, :)      	= repmat((10.5), [num_t, 1, num_fluxes_sinking]);   % SSS; sinking speed; (m/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)
% % % % SinkingSpeed(:, rc_bnth_detritus, :)       	= repmat((25),   [num_t, 1, num_fluxes_sinking]);   % SSS; sinking speed; (m/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)
% % % % SinkingSpeed(:, rc_fishery_offal, :)       	= repmat((40),   [num_t, 1, num_fluxes_sinking]);   % SSS; sinking speed; (m/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)
% % % disp('NOTICE: for vertically-integrated food-web technique IN 2D & 3D ROMS SHELF SETTINGS, sinking is applied ONLY FOR PELAGIC DETRITUS (change manually in code)')
% % % 
% % % % MichaelisMenten_w                               = [0.6 1.0];                                    % SSS; sinking speed; (m/d); [Sm Phytoplankton, Lg Phytoplankton]
% % % % SinkingSpeed(:, looky_ANYPrimaryProducer, :)	= repmat(MichaelisMenten_w, [num_t, 1, num_fluxes_sinking]);	% sinking speed; (m/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)
% % % % -------------------------------------------------------------------------
% % % 
% % % 
% % % %% step 6g: define DIEL VERTICAL MIGRATION (DVM) speeds & duration for each group (m/d) --------
% % % %           NOTE: This is not used in CAFA Northern California Current project, but null variables are still required to run the code
% % % %           NOTE: dusk phase are negative (in "standard" DVM)
% % % DVMinput.MigrationArea_compact	= MigrationArea_compact; % migration as boundary area between boxes; (m2); (3D matrix: num_t X num_fluxes_migration)
% % % DVMinput.MigrationSpeed       	= zeros(num_grps, (num_boxes+1), (num_boxes+1)); % initialize; (m/d); (3D matrix: num_grps, source (num_boxes+1) X destiny (num_boxes+1))
% % % 
% % % DVMspeed_meso2epi            	= zeros(num_grps, 2); % initialize; DVM speed MESO<<-->>EPI; (2D matrix: num_grps X 2 [dusk dawn]);
% % % DVMspeed_bathy2meso          	= zeros(num_grps, 2); % initialize; DVM speed BATHY<<-->>MESO; (2D matrix: num_grps X 2 [dusk dawn]);
% % % 
% % % DVM                             = f_DVMsinusoid_08122021(ECOTRANphysics, ODEinput, DVMinput, t_grid); % calculate Diel Vertical Migration terms; (structure)
% % % % -------------------------------------------------------------------------
% % % 
% % % 
% % % % step 6h: define production loss fractions -------------------------------
% % % %          NOTE: this used for initial conditions only (for dynamic models)
% % % %          NOTE: the ODE accounts for physical removal (or addition) at each time-step; while the static model accounts for physical loss as a reduction in Transfer Efficiency
% % % %                this is the way it was handled in John's original code outline and in the pubs
% % % PhysicalLossFraction        = repmat(ProductionLossScaler', [num_t, 1]);	% (2D matrix: num_t X num_grps); NOTE transpose
% % % PhysicalLossFraction        = PhysicalLossFraction * 0;                     % SSS set to 0 for time-dynamic runs
% % % % -------------------------------------------------------------------------
% % % 
% % % 
% % % %% step 6i: pack physics values for ODE into ODEinput ----------------------
% % % ODEinput.num_boxes                    	= num_boxes;
% % % % ODEinput.BoxLength                    	= BoxLength;                         % (m);  (2D matrix: time X num_boxes)
% % % % ODEinput.BoxHeight                    	= BoxHeight;                         % (m);  (2D matrix: time X num_boxes)
% % % % ODEinput.BoxWidth                      	= BoxWidth;                          % (m);  (2D matrix: time X num_boxes)
% % % ODEinput.BoxVolume                  	= BoxVolume;                         % (m3); (2D matrix: time X num_boxes)
% % % ODEinput.t_grid                      	= t_grid;
% % % ODEinput.dt                             = dt;
% % % ODEinput.num_t                          = num_t;
% % % ODEinput.MLD                        	= MLD;                               % (vertical vector: length = time)
% % % ODEinput.Io                           	= Io;
% % % ODEinput.Kw                           	= Kw;
% % % ODEinput.Kp                           	= Kp;
% % % ODEinput.RetentionScaler              	= repmat(RetentionScaler, [1, (num_boxes)]);	% (scaler 0-1); (2D matrix: num_grps X num_boxes DESTINATION);
% % % ODEinput.SinkingSpeed               	= SinkingSpeed;                      % sinking speed; (m/d); (3D matrix: num_t X num_grps X num_boxes)
% % % ODEinput.flux_domain_import_t        	= zeros(num_grps, (num_boxes+1));	 % initialized variable used in intraODE; (2D matrix: num_grps X num_boxes+1)
% % % ODEinput.flux_domain_export_t       	= zeros(num_grps, (num_boxes+1));	 % initialized variable used in intraODE; (2D matrix: num_grps X num_boxes+1)
% % % ODEinput.flux_domain_import_driver_t	= zeros(num_grps, num_boxes);        % initialized variable used in intraODE; (2D matrix: num_grps X num_boxes)
% % % ODEinput.biomass_plus1                  = zeros(num_grps, 1, (num_boxes+1)); % initialized as all zeros; add 1 clm for boundary fluxes; (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1)
% % % ODEinput.PhysicalLossFraction        	= PhysicalLossFraction;              % (2D matrix: num_t X num_grps)
% % % ODEinput.SINKING_compact             	= SinkingArea_compact .* SinkingSpeed;  % sinking fluxes; (m3/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)                    % box floor area between sinking source box and destination box; (m2); (2D matrix: num_t X num_fluxes_sinking)
% % % ODEinput.MIGRATION_compact              = DVM.MIGRATION_compact;             	% migration fluxes; (m3/d); (3D matrix: num_t X num_grps X num_fluxes_migration)
% % % % *************************************************************************


% *************************************************************************
% OPTIONAL STEP 6 version: load pre-aggregated ROMS output

% SSS -->define directories holding ROMS files
PHYSICSinput.ROMSfile_directory     = '/Volumes/Fortress/Build_Package/For_GitHub/ROMS_driverFiles/'; % directory to ROMS time-series

% load GFDL driver
load([PHYSICSinput.ROMSfile_directory 'PrePrep_GFDL8050_12102022.mat'], 'ODEinput', 'ECOTRANphysics', 'ECOTRANmigration', 'DVM')

% % load IPSL driver
% load([PHYSICSinput.ROMSfile_directory 'PrePrep_IPSL_1980_2050_12_13_2022.mat'], 'ODEinput', 'ECOTRANphysics', 'ECOTRANmigration', 'DVM')

% % load HAD driver
% load([PHYSICSinput.ROMSfile_directory 'PrePrep_HAD_1980_2050_12_13_2022.mat'], 'ODEinput', 'ECOTRANphysics', 'ECOTRANmigration', 'DVM')

num_t                           = ODEinput.num_t;
num_boxes                       = ECOTRANphysics.num_boxes;
num_domains                     = ECOTRANphysics.num_domains;	% number of geographic domains (does not consider depth layers)
num_z                           = ECOTRANphysics.num_z;         % number of depth layers

spatial_BiomassScalers         	= ones(1, num_boxes);	% NCC scalers for estimating initial (or mean) primary producer biomasses across model domain; NOTE: these values are assumed; NOTE: x2 in Box I used to compensate for 30m depth relative to 15 m depths in Boxes II & IV; FFF apply NPZD scalers here

ROMS_temperature_initial     	= ECOTRANphysics.ROMS_temperature_initial;          % (deg C); (horizontal vector: 1 X num_boxes)
ROMS_diatom                     = ECOTRANphysics.diatom_timeseries;                 % (mmole N/m3); (2D matrix: num_t X num_boxes)
ROMS_nanophytoplankton          = ECOTRANphysics.nanophytoplankton_timeseries;      % (mmole N/m3); (2D matrix: num_t X num_boxes)
ROMS_diatom_initial          	= ECOTRANphysics.ROMS_diatom_initial;               % (mmoles N/m3); (horizontal vector: 1 X num_boxes)
ROMS_nanophytoplankton_initial	= ECOTRANphysics.ROMS_nanophytoplankton_initial;	% (mmoles N/m3); (horizontal vector: 1 X num_boxes)

BoxVolume                    	= ECOTRANphysics.BoxVolume;     	% (m3); (2D matrix: num_t X num_boxes)

Io                            	= ECOTRANphysics.Io;             	% time-series of surface PAR light intensity; daily integrated solar raditation at ocean surface averaged across 24 hours; (W m^-2 h^-1); (vertical vector: num_t X 1)
current_light                   = ECOTRANphysics.current_light;     % time-series of surface solar raditation at current time (W/m2); NOTE: current time = midnight when dt = 1 day; (vertical vector: num_t X 1)
sunrise                         = ECOTRANphysics.sunrise;           % time of sunrise (time in h from midnight; 12 = noon, set as a default)
sunset                          = ECOTRANphysics.sunset;            % time of sunset  (time in h from midnight; 12 = noon, set as a default)
Kw                            	= ECOTRANphysics.Kw;              	% Light attenuation_seawater; (scalar)
Kp                             	= ECOTRANphysics.Kp;              	% Light attenuation_phytoplankton (Newberger et al., 2003); (m2/mmol N); (scalar)
MLD                           	= ECOTRANphysics.MLD;            	% mixed-layer depth; (m); (vertical vector: num_t X 1)
EuphoticDepth                  	= repmat(MLD, [1 num_boxes]);     	% FFF (eventually move to f_physics code); depth of euphotic zone, used when converting the vertically-integrated EwE primary producer biomass to biomass/volume; depth; (m); (2D matrix: num_t X num_boxes)

NO3timeseries_conc              = ECOTRANphysics.NO3timeseries_conc;    % NO3 + NO2 concentration interpolated from monthly means; (mmole N/m3); (2D matrix: num_t X num_boxes)
NO3initial_rate                 = ECOTRANphysics.NO3initial_rate;       % initial NO3 input rate; (mmoles N/m3/d); (horizontal vector: 1 X DESTINY (num_boxes))

WWT_to_C                      	= ECOTRANphysics.WWT_to_C;              % (scalar)
atomic_mass_C                  	= ECOTRANphysics.atomic_mass_C;         % (scalar)
C_to_N_phytoplankton          	= ECOTRANphysics.C_to_N_phytoplankton;  % (scalar)

grp_row                      	= 1:num_grps;

calendar_time                   = datevec(ECOTRANphysics.datestart + ECOTRANphysics.t_grid - 1); % the actual calendar time for t_grid

CompactFlux_ADVECTION.fname_CompactFlux = 'pre-aggregated ROMS output';            % file name of f_CompactFluxTimeSeries function
% *************************************************************************





%% *************************************************************************
% STEP 8: prepare external_driver and externalForcing time-series----------
%         NOTE: DEFAULT: The reflective boundary assumption is that biomass of non-external_driver 
%               groups are the same on either side of model domain outer boundaries.

% step 8a: define external_driver and externalForcing time-series ---------
external_driver             = zeros(num_t, 1, (num_boxes+1));	% initialize; (mmole N/m3);   (3D matrix: num_t X 1 X num_boxes+1); NOTE: added layer for external boundary driver (e.g. NO3) biomass for each box
looky_driver                = looky_NO3;                        % identify driver group(s); intitialize with NO3 even if case ExternalDriver_BoundaryConcentration not used
num_drivers                 = length(looky_driver);

externalForcing             = zeros(num_t, 1, num_boxes);       % initialize; (mmole N/m3/d); (3D matrix: num_t X 1 X num_boxes)
looky_externalForcing       = [];                               % identify externalForcing driver group(s)
num_externalForcing_grps	= length(looky_externalForcing);    % number of externally forced groups


% step 8b: Define external forcing rate time-series from 3D ROMS-BGC settings
looky_externalForcing               = [rc_lrg_phyto rc_sml_phyto]; % identify externalForcing driver group(s);	% row address(es) of externally forced input group(s) (e.g., NO3, phytoplankton, juvenile salmon)
num_externalForcing_grps         	= length(looky_externalForcing); % number of externally forced groups

ROMS_diatom                         = ROMS_diatom            .* (qb(rc_lrg_phyto) / 365); % convert to q rate; (mmole N/m3/d); (2D matrix: num_t X num_boxes)
ROMS_nanophytoplankton            	= ROMS_nanophytoplankton .* (qb(rc_sml_phyto) / 365); % convert to q rate; (mmole N/m3/d); (2D matrix: num_t X num_boxes)

ROMS_diatom                         = reshape(ROMS_diatom,            [num_t, 1, num_boxes]);	% (mmole N/m3/d); (3D matrix: num_t X 1 X num_boxes)
ROMS_nanophytoplankton            	= reshape(ROMS_nanophytoplankton, [num_t, 1, num_boxes]);	% (mmole N/m3/d); (3D matrix: num_t X 1 X num_boxes)

externalForcing(:, 1, :)            = ROMS_diatom; % forced external input; (mmole N/m3/d); (3D matrix: num_t X num_externalForcing_grps X num_boxes)
externalForcing(:, 2, :)            = ROMS_nanophytoplankton; % forced external input; (mmole N/m3/d); (3D matrix: num_t X num_externalForcing_grps X num_boxes)

% deactivate nutrient uptake by phytoplankton when driving model with BGC model output
%   NOTE: this change means that initial conditions MUST be defined by the BGC externalForcing driver
disp('USING EXTERNAL DRIVER (BGC model output) -- nutrient recycling & uptake by phytoplankton DEACTIVATED')
EnergyBudget_MC(:, 1:3, :)          = 0; % nutrients NOT used in food web when driving model with BGC model output
ConsumptionBudget_MC(2, 1:3, :)     = 0; % nitrification OFF
ConsumptionBudget_MC(4, 1:3, :)     = 0; % nutrient uptake OFF
ConsumptionBudget_MC(7, 1:3, :)     = 1; % turn on nutrient emigration to 100% to prevent build-up of nutrients in system

% QQQ 12/2/2022 temp patch for deep boxes that don't exist in
% shallow seas AND to catch negative values that MIGHT be from
% original ROMS BGC output
looky_NaN = find(isnan(externalForcing));
externalForcing(looky_NaN) = 0;
looky_negative = find(externalForcing < 0);
externalForcing(looky_negative) = 0;
% QQQ -----------------
        
% -------------------------------------------------------------------------


% step 8b: define external forcing by geographic migrator group(s) --------
% FFF coming soon
% -------------------------------------------------------------------------


% step 8c: pack external_driver, externalForcing, and biomass_boundary for ODE-solver ----------
ODEinput.looky_driver               = looky_driver;             % driver group address(es)
ODEinput.num_drivers                = num_drivers;
ODEinput.external_driver            = external_driver;          % (mmole N/m3); (3D matrix: num_t X num_drivers X num_boxes+1)
% ODEinput.biomass_boundary           = biomass_boundary;	  % NOTE: use for defined boundary conditions (OPTION 2); (mmole N/m3); (3D matrix: num_t X num_grps X num_boxes+1)

ODEinput.looky_externalForcing      = looky_externalForcing;	% row address(es) of externally forced input group(s) (e.g., NO3, juvenile salmon)
ODEinput.num_externalForcing_grps	= num_externalForcing_grps;
ODEinput.externalForcing            = externalForcing;          % forced external input; (mmole N/m3/d); (3D matrix: num_t X num_externalForcing_grps X num_boxes)
% *************************************************************************





%% *************************************************************************
% STEP 9: define functional predator-prey relations------------------------

% step 9a: define Michaelis-Menten functional predator-prey relations------
%         NOTE: for primary producers

% initialize Michaelis-Menten parameters
MichaelisMenten_Vmax     = zeros(num_ANYPrimaryProd, 1, num_MC); % Vmax  = maximum nutrient uptake rate;          (1/d);        (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
MichaelisMenten_KNO3     = zeros(num_ANYPrimaryProd, 1, num_MC); % K     = NO3 half-saturation constant;          (mmol N/m3);  (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
MichaelisMenten_KNH4     = zeros(num_ANYPrimaryProd, 1, num_MC); % K     = NH4 half-saturation constant;          (mmol N/m3);  (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
MichaelisMenten_alpha    = zeros(num_ANYPrimaryProd, 1, num_MC); % alpha = initial slope of light response curve; (m2/W/d);     (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
MichaelisMenten_psi      = zeros(num_ANYPrimaryProd, 1, num_MC); % psi   = NO3 uptake inhibition by NH4;          (m3/mmole N); (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
MichaelisMenten_w        = zeros(num_ANYPrimaryProd, 1, num_MC); % w     = sinking rate;                          (m/d);        (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
MichaelisMenten_eta      = zeros(num_ANYPrimaryProd, 1, num_MC); % eta   = non-grazing mortality;                 (1/d);        (3D matrix: num_ANYPrimaryProd X 1 X num_MC)
% -------------------------------------------------------------------------


% step 9b: select a "linear" or "ECOSIM" response for consumers------------
%         NOTE: FunctionalResponseParams is a function of the CONSUMER (not the producer) and needs to be aligned with the consumer ROWS in ECOTRAN (not producer clms)
%         NOTE: to change one half-sat constant for individual groups, change FunctionalResponseParams([looky_grp])
%                FunctionalResponseParams = 0 is "constant donor-driven"
%                FunctionalResponseParams = 1 is "non-linear" ECOSIM default
FunctionalResponse_CV       = ones(num_grps, 4) * 0;        % SSS --> set uncertainty level for functional response terms; (2D matrix: num_grps->CONSUMERS X 4)
[FunctionalResponseParams, fname_FunctionalResponse_MonteCarlo]	= f_FunctionalResponse_MonteCarlo_05132021(ECOTRAN, ODEinput, FunctionalResponse_CV); % producer "vulnerabilities", m_p; (3D matrix: CONSUMERS X prey group X num_MC) replicated across clms (= producers)

switch switch_FunctionalResponse

    case 'Linear' % constant (predation independent of predator biomass)
        FunctionalResponseParams	= FunctionalResponseParams * 0;     % SSS force to linear for testing & default; (3D matrix: CONSUMERS X prey group X num_MC) replicated across clms (= producers)
        disp('NOTE: DONOR-DRIVEN FUNCTIONAL RESPONSE')
	% end (case 'Linear') --------------
    
    case 'NonLinear_default' % non-linear default
        disp('NOTE: DEFAULT NON-LINEAR FUNCTIONAL RESPONSE')
	% end (case 'NonLinear_default') --------------
    
    case 'NonLinear_alt' % constant (predation independent of predator biomass)
        FunctionalResponseParams	= FunctionalResponseParams * 2;     % SSS force to linear for testing & default; (3D matrix: CONSUMERS X prey group X num_MC) replicated across clms (= producers)
        disp('NOTE: ALTERNATE NON-LINEAR FUNCTIONAL RESPONSE')
	% end (case 'NonLinear_alt') --------------
    
end % switch (switch_FunctionalResponse)-----------------------------------
% *************************************************************************





%% *************************************************************************
% STEP 10: calculate temperature response scalers--------------------------

% step 10a: prep temperature time-series and reference temperatures -------
% redefine temperature time-series for sub-surface groups; (3D matrix: num_t X num_grps X num_boxes)
temperature_timeseries          = ECOTRANphysics.temperature_timeseries;        % (2D matrix: num_t X num_boxes)
VerticalConnectivity_stack      = ECOTRANphysics.VerticalConnectivity_stack;    % (2D matrix: num_domains X num_agg_z); NOTE: num_domains refers to the geographic domains when looking at a 2D map
        
% account for shallow geographic domains where deeper depth layers are not represented (NaN)
%   fill null depth layers (NaN columns) with the deepest non-null depth layer in that geographic domain
for domain_loop = 1:num_domains
    current_domain_boxes	= VerticalConnectivity_stack(domain_loop, :);
    evaluate_timeseries     = temperature_timeseries(1, current_domain_boxes);
    looky_NaN               = min(find(isnan(evaluate_timeseries)));
    if ~isempty(looky_NaN)
        looky_bad_boxes     = VerticalConnectivity_stack(domain_loop, looky_NaN:num_z);
        num_bad_boxes       = length(looky_bad_boxes);
        looky_good_box      = VerticalConnectivity_stack(domain_loop, max([1, (looky_NaN-1)]));
        temperature_timeseries(:, looky_bad_boxes) = temperature_timeseries(:, repmat(looky_good_box, [num_bad_boxes, 1]));
    end % (~isempty(looky_NaN))
end % (domain_loop)
        
temperature_timeseries          = reshape(temperature_timeseries, [num_t, 1, num_boxes]); % (3D matrix: num_t X 1 X num_boxes)
temperature_timeseries          = repmat(temperature_timeseries, [1, num_grps, 1]);       % (3D matrix: num_t X num_grps X num_boxes)        

% re-assign temperature timeseries based on life-history depth distributions
looky_grpType0                  = find(depthdist == 0); % make NO change to temperature timeseries
looky_grpType1                  = find(depthdist == 1); % use only surface depth zone temperature timeseries
looky_grpType2                  = find(depthdist == 2); % use sub-surface depth zone temperature timeseries
looky_grpType3                  = find(depthdist == 3); % use bottom depth zone temperature timeseries

for domain_loop = 1:num_domains
    temperature_timeseries(:, looky_grpType2, VerticalConnectivity_stack(domain_loop, 1)) = mean(temperature_timeseries(:, looky_grpType2, VerticalConnectivity_stack(domain_loop, 2:3)), 3);
    temperature_timeseries(:, looky_grpType3, VerticalConnectivity_stack(domain_loop, 1)) = temperature_timeseries(:, looky_grpType3, VerticalConnectivity_stack(domain_loop, num_z));
end % (domain_loop)

% Just use to calculate reference temperature by group between 1980 and 2010 as reference time period
looky_refPeriod_start   = find(calendar_time(:, 1) == 1980 & calendar_time(:, 2) == 1 & calendar_time(:, 3) == 1);
looky_refPeriod_end     = find(calendar_time(:, 1) == 2010 & calendar_time(:, 2) == 1 & calendar_time(:, 3) == 1);

temp_ref             	= mean(temperature_timeseries(1:10950, :, 1:num_domains), 1, 'omitnan'); % mean across desired reference time period; (3D matrix: 1 X num_grps X num_domains); NOTE: we only need to consider domain boxes 1-15 now
temp_ref_std          	= std(temperature_timeseries(1:10950, :, 1:num_domains), 1, 'omitnan'); % standard deviation across desired reference time period; (3D matrix: 1 X num_grps X num_domains); NOTE: we only need to consider domain boxes 1-15 now

% Reference temperature specific to each ESM defined in switch above; calculated in above loop for reference time period (1980-2010) averaged for each group.
% temperature_reference DEFINED FOR FUNCTIONAL GROUPS and depend on which depth layer the group lives in; (horizontal vector: 1 X num_grps)
temperature_reference           = mean(temp_ref, 3, 'omitnan'); % mean across domains; (horizontal vector: 1 X num_grps)  
% -------------------------------------------------------------------------


% step 10b: Q10 scaler for metabolic costs --------------------------------
% define base Q10 parameters with zoopl and fish set with default value of 2.36 (from Clarke 2004)
Q10base     = [1 1 1 1 1 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 ...
                2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 ...
                2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 2.36 1 1 1 1 ...
                1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];

if length(Q10base) ~= num_grps
    warning('The number of Q10 parameters does not equal the number of groups')
end

% SSS provide optional scaling to Q10 parameter here for sensitivity analyses
Q10_parameter = Q10base;

switch switch_Q10
    
    case 'Q10_ON'
        disp('Q10 metabolic rate scaling is ON')
        % calculate Q10_scaler
        %   takes:
        %       Q10_parameter           (horizontal vector: 1 X num_grps)
        %       temperature_reference	(horizontal vector: 1 X num_grps)
        %       temperature_timeseries  (3D matrix: num_t X num_grps X num_boxes)
        Q10_scaler      = f_physiology_Q10_12082022(Q10_parameter, temperature_reference, temperature_timeseries); % (3D matrix: num_t X num_grps X num_boxes)
	% end (case 'Q10_ON') --------------
        
    case 'Q10_OFF'
        disp('Q10 metabolic rate scaling is OFF')
        Q10_scaler      = []; % QQQ could set to all ones (maybe eliminate a switch below in step 13e?)
	% end (case 'Q10_OFF') --------------
        
end % (switch switch_Q10) ------------------------------------------
% -------------------------------------------------------------------------


% step 10c: Thornton & Lessem scaling for consumption rates ---------------
%           consumption scaling factor (relative to q_max); (3D matrix: num_t X num_grps X num_boxes)
switch switch_ThorntonLessem
    
    case 'ThorntonLessem_ON'
        
        disp('Thornton-Lessem consumption rate scaling is ON')
        
        % Prepare Thornton-Lessem parameters for each species
        ThorntonLessem_parameters                   = f_TLparameterization_08182023(temp_ref, temp_ref_std, looky_grpType1, looky_grpType2, looky_grpType3); % (2D matrix: 8 X num_grps)

        % calculate ThorntonLessem consumption scaler
        [q_TemperatureScaler, fname_ThorntonLessem]	= f_ThorntonLessem_12012022(ThorntonLessem_parameters, temperature_timeseries); % value between 0 and 1; (3D matrix: num_t X num_grps X num_boxes)
        
	% end (case 'ThorntonLessem_ON') --------------

    case 'ThorntonLessem_OFF'
        disp('Thornton-Lessem consumption rate scaling is OFF')
        fname_ThorntonLessem	    = 'f_ThorntonLessem not used';
        ThorntonLessem_parameters   = ones(8, num_grps); %  initialize (2D matrix: 8 X num_grps)
        q_TemperatureScaler         = ones(num_t, num_grps, num_boxes); % QQQ 6/14/2022 Here is were the main error lay. You gave the ThortonLessem parameters values of 1 but you really needed to give the q_TemperatureScaler values of 1 (when ThortonLessem parameters = 1, the q_TemperatureScaler does not = 1)
	% end (case 'ThorntonLessem_OFF') --------------

end % (switch switch_ThorntonLessem) --------------------------------------
% *************************************************************************





%% *************************************************************************
% STEP 11: adjustments to ConsumptionBudget--------------------------------
%           ConsumptionBudget_MC:
%                               1) feces
%                               2) metabolism
%                               3) eggs (reproduction)
%                               4) predation
%                               5) senescence
%                               6) ba (biomass accumulation)
%                               7) em (emigration); NOTE: negative for immigration

% step 11a: Adjust terminal benthic detritus in ConsumptionBudget
%           Removal (sequestration) of terminal benthic detritus is accounted for via emigration (em)
%           Add senescence term to em & set senescence term to 0
%           NOTE: ConsumptionBudget_MC dimensions (3D matrix: 7 X num_grps X num_MC)
ConsumptionBudget_MC(7, looky_terminalBNTHdetritus, :)	= ConsumptionBudget_MC(7, looky_terminalBNTHdetritus, :) + ConsumptionBudget_MC(5, looky_terminalBNTHdetritus, :);
ConsumptionBudget_MC(5, looky_terminalBNTHdetritus, :)	= 0;
% -------------------------------------------------------------------------


% % step 11b: apply and/or test alternate ConsumptionBudget terms------------
% ConsumptionBudget_MC(4, InertTracer2_RC, :)	= 0;	% CB_predation QQQ
% ConsumptionBudget_MC(5, [InertTracer2_RC InertTracer3_RC], :)	= 0;    % CB_senescence QQQ
% ConsumptionBudget_MC(1, InertTracer2_RC, :)	= 0.1;    % CB_feces
% ConsumptionBudget_MC(2, InertTracer2_RC, :)	= 0.1;    % CB_metabolism
% ConsumptionBudget_MC(3, InertTracer2_RC, :)	= 0;    % CB_eggs
% ConsumptionBudget_MC(4, InertTracer2_RC, :)	= 0.1;	% CB_predation
% ConsumptionBudget_MC(5, InertTracer2_RC, :)	= 0.1;    % CB_senescence
% ConsumptionBudget_MC(6, InertTracer2_RC, :)	= 0;	% CB_ba
% ConsumptionBudget_MC(7, InertTracer2_RC, :)	= -0.4;    % CB_em
% % -------------------------------------------------------------------------
% *************************************************************************





%% *************************************************************************
% STEP 12: define individual sub-region box types--------------------------
%         NOTE: this can be replaced by the deliberate definition of
%               individual sub-regional food webs
%         NOTE: Individual submodel food webs are defined below in step 13b)

% step 12a: define specific box types -------------------------------------
%           NOTE: "BoxType" is used to activate/deactive trophic relationships dpending on the location of the spatial grid cell
%           NOTE: for CAFA ecophysiology study, trophic relationships are active in all boxes
%          ConsumptionBudget_BoxType
%                               1) feces
%                               2) metabolism
%                               3) eggs (reproduction)
%                               4) predation
%                               5) senescence
%                               6) ba (biomass accumulation)
%                               7) em (emigration); NOTE: negative for immigration
% use this for vertically-resolved models and 3D models
EnergyBudget_BoxType        = ones(num_grps, num_grps, num_boxes);	% (3D matrix: num_grps X num_grps X num_boxes)
ConsumptionBudget_BoxType	= ones(7, num_grps, num_boxes);         % (3D matrix: 7 X num_grps X num_boxes)
% -------------------------------------------------------------------------


% step 12b: pack BoxType definitions into ODEinput ------------------------
ODEinput.EnergyBudget_BoxType         	= EnergyBudget_BoxType;         % distinguish between surface & sub-surface EnergyBudget; (3D matrix: num_grps X num_grps X num_boxes)
ODEinput.ConsumptionBudget_BoxType      = ConsumptionBudget_BoxType;	% distinguish between surface & sub-surface ConsumptionBudget; (3D matrix: 7 X num_grps X num_boxes)
% *************************************************************************





%% *************************************************************************
% STEP 13: run ODE for each MonteCarlo model-------------------------------

% step 13a: loop through each MonteCarlo model or treatment ---------------
% % QQQ turn off MonteCarlo_loop to run only the "type" model or for debugging 
MonteCarlo_loop     = 1; % MonteCarlo_loop is off for debugging
saveFile            = [SaveFile_directory SaveFile_label '_' date '.mat'];

% for MonteCarlo_loop = 1:num_MC
%     display(['MonteCarlo run ' num2str(MonteCarlo_loop) ' of ' num2str(num_MC)])
%     saveFile                        = [SaveFile_directory SaveFile_label '_' num2str(MonteCarlo_loop) '.mat'];
    
    % step 13b: pick current MonteCarlo model -----------------------------
    current_biomass                	= biomass(:, 1);	% (t WWT/km2); NOTE: these are INITIAL biomass conditions; (vertical vector: num_grps X 1); NOTE: nutrients are zeros
    
    current_pb                      = pb(:, 1)';      	% (1/y); (horizontal vector: 1 X num_grps); NOTE transpose
    current_pb                      = current_pb / 365;	% (1/d); specific growth rate per day; (horizontal vector: 1 X num_grps)
    current_pb([looky_nutrients; looky_eggs; looky_ANYdetritus])	= 1;	% pb for nutrients & detritus are always 1 regardless of time-frame; (horizontal vector: 1 X num_grps)
    
    current_qb                      = qb(:, 1)';      	% (1/y); (horizontal vector: 1 X num_grps); NOTE transpose
    current_qb                      = current_qb / 365;	% (1/d); specific consumption rate per day; (horizontal vector: 1 X num_grps)
    current_qb([looky_nutrients; looky_eggs; looky_ANYdetritus; looky_fleets])	= 1;	% qb for nutrients, detritus, & fleets are always 1 regardless of time-frame; (horizontal vector: 1 X num_grps)
    
    current_EnergyBudget            = EnergyBudget_MC(:, :, MonteCarlo_loop);     	% (2D matrix: num_grps X num_grps)
    current_ConsumptionBudget       = ConsumptionBudget_MC(:, :, MonteCarlo_loop);	% (2D matrix: 7 X num_grps)
    
    current_fate_feces              = fate_feces(:, :, 1);      % (2D matrix: num_ANYdetritus X num_grps); FFF: fates are constant across Monte Carlo models as of now QQQ NO! I have new fates in ECOTRAN_MC, use them now???
	current_fate_metabolism       	= fate_metabolism(:, :, 1);	% (2D matrix: num_nutrients X num_grps); FFF: fates are constant across Monte Carlo models as of now    
    current_fate_eggs              	= fate_eggs(:, :, 1);    	% (2D matrix: num_eggs X num_grps); FFF: fates are constant across Monte Carlo models as of now
    % NOTE fate_predation is calculated below
    current_fate_senescence         = fate_senescence(:, :, 1);	% (2D matrix: num_ANYdetritus X num_grps); FFF: fates are constant across Monte Carlo models as of now
    
    current_TransferEfficiency     	= TransferEfficiency(MonteCarlo_loop, 1:num_grps);	% (horizontal vector: 1 X num_grps)
    
%     current_FunctionalResponseParams      = FunctionalResponseParams(:, :, MonteCarlo_loop); % SSS ACTIVATE IF  USING MONTE CARO FOR FUNCTIONAL RESPONSE; (2D matrix: CONSUMERS (num_grps) X prey group (num_grps)) replicated across clms (= producers)
	current_FunctionalResponseParams      = FunctionalResponseParams(:, :, 1); % QQQ USE IF NOT USING MONTE CARLO FOR FUNCTIONAL RESPONSE; (2D matrix: CONSUMERS (num_grps) X prey group (num_grps)) replicated across clms (= producers)

    current_MichaelisMenten_Vmax	= MichaelisMenten_Vmax(:,  1, MonteCarlo_loop);	% Vmax = maximum nutrient uptake rate; (1/d); (vertical vector: num_ANYPrimaryProd X 1)
    current_MichaelisMenten_KNO3	= MichaelisMenten_KNO3(:,  1, MonteCarlo_loop);	% KNO3 = NO3 half-saturation constant; (mmol N/m3); (vertical vector: num_ANYPrimaryProd X 1)
    current_MichaelisMenten_KNH4	= MichaelisMenten_KNH4(:,  1, MonteCarlo_loop);	% KNH4 = NH4 half-saturation constant; (mmol N/m3);  (vertical vector: num_ANYPrimaryProd X 1)
    current_MichaelisMenten_alpha	= MichaelisMenten_alpha(:, 1, MonteCarlo_loop);	% alpha = initial slope of light response curve; (m2/W/d); (vertical vector: num_ANYPrimaryProd X 1)
    current_MichaelisMenten_psi     = MichaelisMenten_psi(:,   1, MonteCarlo_loop);	% psi = NO3 uptake inhibition by NH4; (m3/mmole N); (vertical vector: num_ANYPrimaryProd X 1)
    current_MichaelisMenten_w       = MichaelisMenten_w(:,     1, MonteCarlo_loop);	% w	= sinking rate; (m/d); (vertical vector: num_ANYPrimaryProd X 1)
	current_MichaelisMenten_eta     = MichaelisMenten_eta(:,   1, MonteCarlo_loop);	% eta = non-grazing mortality; (1/d); (vertical vector: num_ANYPrimaryProd X 1)
    % ---------------------------------------------------------------------

    
    % step 13b: build-up spatial sub-models -------------------------------
    %           NOTE: at this step, variable layers define spatial boxes and no longer define Monte Carlo models
    switch switch_SubModel

        case 'identical_SubModel'  % OPTION 1: use the same food web across the shelf ----------------
    
            disp('NOTE: IDENTICAL regional biological models')
            
            current_biomass                 = repmat(current_biomass, [1, num_boxes]);                          % (t WWT/km2); (2D matrix: num_grps X num_boxes)
            current_biomass                 = current_biomass .* repmat(spatial_BiomassScalers, [num_grps 1]);	% (t WWT/km2); (2D matrix: num_grps X num_boxes)

            current_pb                      = repmat(current_pb,                    [1, 1, num_boxes]);	% specific growth rate; (1/d); (3D matrix: 1 X num_grps X num_boxes)
            current_qb                      = repmat(current_qb,                    [1, 1, num_boxes]);	% specific consumption rate; (1/d); (3D matrix: 1 X num_grps X num_boxes)

            current_EnergyBudget            = repmat(current_EnergyBudget,          [1, 1, num_boxes]); % use for a single spatial definition of the EnergyBudget
            current_EnergyBudget            = current_EnergyBudget .* EnergyBudget_BoxType;             % adjust sub-surface EnergyBudget matrices; allow only certain recycling processes in sub-surface boxes; (fraction); (3D matrix: num_grps X num_grps X num_boxes)

            current_ConsumptionBudget       = repmat(current_ConsumptionBudget,     [1, 1, num_boxes]); % use for a single spatial definition of the ConsumptionBudget
            current_ConsumptionBudget       = current_ConsumptionBudget .* ConsumptionBudget_BoxType;	% adjust sub-surface ConsumptionBudget matrices; allow only certain recycling processes in sub-surface boxes; (fraction); (3D matrix: 7 X num_grps X num_boxes)

            current_fate_feces              = repmat(current_fate_feces,            [1, 1, num_boxes]); % (3D matrix: num_ANYdetritus X num_grps X num_boxes)
            current_fate_metabolism         = repmat(current_fate_metabolism,       [1, 1, num_boxes]); % (3D matrix: num_nutrients X num_grps X num_boxes)
            current_fate_eggs               = repmat(current_fate_eggs,             [1, 1, num_boxes]); % (3D matrix: num_eggs X num_grps X num_boxes)
            % NOTE fate_predation is calculated below
            current_fate_senescence         = repmat(current_fate_senescence,       [1, 1, num_boxes]); % (3D matrix: num_ANYdetritus X num_grps X num_boxes)

            current_TransferEfficiency      = repmat(current_TransferEfficiency,    [1, 1, num_boxes]); % (3D matrix: 1 X num_grps X num_boxes)        
            current_FunctionalResponseParams      = repmat(current_FunctionalResponseParams,    [1, 1, num_boxes]);	% (3D matrix: CONSUMERS (num_grps) X prey group (num_grps) X num_boxes) replicated across clms (= producers)

            current_MichaelisMenten_Vmax	= repmat(current_MichaelisMenten_Vmax,  [1, 1, num_boxes]);	% Vmax = maximum nutrient uptake rate; (1/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
            current_MichaelisMenten_KNO3	= repmat(current_MichaelisMenten_KNO3,  [1, 1, num_boxes]);	% KNO3 = NO3 half-saturation constant; (mmol N/m3); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
            current_MichaelisMenten_KNH4	= repmat(current_MichaelisMenten_KNH4,  [1, 1, num_boxes]);	% KNH4 = NH4 half-saturation constant; (mmol N/m3);  (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
            current_MichaelisMenten_alpha	= repmat(current_MichaelisMenten_alpha, [1, 1, num_boxes]);	% alpha = initial slope of light response curve; (m2/W/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
            current_MichaelisMenten_psi     = repmat(current_MichaelisMenten_psi,   [1, 1, num_boxes]);	% psi = NO3 uptake inhibition by NH4; (m3/mmole N); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
            current_MichaelisMenten_w       = repmat(current_MichaelisMenten_w,     [1, 1, num_boxes]); % w	= sinking rate; (m/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
            current_MichaelisMenten_eta     = repmat(current_MichaelisMenten_eta,   [1, 1, num_boxes]);	% eta = non-grazing mortality; (1/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
        % end (case 'identical_SubModel') --------------
            
        case 'independent_SubModel'  % OPTION 2: use independently defined webs for each shelf zone ----
                                     %           use for multiple definitions of the EnergyBudget
                                     % FFF future code could easily allow for Monte Carlo models for each sub-region

            disp('NOTE: LOADING independently defined regional biological models')
                                      
            % SSS define and load invidividual regional models from files
             
        % end (case 'independent_SubModel') --------------
            
    end % (switch switch_SubModel) ----------------------------------------
	% ---------------------------------------------------------------------
        

	% step 13c: calculate fate_predation ----------------------------------
    %           NOTE: multiple MC versions must replace the single "type" fate_predation from the main ECOTRAN code; 
    sum_predation                                           = sum(current_EnergyBudget(looky_livingANDfleets, :, :)); % (3D matrix: 1 X num_grps X num_boxes)
    current_fate_predation                                  = current_EnergyBudget(looky_livingANDfleets, :, :) ./ repmat(sum_predation, [num_livingANDfleets, 1, 1]); % (3D matrix: num_livingANDfleets X num_grps X num_boxes)
    current_fate_predation(isnan(current_fate_predation))	= 0; % correct div/0 errors
    % ---------------------------------------------------------------------
    
    
    % step 13d: build time-series of varying physiologies -----------------
    %      (rows = time, clms = ECOTRAN groups, layers = spatial boxes)
    %      NOTE: FFF at this point in the code, we can read in time-series changes for each of these terms
    %            this would allow for seasonal reproduction differences (ConsumptionBudget_eggs), & 
    %            seasonal migration changes (ConsumptionBudget_em)
    %      NOTE: the need to break out separate rows of ConsumptionBudget
    %            is to prevent having to deal with 4D matrices:
    %                                                   row 1: feces
    %                                                   row 2: metabolism
    %                                                   row 3: eggs
    %                                                   row 4: predation
    %                                                   row 5: senescence
    %                                                   row 6: ba
    %                                                   row 7: em

	current_pb                      = repmat(current_pb, [num_t, 1, 1]);	% specific growth rate per day; (3D matrix: num_t X num_grps X num_boxes)
	current_qb                      = repmat(current_qb, [num_t, 1, 1]);	% specific consumption rate per day; (3D matrix: num_t X num_grps X num_boxes)
    
    ConsumptionBudget_feces      	= repmat(current_ConsumptionBudget(1, :, :), [num_t, 1, 1]); % (3D matrix: num_t X num_grps X num_boxes)
    ConsumptionBudget_metabolism	= repmat(current_ConsumptionBudget(2, :, :), [num_t, 1, 1]); % (3D matrix: num_t X num_grps X num_boxes)
    ConsumptionBudget_eggs        	= repmat(current_ConsumptionBudget(3, :, :), [num_t, 1, 1]); % (3D matrix: num_t X num_grps X num_boxes)
    ConsumptionBudget_predation     = repmat(current_ConsumptionBudget(4, :, :), [num_t, 1, 1]); % (3D matrix: num_t X num_grps X num_boxes)
    ConsumptionBudget_senescence	= repmat(current_ConsumptionBudget(5, :, :), [num_t, 1, 1]); % (3D matrix: num_t X num_grps X num_boxes)
    ConsumptionBudget_ba            = repmat(current_ConsumptionBudget(6, :, :), [num_t, 1, 1]); % (3D matrix: num_t X num_grps X num_boxes)
    ConsumptionBudget_em            = repmat(current_ConsumptionBudget(7, :, :), [num_t, 1, 1]); % (3D matrix: num_t X num_grps X num_boxes)
    
    
    % rescale ConsumptionBudget_metabolism according to Q10 temperature response
    switch switch_Q10
        case 'Q10_ON'
            
            disp('NOTE: ConsumptionBudget_metabolism Q10 response to temperature change is ON')
            disp('NOTE: ConsumptionBudget components no longer need to sum to 1')

            ConsumptionBudget_metabolism_new	= ConsumptionBudget_metabolism .* Q10_scaler; % (3D matrix: num_t X num_grps X num_boxes)
            ConsumptionBudget_metabolism        = ConsumptionBudget_metabolism_new; % (3D matrix: num_t X num_grps X num_boxes)

        % end (case 'Q10_ON') --------------
            
        case 'Q10_OFF'
            disp('NOTE: ConsumptionBudget_metabolism Q10 response to temperature change is OFF')
        % end (case 'Q10_OFF') --------------
        
    end % (switch switch_Q10) --------------------------------------
    % ---------------------------------------------------------------------
    
    
% 	% step 13e: make scenario changes to ConsumptionBudget ----------------
%     %           NOTE: just a test scenario
%     ConsumptionBudget_ba(950:1000, rc_anchovy, 1)           = -0.5; % QQQ NCC: high mortality for 1 month
%     ConsumptionBudget_senescence(950:1000, rc_anchovy, 1)	= 0.1307; % QQQ NCC: high mortality for 1 month
%     ConsumptionBudget_predation(950:1000, rc_anchovy, 1)    = 0.1136; % QQQ NCC: high mortality for 1 month
%     
%     ConsumptionBudget_ba(950:1000, rc_anchovy, 1)           = 0.5; % QQQ NCC: high mortality for 1 month
%     ConsumptionBudget_ba(950:1000, str2num(SPECIES), 1)   	= 0.5; % QQQ NCC: high mortality for 1 month
%     % --------------------------------------------------------------------- 


    %% step 13f: pack variables needed for ODE -----------------------------
    ODEinput.pb                             = current_pb;                       % (1/d); (3D matrix: num_t X num_grps X num_boxes)
    ODEinput.qb                             = current_qb;                       % (1/d); (3D matrix: num_t X num_grps X num_boxes)
    
    ODEinput.EnergyBudget                   = current_EnergyBudget;             % (proportions); (3D matrix: num_grps X num_grps X num_boxes)
    
    ODEinput.ConsumptionBudget_feces        = ConsumptionBudget_feces;          % (proportions); (3D matrix: num_t X num_grps X num_boxes)
    ODEinput.ConsumptionBudget_metabolism	= ConsumptionBudget_metabolism;     % (proportions); (3D matrix: num_t X num_grps X num_boxes)
    ODEinput.ConsumptionBudget_eggs         = ConsumptionBudget_eggs;           % (proportions); (3D matrix: num_t X num_grps X num_boxes)
    ODEinput.ConsumptionBudget_predation	= ConsumptionBudget_predation;      % (proportions); (3D matrix: num_t X num_grps X num_boxes)
    ODEinput.ConsumptionBudget_senescence	= ConsumptionBudget_senescence;     % (proportions); (3D matrix: num_t X num_grps X num_boxes)
    ODEinput.ConsumptionBudget_ba           = ConsumptionBudget_ba;             % (proportions); (3D matrix: num_t X num_grps X num_boxes)
    ODEinput.ConsumptionBudget_em           = ConsumptionBudget_em;             % (proportions); (3D matrix: num_t X num_grps X num_boxes)

    ODEinput.fate_feces                     = current_fate_feces;               % (proportions); (3D matrix: num_ANYdetritus X num_grps X num_boxes)
    ODEinput.fate_metabolism                = current_fate_metabolism;          % (proportions); (3D matrix: num_nutrients X num_grps X num_boxes)
    ODEinput.fate_eggs                  	= current_fate_eggs;                % (proportions); (3D matrix: num_eggs X num_grps X num_boxes)
    ODEinput.fate_predation                 = current_fate_predation;           % (proportions); (3D matrix: num_livingANDfleets X num_grps X num_boxes)
    ODEinput.fate_senescence                = current_fate_senescence;          % (proportions); (3D matrix: num_ANYdetritus X num_grps X num_boxes)
    
    ODEinput.TransferEfficiency           	= current_TransferEfficiency;       % (3D matrix: 1 X num_grps X num_boxes)
    ODEinput.FunctionalResponseParams     	= current_FunctionalResponseParams;       % (vulnerability of producer, m_p); (3D matrix: CONSUMERS (num_grps) X prey group (num_grps) X num_boxes) replicated across clms (= producers)

    ODEinput.MichaelisMenten_Vmax           = current_MichaelisMenten_Vmax;     % Vmax = maximum nutrient uptake rate; (1/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
    ODEinput.MichaelisMenten_KNO3           = current_MichaelisMenten_KNO3;     % KNO3 = NO3 half-saturation constant; (mmol N/m3); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
    ODEinput.MichaelisMenten_KNH4           = current_MichaelisMenten_KNH4;     % KNH4 = NH4 half-saturation constant; (mmol N/m3);  (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
    ODEinput.MichaelisMenten_alpha          = current_MichaelisMenten_alpha;	% alpha = initial slope of light response curve; (m2/W/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
    ODEinput.MichaelisMenten_psi            = current_MichaelisMenten_psi;      % psi = NO3 uptake inhibition by NH4; (m3/mmole N); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
    ODEinput.MichaelisMenten_w              = current_MichaelisMenten_w;        % w	= sinking rate; (m/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
    ODEinput.MichaelisMenten_eta            = current_MichaelisMenten_eta;      % eta = non-grazing mortality; (1/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes)
    
    ODEinput.q_TemperatureScaler            = q_TemperatureScaler;              % Thornton-Lessem temperature adjustment to consumption rate; value between 0 and 1; (3D matrix: num_t X num_grps X num_boxes)
    % *********************************************************************


    
    
    
    %% *********************************************************************
	% STEP 14: calculate INITIAL production rate conditions----------------
    %          Calculate rates for all groups when model is driven by input
    %            of mean annual ROMS-BioGeoChemichal Primary Production rates
    %          INITIAL production = ingestion inflow

    t_initial                               = 1;
    
    production_initial_driver                           = zeros(1, num_grps, num_boxes);        % initialize DriverProductionVector; (3D matrix: 1 X num_grps X num_boxes)

    production_initial_driver(1, rc_lrg_phyto, :)    	= mean(externalForcing(:, 1, :), 1);	% plug in mean ROMS_diatom input rate; average over entire time-series; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
    production_initial_driver(1, rc_sml_phyto, :)     	= mean(externalForcing(:, 2, :), 1);	% plug in mean ROMS_nanophytoplankton input rate; average over entire time-series; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)

    % finalize initial condition calculations
    [production_initial, fname_InitialProductionRates]	= f_InitialProductionRates_02012022(ODEinput, production_initial_driver, t_initial);  % QQQ NH4 uptake turned OFF; initial or mean production rates (actually consumption inflow); (mmole N/m3/d); (2D matrix: num_grps X num_boxes); NOTE: code does NOT make transfer of bnthNH4 from surface to sub-surface boxes

    % paste in initial NO3 & NH4 input rates
    production_initial(looky_NO3, :)                    = 0;   % initial NO3 are 0 when driving with ROMS-BGC output; (mmole N/m3/d); (2D matrix: num_grps X num_boxes); NOTE the inconsistency in row definitions, nutrient values are concentrations but all other values are production rates!!! QQQ change this to NO3 input rates?
    production_initial(looky_plgcNH4, :)                = 0;   % initial NH4 are 0 when driving with ROMS-BGC output; (mmole N/m3/d); (2D matrix: num_grps X num_boxes); NOTE the inconsistency in row definitions, nutrient values are concentrations but all other values are production rates!!! QQQ change this to NO3 input rates?

        % QQQ 12/2/2022 temp patch to remove NaNs (for deep-layer boxes that don't really exist in shallow seas)
        looky_NaN = find(isnan(production_initial));
        production_initial(looky_NaN) = 0;
        % QQQ -----

        % CCC 1/18/2023 Change initial conditions for those species
        % going extinct in base model
         production_initial(102, :)                 = production_initial(102, :) * 0.25;
         production_initial([49, 50, 53, 104], :)   = production_initial([49, 50, 53 ,104],:) * 0.5;
        % QQQ -----


    production_initial                  	= reshape(production_initial, [num_grps, 1, num_boxes]);	% reshape boxes as layers; (3D matrix: num_grps X 1 X num_boxes); (mmole N/m3/d)
    productionC_initial_repmat            	= repmat(production_initial,   [1, num_grps, 1]);           % replicate groups across columns; (mmole N/m3/d); (3D matrix: num_grps X num_grps X num_boxes); NOTE: this is needed for the ODE
    % ---------------------------------------------------------------------
    
    
    % step 14b: pack initial conditions for ODE ---------------------------
    ODEinput.production_initial            	= production_initial;               % production rates to use as initial conditions; (mmole N/m3/d); (3D matrix: num_grps X 1 X num_boxes); NOTE: used in ODE only for special cases of step-thru debugging or for quadratic functional responses
    ODEinput.productionC_initial_repmat   	= productionC_initial_repmat;       % initial conditions reshaped; (mmole N/m3/d); (3D matrix: num_grps X num_grps X num_boxes); NOTE: replicated vertical vectors across columns
    % *********************************************************************

    
    
    
    
    %% *********************************************************************
    % STEP 15: solve the dynamic model-------------------------------------
    %          NOTE: >> mex mex_ECOTRANode_11272020.cpp -I/usr/local/include/ % compile mex function in matlab
	disp('NOTE: using C++ ODE solver')

    % step 15a: prepare ECOTRAN variables for using C++ ODE solver mex function
	%           Pack parameters & drivers along proper dimensions.
	disp('prepare variables for C++ ODE solver using ROMS BGC drivers...')
	[AddressStruct, TrophicStruct, FuncRespStruct, PhysicsStruct] = f_PrepMexODE_09192022(ODEinput); % pack variables to send to ODE solver; WITH Thornton-Lessem
	% ---------------------------------------------------------------------

	% step 15b: run the model in C++ --------------------------------------
	fname_ECOTRANode         	= 'mex_ECOTRANode_09182022';
	disp(['Running C++ solver using ROMS BGC drivers: ' fname_ECOTRANode])
	tic            
	[output_Cpp]                = mex_ECOTRANode_09182022(AddressStruct, TrophicStruct, FuncRespStruct, PhysicsStruct); % WITH Thornton-Lessem
	time_ODE                    = toc
	% ---------------------------------------------------------------------
       
    
    % step 15c: unstack result (store_ProductionRates) to retrieve spatial boxes
	store_T                     = output_Cpp.mat_obs_times; % (d); (vertical vector: num_t X 1)
	store_ProductionRates       = output_Cpp.mat_obs_states; % (mmole N/m3/d); (2D matrix: num_t X (num_grps*num_boxes))
	re_Y                        = reshape(store_ProductionRates, [num_t, num_grps, num_boxes]); % (mmole N/m3/d); (3D matrix: time X groups X num_boxes)
    % *********************************************************************
    
    
    
    
    
    %% *********************************************************************
    % STEP 16: save run results--------------------------------------------
    
    % step 16a: RUNlog of called functions --------------------------------
    RUNlog.fname_ECOTRANdynamic             = fname_ECOTRANdynamic;
    RUNlog.BiologicalModel                  = BiologicalModel_name;
    RUNlog.PhysicalModel                    = ECOTRANphysics.fname_PhysicalModel;
    RUNlog.MigrationModel                   = ECOTRANmigration.fname_ECOTRANmigration;
    RUNlog.DVMmodel                         = DVM.fname_DVMmodel;
    RUNlog.fname_ReadEwE                    = dat.fname_ReadEwE;
    RUNlog.fname_AggregateBiologicalModel	= EwEResult.fname_AggregateBiologicalModel;
    RUNlog.fname_ECOTRANheart               = ECOTRAN.fname_ECOTRANheart;
    RUNlog.fname_ECOfunction                = ECOTRAN.fname_ECOfunction;                % name of this version of f_ECOfunction
%     RUNlog.fname_RedistributeCannibalism	= ECOTRAN.fname_RedistributeCannibalism;	% FFF this function does not yet pass its name along; name of this version of f_RedistributeCannibalism
    RUNlog.fname_calcEE                     = ECOTRAN.fname_calcEE;                     % file name of this f_calcEEsub-function
%     RUNlog.fname_CalcPredationBudget        = ECOTRAN.fname_CalcPredationBudget;        % FFF this function does not yet pass its name along; file name of the sub-function f_CalcPredationBudget
    RUNlog.fname_E2Epedigree                = ECOTRAN_PEDIGREE.fname_E2Epedigree;       % name of this f_E2Epedigree function
    RUNlog.fname_E2E_MonteCarlo             = ECOTRAN_MC.fname_E2E_MonteCarlo;        	% name of this f_E2E_MonteCarlo function
    RUNlog.fname_LightIntensity             = ECOTRANphysics.fname_LightIntensity;      % name of physical model sub-function
    RUNlog.fname_EvaluateFluxBalance        = ECOTRANphysics.fname_EvaluateFluxBalance;	% name of this f_EvaluateFluxBalance sub-function
    RUNlog.fname_UnCompactFluxTimeSeries	= ECOTRANphysics.fname_UnCompactFluxTimeSeries; % name of this f_UnCompactFluxTimeSeries sub-function
    RUNlog.fname_CalcNetFlux                = ECOTRANphysics.fname_CalcNetFlux;          % name of this f_CalcNetFlux function
    RUNlog.fname_CompactFlux                = CompactFlux_ADVECTION.fname_CompactFlux;            % file name of f_CompactFluxTimeSeries function
    RUNlog.fname_FunctionalResponse_MonteCarlo = fname_FunctionalResponse_MonteCarlo;	% name of this f_FunctionalResponse_MonteCarlo function
    RUNlog.fname_InitialProductionRates     = fname_InitialProductionRates;             % name of f_InitialProductionRates function
    RUNlog.fname_WebProductivity            = 'f_WebProductivity_03272019';             % SSS name of f_WebProductivity function
    RUNlog.fname_physiology_Q10             = '_physiology_Q10_12082022';                            % name of fname_physiology_Q10 function
    RUNlog.fname_ThorntonLessem             = fname_ThorntonLessem;                     % name of f_ThorntonLessem function
    RUNlog.fname_ECOTRANode                 = fname_ECOTRANode;
	RUNlog.fname_VarianceDivision           = 'f_VarianceDivision_12132018';            % SSS name of f_VarianceDivision function
	RUNlog.fname_VarianceMultiplication     = 'f_VarianceMultiplication_12132018';      % SSS name of f_VarianceMultiplication function
    RUNlog.time_ODE                         = time_ODE;                                 % time to run ODE (seconds)
    % ---------------------------------------------------------------------
    
    
    % step 16b: save run results ------------------------------------------
%     save(saveFile, 'RUNlog', 're_Y', 'store_T', 'PHYSICSinput', 'ECOTRANphysics', 'ODEinput'); % save model run
%     save(saveFile, 'RUNlog', 're_Y', 'store_T', 'PHYSICSinput', 'ECOTRANphysics', 'ODEinput','-v7.3'); % save model run
    save(saveFile, 'RUNlog', 're_Y', 'store_T', 'ECOTRANphysics', 'ODEinput', '-v7.3'); % save model run

    % *********************************************************************

% end % (MonteCarlo_loop) -------------------------------------------------


%  end m-file***************************************************************
