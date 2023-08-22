function ECOTRAN_MC = f_E2E_MonteCarlo_08042020(MonteCarloConditions, ECOTRAN, ECOTRAN_PEDIGREE)
% Generate a stack of random, balanced ECOTRAN models.
% by Jim Ruzicka
%   Generates MonteCarlo Production Matrices and production estimates by drawing from a normal or a uniform distribution 
%   about each element of the EnergyBudget.
%
% calls:
%       none
%
% takes:
%       MonteCarloConditions
%           num_MC                      number of random Monte Carlo models to generate
%           DistributionType            distribution to draw random values from ('normal' or 'uniform'); NOTE: code not fully proofed for uniform
%       ECOTRAN
%           EnergyBudget                main model      (2D matrix: num_grps X num_grps); (NOTE: benthic detritus flowing to benthic detritus should be 0)
%           ConsumptionBudget           consumption budget (as fraction of total consumption); (2D matrix: 7 X num_grps)
%                                               row 1: feces
%                                               row 2: metabolism
%                                               row 3: eggs
%                                               row 4: predation
%                                               row 5: senescence
%                                               row 6: ba
%                                               row 7: em
%           LANDINGS                    landings (t/km2/y)                          (2D matrix: num_grps X num_fleets)
%           DISCARDS                    discards (t/km2/y)                          (2D matrix: num_grps X num_fleets)
%           CATCH                       catch (t/km2/y)                             (2D matrix: num_grps X num_fleets)
%           landings_TotalByFleet       total landings by fleet (t/km2/y)           (2D matrix: num_grps X num_fleets)
%           discards_TotalByFleet       total discards by fleet (t/km2/y)           (2D matrix: num_grps X num_fleets)
%           catch_TotalByFleet          total catch by fleet t/km2/y)               (2D matrix: num_grps X num_fleets)
%           DiscardFraction             fraction of catch discarded for each group	(3D matrix: num_grps X num_fleets)
%           GroupType
%           GroupTypeDef_               (several variables)
%           fate_metabolism                                                         (3D matrix: num_nutrients X num_grps X num_MC)
%           fate_eggs                                                               (3D matrix: num_eggs X num_grps X num_MC)
%           fate_feces                                                              (3D matrix: num_ANYdetritus X num_grps X num_MC)
%           fate_senescence                                                         (3D matrix: num_ANYdetritus X num_grps X num_MC)
%           fate_predation                                                          (2D matrix: num_livingANDfleets X num_grps)
%       ECOTRAN_PEDIGREE                (from f_E2Epedigree_12022019)
%           EnergyBudget_CV             CV values for each element of the ECOTRAN EnergyBudget (2D matrix: num_grps X num_grps)
%           ConsumptionBudget_CV        CV values for ConsumptionBudget terms, by row:
%                                               row 1: feces
%                                               row 2: metabolism
%                                               row 3: eggs
%                                               row 4: predation
%                                               row 5: senescence
%                                               row 6: ba
%                                               row 7: em
%           LANDINGS_CV                 landings; (CV); (2D matrix: num_grps X num_fleets)
%           DISCARDS_CV                 discards; (CV); (2D matrix: num_grps X num_fleets)
%           CATCH_CV                    catch; (CV); (2D matrix: num_grps X num_fleets)
%           landings_TotalByFleet_CV	total landings by fleet; (CV); (2D matrix: num_grps X num_fleets)
%           discards_TotalByFleet_CV	total discards by fleet; (CV); (2D matrix: num_grps X num_fleets)
%           catch_TotalByFleet_CV       total catch by fleet; (CV); (2D matrix: num_grps X num_fleets)
%           DiscardFraction_CV          discard fraction of each group by each fleet; (CV); (2D matrix: num_grps X num_fleets)
%
% returns:
%       ECOTRAN_MC
%           EnergyBudget_MC             all Monte Carlo production arrays           (3D matrix: num_grps X num_grps X num_MC)
%           ConsumptionBudget_MC        all Monte Carlo Bioenergetic Budget arrays	(3D matrix: 7 X num_grps X num_MC)
%                                               row 1: feces
%                                               row 2: metabolism
%                                               row 3: eggs
%                                               row 4: predation
%                                               row 5: senescence
%                                               row 6: ba
%                                               row 7: em
%           DiscardFraction_MC          fraction of catch discarded for each group	(3D matrix: num_grps X num_fleets X num_MC)
%           fate_metabolism                                                       	(3D matrix: num_nutrients X num_grps X num_MC)
%           fate_eggs                                                             	(3D matrix: num_eggs X num_grps X num_MC)
%           fate_feces                                                            	(3D matrix: num_ANYdetritus X num_grps X num_MC)
%           fate_senescence                                                        	(3D matrix: num_ANYdetritus X num_grps X num_MC)
%           fate_predation                                                        	(3D matrix: num_livingANDfleets X num_grps X num_MC)
%           fname_E2E_MonteCarlo        name of this f_E2E_MonteCarlo function
%
% Part 1. Generate a random ConsumptionBudget matrix that defines, for each functional group, the
%           fraction of total consumption directed towards each of seven fates: feces, metabolism, 
%           eggs, predation, senescence, biomass accumulation, and emigration.
%         The distribution of uncertainties for each of the 7 terms is defined as variable 
%           ECOTRAN_PEDIGREE.ConsumptionBudget_CV in code f_E2Epedigree_12022019.
%         The sum of each column of the ConsumptionBudget matrix (ConsumptionBudget_MC)
%           sums to 1 to account for 100% of consumption. Feces, metabolism, and egg production 
%           are limited to values between 0 and 1 (actually held within a range slightly > 0 and 
%           slightly < 1 in order to not reduce predation and senescence to 0). Predation and senescence
%           may be any positive value >= 0. Biomass accumulation and emigration may be any value 
%           positive or negative. Feces values for fleets represents discards. Emigration 
%           values for fleets represent landings removed from the system so must be limited 
%           to values between 0 and 1.
%               ConsumptionBudget_MC possible value ranges:
%                   row 1: feces        (0 <--> 1 only)
%                   row 2: metabolism   (0 <--> 1 only)
%                   row 3: eggs         (0 <--> 1 only)
%                   row 4: predation    (any positive value 0+, can be >1 because of ba & em)
%                   row 5: senescence   (any positive value 0+, can be >1 because of ba & em)
%                   row 6: ba           (any value)
%                   row 7: em           (any value)     (fleets: all cells 0 <--> 1 only)
%         Once a random ConsumptionBudget matrix is drawn, differences from 100% consumption are 
%           resolved:
%               1. The sum of feces, metabolism, and eggs kept within an overall maximum value.
%               2. Somatic production (predation and senescence) is adjusted to bring sum of all 
%                   intra-domain terms (feces, metabolism, eggs, predation, senescence) to 1. Thus, 
%                   total predation and senescence are limited by group physiology (feces, metabolism, 
%                   eggs) rather than by randomly drawing from an assumed level of uncertainty. However, 
%                   the scale of predation and senescence relative to each other are set by the
%                   assumed level of uncertainty in the pedigree.
%               3. If net growth (ba) and/or emigration (em) terms are not zero, then somatic production 
%                   (predation and senescence) is raised or lowered to bring the sum of each consumption 
%                   budget column to 1. Predation and senescence are kept at or above the minimum defined values. 
%                   Any imbalance to the system at this point (too much predation) is subtracted from 
%                   the ba term.
% Part 2. Generate a random ECOTRAN EnergyBudget matrix. 
%       - First, all terms within the random matrix are truncated to positive values only. 
%       - Second, the sum of ammonium terms are scaled to equal the metabolism in the Consumption Budget.
%       - Third, egg terms are scaled to equal the egg term in the Consumption Budget.
%       - Fourth, the sum of detritus terms are scaled to equal feces + senescence in the Consumption Budget.
%       - Fifth, the sum of predation terms are scaled to equal predation in the Consumption Budget
%       - Sixth, the defining "type" model is model 1 in the stack of Monte Carlo models
% NOTES:
%	- Layer 1 of EnergyBudget_MC is the defining "type" model.
%	- This code DOES generate alternate models of discard rates for individual 
%       functional groups and fleets. However, discard fractions for all 
%       individual rates must later be adjusted so that the summed landing
%       and discard rates for individual groups match the ConsumptionBudget_MC
%       landings (em) and discard (feces) rates of each fleet.
%   - only predation fates change, all other fates do not change
%   - FFF can consider fate uncertainties in the future
%
% revision date: 8-4-2020


% *************************************************************************
% STEP 1: unpack variables (leave original variable structures unchanged)--

% step 1a: save name of this version of f_E2E_MonteCarlo function ---------
fname_E2E_MonteCarlo            = mfilename; % name of this f_E2E_MonteCarlo function
display(['Running: ' fname_E2E_MonteCarlo])


% step 1b: MonteCarloConditions -------------------------------------------
num_MC                          = MonteCarloConditions.num_MC;              % number of random Monte Carlo models to generate
DistributionType                = MonteCarloConditions.DistributionType;  	% 'normal' or 'uniform'


% step 1c: ECOTRAN budgets ------------------------------------------------
EnergyBudget                    = ECOTRAN.EnergyBudget(:, :, 1);            % this must be a single 2D production matrix, choose top matrix just in case this is a stack of MonteCarlo models
ConsumptionBudget               = ECOTRAN.ConsumptionBudget;                % rows: 1=feces; 2=metabolism; 3=eggs; 4=predation; 5=senescence; 6=ba; 7=em


% step 1d: ECOTRAN fleet info ---------------------------------------------
LANDINGS                        = ECOTRAN.LANDINGS;                         % landings; (t/km2/y); (2D matrix: num_grps X num_fleets)
DISCARDS                        = ECOTRAN.DISCARDS;                         % discards; (t/km2/y); (2D matrix: num_grps X num_fleets)
CATCH                           = ECOTRAN.CATCH;                            % catch; (t/km2/y); (2D matrix: num_grps X num_fleets)
landings_TotalByFleet           = ECOTRAN.landings_TotalByFleet;            % total landings by fleet; (t/km2/y); (2D matrix: num_grps X num_fleets)
discards_TotalByFleet           = ECOTRAN.discards_TotalByFleet;            % total discards by fleet; (t/km2/y); (2D matrix: num_grps X num_fleets)
catch_TotalByFleet              = ECOTRAN.catch_TotalByFleet;               % total catch by fleet; (t/km2/y); (2D matrix: num_grps X num_fleets)
DiscardFraction                 = ECOTRAN.DiscardFraction;                  % fraction of catch discarded for each group; (3D matrix: num_grps X num_fleets)


% step 1e: ECOTRAN GroupTypes & GroupTypeDefs -----------------------------
GroupType                      	= ECOTRAN.GroupType;	% includes nutrients & fleets; (vertical vector: num_grps X 1)
GroupTypeDef_ANYNitroNutr    	= ECOTRAN.GroupTypeDef_ANYNitroNutr;
GroupTypeDef_ANYPrimaryProd    	= ECOTRAN.GroupTypeDef_ANYPrimaryProd;
GroupTypeDef_ANYConsumer      	= ECOTRAN.GroupTypeDef_ANYConsumer;
GroupTypeDef_fleet           	= ECOTRAN.GroupTypeDef_fleet;
GroupTypeDef_eggs              	= ECOTRAN.GroupTypeDef_eggs;
GroupTypeDef_terminalBnthDetr	= ECOTRAN.GroupTypeDef_terminalBnthDetr;
GroupTypeDef_ANYDetritus      	= ECOTRAN.GroupTypeDef_ANYDetritus;
GroupTypeDef_micrograzers      	= ECOTRAN.GroupTypeDef_micrograzers;
GroupTypeDef_bacteria     	    = ECOTRAN.GroupTypeDef_bacteria;


% step 1f: ECOTRAN fates ----------------------------------------------
%          FFF can add fate uncertainties here in the future
fate_metabolism                 = ECOTRAN.fate_metabolism;	% (3D matrix: num_nutrients X num_grps X num_MC)
fate_eggs                       = ECOTRAN.fate_eggs;        % (3D matrix: num_eggs X num_grps X num_MC)
fate_feces                      = ECOTRAN.fate_feces;       % (3D matrix: num_ANYdetritus X num_grps X num_MC)
fate_senescence                 = ECOTRAN.fate_senescence;	% (3D matrix: num_ANYdetritus X num_grps X num_MC)
fate_predation                  = ECOTRAN.fate_predation;	% (2D matrix: num_livingANDfleets X num_grps)


% step 1g: ECOTRAN_PEDIGREE coefficients of variation ---------------------
EnergyBudget_CV                 = ECOTRAN_PEDIGREE.EnergyBudget_CV;             % (2D matrix: num_grps X num_grps)
ConsumptionBudget_CV            = ECOTRAN_PEDIGREE.ConsumptionBudget_CV;        % rows: 1=feces; 2=metabolism; 3=eggs; 4=predation; 5=senescence; 6=ba; 7=em
LANDINGS_CV                     = ECOTRAN_PEDIGREE.LANDINGS_CV;                 % landings; (CV); (2D matrix: num_grps X num_fleets)
DISCARDS_CV                     = ECOTRAN_PEDIGREE.DISCARDS_CV;                 % discards; (CV); (2D matrix: num_grps X num_fleets)
CATCH_CV                        = ECOTRAN_PEDIGREE.CATCH_CV;                    % catch; (CV); (2D matrix: num_grps X num_fleets)
landings_TotalByFleet_CV        = ECOTRAN_PEDIGREE.landings_TotalByFleet_CV;	% total landings by fleet; (CV); (2D matrix: num_grps X num_fleets)
discards_TotalByFleet_CV        = ECOTRAN_PEDIGREE.discards_TotalByFleet_CV;	% total discards by fleet; (CV); (2D matrix: num_grps X num_fleets)
catch_TotalByFleet_CV           = ECOTRAN_PEDIGREE.catch_TotalByFleet_CV;       % total catch by fleet; (CV); (2D matrix: num_grps X num_fleets)
DiscardFraction_CV              = ECOTRAN_PEDIGREE.DiscardFraction_CV;          % discard fraction of each group by each fleet; (CV); (2D matrix: num_grps X num_fleets)
% *************************************************************************





% *************************************************************************
% STEP 2: prepare to create MonteCarlo arrays------------------------------

% step 2a: ECOTRAN group addresses & numbers ------------------------------
looky_nutrients                	= find(floor(GroupType) == GroupTypeDef_ANYNitroNutr); % row addresses of nutrients
looky_ANYPrimaryProd          	= find(floor(GroupType) == GroupTypeDef_ANYPrimaryProd);
looky_ANYconsumer            	= find(floor(GroupType) == GroupTypeDef_ANYConsumer);
looky_bacteria             	    = find(floor(GroupType) == GroupTypeDef_bacteria);
looky_fleets                  	= find(floor(GroupType) == GroupTypeDef_fleet);
looky_livingANDfleets         	= [looky_ANYPrimaryProd; looky_ANYconsumer; looky_bacteria; looky_fleets]; % includes primary producers & bacteria
looky_eggs                   	= find(GroupType == GroupTypeDef_eggs);
looky_ANYdetritus              	= find(floor(GroupType) == GroupTypeDef_ANYDetritus);
looky_terminalBNTHdetritus   	= find(GroupType == GroupTypeDef_terminalBnthDetr);

num_nutrients                   = length(looky_nutrients);
num_eggs                        = length(looky_eggs);
num_detritus                    = length(looky_ANYdetritus);
num_fleets                      = length(looky_fleets);
num_ANYconsumer                	= length(looky_ANYconsumer);
num_livingANDfleets             = length(looky_livingANDfleets);
[num_grps, toss]                = size(EnergyBudget);


% step 2b: initialize result arrays ---------------------------------------
EnergyBudget_MC             = zeros(num_grps, num_grps, num_MC);
ConsumptionBudget_MC        = zeros(7, num_grps, num_MC);
DiscardFraction_MC          = zeros(num_grps, num_fleets, num_MC);


% step 2c: Zero-cell addresses of EnergyBudget & ConsumptionBudget --------
%          NOTE: entries are 0 where consumption budget in model are intentionally defined to be 0
EnergyBudget_ZeroScaler                                         = EnergyBudget;
EnergyBudget_ZeroScaler(EnergyBudget_ZeroScaler>0)              = 1;                   	% (NOTE: all EnergyBudget cells are >=0)
ConsumptionBudget_ZeroScaler                                	= ConsumptionBudget;	% (2D matrix: 7 X num_grps)
ConsumptionBudget_ZeroScaler(ConsumptionBudget_ZeroScaler>0)   	= 1;                  	% (2D matrix: 7 X num_grps)
ConsumptionBudget_ZeroScaler(ConsumptionBudget_ZeroScaler<0)	= 1;                    % NOTE: some ConsumptionBudget cells (ba & em) are <0
ConsumptionBudget_ZeroScaler([1 7], looky_fleets)               = 1;                 	% allow discards (feces) and landings (em) for all fleets; (2D matrix: 7 X num_grps)


% step 2d: assign minimum & maximum terms for ConsumptionBudget -----------
%          NOTE: SSS these are assumed values purposely set to unrealistically wide ranges
min_FECES                       = 0.05;     % absolute minimum feces budget
min_MTBLSM                      = 0.05;     % absolute minimum respiration budget
min_EGGS                        = 0.0001;   % absolute minimum egg budget
min_PREDATION                   = 0.0001; 	% absolute minimum predation budget
min_SENESCENCE                  = 0.0001; 	% absolute minimum senescence budget

min_relative_ba                 = 0.05;     % minimum ba RELATIVE to "type" model ConsumptionBudget ba values
min_relative_em                 = 0.05;     % minimum ba RELATIVE to "type" model ConsumptionBudget em values

min_relative_LANDINGS           = 0.05;     % minimum landings RELATIVE to "type" model ConsumptionBudget landings (em) values; NOTE: FFF could make this negative to simulate fish stocking programs
min_relative_DISCARDS           = 0.05;     % minimum discards RELATIVE to "type" model ConsumptionBudget discards (feces) values;
min_DiscardFraction             = 0;        % absolute minimum fleet DiscardFraction for each group; NOTE: could use minimum and maximum values relative to "type" model DiscardFraction values

min_relative_ba                 = repmat(min_relative_ba, [1, num_grps]);           % minimum ba RELATIVE to "type" model ConsumptionBudget ba values; (horizontal vector: 1 X num_grps)
min_ba                          = min_relative_ba .* ConsumptionBudget(6, :);       % minimum ba ABSOLUTE values; (horizontal vector: 1 X num_grps)
min_relative_em                 = repmat(min_relative_em, [1, num_grps]);           % minimum em RELATIVE to "type" model ConsumptionBudget ba values; (horizontal vector: 1 X num_grps)
min_em                          = min_relative_em .* ConsumptionBudget(7, :);       % minimum em ABSOLUTE values; (horizontal vector: 1 X num_grps)
min_relative_LANDINGS           = repmat(min_relative_LANDINGS, [1, num_fleets]);	% minimum landings RELATIVE to "type" model ConsumptionBudget em values; (horizontal vector: 1 X num_fleets)
min_LANDINGS                    = min_relative_LANDINGS .* ConsumptionBudget(7, looky_fleets);	% minimum landings ABSOLUTE values; (horizontal vector: 1 X num_fleets)
min_relative_DISCARDS           = repmat(min_relative_DISCARDS, [1, num_fleets]);	% minimum discards RELATIVE to "type" model ConsumptionBudget discards (feces) values; (horizontal vector: 1 X num_fleets)
min_DISCARDS                    = min_relative_DISCARDS .* ConsumptionBudget(1, looky_fleets);	% minimum discards ABSOLUTE values; (horizontal vector: 1 X num_fleets)

max_FECES                       = 0.95;                             % absolute maximum feces budget
max_MTBLSM                      = 0.95;                             % absolute maximum respiration budget
max_EGGS                        = 0.95;                             % absolute maximum egg budget

max_LANDINGS                    = 1 - min_DISCARDS;                 % absolute maximum fleet landings (em) budget; (horizontal vector: 1 X num_fleets)
max_DISCARDS                    = 1 - min_LANDINGS;                 % absolute maximum fleet discards (feces) budget; (horizontal vector: 1 X num_fleets)
max_DiscardFraction             = 1;                                % absolute maximum fleet DiscardFraction for each group; NOTE: could use minimum and maximum values relative to "type" model DiscardFraction values

max_FME                         = 0.95;                             % absolute maximum budget for FME (FME = feces + metabolism + eggs)
max_FME                         = repmat(max_FME, [1, num_grps]);	% (horizontal vector: 1 X num_grps)
max_FME(looky_fleets)           = max_DISCARDS;                   	% apply fleet discard maximum; (horizontal vector: 1 X num_grps)
max_FME(looky_ANYdetritus)      = 1;                              	% apply detritus metabolism max; (horizontal vector: 1 X num_grps)

min_ConsumptionBudget                       = [min_FECES; min_MTBLSM; min_EGGS; min_PREDATION; min_SENESCENCE; 1; 1]; % (vertical vector: 7 X 1)
min_ConsumptionBudget                       = repmat(min_ConsumptionBudget, [1, num_grps]);	% (2D matrix: 7 X num_grps)
min_ConsumptionBudget([6 7], :)             = [min_ba; min_em];                                         % set minimum ba & em terms in min_ConsumptionBudget
min_ConsumptionBudget(1, looky_fleets)      = min_DISCARDS;                                             % apply fleet discard (feces) minimum
min_ConsumptionBudget(7, looky_fleets)      = min_LANDINGS;                                             % apply fleet landings (em) minimum
min_ConsumptionBudget                       = min_ConsumptionBudget .* ConsumptionBudget_ZeroScaler;	% zero-out elements that have 0 value in the "type" model; (2D matrix: 7 X num_grps)

max_ConsumptionBudget                       = [max_FECES; max_MTBLSM; max_EGGS; 1; 1; 1; 1];            % (vertical vector: 7 X 1)
max_ConsumptionBudget                       = repmat(max_ConsumptionBudget, [1, num_grps]);             % (2D matrix: 7 X num_grps)
max_ConsumptionBudget(1, looky_fleets)      = max_DISCARDS;                                             % apply fleet discard (feces) maximum
max_ConsumptionBudget(7, looky_fleets)      = max_LANDINGS;                                             % apply fleet landings (em) maximum
max_ConsumptionBudget(2, looky_ANYdetritus)	= 1;                                                        % apply detritus metabolism maximum
max_ConsumptionBudget                       = max_ConsumptionBudget .* ConsumptionBudget_ZeroScaler;	% (2D matrix: 7 X num_grps)


% step 2e: minimum-cell addresses for MonteCarlo EnergyBudget -------------
min_EnergyBudget                = 0.001 * EnergyBudget;     % SSS assumed minimum EnergyBudget values to prevent trace values from going to 0


% step 2f: assign minimum & maximum values to fleet DiscardFraction -------
min_DiscardFraction             = repmat(min_DiscardFraction, [num_grps, num_fleets]);
max_DiscardFraction             = repmat(max_DiscardFraction, [num_grps, num_fleets]);
% *************************************************************************





% *************************************************************************
% STEP 3: prepare for uniform distributions--------------------------------
%         ConsumptionBudget, EnergyBudget, & DiscardFraction
if strcmp(DistributionType, 'uniform')
    
    % step 3a: interval over which to draw random ConsumptionBudget -------
    interval_low_ConsumptionBudget      = ConsumptionBudget - (ConsumptionBudget .* ConsumptionBudget_CV);
    interval_high_ConsumptionBudget     = ConsumptionBudget + (ConsumptionBudget .* ConsumptionBudget_CV);
    % ---------------------------------------------------------------------
    
    
    % step 3b: filter-out improper ConsumptionBudget interval values ------
    %          row 1: feces        (0 <--> 1 only)
    %          row 2: metabolism   (0 <--> 1 only)
    %          row 3: eggs         (0 <--> 1 only)
    %          row 4: predation    (any positive value 0+, can be >1 because of ba & em)
    %          row 5: senescence   (any positive value 0+, can be >1 because of ba & em)
    %          row 6: ba           (any value)
    %          row 7: em           (any value)     (fleets: all cells 0 <--> 1 only)

    % filter out all negative values (except for ba & em rows)
    CheckInterval_low                                    = interval_low_ConsumptionBudget(1:5, :);
    CheckInterval_high                                   = interval_high_ConsumptionBudget(1:5, :);
    CheckInterval_low(CheckInterval_low < 0)             = 0;
    CheckInterval_high(CheckInterval_high < 0)           = 0;
    interval_low_ConsumptionBudget(1:5, :)               = CheckInterval_low;
    interval_high_ConsumptionBudget(1:5, :)              = CheckInterval_high;

    % feces, metabolism, & eggs filter (no values > 1) (assumption that ba & em terms are compensated for by ecology of individuals and not by changes to physiology)
    CheckInterval_low                                    = interval_low_ConsumptionBudget(1:3, :);
    CheckInterval_high                                   = interval_high_ConsumptionBudget(1:3, :);
    CheckInterval_low(CheckInterval_low > 1)             = 1;
    CheckInterval_high(CheckInterval_high > 1)           = 1;
    interval_low_ConsumptionBudget(1:3, :)               = CheckInterval_low;
    interval_high_ConsumptionBudget(1:3, :)              = CheckInterval_high;

    % % FFF set min and max parameter range--- metabolism (assumption that ba & em terms are compensated for by ecology of individuals and not by changes to physiology)
    % CheckInterval_low                                   = interval_low_ConsumptionBudget(2, :);
    % CheckInterval_high                                  = interval_high_ConsumptionBudget(2, :);
    % CheckInterval_low(CheckInterval_low < min_MTBLSM)   = min_MTBLSM;
    % CheckInterval_high(CheckInterval_high < min_MTBLSM) = min_MTBLSM;
    % CheckInterval_low(CheckInterval_low > max_MTBLSM)   = max_MTBLSM;
    % CheckInterval_high(CheckInterval_high > max_MTBLSM) = max_MTBLSM;
    % interval_low_ConsumptionBudget(2, :)                = CheckInterval_low;
    % interval_high_ConsumptionBudget(2, :)               = CheckInterval_high;

    % fleet group filter (em = landings are always between 0 and 1)
    CheckInterval_low                                   = interval_low_ConsumptionBudget(7, looky_fleets);
    CheckInterval_high                                  = interval_high_ConsumptionBudget(7, looky_fleets);
    CheckInterval_low(CheckInterval_low < 0)            = 0;
    CheckInterval_high(CheckInterval_high < 0)          = 0;
    CheckInterval_low(CheckInterval_low > 1)            = 1;
    CheckInterval_high(CheckInterval_high > 1)          = 1;
    interval_low_ConsumptionBudget(7, looky_fleets)     = CheckInterval_low;
    interval_high_ConsumptionBudget(7, looky_fleets)	= CheckInterval_high;
    % ---------------------------------------------------------------------

    
    % step 3c: interval over which to draw random EnergyBudget ------------
    interval_low_EnergyBudget                           = EnergyBudget - (EnergyBudget .* EnergyBudget_CV);
    interval_high_EnergyBudget                          = EnergyBudget + (EnergyBudget .* EnergyBudget_CV);
    % ---------------------------------------------------------------------

    
    % step 3d: filter out improper EnergyBudget interval values -----------
    %          metabolism                      (0 <--> 1 only)
    %          eggs                            (0 <--> 1 only)
    %          predation                       (any positive value 0+, can be >1 because of ba & em)
    %          detritus (feces + senescence)   (any positive value 0+, can be >1 because of ba & em)
    %          any fleet cell                  (0 <--> 1 only) 

    % filter out all negative values (note that EnergyBudget does not have rows for ba & em)
    interval_low_EnergyBudget(interval_low_EnergyBudget < 0)        = 0;
    interval_high_EnergyBudget(interval_high_EnergyBudget < 0)      = 0;

    % filter out metabolism & eggs >1
    CheckInterval_low                                               = interval_low_EnergyBudget([looky_nutrients; looky_eggs], :);
    CheckInterval_high                                              = interval_high_EnergyBudget([looky_nutrients; looky_eggs], :);
    CheckInterval_low(CheckInterval_low > 1)                        = 1;
    CheckInterval_high(CheckInterval_high > 1)                      = 1;
    interval_low_EnergyBudget([looky_nutrients; looky_eggs], :)     = CheckInterval_low;
    interval_high_EnergyBudget([looky_nutrients; looky_eggs], :)	= CheckInterval_high;

    % fleet group filter, no values >1
    CheckInterval_low                                               = interval_low_EnergyBudget(:, looky_fleets);
    CheckInterval_high                                              = interval_high_EnergyBudget(:, looky_fleets);
    CheckInterval_low(CheckInterval_low > 1)                        = 1;
    CheckInterval_high(CheckInterval_high > 1)                      = 1;
    interval_low_EnergyBudget(:, looky_fleets)                      = CheckInterval_low;
    interval_high_EnergyBudget(:, looky_fleets)                     = CheckInterval_high;
    % ---------------------------------------------------------------------
    
    
    % step 3e: interval over which to draw random DiscardFraction ---------
    interval_low_DiscardFraction                = DiscardFraction - (DiscardFraction .* DiscardFraction_CV); % (scaler: 0 to 1); (2D matrix: num_grps X num_fleets)
    interval_high_DiscardFraction               = DiscardFraction + (DiscardFraction .* DiscardFraction_CV); % (scaler: 0 to 1); (2D matrix: num_grps X num_fleets)
    % ---------------------------------------------------------------------
    
    
    % step 3f: filter-out improper DiscardFraction interval values --------
    % filter out all negative values
    CheckInterval_low                           = interval_low_DiscardFraction;
    CheckInterval_high                          = interval_high_DiscardFraction;
    CheckInterval_low(CheckInterval_low < 0)	= 0;
    CheckInterval_high(CheckInterval_high < 0)	= 0;
    interval_low_DiscardFraction                = CheckInterval_low;
    interval_high_DiscardFraction               = CheckInterval_high;

    % filter out values > 1
    CheckInterval_low                           = interval_low_DiscardFraction;
    CheckInterval_high                          = interval_high_DiscardFraction;
    CheckInterval_low(CheckInterval_low > 1)	= 1;
    CheckInterval_high(CheckInterval_high > 1)	= 1;
    interval_low_DiscardFraction                = CheckInterval_low;
    interval_high_DiscardFraction               = CheckInterval_high;    
    % ---------------------------------------------------------------------
    
end % if strcmp(DistributionType, 'uniform')
% *************************************************************************




% *************************************************************************
% STEP 4: generate random models and production vectors--------------------
for trial_loop = 2:num_MC
    
    % step 4a: generate a random ConsumptionBudget ------------------------
    %          row 1: feces
    %          row 2: metabolism
    %          row 3: eggs
    %          row 4: predation
    %          row 5: senescence
    %          row 6: ba
    %          row 7: em           (NOTE: a negative value is an input into the domain)
    if strcmp(DistributionType, 'uniform')
        % for uniform distribution
        random_matrix           = rand(7, num_grps); % draw from UNIFORM distribution
        MonteCarloArray     	= interval_low_ConsumptionBudget + ((interval_high_ConsumptionBudget - interval_low_ConsumptionBudget) .* random_matrix); % produce the new random ConsumptionBudget
    else
        % for normal distribution
        random_matrix           = randn(7, num_grps); % draw from NORMAL distribution
        MonteCarloArray      	= ConsumptionBudget + ((ConsumptionBudget .* ConsumptionBudget_CV) .* random_matrix);
    end % if strcmp(DistributionType, 'uniform')

    % A) truncate all below-minimum values (except for ba & em rows)
      CheckInterval                             = MonteCarloArray(1:5, :);
      MinInterval                               = min_ConsumptionBudget(1:5, :);
      CompInterval                           	= CheckInterval - MinInterval;
      looky_TooSmall                            = find(CompInterval < 0);
      CheckInterval(looky_TooSmall)             = MinInterval(looky_TooSmall);
      MonteCarloArray(1:5, :)                   = CheckInterval;
      
    % B) truncate all above-maximum feces, metabolism, & egg values
      CheckInterval                             = MonteCarloArray(1:3, :);
      MaxInterval                               = max_ConsumptionBudget(1:3, :);
      CompInterval                           	= MaxInterval - CheckInterval;
      looky_TooBig                              = find(CompInterval < 0);
      CheckInterval(looky_TooBig)               = MaxInterval(looky_TooBig);
      MonteCarloArray(1:3, :)                   = CheckInterval;
      
    % C) truncate FME (feces + metabolism + eggs); find FME > max_FME and renormalize FME to max_FME but obey minimum standards
      CheckInterval                             = MonteCarloArray(1:3, :);          % feces, metabolism, and egg terms; (2D matrix: 3 X num_grps)
      MinInterval                               = min_ConsumptionBudget(1:3, :);    % minimum allowed value for each term (feces, metabolism, eggs)
      dfrInternal                               = CheckInterval - MinInterval;      % amount that each term (feces, metabolism, eggs) is above minimum allowed value and can be changed
      total_FME                                 = sum(CheckInterval);               % (horizontal vector: 1 X num_grps)
      dfr_FME                                   = (max_FME - total_FME);            % amount above max_FME that must be reduced; positive values set to 0; (horizontal vector: 1 X num_grps)
      looky_positive                            = find(dfr_FME > 0);
      dfr_FME(looky_positive)                   = 0;
      total_FME                                 = sum(dfrInternal);                 % recalculate using dfrInternal; (horizontal vector: 1 X num_grps)
      total_FME                                 = repmat(total_FME, [3, 1]);        % (2D matrix: 3 X num_grps)
      dfr_FME                                   = repmat(dfr_FME, [3, 1]);          % (2D matrix: 3 X num_grps)
      adjusted_FME                              = CheckInterval + (dfrInternal ./ total_FME) .* dfr_FME; % (2D matrix: 3 X num_grps)
      adjusted_FME(isnan(adjusted_FME))         = 0;                                % set div/0 NaNs to 0; (2D matrix: 3 X num_grps)
      
      % correct for cases where all FME terms are below minimum values
      total_FME                                 = sum(adjusted_FME);                % (horizontal vector: 1 X num_grps)
      looky_zero                                = find(total_FME == 0);             % these are the clms where either all FME terms are below minimu values or that FME terms are not defined (e.g., phytoplankton and fleets)
      FME_fix                                   = zeros(1, num_grps);               % (horizontal vector: 1 X num_grps)
      FME_fix(looky_zero)                       = 1;
      FME_fix                                   = repmat(FME_fix, [3, 1]);          % (2D matrix: 3 X num_grps)
      FME_fix                                   = FME_fix .* ConsumptionBudget_ZeroScaler(1:3, :);	% zero-out cases where feces, metabolism, or eggs are defined as zero in the model
      FME_fix                                   = FME_fix .* MinInterval;           % replace any missing values with minimum values
      adjusted_FME                              = adjusted_FME + FME_fix;
      MonteCarloArray(1:3, :)                   = adjusted_FME;
      
    % D) adjust somatic production (predation & senescence) to bring sum of intra-domain parameters (feces, metabolism, eggs, predation, senescence) to 1
      CheckInterval                             = MonteCarloArray(4:5, :);          % predation and senescence terms; (2D matrix: 2 X num_grps)
      total_PS                                  = sum(CheckInterval);               % current sum of predation and senescence terms; (horizontal vector: 1 X num_grps)
      total_PS                                  = repmat(total_PS, [2, 1]);         % (2D matrix: 2 X num_grps)
      total_FME                                 = sum(MonteCarloArray(1:3, :));     % sum of FME terms (feces, metabolism, eggs); (horizontal vector: 1 X num_grps)
      dfr_PS                                    = 1 - total_FME;                    % value for sum of PS
      dfr_PS                                    = repmat(dfr_PS, [2, 1]);           % (2D matrix: 2 X num_grps)
      adjusted_PS                               = (CheckInterval ./ total_PS) .* dfr_PS;	% (2D matrix: 2 X num_grps)
      adjusted_PS(isnan(adjusted_PS))           = 0;                                % set div/0 NaNs to 0 (e.g., fleets); (2D matrix: 2 X num_grps)
      MonteCarloArray(4:5, :)                   = adjusted_PS;

    % E) balance net growth and emigration (ba + em) with somatic production (predation and senescence)
      total_PS                                  = sum(MonteCarloArray(4:5, :));  	% total somatic production; (horizontal vector: 1 X num_grps);
      total_BAEM                                = sum(MonteCarloArray(6:7, :));    	% total ba & em (biomass accumulation and emigration); (horizontal vector: 1 X num_grps); 
      pred_fraction                             = MonteCarloArray(4, :) ./ total_PS;	% fraction of somatic production going to predation
      pred_fraction(isnan(pred_fraction))       = 0;                             	% set div/0 NaNs to 0
      scen_fraction                             = 1 - pred_fraction;               	% fraction of somatic production going to senescence
      new_PS                                    = total_PS - total_BAEM;            % subtract net BAEM from total somatic production (predation and senescence); (horizontal vector: 1 X num_grps)
      
      % filter out instances of somatic production (= predation + senescence) being reduced below a minimum assumed value
      min_PS                                    = 0.1 * total_PS;                	% SSS assumed minimum somatic production can only be cut by 90% to accomodate ba + em; (horizontal vector: 1 X num_grps)
      looky_TooSmall                            = find(new_PS < min_PS);
      new_PS(looky_TooSmall)                    = min_PS(looky_TooSmall);
      
      % distribute new_PS between predation and senescence and enforce defined minimum P & S values
      adjusted_P                                = new_PS .* pred_fraction;          % (horizontal vector: 1 X num_grps)
      adjusted_S                                = new_PS .* scen_fraction;          % (horizontal vector: 1 X num_grps)
      min_P                                     = min_ConsumptionBudget(4, :);      % minimum predation term with defined zero terms; (horizontal vector: 1 X num_grps)
      CompInterval_P                            = adjusted_P - min_P;               % (horizontal vector: 1 X num_grps)
      looky_TooSmall                            = find(CompInterval_P < 0);         % clm addresses of predation terms that are less than the defined minimum
      adjusted_P(looky_TooSmall)                = min_P(looky_TooSmall);
      adjusted_S(looky_TooSmall)                = adjusted_S(looky_TooSmall) + CompInterval_P(looky_TooSmall); % reduce senescence by amount predation was raised
      min_S                                     = min_ConsumptionBudget(5, :);      % minimum senescence term with defined zero terms; (horizontal vector: 1 X num_grps)
      CompInterval_S                            = adjusted_S - min_S;               % (horizontal vector: 1 X num_grps)
      looky_TooSmall                            = find(CompInterval_S < 0);         % clm addresses of predation terms that are less than the defined minimum
      adjusted_S(looky_TooSmall)                = min_S(looky_TooSmall);
      adjusted_P(looky_TooSmall)                = adjusted_P(looky_TooSmall) + CompInterval_S(looky_TooSmall); % reduce predation by amount senescence was raised
      MonteCarloArray(4, :)                     = adjusted_P;
      MonteCarloArray(5, :)                     = adjusted_S;
      
    % F) put remaining imbalance onto ba
      sum_MonteCarloArray                       = sum(MonteCarloArray);           	% (horizontal vector: 1 X num_grps);
      imbalance                                 = sum_MonteCarloArray - 1;       	% positive imbalance means too much export, take from ba
      MonteCarloArray(6, :)                     = MonteCarloArray(6, :) - imbalance;
      
    % G) fleet landings (discards were brought into range above)
      MonteCarloArray(6, looky_fleets)          = 0;                                    % set ba for fleets to 0
      MonteCarloArray(7, looky_fleets)          = 1 - MonteCarloArray(1, looky_fleets);	% landings (em) are defined by the random draw of discards
      
	% H) apply ConsumptionBudget_ZeroScaler to clean up tiny rounding error in ba terms and makes sure we didn't create ba & em terms where none originally existed
      MonteCarloArray                           = MonteCarloArray .* ConsumptionBudget_ZeroScaler;	% this step cleans up tiny rounding error in ba terms and makes sure we didn't create ba & em terms where none originally existed
    
	% I) final normalization step to make sure MonteCarlo_ConsumptionBudget clms all sum to 1
    %    NOTE: may bring some terms slightly out of the defined min & max range
      total_MonteCarloArray                     = sum(MonteCarloArray);
      total_MonteCarloArray                     = repmat(total_MonteCarloArray, [7, 1]);
      MonteCarloArray                           = MonteCarloArray ./ total_MonteCarloArray;        	
      MonteCarloArray(isnan(MonteCarloArray))   = 0;	% set div/0 NaNs to 0
      MonteCarlo_ConsumptionBudget              = MonteCarloArray; % assign MonteCarloArray to MonteCarlo_ConsumptionBudget
    % ---------------------------------------------------------------------
    
    
    
    % step 4b: generate a random EnergyBudget -----------------------------
    if strcmp(DistributionType, 'uniform')
        % for uniform distribution
        random_matrix        	= rand(num_grps, num_grps); % draw from a UNIFORM distribution
        MonteCarloArray       	= interval_low_EnergyBudget + ((interval_high_EnergyBudget - interval_low_EnergyBudget) .* random_matrix); % produce the new random EnergyBudget    
    else
        % for normal distribution
        random_matrix        	= randn(num_grps, num_grps); % draw from a NORMAL distribution
        MonteCarloArray       	= EnergyBudget + ((EnergyBudget .* EnergyBudget_CV) .* random_matrix);
    end % if strcmp(DistributionType, 'uniform')
    
    % A) re-zero cells that should be 0 (this line should not be needed)
      MonteCarloArray                           = MonteCarloArray .* EnergyBudget_ZeroScaler;
    
    % B) filter out negative values
	  looky_TooSmall                            = find(MonteCarloArray < 0);
	  MonteCarloArray(looky_TooSmall)           = min_EnergyBudget(looky_TooSmall);
    
	% C) set feces, metabolism, eggs, & senescence to MonteCarlo_ConsumptionBudget values
      MonteCarloArray(looky_nutrients, :)       = repmat(MonteCarlo_ConsumptionBudget(2, :), [num_nutrients 1])	.* fate_metabolism;
      MonteCarloArray(looky_eggs, :)            = repmat(MonteCarlo_ConsumptionBudget(3, :), [num_eggs 1])      .* fate_eggs;
      MonteCarlo_feces                          = repmat(MonteCarlo_ConsumptionBudget(1, :), [num_detritus 1])	.* fate_feces;
	  MonteCarlo_senescence                     = repmat(MonteCarlo_ConsumptionBudget(5, :), [num_detritus 1])	.* fate_senescence;
	  MonteCarlo_detritus                       = MonteCarlo_feces + MonteCarlo_senescence;
      MonteCarloArray(looky_ANYdetritus, :)     = MonteCarlo_detritus;      % paste in flows to detritus pools into the MonteCarloArray (EnergyBudget)
      
    % D) normalize predation terms to ConsumptionBudget value
    %    NOTE: special case for terminal benthic detritus group (no senescence in EnergyBudget), increase predation to compensate
      tempMonteCarlo_ConsumptionBudget          = MonteCarlo_ConsumptionBudget;
      tempMonteCarlo_ConsumptionBudget(4, looky_terminalBNTHdetritus) = sum(MonteCarlo_ConsumptionBudget([4 5], looky_terminalBNTHdetritus)); % add senescence to predation for terminal benthic detritus group
      tempMonteCarlo_ConsumptionBudget(5, looky_terminalBNTHdetritus) = 0;
      MonteCarlo_predation                      = sum(MonteCarloArray(looky_livingANDfleets, :));
      scaler_predation                          = tempMonteCarlo_ConsumptionBudget(4, :) ./ MonteCarlo_predation;
      scaler_predation(isnan(scaler_predation))	= 0; % filter div/0 NaNs
      MonteCarloArray(looky_livingANDfleets, :)	= MonteCarloArray(looky_livingANDfleets, :) .* repmat(scaler_predation, [num_livingANDfleets 1]); % (2D matrix: num_grps X num_grps)
      
    % E) assign MonteCarloArray to MonteCarlo_EnergyBudget
      MonteCarlo_EnergyBudget                   = MonteCarloArray;
    % ---------------------------------------------------------------------
    
    
    
    % step 4c: generate a random DiscardFraction --------------------------
    %          NOTE: DiscardFraction defines discard fractions for the
    %          catches of individual functional groups and fleets. This is 
    %          NOT NECESSARILY in agreement with the ConsumptionBudget terms
    %          for each fleet. But, we cannot make adjustments to absolute
    %          landings and discard rates for any group or fleet until we
    %          know (run model to calculate) the production rates of each
    %          group. At that time, in f_StaticScenario (QQQ and
    %          DynamicScenario), the DiscardFraction is readjusted so that
    %          summed landings & discards for each fleet match the fleet
    %          feces & em terms in the ConsumptionBudget
    if strcmp(DistributionType, 'uniform')
        % for uniform distribution
        random_matrix        	= rand(num_grps, num_fleets); % draw from a UNIFORM distribution
        MonteCarloArray       	= interval_low_DiscardFraction + ((interval_high_DiscardFraction - interval_low_DiscardFraction) .* random_matrix); % produce the new random DiscardFraction    
    else
        % for normal distribution
        random_matrix        	= randn(num_grps, num_fleets); % draw from a NORMAL distribution
        MonteCarloArray       	= DiscardFraction + ((DiscardFraction .* DiscardFraction_CV) .* random_matrix); % produce the new random DiscardFraction
    end % if strcmp(DistributionType, 'uniform')

    % A) truncate all below-minimum values
      CheckInterval                             = MonteCarloArray;
      MinInterval                               = min_DiscardFraction;
      CompInterval                           	= CheckInterval - MinInterval;
      looky_TooSmall                            = find(CompInterval < 0);
      CheckInterval(looky_TooSmall)             = MinInterval(looky_TooSmall);
      MonteCarloArray                           = CheckInterval;

    % B) truncate all above-maximum values
      CheckInterval                             = MonteCarloArray;
      MaxInterval                               = max_DiscardFraction;
      CompInterval                           	= MaxInterval - CheckInterval;
      looky_TooBig                              = find(CompInterval < 0);
      CheckInterval(looky_TooBig)               = MaxInterval(looky_TooBig);
      MonteCarloArray                           = CheckInterval;    
    
	% C) assign MonteCarloArray to MonteCarlo_DiscardFraction
      MonteCarlo_DiscardFraction                = MonteCarloArray; % (2D matrix: num_grps X num_fleets)
    % ---------------------------------------------------------------------
    
    
    
    % step 4d: store results ----------------------------------------------
	EnergyBudget_MC(:, :, trial_loop)           = MonteCarlo_EnergyBudget;
    ConsumptionBudget_MC(1:7, :, trial_loop)	= MonteCarlo_ConsumptionBudget;
    DiscardFraction_MC(:, :, trial_loop)        = MonteCarlo_DiscardFraction;
    % ---------------------------------------------------------------------
    
end % (trial_loop)
% *************************************************************************





% *************************************************************************
% STEP 5: store defining "type" model in layer 1---------------------------
EnergyBudget_MC(:, :, 1)            = EnergyBudget;         % store defining "type" model in layer 1
ConsumptionBudget_MC(:, :, 1)       = ConsumptionBudget;	% store defining "type" model in layer 1
DiscardFraction_MC(:, :, 1)         = DiscardFraction;      % store defining "type" model in layer 1
% *************************************************************************





% *************************************************************************
% STEP 6: calculate Monte Carlo fates--------------------------------------
% step 6a: replicate fates of feces, metabolism, eggs, & senescence -------
%          FFF remove this step when fate uncertainties are considered
fate_metabolism                         = repmat(fate_metabolism, [1, 1, num_MC]); 	% (3D matrix: num_nutrients X num_grps X num_MC)
fate_eggs                               = repmat(fate_eggs, [1, 1, num_MC]);      	% (3D matrix: num_eggs X num_grps X num_MC)
fate_feces                              = repmat(fate_feces, [1, 1, num_MC]);    	% (3D matrix: num_ANYdetritus X num_grps X num_MC)
fate_senescence                         = repmat(fate_senescence, [1, 1, num_MC]);	% (3D matrix: num_ANYdetritus X num_grps X num_MC)


% step 6b: calculate fate_predation----------------------------------------
%          this is the distribution of the total predation of each producer among all its consumers
sum_predation                         	= sum(EnergyBudget_MC(looky_livingANDfleets, :, :), 1); % (3D matrix: 1 X num_grps X num_MC)
fate_predation                          = EnergyBudget_MC(looky_livingANDfleets, :, :) ./ repmat(sum_predation, [num_livingANDfleets, 1, 1]); % (3D matrix: num_livingANDfleets X num_grps X num_MC)
fate_predation(isnan(fate_predation))	= 0; % correct div/0 errors
% *************************************************************************





% *************************************************************************
% STEP 7: pack results-----------------------------------------------------
ECOTRAN_MC.EnergyBudget_MC              = EnergyBudget_MC;              % (3D matrix: num_grps X num_grps X num_MC)
ECOTRAN_MC.ConsumptionBudget_MC         = ConsumptionBudget_MC;         % (3D matrix: 7 X num_grps X num_MC)
ECOTRAN_MC.DiscardFraction_MC           = DiscardFraction_MC;           % (3D matrix: num_grps X num_fleets X num_MC)

ECOTRAN_MC.fate_metabolism              = fate_metabolism;              % (3D matrix: num_nutrients X num_grps X num_MC)
ECOTRAN_MC.fate_eggs                    = fate_eggs;                    % (3D matrix: num_eggs X num_grps X num_MC)
ECOTRAN_MC.fate_feces                   = fate_feces;                   % (3D matrix: num_ANYdetritus X num_grps X num_MC)
ECOTRAN_MC.fate_senescence              = fate_senescence;              % (3D matrix: num_ANYdetritus X num_grps X num_MC)
ECOTRAN_MC.fate_predation               = fate_predation;               % (3D matrix: num_livingANDfleets X num_grps X num_MC)

ECOTRAN_MC.fname_E2E_MonteCarlo         = fname_E2E_MonteCarlo;         % name of this f_E2E_MonteCarlo function
% *************************************************************************


% end m-file***************************************************************