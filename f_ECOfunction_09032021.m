function [SINGLE_ECO] = f_ECOfunction_09032021(ModelDefinitions, CurrentModel)
% Generate a single ECOTRAN model for 1 "type" EwE model or 1 MonteCarlo EwE model.
%   NOTE: Cannibalism terms are removed from the EnergyBudget matrix. To account
%         for this, feces & metabolism terms are increased by the amount of
%         cannibalism in the diet. This effectively reduces Transfer Efficiency by
%         the correct amount.
%   NOTE: A correction is made for egg production in food web models where eggs are
%         parameterized as a "detritus" production term (see STEP 4)
%
% by Jim Ruzicka
%
% calls:
%       f_RedistributeCannibalism_11202019      remove cannibalism terms on diagonal of matrix EwE_diet
%       f_calcEE_12292020                       calculate Ecotrophic Efficiency 
%       f_calcPredationBudget_12102019      	for each producer, p, and consumer, c: ((b_pc * Q_c) / M2_p)) = the fraction of total predation on each producer, p, going to each consumer, c; (2D matrix: num_grps X num_grps; consumers X producers)
%
% takes:
%       ModelDefinitions
%       CurrentModel
%
% returns:
%       four budget matrices:
%           SINGLE_ECO.BioenergeticBudget	(2D matrix: 3 X num_grps)
%                                               1) feces
%                                               2) metabolism
%                                               3) production
%           SINGLE_ECO.ProductionBudget     (2D matrix: 5 X num_grps)
%                                               1) eggs (reproduction)
%                                               2) predation
%                                               3) senescence
%                                               4) ba (biomass accumulation)
%                                               5) em (emigration); NOTE: negative for immigration
%           SINGLE_ECO.ConsumptionBudget	(2D matrix: 7 X num_grps)
%                                               1) feces
%                                               2) metabolism
%                                               3) eggs (reproduction)
%                                               4) predation
%                                               5) senescence
%                                               6) ba (biomass accumulation)
%                                               7) em (emigration); NOTE: negative for immigration
%           SINGLE_ECO.EnergyBudget         (2D matrix: num_grps X num_grps); main ECOTRAN matrix defining flow from each producer to each consumer; (aka "production matrix", "energy budget"); (no ba & em rows)
%       metabolism, egg, feces, & senescence fates:
%           SINGLE_ECO.fate_metabolism      (2D matrix: num_nutrients X num_grps)
%           SINGLE_ECO.fate_eggs           	(2D matrix: num_eggs X num_grps)
%           SINGLE_ECO.fate_feces        	(2D matrix: num_ANYdetritus X num_grps)
%           SINGLE_ECO.fate_senescence     	(2D matrix: num_ANYdetritus X num_grps)
%           SINGLE_ECO.fate_predation      	(2D matrix: num_livingANDfleets X num_grps)
%       ecotrophic efficiency terms:
%           SINGLE_ECO.ee_eggs            	(vertical vector: num_grps X 1)
%           SINGLE_ECO.ee_predation      	(vertical vector: num_grps X 1)
%           SINGLE_ECO.ee_ba              	(vertical vector: num_grps X 1)
%           SINGLE_ECO.ee_em             	(vertical vector: num_grps X 1)
%           SINGLE_ECO.ee                  	(vertical vector: num_grps X 1)
%           SINGLE_ECO.TransferEfficiency 	(vertical vector: num_grps X 1)
%       physiology parameters:
%           SINGLE_ECO.ae                 	(vertical vector: num_grps X 1)
%           SINGLE_ECO.pq                 	(vertical vector: num_grps X 1)
%           SINGLE_ECO.pb                 	(vertical vector: num_grps X 1)
%       diet & consumption terms: (NOTE: I don't think it is necessary to pass along these terms)
%           SINGLE_ECO.DIET_NoCannibalism                   (2D matrix: num_grps X num_grps)
%           SINGLE_ECO.diet_cannibalism                    	(horizontal vector: 1 X num_grps)
%           SINGLE_ECO.consumption_domestic_cannibalism   	cannibalism rates; (t/km2/y); (vertical vector: num_grps X 1)
%           SINGLE_ECO.CONSUMPTION_NoCannibalism          	(2D matrix: num_grps X num_grps)
%           SINGLE_ECO.predation_total                   	total predation ON each producer group, p; includes cannibalism; (M2_p); (t/km2/yr); (vertical vector: num_grps X 1)
%           SINGLE_ECO.predation_total_NoCannibalism       	total predation ON each producer group, p; does NOT include cannibalism; (M2_p); (t/km2/yr); (vertical vector: num_grps X 1)
%           SINGLE_ECO.consumption_domestic_NOcannibalism   total consumption BY each consumer c; does NOT include cannibalism; (t/km2/yr); (horizontal vector: num_grps X 1)
%           SINGLE_ECO.consumption_domestic_cannibalism  	total cannibalism BY each consumer c; (t/km2/yr); (horizontal vector: num_grps X 1); NOTE transpose
%           SINGLE_ECO.ECOfunction_name                     name of this version of f_ECOfunction
%
% revision date: 9-3-2021
%       JR 12/29/2020 adjusting ConsumptionBudget & EnergyBudget to account for cannibalism
%       JR 2/6/2021 setting detritus & egg pb, qb, pq values = 1
%       9/3/2021 fixed cannibalism correcton and tested output
% FFF in future allow forced metabolism & forced egg production?


% *************************************************************************
% STEP 1: unpack terms for the current model ------------------------------
fname_ECOfunction                = mfilename; % save name of this m-file to keep in saved model results
disp(['   Running: ' fname_ECOfunction])

% step 1a: GroupType ------------------------------------------------------
GroupType                           = ModelDefinitions.GroupType;	% includes nutrients & fleets; (vertical vector: num_grps X 1)


% step 1b: group type definitions -----------------------------------------
GroupTypeDef_ANYNitroNutr           = ModelDefinitions.GroupTypeDef_ANYNitroNutr;
GroupTypeDef_NO3                    = ModelDefinitions.GroupTypeDef_NO3;
GroupTypeDef_plgcNH4                = ModelDefinitions.GroupTypeDef_plgcNH4;
GroupTypeDef_bnthNH4                = ModelDefinitions.GroupTypeDef_bnthNH4;
GroupTypeDef_ANYPrimaryProd         = ModelDefinitions.GroupTypeDef_ANYPrimaryProd;
GroupTypeDef_ANYConsumer            = ModelDefinitions.GroupTypeDef_ANYConsumer;
GroupTypeDef_fleet                  = ModelDefinitions.GroupTypeDef_fleet;
GroupTypeDef_eggs                   = ModelDefinitions.GroupTypeDef_eggs;
GroupTypeDef_terminalPlgcDetr       = ModelDefinitions.GroupTypeDef_terminalPlgcDetr;
GroupTypeDef_terminalBnthDetr       = ModelDefinitions.GroupTypeDef_terminalBnthDetr;
GroupTypeDef_ANYDetritus            = ModelDefinitions.GroupTypeDef_ANYDetritus;
% GroupTypeDef_micrograzers           = ModelDefinitions.GroupTypeDef_micrograzers;
GroupTypeDef_bacteria               = ModelDefinitions.GroupTypeDef_bacteria;


% step 1c: find group addresses -------------------------------------------
looky_NO3                       	= find(GroupType == GroupTypeDef_NO3);
looky_plgcNH4                       = find(GroupType == GroupTypeDef_plgcNH4);
looky_bnthNH4                       = find(GroupType == GroupTypeDef_bnthNH4);
looky_NH4                           = find(GroupType == GroupTypeDef_plgcNH4 | GroupType == GroupTypeDef_bnthNH4);
looky_nutrients                     = find(floor(GroupType) == GroupTypeDef_ANYNitroNutr);      % row addresses of E2E nutrients
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
looky_TermDetINdetritus          	= find(GroupType(looky_ANYdetritus) == GroupTypeDef_terminalPlgcDetr | GroupType(looky_ANYdetritus) == GroupTypeDef_terminalBnthDetr); % row addresses of terminal detritus groups WITHIN detritus sub-set of groups


% step 1d: numbers of group types -----------------------------------------
num_grps                            = ModelDefinitions.num_grps;              	% number of aggregated groups, including nutrients & fleets
num_NO3                             = ModelDefinitions.num_NO3;
num_plgcNH4                         = ModelDefinitions.num_plgcNH4;
num_bnthNH4                         = ModelDefinitions.num_bnthNH4;
num_NH4                             = ModelDefinitions.num_NH4;
num_nutrients                       = ModelDefinitions.num_nutrients;
num_eggs                            = ModelDefinitions.num_eggs;
num_ANYdetritus                     = ModelDefinitions.num_ANYdetritus;         % ANY detritus groups (w/o eggs)
num_terminalANYdetritus             = ModelDefinitions.num_terminalANYdetritus;
num_livingANDfleets                 = ModelDefinitions.num_livingANDfleets;


% step 1e: nutrient fates -------------------------------------------------
Oxidation_NH4                       = ModelDefinitions.Oxidation_NH4;    	% fraction of NH4 produced oxidized directly back to NO3 abiologically; this should take precedence over phytoplankton uptake of NH4 in code below; (vertical vector: num_NH4 X 1)
PhytoUptake_NO3                     = ModelDefinitions.PhytoUptake_NO3;  	% (vertical vector: num_PrimaryProducers X 1)
PhytoUptake_NH4                     = ModelDefinitions.PhytoUptake_NH4;     % (vertical vector: num_PrimaryProducers X 1)


% step 1f: implicit bacterial metabolism of detritus ----------------------
%          NOTE: applies to TERMINAL detritus pools
PelagicBacterialReduction           = ModelDefinitions.PelagicBacterialReduction; % fraction of terminal pelagic detritus oxidized to pelagic NH4
BenthicBacterialReduction           = ModelDefinitions.BenthicBacterialReduction; % fraction of terminal benthic detritus oxidized to benthic NH4


% step 1g: parameters -----------------------------------------------------
production                          = CurrentModel.production;	% (t/km2/y); (vertical vector: num_grps X 1)
pb                                  = CurrentModel.pb;       	% (1/y); (vertical vector: num_grps X 1) QQQ turn off??? recalculate later???
qb                                  = CurrentModel.qb;       	% (1/y); (vertical vector: num_grps X 1)
pq                                  = CurrentModel.pq;       	% (dimensionless); (vertical vector: num_grps X 1) QQQ recalculate later???
ae                                  = CurrentModel.ae;       	% (dimensionless); (vertical vector: num_grps X 1)
ba                                  = CurrentModel.ba;       	% (t/km2/y); (vertical vector: num_grps X 1); NOTE: still in absolute production rate units at this point
em                                  = CurrentModel.em;       	% (t/km2/y); (vertical vector: num_grps X 1); NOTE: still in absolute production rate units at this point
% EggProduction_forced                = CurrentModel.EggProduction_forced;	% (t/km2/y); (vertical vector: num_grps X 1); FFF for future


% step 1h: diet & consumption ---------------------------------------------
CONSUMPTION                         = CurrentModel.CONSUMPTION;         % (t/km2/y); (2D matrix: num_grps X num_grps); NOTE: there ARE non-zero entries for eggs, detritus, & fleets but all zero entries for nutrients & primary producers
% consumption_total                   = CurrentModel.consumption_total;	% (t/km2/y); (vertical vector: num_grps X 1); NOTE: includes import diet (which is not in CONSUMPTION); NOTE: includes primary producer uptake of nutrients;
DIET                                = CurrentModel.DIET;                % (2D matrix: num_grps X num_grps); NOTE: there ARE non-zero entries for eggs, detritus, & fleets but all zero entries for nutrients & primary producers


% step 1i: metabolism, egg, feces, & senescence fates ---------------------
fate_eggs                           = ModelDefinitions.fate_eggs';          % EwE egg fate; (proportions); (2D matrix: num_eggs X num_grps); NOTE transpose; NOTE: does not sum to 1 (eggs + detritus sum to 1)
fate_EwEdetritus                    = ModelDefinitions.fate_EwEdetritus';	% EwE detritus fate; (proportions); (2D matrix: num_ANYdetritus X num_grps); NOTE transpose; NOTE: does not include eggs & does not sum to 1 when group has eggs (eggs + detritus sum to 1)
fate_metabolism                     = ModelDefinitions.fate_metabolism;   	% ECOTRAN metabolism fate; (proportions); (2D matrix: num_NH4 X num_grps)
fate_feces                          = ModelDefinitions.fate_feces;       	% ECOTRAN feces detritus fate; (proportions); (2D matrix: num_terminalANYdetritus X num_grps)
fate_senescence                     = ModelDefinitions.fate_senescence;  	% ECOTRAN senescence detritus fate; (proportions); (2D matrix: num_terminalANYdetritus X num_grps)
% *************************************************************************





% *************************************************************************
% STEP 2: prep DIET & CONSUMPTION matrices---------------------------------

% step 2a: remove cannibalism from DIET matrix ----------------------------
[DIET_NoCannibalism, diet_cannibalism]          = f_RedistributeCannibalism_11202019(DIET); % remove cannibalism terms on diagonal of DIET matrix
%                                                   DIET_NoCannibalism (2D matrix: num_grps X num_grps)
%                                                   diet_cannibalism (horizontal vector: 1 X num_grps)


% step 2b: remove cannibalism from CONSUMPTION matrix ---------------------
consumption_domestic_cannibalism             	= diag(CONSUMPTION);   	% cannibalism rates; (t/km2/y); (vertical vector: num_grps X 1)
CONSUMPTION_NoCannibalism                       = CONSUMPTION;       	% initialize matrix
CONSUMPTION_NoCannibalism(1:(num_grps+1):end)	= zeros(1, num_grps);	% SSS turn off for testing with cannibalism on; assign zeros along diagonal; (2D matrix: num_grps X num_grps)


% step 2c: calculate total absolute predation on each producer, p; M2_p = sigma(b_pc * C_c) for all predators, c
predation_total                                 = sum(CONSUMPTION(:, looky_livingANDfleets), 2);                  % total predation ON each producer group, p; includes cannibalism; (M2_p); (t/km2/yr); (vertical vector: num_grps X 1)
predation_total_NoCannibalism                   = sum(CONSUMPTION_NoCannibalism(:, looky_livingANDfleets), 2);    % total predation ON each producer group, p; does NOT include cannibalism; (M2_p); (t/km2/yr); (vertical vector: num_grps X 1)


% step 2d: total consumption by each consumer c
%          NOTE: use of consumption_domestic (since we are considering only consumption within domain and not import diet)
consumption_domestic_NOcannibalism              = sum(CONSUMPTION_NoCannibalism, 1); % total consumption BY each consumer c; does NOT include cannibalism; (t/km2/yr); (horizontal vector: num_grps X 1)
consumption_domestic_cannibalism                = consumption_domestic_cannibalism'; % total cannibalism BY each consumer c; (t/km2/yr); (horizontal vector: num_grps X 1); NOTE transpose


% SSS turn cannibalism ON/OFF
%     -- use when cannibalism is present
%          disp('   -->WARNING: cannibalism is left ACTIVE for testing')
%          DIET_NoCannibalism     = DIET;
%          diet_cannibalism       = NaN;
% 
%     -- use when cannibalism is removed
         disp(['      -->WARNING in ' fname_ECOfunction ': cannibalism is REMOVED. Feces & metabolism are increased by the cannibalism fraction of diet and this results in reduced TE.'])
% *************************************************************************





% *************************************************************************
% STEP 3: metabolism & nutrients -- fix for BioenergeticBudget, ProductionBudget, ConsumptionBudget and adjust fates for PredationBudget

% step 3a: reduce pq of NH4 pools to account for oxidation to NO3 ---------
pq(looky_NH4)                   = pq(looky_NH4) - Oxidation_NH4;	% reduce NH4 pq by fraction oxidized to NO3; NOTE: these terms are defined separately for pelagic & benthic NH4 pools


% step 3b: set detritus & egg pb, qb, & pq values to 1 --------------------
% pb([looky_eggs; looky_ANYdetritus]) = 1; % pb is redefined below QQQ leave off ???
qb([looky_eggs; looky_ANYdetritus]) = 1;
pq([looky_eggs; looky_ANYdetritus]) = 1;


% step 3c: apply any implicit bacterial metabolism of detritus ------------
pq(looky_terminalPLGCdetritus)	= (1 - PelagicBacterialReduction);
pq(looky_terminalBNTHdetritus)	= (1 - BenthicBacterialReduction);
% warning if there are explicitly defined bacterial groups AND implicitly defined bacterial reduction of detritus
if ~isempty(looky_bacteria) && sum([PelagicBacterialReduction; BenthicBacterialReduction]) > 0
    disp(['      -->WARNING in ' fname_ECOfunction ': implicit bacterial metabolism AND explicit bacteria groups are BOTH defined'])
end
% *************************************************************************





% *************************************************************************
% STEP 4: calculate Ecotrophic Efficiencies (distribution of production)
%         NOTE: additional modifications made to ee terms due to EwE eggs-as-detritus definitions are made in STEP 5
%         NOTE: Transfer Efficiencies (distribution of consumption, fractions of Q_p) are calculated in STEP 5

% step 4a: calculate initial ee terms -------------------------------------
randomEwE.production            = production;                   % (t/km2/y); (vertical vector: num_grps X 1)
randomEwE.ba                    = ba;                           % (t/km2/y); (vertical vector: num_grps X 1); NOTE: still in absolute production rate units at this point
randomEwE.em                    = em;                           % (t/km2/y); (vertical vector: num_grps X 1); NOTE: still in absolute production rate units at this point
% randomEwE.EggProduction     	= EggProduction_forced;         % (t/km2/y); (vertical vector: num_grps X 1); FFF for future
randomEwE.EggProduction     	= zeros(num_grps, 1);           % (t/km2/y); (vertical vector: num_grps X 1); FFF EggProduction_forced is for future, set to all zeros for now
randomEwE.looky_eggsANDdetritus = looky_eggsANDdetritus;
randomEwE.looky_nutrients       = looky_nutrients;

% ee terms with cannibalism INCLUDED
%	NOTE: use senescence term calculated WITH cannibalism included
randomEwE.CONSUMPTION           = CONSUMPTION;                  % (t/km2/y); (2D matrix: num_grps X num_grps); NOTE: use of CONSUMPTION (includes cannibalism)
EcotrophicEfficiency_cannibal	= f_calcEE_12292020(randomEwE);
ee_predation_cannibal       	= EcotrophicEfficiency_cannibal.ee_predation; % predation as proportion of production budget; (vertical vector: num_grps X 1)
ee_cannibal                 	= EcotrophicEfficiency_cannibal.ee;           % fraction of all production used within food web and NOT senescence; (vertical vector: num_grps X 1); NOTE: in EwE this traditionally did not include em but it does in ECOTRAN
ee_senescence_cannibal          = (1 - ee_cannibal);                          % senescence; (vertical vector: num_grps X 1)
                                                                              %   NOTE: senescence is not really an ee term but I add this as a specific term for convenience when doing egg & reproduction adjustements later
                                                                              %   NOTE: use of ee that includes cannibalism in ee_predation to calculate ee_senescence_cannibal
% ee terms with cannibalism REMOVED
%   NOTE: ee_predation & ee change if cannibalism is removed but other ee terms do not (ee_eggs, ee_ba, ee_em)
randomEwE.CONSUMPTION           = CONSUMPTION_NoCannibalism;    % (t/km2/y); (2D matrix: num_grps X num_grps); NOTE: use of CONSUMPTION_NoCannibalism
EcotrophicEfficiency            = f_calcEE_12292020(randomEwE);
ee_eggs                         = EcotrophicEfficiency.ee_eggs;         % eggs (or other reproduction) as proportion of production budget; (vertical vector: num_grps X 1); NOTE: at this point ee_eggs is defined by forced egg production, eggs as EwE detritus fates are accounted for in ee at STEP 5;
ee_predation                    = EcotrophicEfficiency.ee_predation;    % predation as proportion of production budget; (vertical vector: num_grps X 1)
ee_ba                           = EcotrophicEfficiency.ee_ba;           % biomass accumulation (local production) as proportion of production budget; (vertical vector: num_grps X 1)
ee_em                           = EcotrophicEfficiency.ee_em;           % emigration as proportion of production budget; (vertical vector: num_grps X 1)

% step 4b: adjust nutrient ee terms ---------------------------------------
ee_predation(looky_NO3)             = round2(sum(PhytoUptake_NO3), 7);  % NOTE use of round2 here because of rounding error in the 10th decimal point
ee_predation(looky_NH4)             = round2([sum(PhytoUptake_NH4); sum(PhytoUptake_NH4)], 7); % FFF can treat pelagic & benthic pools separately; NOTE: use of round2 here because of rounding error in the 10th decimal point

ee_predation_cannibal(looky_NO3)	= round2(sum(PhytoUptake_NO3), 7);  % NOTE use of round2 here because of rounding error in the 10th decimal point
ee_predation_cannibal(looky_NH4)	= round2([sum(PhytoUptake_NH4); sum(PhytoUptake_NH4)], 7); % FFF can treat pelagic & benthic pools separately; NOTE: use of round2 here because of rounding error in the 10th decimal point

ee_senescence_cannibal(looky_NO3)	= 1 - (ee_predation_cannibal(looky_NO3) + ee_ba(looky_NO3) + ee_em(looky_NO3));
ee_senescence_cannibal(looky_NH4)	= 1 - (ee_predation_cannibal(looky_NH4) + ee_ba(looky_NH4) + ee_em(looky_NH4));
% *************************************************************************





% *************************************************************************
% STEP 5: egg correction (for EwE eggs-as-detritus assignments) -----------
%         1) remove EwEegg_fraction from feces
%         2) add "feces" eggs to production by increasing ae & pq
%         3) remove EwEegg_fraction from senescence production
%         4) put "senescence" eggs onto egg production (ee_eggs)
%         5) redistribute the new "feces" egg production to ee_eggs and redistribute all other production budget terms (ee_terms) downward

% step 5a: allow for use of only EwE eggs-as-detritus at this time --------
%          FFF: will allow forced egg production in future. This will
%              require adjustment to ProductionBudget, ConsumptionBudget, & EnergyBudget
% warn if eggs are defined as forced egg production
sum_eggs_forced             = sum(ee_eggs); % at this point ee_eggs is defined by forced egg production
if sum_eggs_forced > 0
    disp(['      -->WARNING in ' fname_ECOfunction ': forced egg production is not yet supported. Use EwE eggs-as-detritus OR enter egg production directly into modified ProductionBudget, ConsumptionBudget, & EnergyBudget'])
end
% -------------------------------------------------------------------------


% step 5b: egg fraction of feces & senescence and corrections to production & physiology terms
EwEegg_fraction             = sum(fate_eggs, 1);          	% fraction of senescence & feces reassigned to eggs; (horizontal vector: 1 X num_grps)
feces_EggCorrection         = (1 - ae) .* EwEegg_fraction';	% fraction of consumption going to eggs, taken from feces; NOTE: this will increase production; (vertical vector: num_grps X 1); NOTE transpose
ae                          = ae + feces_EggCorrection;   	% remove eggs from feces by INCREASING ae by egg fraction of feces; (vertical vector: num_grps X 1)
pq                          = pq + feces_EggCorrection;   	% add eggs (& reproduction) to production by INCREASING pq by egg fraction of feces; (vertical vector: num_grps X 1)
pb                          = pq .* qb;                     % recalculate pb to include egg (reproduction) production QQQ pb recalculated again later??
% -------------------------------------------------------------------------


% step 5c: remove egg fraction from senescence, add to egg production -----
ee_eggs                 	= ee_senescence_cannibal .* EwEegg_fraction';       % (vertical vector: num_grps X 1); NOTE transpose; NOTE: this line must come before the ee_senescence_cannibal correction line
ee_senescence_cannibal      = ee_senescence_cannibal .* (1 - EwEegg_fraction)';	% (vertical vector: num_grps X 1); NOTE transpose
% -------------------------------------------------------------------------


% step 5d: redistribute "surplus" production from "feces" eggs ------------
%          redistribute this additional production onto eggs only
%          NOTE: INCREASE ee_eggs (all additional production, formerly going to feces, goes onto egg production)
%          NOTE: DECREASE ee_predation, ee_senescence_cannibal, ee_ba, & ee_em by
%                the fraction of production going to "feces" eggs (distribute
%                this reduction based upon the fraction of production going to
%                each term (ee_predation, ee_senescence_cannibal, ee_ba, or ee_em). The
%                remaining fraction of "feces" egg production (1-ee_eggs) goes to eggs.
%                This algorithm was tested (4/1/2019) and correctly puts all "surplus"
%                "eggs-as-feces" surplus production (created by the revised
%                pq term) onto eggs.
%          NOTE: these are changes to production budget terms, reduction to individual ee_terms does not mean a reduction in their production rates
%          NOTE: term (feces_EggCorrection./pq) = fraction of production that is derived from "feces" eggs
ee_eggs                 	= ee_eggs                + ((feces_EggCorrection./pq) .* (1-ee_eggs)); 	% (vertical vector: num_grps X 1)
ee_predation            	= ee_predation           - ((feces_EggCorrection./pq) .* ee_predation); % (vertical vector: num_grps X 1); cannibalism REMOVED
ee_predation_cannibal       = ee_predation_cannibal	 - ((feces_EggCorrection./pq) .* ee_predation_cannibal); 	% (vertical vector: num_grps X 1); cannibalism INCLUDED
ee_senescence_cannibal      = ee_senescence_cannibal - ((feces_EggCorrection./pq) .* ee_senescence_cannibal);	% (vertical vector: num_grps X 1)
ee_ba                     	= ee_ba                  - ((feces_EggCorrection./pq) .* ee_ba);       	% (vertical vector: num_grps X 1)
% ee_em                      	= ee_em                  - ((feces_EggCorrection./pq) .* ee_em);        % NOTE: Keep this commented off, make NO egg correction to decrease production going to ee_em; (vertical vector: num_grps X 1)

%          NOTE: do not decrease production going to ee_em, redistribute this egg adjustment of ee_em to ee_predation, ee_senescence_cannibal, and ee_ba
em_EggCorrection            = (feces_EggCorrection./pq) .* ee_em;

temp_sum                    = ee_predation          + ee_senescence_cannibal + ee_ba;
temp_sum_cannibal        	= ee_predation_cannibal + ee_senescence_cannibal + ee_ba;

ee_predation                = ee_predation           - (em_EggCorrection .* (ee_predation           ./ temp_sum));
ee_predation_cannibal     	= ee_predation_cannibal	 - (em_EggCorrection .* (ee_predation_cannibal  ./ temp_sum_cannibal));
ee_senescence_cannibal      = ee_senescence_cannibal - (em_EggCorrection .* (ee_senescence_cannibal ./ temp_sum_cannibal));
ee_ba                       = ee_ba                  - (em_EggCorrection .* (ee_ba                  ./ temp_sum_cannibal)); % QQQ change this to temp_sum??

ee_predation(isnan(ee_predation))                     = 0; % set div/0 errors (e.g., fleets) to 0
ee_predation_cannibal(isnan(ee_predation_cannibal))	  = 0; % set div/0 errors (e.g., fleets) to 0
ee_senescence_cannibal(isnan(ee_senescence_cannibal)) = 0; % set div/0 errors (e.g., fleets) to 0
ee_ba(isnan(ee_ba))                                   = 0; % set div/0 errors (e.g., fleets) to 0
% *************************************************************************





% *************************************************************************
% STEP 6: prep ConsumptionBudget terms-------------------------------------
%         NOTE: "cb" prefix stands for ConsumptionBudget 
%               as opposed to "ee" which stands for EcotrophicEfficiency,
%               which are ProductionBudget terms

% step 6a: cb terms (without cannibalism corrections) ---------------------
temp_ee_em                  = ee_em;	% zero-out fleet em within temporary ee_em variable
temp_ee_em(looky_fleets)	= 0;        % Used for fleets in ConsumptionBudget calculations

cb_production_cannibal    	= pq        .* (1 - temp_ee_em); % fraction of consumption going to production; (vertical vector: num_grps X 1)
cb_feces                    = (1 - ae)  .* (1 - temp_ee_em); % fraction of consumption going to feces; (vertical vector: num_grps X 1); NOTE: correction for non-fleet import & export terms
cb_metabolism               = (ae - pq) .* (1 - temp_ee_em); % fraction of consumption going to metabolism (NH4 production); Metabolism = 1 - PQ - feces = (1 - (EwE_PQ + (1-EwE_AE))); (vertical vector: num_grps X 1); NOTE: correction for non-fleet import & export terms

cb_eggs                     = pq .* ee_eggs;               % fraction of consumption going to eggs (& reproduction); (vertical vector: num_grps X 1)
cb_predation                = pq .* ee_predation;          % fraction of consumption going to predation; (vertical vector: num_grps X 1)
cb_predation_cannibal       = pq .* ee_predation_cannibal; % fraction of consumption going to predation; (vertical vector: num_grps X 1)

cb_senescence               = pq .* ee_senescence_cannibal;	 % fraction of consumption going to senescence; (vertical vector: num_grps X 1)
cb_ba                       = pq .* ee_ba;           % fraction of consumption going to ba; (vertical vector: num_grps X 1)
cb_em                       = ee_em;                 % fraction of consumption going to em; (vertical vector: num_grps X 1); NOTE: cb_em = ee_em, the logic being if 10% of diet is import diet, then both 10% of consumption and 10% of production is from import (or fated to export for a positive em), fleets are excluded
cb_em(looky_fleets)         = pq(looky_fleets) .* ee_em(looky_fleets);	% for fleets, landings are the em term, landings are the fraction of catch "emigrating" from the system; (vertical vector: num_grps X 1)
% -------------------------------------------------------------------------


% step 6b: cannibalism corrections for ConsumptionBudget ------------------
cannibal_correction_cb         = cb_predation_cannibal - cb_predation; % portion of cb_predation due to cannibalism

% Increase cb_feces & cb_metabolism to account for cannibalism
feces_fraction              = cb_feces      ./ (cb_feces + cb_metabolism);
metabolism_fraction         = cb_metabolism ./ (cb_feces + cb_metabolism);
feces_fraction(isnan(feces_fraction))           = 0; % correct div/0 error (nutrients & detritus)
metabolism_fraction(isnan(metabolism_fraction))	= 0; % correct div/0 error (nutrients & detritus)
cb_feces                    = cb_feces      + (feces_fraction      .* cannibal_correction_cb);
cb_metabolism               = cb_metabolism + (metabolism_fraction .* cannibal_correction_cb);

% Reduce cb_production to account for cannibalism
cb_production               = cb_production_cannibal - cannibal_correction_cb; % reduced to correct for cannibalism

% Reduce pq & pb to account for cannibalism
pq                          = cb_production; % pq is now modified to account for cannibalism
pb                          = pq .* qb;      % pb is now modified to account for cannibalism
% -------------------------------------------------------------------------


% step 6c: cannibalism corrections for ProductionBudget (ee) terms --------
cannibal_correction_ee      = cb_production_cannibal ./ cb_production;
cannibal_correction_ee(isnan(cannibal_correction_ee)) = 0; % correct div/0 error (probably not necessary) 
ee_eggs                     = ee_eggs                .* cannibal_correction_ee; % eggs & reproduction; (vertical vector: num_grps X 1)
ee_predation                = ee_predation           .* cannibal_correction_ee;	% predation (M2_p / P_p); (vertical vector: num_grps X 1)
ee_senescence               = ee_senescence_cannibal .* cannibal_correction_ee; % senescence (M0_p / P_p); fraction of production going to other mortality (non-consumed production); (vertical vector: num_grps X 1)
ee_ba                       = ee_ba                  .* cannibal_correction_ee; % ba (ba_p / P_p); (vertical vector: num_grps X 1)
ee_em                       = ee_em                  .* cannibal_correction_ee; % em (em_p / P_p); includes physical transport out of model domain; immigration is negative production, emigration is positive production; (vertical vector: num_grps X 1)
% -------------------------------------------------------------------------


% step 6d: finalize ee term -----------------------------------------------
%          NOTE: addition of ee_eggs term (caused by taking eggs from feces & senescence). 
%                This will cause E2E ee to be higher than EwE ee
%          NOTE: use ee_predation_cannibal (predation that includes the cannibalism component) 
%                because senescence is NOT changed by the cannibalism
%                correction. Cannibalism is a predation term. Removing
%                cannibalism from the EnergyBudget does NOT cause a change 
%                in senescence.
%                Removing cannibalism from the EnergyBudget DOES require 
%                accounting for reduced TE due to the cannibalism contribution 
%                to feces & metabolism. Feces & metabolism are increased to 
%                account for cannibalism.
ee                          =  ee_eggs + ee_predation + ee_ba + ee_em; % (vertical vector: num_grps X 1)
% -------------------------------------------------------------------------


% step 6e: calculate TransferEfficiency -----------------------------------
%          NOTE: (M2_p + EM_p + BA_p) / Q_p = fraction of consumption by each group c passed upwards in food web 
%                                             via predation, ba, & em
TransferEfficiency      	= ee .* pq; % (vertical vector: num_grps X 1)
% *************************************************************************





% *************************************************************************
% STEP 7: prep nutrient, egg, & detritus fates-----------------------------

% step 7a: combine EwE & E2E detritus fates into unified fate rules -------
%          NOTE: EwE fates defined for terminal detritus get redistributed into
%                E2E feces & senescence detritus fates
%          NOTE: EwE fates for non-terminal detritus pools are not altered 
%                (these EwE fates are the same whether feces or senescence)
%          NOTE: fate assignments are correct as long as the order of
%                detritus pools has not been changed
fate_EwEdetritus_unified                            = fate_EwEdetritus ./ repmat(sum(fate_EwEdetritus, 1), [num_ANYdetritus, 1]);	% renormalize EwE detritus fates; (2D matrix: num_ANYdetritus X num_grps)
fate_EwEdetritus_unified(isnan(fate_EwEdetritus_unified)) = 0;                                      % fix div/0 errors (nutrients & terminal benthic detritus)
fate_feces_unified                                  = fate_EwEdetritus_unified;                     % initiaize; (2D matrix: num_ANYdetritus X num_grps)
fate_senescence_unified                             = fate_EwEdetritus_unified;                     % initiaize; (2D matrix: num_ANYdetritus X num_grps)
total_EwETerminalDetritus                           = repmat(sum(fate_EwEdetritus_unified(looky_TermDetINdetritus, :), 1), [num_terminalANYdetritus, 1]);	% total terminal detritus in EwE detritus fates; (2D matrix: num_terminalANYdetritus X num_grps)
fate_feces_unified(looky_TermDetINdetritus, :)    	= fate_feces      .* total_EwETerminalDetritus;	% ECOTRAN feces detritus fate; (proportions); (2D matrix: num_ANYdetritus X num_grps)
fate_senescence_unified(looky_TermDetINdetritus, :)	= fate_senescence .* total_EwETerminalDetritus;	% ECOTRAN senescence detritus fate; (proportions); (2D matrix: num_ANYdetritus X num_grps)


% step 7b: normalize egg fates to sum to 1 --------------------------------
fate_eggs                                           = fate_eggs ./ repmat(sum(fate_eggs, 1), [num_eggs, 1]);	% egg fate; (proportions); (2D matrix: num_eggs X num_grps);
fate_eggs(isnan(fate_eggs))                         = 0;                                                        % fix div/0 error (groups without eggs)


% step 7c: add abiotic oxidation of NH4 to NO3 as a metabolism fate -------
fate_metabolism((num_NO3+1):(num_NH4+num_NO3), :)	= fate_metabolism;
fate_metabolism(1:num_NO3, :)                       = 0;
fate_metabolism(1:num_NO3, looky_NH4)               = 1;
% *************************************************************************





% *************************************************************************
% STEP 8: construct energy flow budgets------------------------------------
%         NOTE: There are four energy flow budgets (unitless fractions, 0 to 1):
%               A) BioenergeticBudget	(fraction of consumption by each producer group, p, going to production, feces = non-assimilated consumption, or metabolism)
%               B) ProductionBudget     (fraction of total production of each producer, p, going to predation, senescence, emigration (or immigration, physical transport), or biomass accumulation)
%               C) ConsumptionBudget    (detailed fate of all consumption)
%               D) EnergyBudget         (fraction of total consumption by each producer p going to each consumer c, metabolism, predation, eggs, detritus, or fleets)

% step 8a: initialize variables -------------------------------------------
BioenergeticBudget          = zeros(3, num_grps);
%                               1) feces
%                               2) metabolism
%                               3) production
ProductionBudget            = zeros(5, num_grps);
%                               1) eggs (reproduction)
%                               2) predation
%                               3) senescence
%                               4) ba (biomass accumulation)
%                               5) em (emigration); NOTE: negative for immigration
ConsumptionBudget           = zeros(7, num_grps);
%                               1) feces
%                               2) metabolism
%                               3) eggs (reproduction)
%                               4) predation
%                               5) senescence
%                               6) ba (biomass accumulation)
%                               7) em (emigration); NOTE: negative for immigration


% step 8b: calculate BioenergeticBudget for each producer, p --------------
%          fraction of total consumption (Q) going to metabolism, feces (non-assimilated consumption), or production
% QQQ NOTE: will not sum to 1 if there is import diet or (non-fleet) export; maybe I still wnat this to sum to 1 and handle the import/export terms in the ConsumptionBudget?
BioenergeticBudget(1, :)	= cb_feces';        % fraction of consumption going to feces; (horizontal vector: 1 X num_grps); NOTE transpose
BioenergeticBudget(2, :)	= cb_metabolism';	% fraction of consumption going to metabolism (NH4 production); Metabolism = 1 - PQ - feces = (1 - (EwE_PQ + (1-EwE_AE))); (horizontal vector: 1 X num_grps); NOTE transpose
BioenergeticBudget(3, :)    = cb_production';	% fraction of consumption going to production; (horizontal vector: 1 X num_grps); NOTE transpose


% step 8c: calculate ProductionBudget for each group, p -------------------
%          fraction of production going to eggs & reproduction, predation, senescence, biomass accumulation, or emigration (incl. immigration, physical transport)
ProductionBudget(1, :)      = ee_eggs';         % eggs & reproduction; NOTE: transpose
ProductionBudget(2, :)      = ee_predation';	% predation (M2_p / P_p); NOTE: transpose
ProductionBudget(3, :)      = ee_senescence';   % senescence (M0_p / P_p); fraction of production going to other mortality (non-consumed production); NOTE: cannibalism REMOVED; prior to 12/29/2020 elimination of cannibalism predation resulted in increased senescence; NOTE: transpose
ProductionBudget(4, :)      = ee_ba';           % ba (ba_p / P_p); NOTE: transpose
ProductionBudget(5, :)      = ee_em';           % em (em_p / P_p); includes physical transport out of model domain; immigration is negative production, emigration is positive production; NOTE: transpose


% step 8d: calculate ConsumptionBudget for each group, p ------------------
%          fraction of consumption going to feces, metabolism, eggs, predation, senescence, biomass accumulation, or emigration (incl. immigration, physical transport)
ConsumptionBudget(1, :)     = cb_feces';        % feces; (horizontal vector: 1 X num_grps); NOTE transpose
ConsumptionBudget(2, :)     = cb_metabolism';	% metabolism; (horizontal vector: 1 X num_grps); NOTE transpose
ConsumptionBudget(3, :)     = cb_eggs';         % eggs; (horizontal vector: 1 X num_grps); NOTE transpose
ConsumptionBudget(4, :)     = cb_predation';	% predation; (horizontal vector: 1 X num_grps); NOTE transpose
ConsumptionBudget(5, :)     = cb_senescence';	% senescence; (horizontal vector: 1 X num_grps); NOTE transpose
ConsumptionBudget(6, :)     = cb_ba';           % ba; (horizontal vector: 1 X num_grps); NOTE transpose
ConsumptionBudget(7, :)     = cb_em';           % em; (horizontal vector: 1 X num_grps); NOTE transpose


% step 8e: calculate PredationBudget --------------------------------------
%          NOTE: import consumption & diet purposely not included in these calculations
%          NOTE: there ARE non-zero entries for eggs, detritus, & fleets but 
%                all zero entries for nutrients & primary producers in
%                consumption & diet variables
%          NOTE: We do NOT exclude any flow to nutrients, eggs, or detritus
%                when calculating the PredationBudget. BUT, we DO exclude 
%                any flow nutrients, eggs, or detritus in other functions
%                when we are specifically calculating predation terms in
%                ProductionBudget & ConsumptionBudget
%          NOTE: PredationBudget includes flow to eggs & detritus pools, but
%                this is later accounted for in the calculation of the
%                EnergyBudget by scaling consumer & fleet rows by the predation
%                terms of the ConsumptionBudget
PredationBudget             = f_calcPredationBudget_12102019(consumption_domestic_NOcannibalism, DIET_NoCannibalism); % (2D matrix: num_grps X num_grps; consumers X producers)
% *************************************************************************





% *************************************************************************
% STEP 9: complete the EnergyBudget----------------------------------------
%         NOTE: this is ECOTRAN matrix A_cp (here called BigMatrix)

% step 9a: initialize -----------------------------------------------------
BigMatrix                                       = PredationBudget; % initialize; (2D matrix: num_grps X num_grps)


% step 9b: paste nutrient flow to primary producers into BigMatrix --------
BigMatrix(looky_ANYPrimaryProd, looky_NO3)     	= repmat(PhytoUptake_NO3, [1, num_NO3]);            % distribute NO3 among primary producers
BigMatrix(looky_ANYPrimaryProd, looky_plgcNH4)	= repmat(PhytoUptake_NH4, [1, num_plgcNH4]);        % distribute pelagic NH4 among primary producers; rescaled to allow NH4 oxidation to NO3
BigMatrix(looky_ANYPrimaryProd, looky_bnthNH4)	= repmat(PhytoUptake_NH4, [1, num_bnthNH4]);        % distribute benthic NH4 among primary producers; rescaled to allow NH4 oxidation to NO3


% step 9c: scale predator rows in BigMatrix by predation term in ConsumptionBudget
%          NOTE: special case for terminal benthic detritus group(s) where predator rows are 
%                scaled by predation + senescence so that their column sum(s) add to 1
sum_predation                                   = sum(BigMatrix(looky_livingANDfleets, :), 1);      % (horizontal vector: 1 X num_grps)
temp_livingANDfleets                            = BigMatrix(looky_livingANDfleets, :) ./ repmat(sum_predation, [num_livingANDfleets, 1]); % temporarily renormalize predation rows to 1; (2D matrix: num_livingANDfleets X num_grps)
temp_livingANDfleets(isnan(temp_livingANDfleets)) = 0;                  % fix div/0 errors
temp_cb_predation                               = cb_predation; 
temp_cb_predation(looky_terminalBNTHdetritus)   = cb_predation(looky_terminalBNTHdetritus) + cb_senescence(looky_terminalBNTHdetritus); % special case for terminal benthic detritus group(s)
temp_livingANDfleets                            = temp_livingANDfleets .* repmat(temp_cb_predation', [num_livingANDfleets, 1]); % rescale by predation term of ConsumptionBudget; (2D matrix: num_livingANDfleets X num_grps); NOTE transpose
BigMatrix(looky_livingANDfleets, :)             = temp_livingANDfleets;                             % paste predator rows back into BigMatrix; (2D matrix: num_grps X num_grps)


% step 9d: paste eggs (& reproduction) term from ConsumptionBudget into BigMatrix
temp_eggs                                       = repmat(ConsumptionBudget(3, :), [num_eggs, 1]);	% (2D matrix: num_eggs X num_grps)
temp_eggs                                       = temp_eggs .* fate_eggs;                           % redistribute among egg rows as defined by fate_eggs; (2D matrix: num_eggs X num_grps)
BigMatrix(looky_eggs, :)                        = temp_eggs;                                        % paste egg (& reproduction) rows back into BigMatrix; (2D matrix: num_grps X num_grps)


% step 9e: paste feces & senescence terms from ConsumptionBudget into appropriate BigMatrix detritus rows
temp_feces                                      = repmat(ConsumptionBudget(1, :), [num_ANYdetritus, 1]);	% (2D matrix: num_ANYdetritus X num_grps)
temp_feces                                      = temp_feces .* fate_feces_unified;                         % redistribute feces among detritus rows as defined by fate_feces_unified; (2D matrix: num_ANYdetritus X num_grps)
temp_senescence                                 = repmat(ConsumptionBudget(5, :), [num_ANYdetritus, 1]);    % (2D matrix: num_ANYdetritus X num_grps)
temp_senescence                                 = temp_senescence .* fate_senescence_unified;               % redistribute senescence among detritus rows as defined by fate_senescence_unified; (2D matrix: num_ANYdetritus X num_grps)
BigMatrix(looky_ANYdetritus, :)                 = temp_feces + temp_senescence;                             % paste detritus rows back into BigMatrix; (2D matrix: num_grps X num_grps)


% step 9f: paste metabolism term from ConsumptionBudget into BigMatrix ----
%          NOTE: includes abiotic oxidation of NH4 to NO3
temp_metabolism                                 = repmat(ConsumptionBudget(2, :), [num_nutrients, 1]);      % (2D matrix: num_nutrients X num_grps)
temp_metabolism                                 = temp_metabolism .* fate_metabolism;                       % redistribute metabolism among nutrient rows as defined by fate_metabolism; (2D matrix: num_ANYdetritus X num_grps)
BigMatrix(looky_nutrients, :)                   = temp_metabolism;                                          % paste nutrient rows back into BigMatrix; (2D matrix: num_grps X num_grps)


% step 9g: error-checking of budgets --------------------------------------
BioenergeticBudget_checksum                     = round2(sum(BioenergeticBudget), 7); % (horizontal vector: 1 X num_grps)
looky_error                                     = find(BioenergeticBudget_checksum ~= 1);
if ~isempty(looky_error)
    disp(['      -->WARNING in ' fname_ECOfunction ': BioenergeticBudget does not sum to 1 (confirm non-zero em) at grps: ' num2str(looky_error)])
end

ProductionBudget_checksum                       = round2(sum(ProductionBudget), 7); % (horizontal vector: 1 X num_grps)
looky_error                                     = find(ProductionBudget_checksum ~= 1);
if ~isempty(looky_error)
    disp(['      -->WARNING in ' fname_ECOfunction ': ProductionBudget does not sum to 1 at grps: ' num2str(looky_error)])
end

EnergyBudget_checksum                        	= round2((sum(BigMatrix) + (cb_ba' + cb_em')), 7); % (horizontal vector: 1 X num_grps); NOTE: ba & em terms are not included in EnergyBudget
looky_error                                     = find(EnergyBudget_checksum ~= 1);
if ~isempty(looky_error)
    disp(['      -->WARNING in ' fname_ECOfunction ': EnergyBudget does not sum to 1 at grp(s): ' num2str(looky_error)])
end
% *************************************************************************





% *************************************************************************
% STEP 10: calculate fate_predation----------------------------------------
%          this is the distribution of the total predation of each producer among all its consumers
sum_predation                         	= sum(BigMatrix(looky_livingANDfleets, :), 1); % (horizontal vector: 1 X num_grps)
fate_predation                          = BigMatrix(looky_livingANDfleets, :) ./ repmat(sum_predation, [num_livingANDfleets, 1]); % (2D matrix: num_livingANDfleets X num_grps)
fate_predation(isnan(fate_predation))	= 0; % correct div/0 errors
% *************************************************************************





% *************************************************************************
% STEP 11: package results for export--------------------------------------
%          NOTE: only returning parameters that were created or altered within this function

% step 11a: four budget matrices ------------------------------------------
SINGLE_ECO.BioenergeticBudget       = BioenergeticBudget;   % (2D matrix: 3 X num_grps)
%                                           1) feces
%                                           2) metabolism
%                                           3) production
SINGLE_ECO.ProductionBudget         = ProductionBudget; % (2D matrix: 5 X num_grps)
%                                           1) eggs (reproduction)
%                                           2) predation
%                                           3) senescence
%                                           4) ba (biomass accumulation)
%                                           5) em (emigration); NOTE: negative for immigration
SINGLE_ECO.ConsumptionBudget        = ConsumptionBudget; % (2D matrix: 7 X num_grps)
%                                           1) feces
%                                           2) metabolism
%                                           3) eggs (reproduction)
%                                           4) predation
%                                           5) senescence
%                                           6) ba (biomass accumulation)
%                                           7) em (emigration); NOTE: negative for immigration
SINGLE_ECO.EnergyBudget             = BigMatrix; % (2D matrix: num_grps X num_grps)


% step 11b: metabolism, egg, feces, & senescence fates --------------------
SINGLE_ECO.fate_metabolism          = fate_metabolism;          % (2D matrix: num_nutrients X num_grps)
SINGLE_ECO.fate_eggs                = fate_eggs;                % (2D matrix: num_eggs X num_grps)
SINGLE_ECO.fate_feces               = fate_feces_unified;       % (2D matrix: num_ANYdetritus X num_grps)
SINGLE_ECO.fate_senescence          = fate_senescence_unified;	% (2D matrix: num_ANYdetritus X num_grps)
SINGLE_ECO.fate_predation           = fate_predation;           % (2D matrix: num_livingANDfleets X num_grps)


% step 11c: ecotrophic efficiency terms -----------------------------------
SINGLE_ECO.ee_eggs                  = ee_eggs;              % (vertical vector: num_grps X 1)
SINGLE_ECO.ee_predation             = ee_predation;         % (vertical vector: num_grps X 1)
SINGLE_ECO.ee_ba                    = ee_ba;                % (vertical vector: num_grps X 1)
SINGLE_ECO.ee_em                    = ee_em;                % (vertical vector: num_grps X 1)
SINGLE_ECO.ee                       = ee;                   % (vertical vector: num_grps X 1)
SINGLE_ECO.TransferEfficiency       = TransferEfficiency;	% (vertical vector: num_grps X 1)


% step 11d: physiology parameters -----------------------------------------
SINGLE_ECO.ae                       = ae; % (vertical vector: num_grps X 1)
SINGLE_ECO.pb                       = pb; % (vertical vector: num_grps X 1)
SINGLE_ECO.qb                       = qb; % (vertical vector: num_grps X 1)
SINGLE_ECO.pq                       = pq; % (vertical vector: num_grps X 1)


% step 11e: diet & consumption terms --------------------------------------
%          (NOTE: I don't think it is necessary to pass along these terms)
SINGLE_ECO.DIET_NoCannibalism                   = DIET_NoCannibalism;                 % (2D matrix: num_grps X num_grps)
SINGLE_ECO.diet_cannibalism                     = diet_cannibalism;                   % (horizontal vector: 1 X num_grps)
SINGLE_ECO.consumption_domestic_cannibalism     = consumption_domestic_cannibalism;   % cannibalism rates; (t/km2/y); (vertical vector: num_grps X 1)
SINGLE_ECO.CONSUMPTION_NoCannibalism            = CONSUMPTION_NoCannibalism;          % (2D matrix: num_grps X num_grps)
SINGLE_ECO.predation_total                      = predation_total;                    % total predation ON each producer group, p; includes cannibalism; (M2_p); (t/km2/yr); (vertical vector: num_grps X 1)
SINGLE_ECO.predation_total_NoCannibalism        = predation_total_NoCannibalism;      % total predation ON each producer group, p; does NOT include cannibalism; (M2_p); (t/km2/yr); (vertical vector: num_grps X 1)
SINGLE_ECO.consumption_domestic_NOcannibalism	= consumption_domestic_NOcannibalism; % total consumption BY each consumer c; does NOT include cannibalism; (t/km2/yr); (horizontal vector: num_grps X 1)
SINGLE_ECO.consumption_domestic_cannibalism     = consumption_domestic_cannibalism;   % total cannibalism BY each consumer c; (t/km2/yr); (horizontal vector: num_grps X 1); NOTE transpose
SINGLE_ECO.fname_ECOfunction                    = fname_ECOfunction;          	      % name of this version of f_ECOfunction
% SINGLE_ECO.fname_RedistributeCannibalism        = fname_RedistributeCannibalism;	    % FFF this function does not yet pass its name along; name of this version of f_RedistributeCannibalism
SINGLE_ECO.fname_calcEE                         = EcotrophicEfficiency.fname_calcEE;  % file name of this f_calcEEsub-function
% SINGLE_ECO.fname_CalcPredationBudget            = fname_CalcPredationBudget;          % FFF this function does not yet pass its name along; file name of the sub-function f_CalcPredationBudget
% *************************************************************************


% end of m-file************************************************************