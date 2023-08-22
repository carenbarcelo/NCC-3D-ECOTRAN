function EcotrophicEfficiency = f_calcEE_12292020(randomEwE)
% calculate Ecotrophic Efficiency 
% by Jim Ruzicka
%
% calls:
%       none
%
% takes:
%       randomEwE.
%           CONSUMPTION             consumption matrix from balanced food web; (t/km2/y); (2D matrix: num_grps (producers) X num_grps (consumers))
%       	production              (t/km2/y); (vertical vector: num_grps X 1)
%       	ba
%       	em
%       	EggProduction           (t/km2/y); (vertical vector: num_grps X 1)
%       	looky_eggsANDdetritus
%       	looky_nutrients
%
% returns:
%       EcotrophicEfficiency
%       	ee_eggs
%       	ee_predation        (vertical vector)
%       	ee_ba               (vertical vector)
%       	ee_em               (vertical vector)
%       	ee                  (vertical vector)
%       	fname_calcEE        name of this f_calcEE code file
%
% NOTE: allows for independent definition of egg production
%       (as opposed to EwE eggs-as-detritus definition)
%
% revision date: 12-29-2020
%                   12/29/2020 revision only cleans up comments


% *************************************************************************
% STEP 1: load variable structures (don't alter original values)-----------
fname_calcEE            = mfilename; % save name of this m-file to keep in saved model results
display(['   Running: ' fname_calcEE])

CONSUMPTION           	= randomEwE.CONSUMPTION;            % (t/km2/y); (2D matrix: num_grps X num_grps); NOTE: there ARE non-zero entries for eggs, detritus, & fleets but all zero entries for nutrients & primary producers
production              = randomEwE.production;             % (t/km2/y); (vertical vector: num_grps X 1)
ba                      = randomEwE.ba;                     % (t/km2/y); (vertical vector: num_grps X 1)
em                      = randomEwE.em;                     % (t/km2/y); (vertical vector: num_grps X 1)
EggProduction           = randomEwE.EggProduction;          % (t/km2/y); (vertical vector: num_grps X 1)
looky_eggsANDdetritus	= randomEwE.looky_eggsANDdetritus;
looky_nutrients         = randomEwE.looky_nutrients;
% *************************************************************************



% *************************************************************************
% STEP 2: correct nutrient, egg, detritus, & fleet terms-------------------
CONSUMPTION(:, looky_eggsANDdetritus)	= 0;	% set egg & detritus clms to 0 in CONSUMPTION (don't include flow to eggs & detritus as predation)
CONSUMPTION(:, looky_nutrients)         = 0;   	% set nutrient clms to 0 in CONSUMPTION (don't include flow to nutrients as predation)
em(looky_eggsANDdetritus)               = 0;	% set egg & detritus clms to 0 in em (don't include flow to eggs & detritus as emigration)
em(looky_nutrients)                     = 0;   	% set nutrient clms to 0 in em (don't include flow to nutrients as emigration)
production(looky_nutrients)             = 1;    % this prevents div/0 errors for nutrients; this change is isolated to this function
% *************************************************************************



% *************************************************************************
% STEP 3: calculate total predation & fishing pressure on each group-------
TotalPredation          = sum(CONSUMPTION, 2); % total predation on each grp; (t/km2/y); (2D matrix: num_grps X num_grps)
% *************************************************************************



% *************************************************************************
% STEP 4: Ecotrophic Efficiencies (distribution of production)-------------
ee_eggs                 = EggProduction  ./ production; % (EGG_p / P_p) = fraction of production by each group, p, directed to eggs (or reproduction)
ee_predation            = TotalPredation ./ production; % (M2_p  / P_p) = fraction of production by each group, p, passed upwards in food web via predation; note omission of cannibalism
ee_ba                   = ba             ./ production; % (BA_p  / P_p) = fraction of production by each group, p, directed to biomass accumulation
ee_em                   = em             ./ production; % (EM_p  / P_p) = fraction of production by each group, p, directed to emigration (immigration is negative production)
ee                      = ee_eggs + ee_predation + ee_ba + ee_em; % (EE_p = (EGG_p + M2_p + EM_p + BA_p) / P_p) = Ecotrophic Efficiency; NOTE: EwE definition does not include eggs
% *************************************************************************



% *************************************************************************
% STEP 5: some error checking (alert for NaN errors)-----------------------
looky_NaN               = find(isnan(ee));
if ~isempty(looky_NaN)
    display('   -->WARNING: NaN error found in recalculated Ecotrophic Efficiency (ee).')
end
% *************************************************************************



% *************************************************************************
% STEP 6: pack results for export------------------------------------------
EcotrophicEfficiency.ee_eggs         = ee_eggs;         % eggs (or other reproduction)
EcotrophicEfficiency.ee_predation    = ee_predation;
EcotrophicEfficiency.ee_ba           = ee_ba;           % biomass accumulation (local production)
EcotrophicEfficiency.ee_em           = ee_em;           % emigration
EcotrophicEfficiency.ee              = ee;
EcotrophicEfficiency.fname_calcEE    = fname_calcEE;    % name of this f_calcEE code file
% *************************************************************************


% end m-file***************************************************************