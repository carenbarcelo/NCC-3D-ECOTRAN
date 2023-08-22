function [production_initial, FunctionName] = f_InitialProductionRates_02012022(ODEinput, production_driver, t_current)
% calculate initial or mean production conditions
% by Jim Ruzicka
%
% calls:
%       f_WebProductivity_03272019
%
% takes: 
%       ODEinput
%           EnergyBudget            consumption fates; (proportions); (3D matrix: num_grps X num_grps X num_boxes)
%           PhysicalLossFraction	PhysicalLossFraction @ t; (2D matrix: time X num_grps)
%           TransferEfficiency      (horizontal vector: 1 X num_grps)
%           looky_nutrients         (addresses)
%           num_boxes               (scaler)
%           num_grps                (scaler)
%       production_driver           production driver vector for each box; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
%       t_current                   day of year
%
% returns:
%       production_initial          initial or mean production (actually consumption) rates; (mmole N/m3/d); (2D matrix: num_grps X num_boxes)
%       FunctionName                name of this f_InitialProductionRates function
%
% revision date: 2-1-2022


% *************************************************************************
% STEP 1: unpack givens----------------------------------------------------
% step 1a: log function file name -----------------------------------------
fname_InitialProductionRates	= mfilename; % save name of this m-file to keep in saved model results
display(['Running: ' fname_InitialProductionRates])
FunctionName                    = fname_InitialProductionRates; % name of this f_InitialProductionRates function

% step 1b: unpack givens --------------------------------------------------
EnergyBudget                = ODEinput.EnergyBudget;            % (proportions); (3D matrix: num_grps X num_grps X num_boxes)
PhysicalLossFraction        = ODEinput.PhysicalLossFraction;	% PhysicalLossFraction @ t; (2D matrix: num_t X num_grps) QQQ should define for individual boxes (but usually = 0 in time-dynamic runs)
TransferEfficiency          = ODEinput.TransferEfficiency;      % (3D matrix: 1 X num_grps X num_boxes)
looky_nutrients             = ODEinput.looky_nutrients;
looky_plgcNH4               = ODEinput.looky_plgcNH4;
looky_bnthNH4               = ODEinput.looky_bnthNH4;
num_boxes                   = ODEinput.num_boxes;
num_grps                    = ODEinput.num_grps;
% *************************************************************************



% *************************************************************************
% STEP 2: prep-------------------------------------------------------------
% step 2a: ----------------------------------------------------------------
PhysicalLossFraction_t                  = PhysicalLossFraction(t_current, :);	% (horizontal vector: 1 X num_grps)
production_initial                      = zeros(num_grps, num_boxes);           % initialize; (mmole N/m3/d); (2D matrix: num_grps X num_boxes)

% step 2b: evaluate type of driver (nutrient vs living group) -------------
%          NOTE: TURN OFF NUTRIENT RECYCLING IF DRIVEN BY ANYTHING OTHER THAN NO3
% % % TestInput_base                          = production_driver; % production driver vector for each Box; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
% % % TestInput_base(:, looky_nutrients, :)   = [];
% % % TestInput_base                          = sum(TestInput_base, 3); % flatten TestInput_base; if evaluation conditions hold for 1 box, we say that they hold for all boxes
% % % if max(TestInput_base > 0)
% % % 	% NOTE: this method cannot handle mixed nutrient & non-nutrient drivers in the same run
% % %     disp(['   -->NOTICE in ' fname_InitialProductionRates ': non-NO3 driver. Nutrient uptake in EnergyBudget is DEACTIVED (nutrient columns set to 0)']);
% % %     EnergyBudget(:, looky_nutrients, :) = 0; % shut off recycling; turn off flow of nutrients to any other box (i.e., primary producers); NOTE: there is still flow INTO nutrients but no flow FROM nutrients
% % % else
% % %     disp(['   -->NOTICE in ' fname_InitialProductionRates ': NO3 driver. Nutrient uptake & recycling in EnergyBudget is ACTIVE']); 
% % % end

disp(['   -->NOTICE in ' fname_InitialProductionRates ': Ammonium uptake in EnergyBudget is DEACTIVED (NH4 columns set to 0)']);
EnergyBudget(:, [looky_plgcNH4 looky_bnthNH4], :) = 0; % shut off recycling; turn off flow of NH4 to any other box (i.e., primary producers, nitrification to NO3); NOTE: there is still flow INTO NH4 but no flow FROM NH4
% *************************************************************************



% *************************************************************************
% step 3: calculate initial production rates for each model box domain ----
%         (mmole N/m3/d); (2D matrix: num_grps X num_boxes)
%         NOTE: no ZEROS are allowed else NaN errors else NaN errors from functional response equations
%         NOTE: when using the full BioenergeticBudget matrix, 
%                 this is really total consumption (Q) flowing into each box rather than production (P)
for box_loop = 1:num_boxes
    current_TransferEfficiency      = TransferEfficiency(1, :, box_loop);	% (horizontal vector: 1 X num_grps)
    current_ProductionDriver        = production_driver(1, :, box_loop);	% production driver for current Box; (mmole N/m3/d); (horizontal vector: 1 X num_grps)
    current_EnergyBudget            = EnergyBudget(:, :, box_loop);         % (2D matrix: num_grps X num_grps)
    production_initial(:, box_loop) = f_WebProductivity_03272019(current_TransferEfficiency, current_EnergyBudget, current_ProductionDriver, PhysicalLossFraction_t); % production (actually CONSUMPTION) rates as "initial" or "mean" conditions; (mmole N/m3/d); (vertical vector: num_grps X 1)
end
% *************************************************************************



% *************************************************************************
% STEP 4: transfer benthic NH4 into pelagic NH4 pool-----------------------
%         pool pelagic NH4 & benthic NH4 in all boxes
production_initial(looky_plgcNH4, :)        = production_initial(looky_plgcNH4, :) + production_initial(looky_bnthNH4, :); % (mmole N/m3/d); (2D matrix: num_E2Egroups X num_boxes)
production_initial(looky_bnthNH4, :)        = 0;                                    % zero out benthic NH4 pool after combining with the pelagic NH4 pool; (mmole N/m3/d); (2D matrix: num_E2Egroups X num_boxes)
% *************************************************************************


% m-file*******************************************************************