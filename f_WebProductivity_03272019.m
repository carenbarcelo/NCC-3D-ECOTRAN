function production_rate = f_WebProductivity_03272019(TransferEfficiency, EnergyBudget, InputProductionVector, PhysicalLossFraction)
% based on John Steele code
% calculate web production rates given production input at bottom of web
% takes:
%       TransferEfficiency          dimensionless                                       (horizontal vector: 1 X num_grps)
%       EnergyBudget                dimensionless; A_cp                                 (2D matrix: num_grps X num_grps)
%       InputProductionVector       biomass density / time                              (horizontal vector: 1 X num_grps)
%       PhysicalLossFraction        loss fraction of production; export is positive     (horizontal vector: 1 X num_grps)
% returns:
%       WebProductivity             biomass density / time                              (vertical vector: num_grps X 1)
% revision date: 3-27-2019


% *************************************************************************
effective_TransferEfficiency	= TransferEfficiency ./ (1 + PhysicalLossFraction);
inverted_TransferEfficiency     = 1 ./ effective_TransferEfficiency;
K                               = diag(inverted_TransferEfficiency, 0);
production_rate                 = inv(K - EnergyBudget) * InputProductionVector';     % production rates; (vertical vector: num_grps X 1)


% end function*************************************************************