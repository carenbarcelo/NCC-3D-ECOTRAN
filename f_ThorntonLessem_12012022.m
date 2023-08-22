function [q_TemperatureScaler, fname_ThorntonLessem] = f_ThorntonLessem_12012022(ThorntonLessem_parameters, T_timeseries)
% Thornton & Lessem temperature adjustment model for consumption
%
% Written by Jim Ruzicka & Caren Barcel—
%
% [Based on: from Megrey, Rose, Werner, Klumb, Hay (2000). A generalized fish bioenergetics/biomass model with an application to Pacific herring [https://www.pices.int/publications/scientific_reports/report20/rep20_model_rex_megrey.pdf]]
%
% calls: none
%
% takes:
%   ECOTRANphysics
%   ThorntonLessem_parameters                                               (2D matrix: 8 X num_grps)
%       row 1 = te1         	Temperature for xk1 (in ºC)
%       row 2 = te2           	Temperature for xk2 (in ºC)
%       row 3 = te3          	Temperature for xk3 (in ºC)
%       row 4 = te4         	Temperature for xk4 (in ºC)
%       row 5 = xk1          	Proportion of q_max at te1      
%       row 6 = xk2           	Proportion of q_max at te2      
%       row 7 = xk3           	Proportion of q_max at te3      
%       row 8 = xk4           	Proportion of q_max at te4
%	T_timeseries                                                            (3D matrix: num_t X num_grps X num_boxes)
%
% returns:
%   q_TemperatureScaler     consumption scaling factor (relative to q_max)	(3D matrix: num_t X num_grps X num_boxes)
%
% revision date: 12/1/2022
%   12/1/2022 Changed temperature timeseries to NOT be bundled into ECOTRANphysics; updated temperature time-series to allow temperatures from different depth zones to be applied to different groups depending on their life-history (useful for vertically-integrated food web technique) 


% *************************************************************************
% STEP 1: unpack parameters------------------------------------------------
fname_ThorntonLessem	= mfilename; % save name of this m-file to keep in saved model results
display(['Running: ' fname_ThorntonLessem])

[num_t, num_grps, num_boxes] = size(T_timeseries);

te1                     = ThorntonLessem_parameters(1, :); %  Temperature for xk1 (in ºC); (horizontal vector: 1 X num_grps)
te2                     = ThorntonLessem_parameters(2, :); %  Temperature for xk2 (in ºC); (horizontal vector: 1 X num_grps)
te3                     = ThorntonLessem_parameters(3, :); %  Temperature for xk3 (in ºC); (horizontal vector: 1 X num_grps)
te4                     = ThorntonLessem_parameters(4, :); %  Temperature for xk4 (in ºC); (horizontal vector: 1 X num_grps)
xk1                     = ThorntonLessem_parameters(5, :); %  Proportion of q_max at te1; (horizontal vector: 1 X num_grps)
xk2                     = ThorntonLessem_parameters(6, :); %  Proportion of q_max at te2; (horizontal vector: 1 X num_grps)
xk3                     = ThorntonLessem_parameters(7, :); %  Proportion of q_max at te3; (horizontal vector: 1 X num_grps)
xk4                     = ThorntonLessem_parameters(8, :); %  Proportion of q_max at te4; (horizontal vector: 1 X num_grps)
% *************************************************************************



% *************************************************************************
% STEP 2: replicate parameters by time and by subregions and temperature by groups
check_num_grps      	= length(te1);

if check_num_grps ~= num_grps
    error('ERROR: number of Thornton-Lessem parameters DOES NOT EQUAL number of groups')
end

te1                     = repmat(te1, [num_t, 1, num_boxes]); % (3D matrix: num_t X num_grps X num_boxes)
te2                     = repmat(te2, [num_t, 1, num_boxes]); % (3D matrix: num_t X num_grps X num_boxes)
te3                     = repmat(te3, [num_t, 1, num_boxes]); % (3D matrix: num_t X num_grps X num_boxes)
te4                     = repmat(te4, [num_t, 1, num_boxes]); % (3D matrix: num_t X num_grps X num_boxes)
xk1                     = repmat(xk1, [num_t, 1, num_boxes]); % (3D matrix: num_t X num_grps X num_boxes)
xk2                     = repmat(xk2, [num_t, 1, num_boxes]); % (3D matrix: num_t X num_grps X num_boxes)
xk3                     = repmat(xk3, [num_t, 1, num_boxes]); % (3D matrix: num_t X num_grps X num_boxes)
xk4                     = repmat(xk4, [num_t, 1, num_boxes]); % (3D matrix: num_t X num_grps X num_boxes)
% *************************************************************************



% *************************************************************************
% STEP 3: calculate consumption scaling factor-----------------------------
tt5                     = 1 ./ (te2 - te1);                                 % (3D matrix: num_t X num_grps X num_boxes)
t5                      = tt5 .* log(xk2 .* (1.0 - xk1) ./ (0.02 * xk1));	% (3D matrix: num_t X num_grps X num_boxes)
t4                      = exp(t5 .* (T_timeseries - te1));        	% (3D matrix: num_t X num_grps X num_boxes)

tt7                     = 1 ./ (te4 - te3);                                 % (3D matrix: num_t X num_grps X num_boxes)
t7                      = tt7 .* log(xk3 .* (1.0 - xk4) ./ (0.02 * xk4));	% (3D matrix: num_t X num_grps X num_boxes)
t6                      = exp(t7 .* (te4 - T_timeseries));         	% (3D matrix: num_t X num_grps X num_boxes)

gcta                    = (xk1 .* t4) ./ (1.0 + xk1 .* (t4 - 1.0));         % (3D matrix: num_t X num_grps X num_boxes)
gctb                    = (xk4 .* t6) ./ (1.0 + xk4 .* (t6 - 1.0));         % (3D matrix: num_t X num_grps X num_boxes)

q_TemperatureScaler     = gcta .* gctb;                                     % (3D matrix: num_t X num_grps X num_boxes)

q_TemperatureScaler(isnan(q_TemperatureScaler)) = 1; % set scalers for un-parameterized groups (NaNs) to 1;
q_TemperatureScaler(isinf(q_TemperatureScaler)) = 1; % set scalers for un-parameterized groups (NaNs) to 1; errors probably all show up as NaNs
% *************************************************************************


% end m-file***************************************************************
