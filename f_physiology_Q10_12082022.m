function [Q10_scaler, fname_physiology_Q10] = f_physiology_Q10_12082022(Q10_parameter, T_reference, T_timeseries)
%
% Written by Caren Barceló & Jim Ruzicka
%
% Calculate Q10 metabolic rate scaler
%
% Q10 = van't Hoff's equation: Q10 = (k2/k1)^10/(T-Treference) 
%       where k2 is reaction rate at temperature T, and
%             k1 is reaction rate at temperature Treference.
%
% Rearrange to get:
%   Q10_scaler = (k2/k1) = Q10 ^ ((T-Treference) / 10)
%
% calls:
%   none
%
% takes:
%	Q10_parameter                   (horizontal vector: 1 X num_grps)
%   T_timeseries                    (3D matrix: num_t X num_grps X num_boxes)
%   T_reference                     (horizontal vector: 1 X num_grps)
%
% returns:
%   Q10_scaler                      (3D matrix: num_t X num_grps X num_boxes)
%
% revision date: 12/8/2022
%       6/17/2022 added use of different reference temperatures for different groups and for different sub-region boxes
%       11/28/2022 only cleaned up comments
%       11/29/2022 Changed reference temperatures to NOT be bundled into ECOTRANphysics
%       11/30/2022 Changed temperature timeseries to NOT be bundled into ECOTRANphysics
%       12/6/2022 corrected error by removing unneccesary repmat step

% *************************************************************************
% STEP 1: set operating conditions & parameters
fname_physiology_Q10	= mfilename; % save name of this m-file to keep in saved model results

% step 1a: unpack variables -----------------------------------------------
% temperature_timeseries	= ECOTRANphysics.temperature_timeseries; % (2D matrix: num_t X num_boxes)
% temperature_dates         = ECOTRANphysics.temperature_dates; % (vertical vector: num_t X 1); not used
% num_t                     = ECOTRANphysics.num_t;
% num_boxes                 = ECOTRANphysics.num_boxes;

num_grps_parameter            	= length(Q10_parameter);

[num_t, num_grps, num_boxes]	= size(T_timeseries);

if num_grps_parameter ~= num_grps
    error('ERROR: number of groups in temperature timeseries does not equal number of Q10 parameters')
end

% *************************************************************************
% STEP 2: calculate Q10 scaling factor-------------------------------------
T_reference         = repmat(T_reference, [num_t, 1]);                   % (2D matrix: num_t X num_grps);
T_reference         = repmat(T_reference, [1, 1, num_boxes]);            % (3D matrix: num_t X num_grps X num_boxes)

Q10_repmat         	= repmat(Q10_parameter, [num_t, 1, num_boxes]);      % (3D matrix: num_t X num_grps X num_boxes)

Q10_scaler        	= Q10_repmat .^ ((T_timeseries - T_reference) / 10); % (3D matrix: num_t X num_grps X num_boxes)

% end m-file***************************************************************
