function SurfaceRadiation = f_LightIntensity_12112020(LightParams)
% calculate light intensity
% Solar radiation, based on Brock (1981, Ecol. Modelling 14:1-19)
% NOTE: t_grid must be defined so that t_grid = 1 is Jan. 1
%
% calls:
%       none
%
% takes:
%	LightParams
%       t_grid          ordinal, decimal day	(vertical vector)
%       latitude        latitude (decimal degress north)	(scalar)
%       SolarConstant	Solar constant (W/m2) recent average = 1366 (Froehlich & Lean (1998, Geographical Research Letters 25(23):4377-4380)) (scalar)
%       par_frac        fraction of 40(spectrum range???) that is PAR, clear atmos. (Fasham, Ducklow, & McKelvie (1990, J. Mar. Res. 48:591-639), p. 604)	(scalar)
%       surf_trans      1-reflectance (Parsons, Takahachi, & Hargrave (1984, Biological Oceanographic Processes, 3rd ed.), p. 69)	(scalar)
%       cloud_trans     1.0 - (0.7 * cloud_cover); (Evans and Parslow (1985, Biol. Oceanogr. 3(3):327-347), p. 344); (scalar)
%
% returns:
%	SurfaceRadiation
%       sunrise                 time of sunrise (time in h from midnight; 12 = noon, set as a default)
%       sunset                  time of sunset  (time in h from midnight; 12 = noon, set as a default)
%       day_length              hours of daylight
%       current_light           surface solar raditation at current time (W/m2) (vertical vector)
%       daily_integrated_light  daily integrated solar raditation at ocean surface (W m^-2 d^-1) (vertical vector)
%       daily_average_light     daily mean solar raditation at ocean surface averaged across 24 hours (W m^-2 h^-1) (vertical vector)
%       fname_LightIntensity	file name of this version of f_LightIntensity
%
% revision date: 12-11-2020
%           12/11/2020 removed "round" function from local_time calculation


% *************************************************************************
% STEP 1: unpack operating parameters--------------------------------------
fname_LightIntensity	= mfilename; % save name of this f_LightIntensity function
display(['   Running: ' fname_LightIntensity])

t_grid                  = LightParams.t_grid;   % ordinal, decimal day; (vertical vector)
latitude                = LightParams.latitude; % latitude (decimal degress north); (scalar)

SolarConstant           = LightParams.SolarConstant;	% Solar constant (W/m2) recent average = 1366 (Froehlich & Lean (1998, Geographical Research Letters 25(23): 4377-4380)); (scalar)
par_frac                = LightParams.par_frac;         % fraction of 40 that is PAR, clear atmos. (Fasham, Ducklow, & McKelvie (1990, J. Mar. Res. 48:591-639), p. 604); (scalar)
surf_trans              = LightParams.surf_trans;       % 1-reflectance (Parsons, Takahachi, & Hargrave (1984, Biological Oceanographic Processes, 3rd ed.), p. 69); (scalar)
cloud_trans             = LightParams.cloud_trans;      % 1.0 - (0.7 * cloud_cover); (Evans and Parslow (1985, Biol. Oceanogr. 3(3):327-347), p. 344); (scalar)
% *************************************************************************



% *************************************************************************
% STEP 2: calculate working parameters-------------------------------------
%         NOTE: trigonmetric functions are set for degrees
julian_day	= floor(t_grid);                                % (vertical vector); NOTE: not necessary to correct t for multiple years
D1          = 23.45 * sind(360 * (284 + julian_day) / 365);	% declination, angular distance between sun and equator at solar noon (degrees); eq. 1; (vertical vector)
W1          = acosd(-(tand(D1) * tand(latitude)));          % sunset hour-angle, angle between setting sun and the south point (degrees); eq. 3
day_length	= (W1 / 15) * 2;                                % day length (h); eq. 4; (vertical vector)
sunrise     = 12 - 1/2 * day_length;                        % time of sunrise; (vertical vector)
sunset      = 12 + 1/2 * day_length;                        % time of sunset; (vertical vector)
R1          = 1 ./ (1 + (0.033 * cosd(360 * julian_day / 365))).^0.5; % radius vector, used to correct for ellipticity of Earth's orbit (dimensionless); eq. 2; (vertical vector)
% *************************************************************************



% *************************************************************************
% STEP 3: current light intensity at local_time----------------------------
%         local_time is time in h from midnight (12 = noon, set as a default)
local_time      = ((t_grid - julian_day) * 24);                 % local_time is time in h from midnight (12 = noon, set as a default); (vertical vector)
W2              = (local_time - 12) * 15;                       % hour-angle (degrees); eq. 5; (vertical vector)
I1              = (SolarConstant ./ R1.^2) .* (sind(D1) .* sind(latitude) + cosd(D1) .* cosd(latitude) .* cosd(W2)); % top-of-atmosphere solar radiation at local_time (W/m2); eq. 7; (vertical vector)
I1(I1 < 0)      = 0;                                            % filter out negative values (caused by after sunset times); (vertical vector)
current_light	= I1 * (par_frac * surf_trans * cloud_trans);	% surface solar raditation at current time (W/m2); (vertical vector)
% *************************************************************************



% *************************************************************************
% STEP 4: daily integrated light intensity on julian_day-------------------
I6                      = (24 / pi) * (SolarConstant ./ R1.^2) .* (W1 .* (pi / 180) .* sind(latitude) .* sind(D1) + sind(W1) .* cosd(latitude) .* cosd(D1));  % integrated daily top-of-atmosphere (W m^-2 d^-1); (vertical vector)
I6_average              = (I6 / 24) * par_frac * surf_trans * cloud_trans;	% daily integrated solar raditation at ocean surface averaged across 24 hours (W m^-2 h^-1); (vertical vector)
daily_integrated_light	= I6 * (par_frac * surf_trans * cloud_trans);       % daily integrated solar raditation at ocean surface (W m^-2 d^-1)
daily_average_light     = I6_average;                                       % daily mean solar raditation at ocean surface averaged across 24 hours (W m^-2 h^-1)
% *************************************************************************



% *************************************************************************
% STEP 5: pack results for export------------------------------------------
SurfaceRadiation.sunrise                    = sunrise;                  % time of sunrise (time in h from midnight; 12 = noon, set as a default)
SurfaceRadiation.sunset                     = sunset;                   % time of sunset  (time in h from midnight; 12 = noon, set as a default)
SurfaceRadiation.day_length                 = day_length;               % hours of daylight

SurfaceRadiation.local_time                 = local_time;               % time in h from midnight (12 = noon, set as a default)
SurfaceRadiation.current_light              = current_light;            % surface solar raditation at current time (W/m2)

SurfaceRadiation.daily_integrated_light     = daily_integrated_light;	% daily integrated solar raditation at ocean surface (W m^-2 d^-1)
SurfaceRadiation.daily_average_light        = daily_average_light;      % daily mean solar raditation at ocean surface averaged across 24 hours (W m^-2 h^-1)

SurfaceRadiation.fname_LightIntensity       = fname_LightIntensity;     % file name of this version of f_LightIntensity
% *************************************************************************


% end m-file***************************************************************