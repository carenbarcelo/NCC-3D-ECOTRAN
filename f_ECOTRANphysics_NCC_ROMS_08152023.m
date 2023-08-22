function [ECOTRANphysics] = f_ECOTRANphysics_NCC_ROMS_08152023(PHYSICSinput)
% returns ROMS Advection exchanges between boxes
% use for NCC upwelling setting
%
% calls:
%       f_ROMS_GridPrep_NCC_08152023            Map ROMS grid to ECOTRAN grid
%       f_ROMS_FluxPrep_NCC_08152023            express ROMS fluxes in DESTINY<--SOURCE format on ECOTRAN grid
%             f_CompactFluxTimeSeries_11182019      compact physical flux time-series, arranged as 3D matrix (time X source box X destiny box), eliminate source & destiny infomration for any subregions that never communicate
%             f_UnCompactFluxTimeSeries_12112019    UnCompact a previously compacted physical flux time-series to provide IMPORT & EXPORT fluxes for each box and the domain as a whole
%                 f_calcNetFlux_12112019                calculate net flux into and net flux out of each model box and across outer domain boundaries
%             f_EvaluateFluxBalance_11262021        examine for flux time-series for imbalances IN & OUT of individual boxes and IN & OUT of the overall domain
%       f_LightIntensity_12112020           instantaneous (W/m2), daily mean averaged across 24 h (W m^-2 h^-1), & daily integrated (W m^-2 d^-1) solar raditation at ocean surface; (vertical vector: num_t X 1)
%
% takes:
%       internally defined time-series of advection & mixing rates
%
% returns ECOTRANphysics:
%       AdvectionMatrix         advection rate = volume transported per day; (m3/d)
%                               (3D matrix: time X  source box X destination box)
%                               (values can be positive or negative)
%       VerticalMixMatrix       vertical mixing rate = volume mixed per day; (m3/d)
%                               (3D matrix: time X  source box X destination box)
%                               (values can be positive or negative; sign defines direction of net mixing)
%
% NOTE: What happens in cases when there is a complete advective washout of a spatial box?)
%       I added cap to advection rate to prevent >100% box washout (but do not account for river flux)
%       Added warning for days when AdvectionMatrix > BoxVolume; see ECOTRANphysics.BoxWashoutTime
%
% revision date: 5-11-2023
%       5/30/2022 stack multiple ROMS years
%       8/29/2022 calls revised f_ROMS_FluxPrep_NCC_08292022
%       9/2/2022 retained uncorrected flux error and total horizontal flux for comparison
%       9/5/2022 cleaned up code
%       9/13/2022 processing temperature info
%       10/16/2022 processing BGC info
%       11/28/2022 added VerticalConnectivity_stack in f_ROMS_GridPrep telling which boxes lay directly under each surface box; clm 1 is suface, clm 2 is depth layer 2, final column is the bottom
%       11/30/2022 corrected f_ROMS_FluxPrep_ for weighted averaging of aggregated BGC info
%       12/02/2022 corrected BGC interpolation error (interp1 doesn't like NaNs)
%       5/11/2022 Moved all directory info into the PHYSICSinput structure so these can be defined in ECOTRANdynamic script


% *************************************************************************
% STEP 1: basic domain parameters (NCC upwelling)--------------------------
fname_PhysicalModel      = mfilename; % save name of this m-file to keep in saved model results
display(['Running: ' fname_PhysicalModel])

% files & directories
ROMSfile_directory      = PHYSICSinput.ROMSfile_directory;      % directory with ROMS netcdf files
filename_list           = PHYSICSinput.filename_list;           % list of ROMS output files to read (one file for each year)
% -------------------------------------------------------------------------


% step 1a: physical parameters --------------------------------------------
latitude                = 44.6516667;              % latitude of NH line (decimal degress north) based on LTOP transects
rho_water               = 1000;                    % seawater density; (kg/m3)
rho_air                 = 1.220;                   % air density (kg/m3) (Schwing et al. 1996)
coriollis               = 1.025 * 10^(-4);         % coriollis factor (1/s); (NH-Line NNPPZD model: Ruzicka et al. 2011)
midmonth_doy            = [15 45 74 105 135 166 196 227 258 288 319 349]; % mid-month day-of-year

% step 1b: time domain and vectors ----------------------------------------
datestart               = PHYSICSinput.datestart;  % enter starting date
dateend                 = PHYSICSinput.dateend;    % enter ending date
dt                      = PHYSICSinput.dt;         % t-step (d)
t_grid                  = PHYSICSinput.t_grid;     % t_grid; (days, starting from 1); (vertical vector)
num_t                   = length(t_grid);          % length of t_grid; (scaler)
min_t                   = min(t_grid);
max_t                   = max(t_grid);
% num_years               = ceil(max_t/365);

datevec_start           = datevec(datestart);
datevec_end             = datevec(dateend);
yearvec                 = datevec_start:1:datevec_end;
num_years               = length(unique(yearvec));

ROMS_refTime_matlab     = datenum([1900, 1, 1, 0, 0, 0]);
t_grid_real             = PHYSICSinput.t_grid_real;
% *************************************************************************





% *************************************************************************
% STEP 2: conversion factors-----------------------------------------------
atomic_mass_N             = 14.00674;   % (g N / mole N)
atomic_mass_C             = 12.0107;    % (g C / mole C)

Chla_to_N_diatoms         = 2.19;       % (mg Chla mmoles N^-1) (Dickson & Wheeler 1995 L&O 40(3):533-543)
% Chla_to_N_diatoms         = 2.5;        % 2.49 to 2.67 (mg Chla mmoles N^-1) (Chan 1980 J Phycol 16:428-432)
Chla_to_N_dinoflagellates = 0.84;       % 0.73 to 0.95 (mg Chla mmoles N^-1) (Chan 1980 J Phycol 16:428-432)

% C_to_N_phytoplankton      = 7.3;        % general phytoplankton C:N (moles C / moles N) (Geider & LaRoche 2002 Eur. J. Phycol 37:1-17)
C_to_N_phytoplankton      = (106/16);   % Redfield ratio = C:N:P = 106:16:1 ---> C:N = 106:16 = 6.625 mole C/mole N
C_to_N_zooplankton        = 5;          % general zooplankton C:N (moles C moles N^-1) for mesozooplankton (W. Peterson pers com); C:N=5.16-5.26 [Schneider 1990 (Mar Biol 106:219-225)]; C:N=3.9 [Uye 1982]

WWT_to_C                  = 8.77;       % (g WWT / g C); (Steele et al. 2007), NOTE: I could not find this value in 2007 citation, must have gotten value directly from John, value is close to other lit. values; {Strickland, 1966 #1521} says "From the foregoing, it appears that the 10-12 mg algal weight:mg C is the most suitable factor for estimating phytoplankton biomass as ?wet weight?";
mgWWT_per_mmolesN         = (C_to_N_phytoplankton * atomic_mass_C * WWT_to_C); % (mg WWT/mmoles N)

mgram_to_gram             = 1000;       % (mg/g)
gram_to_ton               = 1000000;	% (g/t)
mgram_to_ton              = mgram_to_gram * gram_to_ton; % (mg/t)

W_to_microE               = 4.6;        % convert (W m^-2) to (microE m^-2 s^-1) at 400-700nm (Brock 1981 Ecol Model 14: 1-19)
% *************************************************************************





% *************************************************************************
% STEP 3: prepare ROMS grid------------------------------------------------
ROMSgrid                    = f_ROMS_GridPrep_NCC_08152023(PHYSICSinput);
num_boxes                   = ROMSgrid.num_boxes;
num_domains                 = ROMSgrid.num_domains; % number of geographic domains (does not consider depth layers)
num_z                       = ROMSgrid.num_agg_z; % number of depth layers

VerticalConnectivity        = ROMSgrid.VerticalConnectivity; % 0 or 1; (2D matrix: DESTINY-->(num_boxes+1) (148) X SOURCE-->(num_boxes+1) (148))
VerticalConnectivity_stack	= ROMSgrid.VerticalConnectivity_stack; % (2D matrix: num_domains X num_agg_z); tells which boxes lay directly under each surface box; clm 1 is suface, clm 2 is depth layer 2, final column is the bottom
looky_BottomBoxes           = ROMSgrid.looky_BottomBoxes; % addresses of bottom (benthic) boxes
VerticalConnectivity        = permute(VerticalConnectivity, [3, 2, 1]); % connectivity between source & destiny boxes; 0 or 1; (3D matrix: 1 X SOURCE-->(num_boxes+1) (148) X DESTINY-->(num_boxes+1) (148))
VerticalConnectivity(1, looky_BottomBoxes, end)	= 1; % sinking from box 1 to benthos; NOTE: not necessary if there is a specifically defined benthic box in the physical model (e.g., GoMexOcn); QQQ need to add benthic box to NCC and remove export sinking out of boundary
% *************************************************************************





% *************************************************************************
% STEP 4: process ROMS fluxes & biogeochemical model for specified years---

% step 4a: process ROMS results -------------------------------------------

% NOTE: row length = time initialization will be short by leap-years, but varaibles will grow within year_loop as necessary
build_ocean_time        = zeros((365*num_years), 1);                              % initialize; (vertical vector: (365*num_years) X 1)
build_BoxVolume         = zeros((365*num_years), num_boxes);                      % initialize; (2D matrix: (365*num_years) X num_boxes)
build_ADVECTION         = zeros((365*num_years), (num_boxes+1), (num_boxes+1));   % initialize; (3D matrix: (365*num_years) X num_boxes+1 X num_boxes+1)
build_HORIZONTALMIXING	= zeros((365*num_years), (num_boxes+1), (num_boxes+1));   % initialize; (3D matrix: (365*num_years) X num_boxes+1 X num_boxes+1)
build_VERTICALMIXING	= zeros((365*num_years), (num_boxes+1), (num_boxes+1));   % initialize; (3D matrix: (365*num_years) X num_boxes+1 X num_boxes+1)
build_SINKING           = zeros((365*num_years), (num_boxes+1), (num_boxes+1));   % initialize; (3D matrix: (365*num_years) X num_boxes+1 X num_boxes+1)
build_error_record      = zeros(num_domains,     2,             (365*num_years)); % initialize; (3D matrix: num_domains     X 2           X (365*num_years)); NOTE: all flux error should be restricted to 1 ROMS depth layer, if there is more than 1 depth layer with error then ROMSflux.error_record will be a 4D matrix and this build varialbe will need to be adjusted to accomodate 4D

build_temperature     	= zeros((365*num_years), num_boxes);   % initialize; (2D matrix: (365*num_years) X num_boxes)
build_diatom            = zeros((365*num_years), num_boxes);   % initialize; (2D matrix: (365*num_years) X num_boxes)
build_nanophytoplankton	= zeros((365*num_years), num_boxes);   % initialize; (2D matrix: (365*num_years) X num_boxes)

pointer_1               = 0; % initialize
num_t_ROMS              = 0; % initialize

for year_loop = 1:num_years
    
    pointer_1           = pointer_1 + 1;
    
    current_year        = yearvec(year_loop);
    looky_file          = find(~cellfun(@isempty, strfind(filename_list, num2str(current_year))));
    if isempty(looky_file)
        error('ERROR: missing ROMS netcdf file')
    end
    
    readFile_FluxYear	= [ROMSfile_directory filename_list{looky_file}];
    disp(['Processing year: ' num2str(current_year) ' ROMS file: ' filename_list{looky_file}])
    
    ROMSflux            = f_ROMS_FluxPrep_NCC_08152023(ROMSgrid, readFile_FluxYear); % NOTE: this code will provide compacted fluxes for the current ROMS year, but compaction step will be repeated AFTER all ROMS years are stacked
        
    current_num_t_ROMS	= ROMSflux.num_t_ROMS; % number of ROMS time-points in current year
    num_t_ROMS          = num_t_ROMS + current_num_t_ROMS;
    pointer_2           = pointer_1 + current_num_t_ROMS - 1;
    
    build_ocean_time(pointer_1:pointer_2, 1)                                        = ROMSflux.ocean_time; % seconds since 1900-01-01 00:00:00; (vertical vector: num_t_ROMS X 1)
    build_BoxVolume(pointer_1:pointer_2, 1:num_boxes)                               = ROMSflux.BoxVolume; % box volume; (m3); (2D matrix: num_t_ROMS (366) X num_boxes (147))
    
    build_ADVECTION(pointer_1:pointer_2, 1:(num_boxes+1), 1:(num_boxes+1))          = ROMSflux.ADVECTION;        % advection rate;                (m3/d); (3D matrix: current_num_t_ROMS X SOURCE (num_boxes + boundary) X DESTINY (num_boxes + boundary))
    build_HORIZONTALMIXING(pointer_1:pointer_2, 1:(num_boxes+1), 1:(num_boxes+1))	= ROMSflux.HORIZONTALMIXING; % horizontal mixing rate;        (m3/d); (3D matrix: current_num_t_ROMS X SOURCE (num_boxes + boundary) X DESTINY (num_boxes + boundary))
    build_VERTICALMIXING(pointer_1:pointer_2, 1:(num_boxes+1), 1:(num_boxes+1))     = ROMSflux.VERTICALMIXING;   % vertical mixing rate;          (m3/d); (3D matrix: current_num_t_ROMS X SOURCE (num_boxes + boundary) X DESTINY (num_boxes + boundary))
    build_SINKING(pointer_1:pointer_2, 1:(num_boxes+1), 1:(num_boxes+1))            = ROMSflux.SINKING;          % sinking source box floor area; (m2);   (3D matrix: current_num_t_ROMS X SOURCE (num_boxes + boundary) X DESTINY (num_boxes + boundary))
    
    build_error_record(:, 1:2, pointer_1:pointer_2, :)                           	= ROMSflux.error_record; % total flux error; (m3/s); (4D matrix: DESTINY-->num_domains+1 X 2 [total horizontal flux; flux error] X num_t_ROMS (366) X length(looky_error_depth)); NOTE: all flux error should be restricted to 1 ROMS depth layer, if there is more than 1 depth layer with error then ROMSflux.error_record will be a 4D matrix and this build varialbe will need to be adjusted to accomodate 4D
    
    % BGC terms
    build_temperature(pointer_1:pointer_2, :)                                       = ROMSflux.ROMS_temperature; % temperature time-series (deg C); (2D matrix: num_t_ROMS X num_boxes)
    build_diatom(pointer_1:pointer_2, :)                                            = ROMSflux.ROMS_diatom; % diatom time-series (mmole N/m3); (2D matrix: num_t_ROMS X num_boxes)
    build_nanophytoplankton(pointer_1:pointer_2, :)                                 = ROMSflux.ROMS_nanophytoplankton; % nanophytoplankton time-series (mmole N/m3); (2D matrix: num_t_ROMS X num_boxes)

    pointer_1           = pointer_2;
    
end % year_loop
% -------------------------------------------------------------------------


% step 4b: adjust variable names ------------------------------------------
ROMS_time                       = build_ocean_time; % seconds since 1900-01-01 00:00:00; (vertical vector: num_t_ROMS X 1)
% % % ROMS_time = [1:86400:(86400*num_t_ROMS)]'; % (vertical vector: num_t_ROMS X 1)

temp_ROMS_time = zeros(pointer_2, 1); % initialize; (vertical vector: pointer_2 X 1)
for time_loop = 1:pointer_2
    temp_ROMS_time(time_loop) = addtodate(ROMS_refTime_matlab, ROMS_time(time_loop), 'second');
end

BoxVolume                       = build_BoxVolume; % box volume; (m3); (2D matrix: num_t_ROMS X num_boxes)
BoxFloorArea                    = ROMSflux.BoxFloorArea; % box floor area; (m2); (horizontal vector: 1 X num_boxes); NOTE: this does not change over time so just use the last ROMS year; FFF move to f_ROMS_GridPrep_NCC_ code
% -------------------------------------------------------------------------


% step 4c: compact fluxes -------------------------------------------------
CompactFlux_ADVECTION           = f_CompactFluxTimeSeries_11182019(build_ADVECTION);
CompactFlux_HORIZONTALMIXING	= f_CompactFluxTimeSeries_11182019(build_HORIZONTALMIXING);
CompactFlux_VERTICALMIXING      = f_CompactFluxTimeSeries_11182019(build_VERTICALMIXING);
CompactFlux_SINKING             = f_CompactFluxTimeSeries_11182019(build_SINKING); % compact SINKING as box floor areas and connectivity information; apply functional group sinking speeds in ECOTRANdynamic_ code
% FFF add migration conncectivity here?
% -------------------------------------------------------------------------


% step 4d: interpolate ROMS_time to t_grid --------------------------------
% temp_ROMS_time                           	= (t_grid(1):1:(ROMS_time(end) - ROMS_time(1))/86400 + 1)'; % time starting from t_grid(1); (vertical vector: num_t_ROMS X 1); NOTE: 86400 = seconds in 1 day

old_ADVECTION_flux                       	= CompactFlux_ADVECTION.compact_flux;
% % CompactFlux_ADVECTION.compact_flux          = interp1(temp_ROMS_time, old_ADVECTION_flux, t_grid, 'pchip'); % (2D matrix: num_t X CompactFlux_ADVECTION.num_fluxes)
CompactFlux_ADVECTION.compact_flux          = interp1(temp_ROMS_time, old_ADVECTION_flux, t_grid_real, 'pchip'); % (2D matrix: num_t X CompactFlux_ADVECTION.num_fluxes)

old_HORIZONTALMIXING_flux                 	= CompactFlux_HORIZONTALMIXING.compact_flux;
CompactFlux_HORIZONTALMIXING.compact_flux	= interp1(temp_ROMS_time, old_HORIZONTALMIXING_flux, t_grid_real, 'pchip'); % (2D matrix: num_t X CompactFlux_HORIZONTALMIXING.num_fluxes)

old_VERTICALMIXING_flux                     = CompactFlux_VERTICALMIXING.compact_flux;
CompactFlux_VERTICALMIXING.compact_flux     = interp1(temp_ROMS_time, old_VERTICALMIXING_flux, t_grid_real, 'pchip'); % (2D matrix: num_t X CompactFlux_VERTICALMIXING.num_fluxes)

old_SINKING_flux                            = CompactFlux_SINKING.compact_flux;
CompactFlux_SINKING.compact_flux            = interp1(temp_ROMS_time, old_SINKING_flux, t_grid_real, 'pchip'); % (2D matrix: num_t X CompactFlux_SINKING.num_fluxes)

old_BoxVolume                               = BoxVolume;
BoxVolume                                   = interp1(temp_ROMS_time, old_BoxVolume, t_grid_real, 'pchip'); % (2D matrix: num_t X num_boxes)

% BGC terms -------
build_temperature(isnan(build_temperature))             = -99999; % convert NaNs, interp1 doesn't like NaNs
build_diatom(isnan(build_diatom))                       = -99999; % convert NaNs, interp1 doesn't like NaNs
build_nanophytoplankton(isnan(build_nanophytoplankton))	= -99999; % convert NaNs, interp1 doesn't like NaNs

ROMS_temperature                         	= interp1(temp_ROMS_time, build_temperature, t_grid_real, 'pchip');         % (2D matrix: num_t X num_boxes)
ROMS_diatom                                 = interp1(temp_ROMS_time, build_diatom, t_grid_real, 'pchip');              % (2D matrix: num_t X num_boxes)
ROMS_nanophytoplankton                      = interp1(temp_ROMS_time, build_nanophytoplankton, t_grid_real, 'pchip');	% (2D matrix: num_t X num_boxes)

ROMS_temperature(ROMS_temperature == -99999)                = NaN; % convert null entries back to NaN
ROMS_diatom(ROMS_diatom == -99999)                          = NaN; % convert null entries back to NaN
ROMS_nanophytoplankton(ROMS_nanophytoplankton == -99999)	= NaN; % convert null entries back to NaN
% -------------------------------------------------------------------------


% step 4e: initial values for BGC driver(s) -------------------------------
%          NOTE: Initial values are taken as mean across the entire time-series
%                Other options are possible such as average across upwelling season(s)
ROMS_temperature_initial        = mean(ROMS_temperature, 1); % (horizontal vector: 1 X num_boxes)
ROMS_diatom_initial             = mean(ROMS_diatom, 1); % (horizontal vector: 1 X num_boxes)
ROMS_nanophytoplankton_initial	= mean(ROMS_nanophytoplankton, 1); % (horizontal vector: 1 X num_boxes)
% -------------------------------------------------------------------------


% % step 4f: visualize error ------------------------------------------------
% figure; hold on
% plot(1:5124, squeeze(build_error_record(:, 1, :)), 'b:')
% plot(1:5124, squeeze(build_error_record(:, 2, :)), 'r')
% *************************************************************************





% *************************************************************************
% STEP 6: light parameters-------------------------------------------------
LightParams.latitude        = latitude;   % latitude of Deep-Water Horizon oil spill (decimal degress north)
LightParams.SolarConstant	= 1366;       % Solar const. (W m^-2) recent ave. (Froehlich & Lean (1998, Geographical Research Letters 25(23): 4377-4380))

LightParams.par_frac        = 0.43;       % fraction of 40 that is PAR, clear atmos. (Fasham, Ducklow, & McKelvie (1990, J. Mar. Res. 48:591-639), p. 604)
LightParams.surf_trans      = 0.85;       % 1-reflectance (Parsons, Takahachi, & Hargrave (1984, Biological Oceanographic Processes, 3rd ed.), p. 69)
LightParams.cloud_cover     = 0.50;       % fitting daily average and weekly-smoothed light model curve to weekly average of HMSC roof data by eye
% LightParams.cloud_cover     = 0.70;       % Ave. 1990-97, COADS Standard 1-deg. for 125.5W, 45.5N
LightParams.cloud_trans     = 1.0 - (0.7 * LightParams.cloud_cover);  % (Evans and Parslow (1985, Biol. Oceanogr. 3(3):327-347), p. 344)
LightParams.t_grid          = t_grid;

SurfaceRadiation            = f_LightIntensity_12112020(LightParams);	% instantaneous (W/m2), daily mean averaged across 24 h (W m^-2 h^-1), & daily integrated (W m^-2 d^-1) solar raditation at ocean surface
daily_average_light         = SurfaceRadiation.daily_average_light;     % daily integrated solar raditation at ocean surface averaged across 24 hours (W m^-2 h^-1); (vertical vector: num_t X 1)
%                                                 NOTE: the daily average is also what was used by Spitz in her NPZD models
%                                                 NOTE: the daily integrated value avergaed over 24 hours is almost identical to the mean of the current instantaneous light intensity averaged over 24 hours (I checked)
daily_integrated_light      = SurfaceRadiation.daily_integrated_light;	% time-series of daily integrated solar raditation at ocean surface (W m^-2 d^-1)<--QQQ check units, integrated value should not be a rate but a TOTAL-->(W m^-2)

current_light               = SurfaceRadiation.current_light;           % time-series of surface solar raditation at current time (W/m2); NOTE: current time = midnight when dt = 1 day
sunrise                     = SurfaceRadiation.sunrise;                  % time of sunrise (time in h from midnight; 12 = noon, set as a default)
sunset                      = SurfaceRadiation.sunset;                   % time of sunset  (time in h from midnight; 12 = noon, set as a default)

Kw                          = 0.067;    % Light attenuation_seawater (Newberger et al., 2003); (1/m)
Kp                          = 0.0095;   % Light attenuation_phytoplankton (Newberger et al., 2003); (m2/mmol N)

% NOTE: Mixed Layer Depth (only used for light intensity calculations, not used for advection rate calculations)
maxMLD                      = 15;                                                       % max MLD; (m) % alt 40
minMLD                      = 15;                                                       % min MLD; (m) % alt 10
rangeMLD                    = maxMLD - minMLD;                                          % set minimum mixed layer depth and annual range of the mixed layer depth; (m)
MLD                         = (cos(((t_grid/365)*pi)*2) * (rangeMLD * 0.5) + (minMLD + rangeMLD * 0.5)); % MLD; (m); (vertical vector: time X 1)
% *************************************************************************





% *************************************************************************
% STEP 8: pack results for export------------------------------------------
% step 8a: save names of this m-file and sub-functions --------------------
ECOTRANphysics.fname_PhysicalModel              = fname_PhysicalModel;

ECOTRANphysics.fname_ROMS_GridPrep              = ROMSgrid.fname_ROMS_GridPrep;
ECOTRANphysics.fname_ROMS_FluxPrep              = ROMSflux.fname_ROMS_FluxPrep;

ECOTRANphysics.fname_CompactFlux                = ROMSflux.fname_CompactFlux;           % file name of f_CompactFluxTimeSeries function
ECOTRANphysics.fname_UnCompactFluxTimeSeries	= ROMSflux.fname_UnCompactFluxTimeSeries; % name of this f_UnCompactFluxTimeSeries sub-function
ECOTRANphysics.fname_CalcNetFlux                = ROMSflux.fname_CalcNetFlux;           % name of this f_CalcNetFlux function
ECOTRANphysics.fname_EvaluateFluxBalance        = ROMSflux.fname_EvaluateFluxBalance;	% name of this f_EvaluateFluxBalance sub-function

ECOTRANphysics.fname_LightIntensity             = SurfaceRadiation.fname_LightIntensity;
% -------------------------------------------------------------------------


% step 8b: time dimensions ------------------------------------------------
ECOTRANphysics.datestart                        = datestart;	% starting date
ECOTRANphysics.dateend                          = dateend;      % ending date
ECOTRANphysics.t_grid                           = t_grid;       % (vertical vector: num_t X 1)
ECOTRANphysics.num_t                            = num_t;
ECOTRANphysics.dt                               = dt;           % time-step (d)

ECOTRANphysics.ROMS_time                        = ROMS_time;    % seconds since 1900-01-01 00:00:00; (vertical vector: 1 X days)
ECOTRANphysics.num_t_ROMS                       = num_t_ROMS;
% -------------------------------------------------------------------------


% step 8c: space dimensions -----------------------------------------------
ECOTRANphysics.num_boxes                        = num_boxes;
ECOTRANphysics.num_domains                      = num_domains;	% number of geographic domains (does not consider depth layers)
ECOTRANphysics.num_z                            = num_z;        % number of depth layers
ECOTRANphysics.BoxVolume                        = BoxVolume;    % box volume; (m3); (2D matrix: num_t_ROMS X num_boxes)
ECOTRANphysics.BoxFloorArea                     = BoxFloorArea;	% box floor area; (m2); (horizontal vector: 1 X num_boxes (147))
% QQQ other box dimensions are available from f_ROMS_FluxPrep
% -------------------------------------------------------------------------


% step 8d: spatial relationships ------------------------------------------
ECOTRANphysics.VerticalConnectivity             = VerticalConnectivity;     % connectivity between source & destiny boxes; 0 or 1; (3D matrix: 1 X SOURCE-->(num_boxes+1) (148) X DESTINY-->(num_boxes+1) (148))
ECOTRANphysics.VerticalConnectivity_stack       = VerticalConnectivity_stack; % (2D matrix: num_domains X num_agg_z); tells which boxes lay directly under each surface box; clm 1 is suface, clm 2 is depth layer 2, final column is the bottom
ECOTRANphysics.looky_BottomBoxes                = looky_BottomBoxes;        % addresses of bottom (benthic) boxes
% -------------------------------------------------------------------------


% step 8e: physical driver time-series ------------------------------------
ECOTRANphysics.ADVECTION                    = build_ADVECTION;          % advection rate                (m3/d); (3D matrix: num_t X SOURCE (num_boxes + boundary) X DESTINY (num_boxes + boundary))
ECOTRANphysics.HORIZONTALMIXING             = build_HORIZONTALMIXING;	% horizontal mixing rate        (m3/d); (3D matrix: num_t X SOURCE (num_boxes + boundary) X DESTINY (num_boxes + boundary))
ECOTRANphysics.VERTICALMIXING               = build_VERTICALMIXING;     % vertical mixing rate          (m3/d); (3D matrix: num_t X SOURCE (num_boxes + boundary) X DESTINY (num_boxes + boundary))
ECOTRANphysics.SINKING                   	= build_SINKING;            % sinking source box floor area (m2);   (3D matrix: num_t_ROMS X SOURCE-->(num_boxes+1) X DESTINY-->(num_boxes+1))

ECOTRANphysics.CompactFlux_ADVECTION        = CompactFlux_ADVECTION;
ECOTRANphysics.CompactFlux_HORIZONTALMIXING	= CompactFlux_HORIZONTALMIXING;
ECOTRANphysics.CompactFlux_VERTICALMIXING	= CompactFlux_VERTICALMIXING;
ECOTRANphysics.CompactFlux_SINKING          = CompactFlux_SINKING;      % compact SINKING as box floor areas and connectivity information; apply functional group sinking speeds in ECOTRANdynamic_ code

ECOTRANphysics.error_record                 = build_error_record;       % total flux error & horizontal flux for comparison; (m3/s); (4D matrix: DESTINY-->num_domains+1 X 2 [total horizontal flux; flux error] X num_t_ROMS (366) X length(looky_error_depth))
% -------------------------------------------------------------------------


% step 8f: light model info -----------------------------------------------
ECOTRANphysics.Io                           = daily_average_light;    	% time-series of surface PAR light intensity; daily mean solar raditation at ocean surface averaged across 24 hours; (W m^-2 h^-1) (Brock, 1981, R. Ozretich unpub. data); (vertical vector: num_t X 1)
ECOTRANphysics.daily_integrated_light       = daily_integrated_light;	% time-series of daily integrated solar raditation at ocean surface (W m^-2 d^-1)<--QQQ check units, integrated value should not be a rate but a TOTAL-->(W m^-2); (vertical vector: num_t X 1)
ECOTRANphysics.current_light                = current_light;            % time-series of surface solar raditation at current time (W/m2); NOTE: current time = midnight when dt = 1 day; (vertical vector: num_t X 1)
ECOTRANphysics.sunrise                      = sunrise;
ECOTRANphysics.sunset                       = sunset;
ECOTRANphysics.Kw                           = Kw;                       % Light attenuation of seawater; (scalar)
ECOTRANphysics.Kp                           = Kp;                     	% Light attenuation of phytoplankton (Newberger et al., 2003); (m2/mmol N); (scalar)
ECOTRANphysics.MLD                          = MLD;                   	% time-series of Mixed Layer Depth (for light intensity calculations, NOT advection calcs); (m); (vertical vector: num_t X 1)
% -------------------------------------------------------------------------


% step 8h: conversion factors ---------------------------------------------
ECOTRANphysics.atomic_mass_C                = atomic_mass_C;             % (g C / mole C)
ECOTRANphysics.atomic_mass_N                = atomic_mass_N;             % (g N / mole N)
ECOTRANphysics.C_to_N_phytoplankton         = C_to_N_phytoplankton;      % Redfield ratio = C:N:P = 106:16:1 ---> C:N = 106:16 = 6.625 mole C/mole N
ECOTRANphysics.C_to_N_zooplankton           = C_to_N_zooplankton;        % general zooplankton C:N (moles C moles N^-1) for mesozooplankton (W. Peterson pers com); C:N=5.16-5.26 [Schneider 1990 (Mar Biol 106:219-225)]; C:N=3.9 [Uye 1982]
ECOTRANphysics.Chla_to_N_diatoms            = Chla_to_N_diatoms;         % (mg Chla mmoles N^-1) (Dickson & Wheeler 1995 L&O 40(3):533-543)
ECOTRANphysics.Chla_to_N_dinoflagellates    = Chla_to_N_dinoflagellates; % 0.73 to 0.95 (mg Chla mmoles N^-1) (Chan 1980 J Phycol 16:428-432)
ECOTRANphysics.mgWWT_per_mmolesN            = mgWWT_per_mmolesN;         % (mg WWT/mmoles N); for phytoplankton
ECOTRANphysics.WWT_to_C                     = WWT_to_C;                  % (g WWT / g C); (Steele et al. 2007); QQQ recheck this for phytoplankton


% step 8i: temperature & BGC info -----------------------------------------
ECOTRANphysics.temperature_timeseries           = ROMS_temperature;         % (deg C); (2D matrix: num_t X num_boxes)
ECOTRANphysics.temperature_dates                = t_grid_real;
ECOTRANphysics.temperature_reference            = repmat(11.02, [1, num_boxes]);

ECOTRANphysics.diatom_timeseries                = ROMS_diatom;              % (mmole N/m3); (2D matrix: num_t X num_boxes)
ECOTRANphysics.nanophytoplankton_timeseries     = ROMS_nanophytoplankton;	% (mmole N/m3); (2D matrix: num_t X num_boxes)

ECOTRANphysics.ROMS_temperature_initial         = ROMS_temperature_initial; % (deg C); (horizontal vector: 1 X num_boxes)
ECOTRANphysics.ROMS_diatom_initial              = ROMS_diatom_initial; % (mmoles N/m3); (horizontal vector: 1 X num_boxes)
ECOTRANphysics.ROMS_nanophytoplankton_initial	= ROMS_nanophytoplankton_initial; % (mmoles N/m3); (horizontal vector: 1 X num_boxes)
% *************************************************************************


% end m-file***************************************************************