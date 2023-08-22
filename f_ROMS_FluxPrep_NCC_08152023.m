function ROMSflux = f_ROMS_FluxPrep_NCC_08152023(ROMSgrid, readFile_FluxYear)
% express ROMS fluxes in DESTINY<--SOURCE format on ECOTRAN grid
%
% takes:
%       ROMSgrid            from f_ROMS_GridPrep_NCC_11012021; structure variable defining how ROMS grid cells map into ECOTRAN grid and their connectivity with each other
%       readFile_FluxYear   ROMS flux & tracer data for one year
%
% returns:
%       ROMSflux structure variable defining fluxes & tracers in DESTINY-->SOURCE format
%
% calls:
%	f_CompactFluxTimeSeries_11182019      compact physical flux time-series, arranged as 3D matrix (time X source box X destiny box), eliminate source & destiny infomration for any subregions that never communicate
%	f_UnCompactFluxTimeSeries_12112019    UnCompact a previously compacted physical flux time-series to provide IMPORT & EXPORT fluxes for each box and the domain as a whole
%       f_calcNetFlux_12112019                calculate net flux into and net flux out of each model box and across outer domain boundaries
%	f_EvaluateFluxBalance_11262021        examine for flux time-series for imbalances IN & OUT of individual boxes and IN & OUT of the overall domain
%
% revision date: 8/15/2023
%       5/30/2022 deactiviated write-over readFile_FluxYear = wc12_avg_2008_trimmed_2.nc
%       8/5/2022 fixing flux error correction for cases when error is disttrbuted to boxes with no horizontal connectivity to any other box at time t
%       8/29/2022 fixed flux error correction
%       9/2/2022 processing biogeochemical model output
%       9/13/2022 processing temeperature (biogeochemical) model output
%       10/16/2022 processing BGC (biogeochemical) model output
%       11/29/2022 f_ROMS_FluxPrep_NCC_10162022_B is for debugging work
%       11/30/2022 corrected f_ROMS_FluxPrep_ for weighted averaging of aggregated BGC info
%       8/15/2023 cleaned up notation


% *************************************************************************
% STEP 1: set operating conditions-----------------------------------------
fname_ROMS_FluxPrep     = mfilename; % save name of this m-file to keep in saved model results
display(['Running: ' fname_ROMS_FluxPrep])
% -------------------------------------------------------------------------


% step 1b: unpack ROMS grid and ECOTRAN aggregation info ------------------
fname_ROMS_GridPrep             = ROMSgrid.fname_ROMS_GridPrep;

z_definition                    = ROMSgrid.z_definition;

mask_v                          = ROMSgrid.mask_v;
mask_u                          = ROMSgrid.mask_u;
mask_rho                        = ROMSgrid.mask_rho;

num_domains                     = ROMSgrid.num_domains;
num_domainsWboundary            = ROMSgrid.num_domainsWboundary;
num_agg_z                       = ROMSgrid.num_agg_z;
num_boxes                       = ROMSgrid.num_boxes;

z_w                             = ROMSgrid.z_w;  % z_w are the depths @ each horizontal depth-layer face at the w point (top & bottom face depths at the cell center); (m); (3D matrix: rho longitude (51) X rho latitude (81) X num_z+1 (43))
z_wv                            = ROMSgrid.z_wv; % z_wv is the depth of each horizontal depth-layer face at the v point
z_wu                            = ROMSgrid.z_wu; % z_wu is the depth of each horizontal depth-layer face at the u point

v_height                        = ROMSgrid.v_height; % height of cell @ v-point; (bottom minus top);   (m); (3D matrix: v longitude (51) X v latitude (80) X num_z (42)); NOTE: JR confirms sum of these heights = bathymetry @ v
u_height                        = ROMSgrid.u_height; % height of cell @ u-point; (bottom minus top);   (m); (3D matrix: u longitude (50) X u latitude (81) X num_z (42))
w_height                        = ROMSgrid.w_height; % height of cell @ center; (bottom minus top);    (m); (3D matrix: rho longitude (51) X rho latitude (81) X num_z (42))
h                               = ROMSgrid.h; % bathymetry at RHO-points; (2D matrix: rho longitude (51 west:east) X rho latitude (81 south:north))

dist_NS_globe                   = ROMSgrid.dist_NS_globe; % N<-->S lengths relative to globe; (m); abs(NORTH minus SOUTH); (2D matrix: u longitude (50) X v latitude-1 (79));
dist_EW_globe                   = ROMSgrid.dist_EW_globe; % W<-->E lengths relative to globe; (m); abs(EAST minus WEST); (2D matrix: u longitude-1 (49) X v latitude (80));
dist_NS_face                    = ROMSgrid.dist_NS_face; % N<-->S face lengths; (m); (2D matrix: u longitude (50) X v latitude-1 (79)); NOTE: x_psi(south) - x_psi(north)
dist_EW_face                    = ROMSgrid.dist_EW_face; % E<-->W face lengths; (m); (2D matrix: u longitude-1 (49) X v latitude (80)); NOTE: y_psi(east) - y_psi(west)

area_NS_globe                   = ROMSgrid.area_NS_globe;    % N<-->S ROMS lateral cell area based on EW & NS distances relative to the globe; (m2); (3D matrix: v longitude-2 (49) X v latitude (80) X num_z (42))
area_EW_globe                   = ROMSgrid.area_EW_globe;    % E<-->W ROMS lateral cell area based on EW & NS distances relative to the globe; (m2); (3D matrix: v longitude-2 (49) X v latitude (80) X num_z (42))
area_floor_globe                = ROMSgrid.area_floor_globe; % cell floor area based on EW & NS distances relative to the globe; (m2); (3D matrix: rho longitude-2 (49) X rho latitude-2 (79), num_w (43)); (includes only values at rho points the inside domain)
area_floor_face                 = ROMSgrid.area_floor_face;	 % cell floor area based on EW & NS distances relative to face lengths; (m2); (3D matrix: rho longitude-2 (49) X rho latitude-2 (79), num_w (43)); (includes only values at rho points the inside domain)

LatitudinalConnectivity         = ROMSgrid.LatitudinalConnectivity; % values = 0 or 1; neighbor = 1; (2D matrix: num_domainsWboundary X num_domainsWboundary)
LatitudeFluxConnectivity_W_v	= ROMSgrid.LatitudeFluxConnectivity_W_v; % western v-point for ALL overlapping ECOTRAN boxes; (2D matrix: num_domainsWboundary X num_domainsWboundary)
LatitudeFluxConnectivity_E_v	= ROMSgrid.LatitudeFluxConnectivity_E_v; % eastern v-point for ALL overlapping ECOTRAN boxes; (2D matrix: num_domainsWboundary X num_domainsWboundary)
grid_addresses_v                = ROMSgrid.grid_addresses_v; % (2D matrix: num_domains X 7); [south north west east (num cells lat) (num cells lon) (num cells total)]

LongitudinalConnectivity        = ROMSgrid.LongitudinalConnectivity; % values = 0 or 1; neighbor = 1; (2D matrix: num_domainsWboundary X num_domainsWboundary)
LongitudeFluxConnectivity_S_u	= ROMSgrid.LongitudeFluxConnectivity_S_u; % southern u-point for ALL overlapping ECOTRAN boxes; (2D matrix: num_domainsWboundary X num_domainsWboundary)
LongitudeFluxConnectivity_N_u	= ROMSgrid.LongitudeFluxConnectivity_N_u; % northern u-point for ALL overlapping ECOTRAN boxes; (2D matrix: num_domainsWboundary X num_domainsWboundary)
grid_addresses_u                = ROMSgrid.grid_addresses_u; % (2D matrix: num_domains X 7); [south north west east (num cells lat) (num cells lon) (num cells total)]

VerticalConnectivity            = ROMSgrid.VerticalConnectivity; % 0 or 1; (2D matrix: DESTINY-->(num_boxes+1) (148) X SOURCE-->(num_boxes+1) (148))

looky_BottomBoxes               = ROMSgrid.looky_BottomBoxes; % addresses of bottom (benthic) boxes

grid_addresses_rho            	= ROMSgrid.grid_addresses_rho; % (2D matrix: num_domains X 7); [south north west east (num cells lat) (num cells lon) (num cells total)]
% -------------------------------------------------------------------------


% step 1c: define constants -----------------------------------------------
south       = 1; % grid_addresses column definitions
north       = 2;
west        = 3;
east        = 4;
southwest	= 5;
southeast	= 6;
northwest	= 7;
northeast	= 8;
% *************************************************************************





% *************************************************************************
% STEP 2: load ROMS variables from NetCDF files----------------------------
% step 2a: load ROMS NCC flow variables -----------------------------------
ncid        = netcdf.open(readFile_FluxYear, 'NC_NOWRITE'); % Open readFile

% values at rho-points (cell centers)
varid       = netcdf.inqVarID(ncid, 'lat_rho'); % Get variable ID of the first variable, given its name
lat_rho  	= netcdf.getVar(ncid, varid);       % 'latitude of RHO-points'; (2D matrix: rho longitude (51 west:east) X rho latitude (81 south:north))
varid       = netcdf.inqVarID(ncid, 'lon_rho'); 
lon_rho   	= netcdf.getVar(ncid, varid);       % 'longitude of RHO-points'; (2D matrix: rho longitude (51 west:east) X rho latitude (81 south:north))

varid       = netcdf.inqVarID(ncid, 'w');
w           = netcdf.getVar(ncid, varid);       % time-averaged vertical momentum component; (m/s); (4D matrix: rho longitude (51 west:east) X rho latitude (81 south:north) X num_z+1 (43; s_w) X num_days)

varid       = netcdf.inqVarID(ncid, 'omega');
omega       = netcdf.getVar(ncid, varid);       % time-averaged S-coordinate vertical momentum component; (m/s); NOTE: comments say units are (m3/s), but this won't allow volume balance while (m/s) DOES WORK; (4D matrix: rho longitude (51 west:east) X rho latitude (81 south:north) X num_z+1 (43; s_w) X num_days)

varid       = netcdf.inqVarID(ncid, 's_rho');
s_rho       = netcdf.getVar(ncid, varid);       % 'S-coordinate at RHO-points'; (0 to -1; up is towards 0); (vertical vector: num_z (42; s_rho) X 1);
varid       = netcdf.inqVarID(ncid, 's_w');
s_w         = netcdf.getVar(ncid, varid);       % 'S-coordinate at W-points'; (0 to -1; up is "+"; up is towards 0); (vertical vector: num_z+1 (43; s_w) X 1);

varid       = netcdf.inqVarID(ncid, 'zeta');
zeta      	= netcdf.getVar(ncid, varid);       % free surface at RHO-points; (m); (3D matrix: rho longitude (51 west:east) X rho latitude (81 south:north) X num_t_ROMS)


% biogeochemical model variables
varid               = netcdf.inqVarID(ncid, 'temp');
temperature         = netcdf.getVar(ncid, varid);       % potential temperature; (Celsius); (4D matrix: rho longitude (51 west:east) X rho latitude (81 south:north) X num_z (42; s_rho) X num_days)

varid               = netcdf.inqVarID(ncid, 'nanophytoplankton');
nanophytoplankton	= netcdf.getVar(ncid, varid); % time-averaged nanophytoplankton biomass; (millimole N /m3); % (4D matrix: rho longitude (51 west:east) X rho latitude (81 south:north) X num_z (42; s_rho) X num_days)
varid               = netcdf.inqVarID(ncid, 'diatom');
diatom              = netcdf.getVar(ncid, varid); % time-averaged diatom biomass; (millimole N /m3); % (4D matrix: rho longitude (51 west:east) X rho latitude (81 south:north) X num_z (42; s_rho) X num_days)


% values at v-points (cell N/S boundaries)
varid       = netcdf.inqVarID(ncid, 'lat_v');
lat_v       = netcdf.getVar(ncid, varid);       % 'latitude of V-points'; (2D matrix: v longitude (51 west:east) X v latitude (80 south:north));
varid       = netcdf.inqVarID(ncid, 'lon_v');
lon_v       = netcdf.getVar(ncid, varid);       % 'longitude of V-points'; (2D matrix: v longitude (51 west:east) X v latitude (80 south:north));
varid       = netcdf.inqVarID(ncid, 'v');
v           = netcdf.getVar(ncid, varid);       % time-averaged v-momentum component; (m/s); (4D matrix: v longitude (51 west:east) X v latitude (80 south:north) X num_z (42; s_rho) X num_days)

% values at u-points (cell E/W boundaries)
varid       = netcdf.inqVarID(ncid, 'lat_u');
lat_u       = netcdf.getVar(ncid, varid);       % 'latitude of U-points'; (2D matrix: u longitude (52 west:east) X u latitude (81 south:north))
varid       = netcdf.inqVarID(ncid, 'lon_u');
lon_u       = netcdf.getVar(ncid, varid);       % 'longitude of U-points'; (2D matrix: u longitude (52 west:east) X u latitude (81 south:north))
varid       = netcdf.inqVarID(ncid, 'u');
u           = netcdf.getVar(ncid, varid);       % time-averaged u-momentum component; (m/s); (4D matrix: u longitude (52 west:east) X u latitude (81 south:north) X num_z (42; s_rho) X num_days)

% other variables
varid       = netcdf.inqVarID(ncid, 'hc');
hc          = netcdf.getVar(ncid, varid);       % 'S-coordinate parameter, critical depth'; (m); scalar

varid       = netcdf.inqVarID(ncid, 'ocean_time');
ocean_time	= netcdf.getVar(ncid, varid);       % seconds since 1900-01-01 00:00:00; (vertical vector: 366 X 1)

% Close the NetCDF file
netcdf.close(ncid)
% -------------------------------------------------------------------------


% step 2b: clear temporary variables --------------------------------------
clear ncid varid
% *************************************************************************





% *************************************************************************
% STEP 3: trim full grid ROMS terms to only those WITHIN the NCC ----------
% step 3a: round lats & lons (needed because of TINY rounding errors) -----
lat_rho             = round(lat_rho, 2);
lon_rho             = round(lon_rho, 2);
lat_v               = round(lat_v, 2);      % (2D matrix: v longitude (51 west:east) X v latitude (80 south:north));
lon_v               = round(lon_v, 2);
lat_u               = round(lat_u, 2);
lon_u               = round(lon_u, 2);
% -------------------------------------------------------------------------


% step 3b: make sure values align to Arakawa-C grid -----------------------
%          NOTE: apply to the NCC sub-setted ROMS flow variables file
% longitude range (outer E/W boundaries defined by rho & v)
WestBoundary        = min(lon_rho(:, 1)); % NOTE: assumes longitudes run down rows; NOTE: assumes location entirely in western hemisphere (negative longtitudes)
EastBoundary        = max(lon_rho(:, 1)); % NOTE: assumes longitudes run down rows; NOTE: assumes location entirely in western hemisphere (negative longtitudes)

% trim u terms to only those WITHIN the rho longitude range
looky_TooWest	= find(lon_u(:, 1) < WestBoundary); % NOTE: assumes longitudes run down rows; probably an empty matrix
if ~isempty(looky_TooWest)
    lon_u(looky_TooWest, :)     = []; % delete rows west of Arakawa-C longitude domain
    lat_u(looky_TooWest, :)     = []; % delete rows west of Arakawa-C longitude domain
    u(looky_TooWest, :, :, :)	= []; % delete rows west of Arakawa-C longitude domain; (m/s)
end
    
looky_TooEast	= find(lon_u(:, 1) > EastBoundary); % NOTE: assumes longitudes run down rows; probably an empty matrix
if ~isempty(looky_TooEast)
    lon_u(looky_TooEast, :)     = []; % delete rows east of Arakawa-C longitude domain
    lat_u(looky_TooEast, :)     = []; % delete rows east of Arakawa-C longitude domain
    u(looky_TooEast, :, :, :)	= []; % delete rows east of Arakawa-C longitude domain; (m/s)
end
% -----

% latitude range (outer N/S boundaries defined by rho & u)
NorthBoundary       = max(lat_rho(1, :)); % NOTE: assumes latitudes run across columns; NOTE: assumes location entirely in northern hemisphere
SouthBoundary       = min(lat_rho(1, :)); % NOTE: assumes latitudes run across columns; NOTE: assumes location entirely in northern hemisphere

% trim v terms to only those WITHIN the rho latitude range
looky_TooNorth = find(lat_v(1, :) >= NorthBoundary); % NOTE: assumes latitudes run across columns; probably an empty matrix
if ~isempty(looky_TooNorth)
    lon_v(:, looky_TooNorth)	= []; % delete rows north of Arakawa-C latitude domain
    lat_v(:, looky_TooNorth)	= []; % delete rows north of Arakawa-C latitude domain
    v(:, looky_TooNorth, :, :)	= []; % delete rows north of Arakawa-C latitude domain; (m/s)
end
    
looky_TooSouth = find(lat_v(1, :) <= SouthBoundary); % NOTE: assumes latitudes run across columns; probably an empty matrix
if ~isempty(looky_TooSouth)
    lon_v(:, looky_TooSouth)	= []; % delete rows south of Arakawa-C latitude domain
    lat_v(:, looky_TooSouth)	= []; % delete rows south of Arakawa-C latitude domain
    v(:, looky_TooSouth, :, :)	= []; % delete rows south of Arakawa-C latitude domain; (m/s)
end

% QQQ could add more error checking to compare lat & long of rho, u, & v
% -------------------------------------------------------------------------


% step 3c: get final matrix sizes -----------------------------------------
[rows_rho, clms_rho]	= size(lat_rho);	% rows = w longitudes (rho=v, 51); clms = w latitudes (rho=u, 81)
[rows_h, clms_h]        = size(h);          % rows = longitudes (h=rho); clms = latitudes (h=rho)
num_z                   = length(s_rho);
num_w                   = length(s_w);      % num_w = num_z + 1;
[rows_v, clms_v]        = size(lat_v);      % rows = v longitudes (v=rho, 51); clms = v latitudes (v=rho-1, u-1, 80)
[rows_u, clms_u]        = size(lat_u);      % rows = u longitudes (u=rho+1, v+1, 50); clms = latitudes (u=rho, 81)
num_t_ROMS             	= length(ocean_time);

[rows_temperature, clms_temperature, ~, ~]	= size(temperature);	% rows = temperature longitudes (rho=v, 51); clms = temperature latitudes (rho=u, 81)
[rows_BGC, clms_BGC, ~, ~]                  = size(diatom);         % rows = diatom longitudes (rho=v, 51); clms = diatom latitudes (rho=u, 81)

if (rows_rho ~= rows_h) || (clms_rho ~= clms_h)
    error('rho & h matrices are different sizes')
end

if (ROMSgrid.rows_v ~= rows_v) || (ROMSgrid.clms_v ~= clms_v)
	error('v matrix is different size from ROMS_GridPrep function')
end
  
if (ROMSgrid.rows_u ~= rows_u) || (ROMSgrid.clms_u ~= clms_u)
	error('u matrix is different size from ROMS_GridPrep function')
end
    
if (ROMSgrid.rows_rho ~= rows_rho) || (ROMSgrid.clms_rho ~= clms_rho)
	error('rho matrix is different size from ROMS_GridPrep function')
end

if (ROMSgrid.rows_h ~= rows_h) || (ROMSgrid.clms_h ~= clms_h)
	error('h matrix is different size from ROMS_GridPrep function')
end

if ROMSgrid.num_z ~= num_z
	error('z vector is different size from ROMS_GridPrep function')
end

if ROMSgrid.num_w ~= num_w
	error('w vector is different size from ROMS_GridPrep function')
end

if (rows_rho ~= rows_temperature) || (clms_rho ~= clms_temperature)
    error('rho & temperature matrices are different sizes')
end

if (rows_rho ~= rows_BGC) || (clms_rho ~= clms_BGC)
    error('rho & BGC matrices are different sizes')
end
% -------------------------------------------------------------------------


% step 3d: apply masks ----------------------------------------------------
v       = v     .* repmat(mask_v, [1, 1, num_z, num_t_ROMS]);   % (m/s); (4D matrix: v longitude (51) X v latitude (80) X num_z (42) X num_t_ROMS (366))
u       = u     .* repmat(mask_u, [1, 1, num_z, num_t_ROMS]);   % (m/s); (4D matrix: u longitude (50) X u latitude (81) X num_z (42) X num_t_ROMS (366))
w       = w     .* repmat(mask_rho, [1, 1, num_w, num_t_ROMS]); % (m/s); (4D matrix: rho longitude (51) X rho latitude (81) X num_w (43) X num_t_ROMS (366))
omega	= omega .* repmat(mask_rho, [1, 1, num_w, num_t_ROMS]); % (m/s); (4D matrix: rho longitude (51) X rho latitude (81) X num_w (43) X num_t_ROMS (366))
zeta	= zeta  .* repmat(mask_rho, [1, 1, num_t_ROMS]);        % (m);   (3D matrix: rho longitude (51) X rho latitude (81) X num_t_ROMS (366))

% biological terms
temperature         = temperature       .* repmat(mask_rho, [1, 1, num_z, num_t_ROMS]);   % (m/s); (4D matrix: rho longitude (51) X rho latitude (81) X num_z (42) X num_t_ROMS (366))
diatom              = diatom            .* repmat(mask_rho, [1, 1, num_z, num_t_ROMS]);   % (m/s); (4D matrix: rho longitude (51) X rho latitude (81) X num_z (42) X num_t_ROMS (366))
nanophytoplankton	= nanophytoplankton	.* repmat(mask_rho, [1, 1, num_z, num_t_ROMS]);   % (m/s); (4D matrix: rho longitude (51) X rho latitude (81) X num_z (42) X num_t_ROMS (366))
% -------------------------------------------------------------------------


% step 3e: clear temporary variables --------------------------------------
clear looky_* mask_u mask_v
% *************************************************************************





% *************************************************************************
% STEP 4: calculate ROMS cell volume fluxes--------------------------------
% step 4a: calculate horizontal volume fluxes through cell faces ----------
%          NOTE: rectangular face
repmat_area_EW_v        = repmat(area_EW_globe, [1, 1, 1, num_t_ROMS]);      % E<-->W ROMS cell area; (m2); (4D matrix: v longitude-2 (49) X v latitude (80) X num_z (42) X num_t_ROMS (366))
flux_NS                 = v(2:(end-1), :, :, :) .* repmat_area_EW_v;	% N<-->S volume flux;  (m3/s); (4D matrix: v longitude-2 (49) X v latitude (80) X num_z (42) X num_t_ROMS (366)); NOTE: excludes far west & far east v terms that are outside of complete ROMS cells

repmat_area_NS_u        = repmat(area_NS_globe, [1, 1, 1, num_t_ROMS]);	    % N<-->S ROMS cell area; (m2); (4D matrix: u longitude (50) X u latitude-2 (79) X num_z (42) X num_t_ROMS (366))
flux_EW                 = u(:, 2:(end-1), :, :) .* repmat_area_NS_u;	% E<-->W volume flux;  (m3/s); (4D matrix: u longitude (50) X u latitude-2 (79) X num_z (42) X num_t_ROMS (366)); NOTE: excludes far north & far south u terms that are outside of complete ROMS cells
% -------------------------------------------------------------------------


% step 4b: calculate vertical volume fluxes through cell faces ------------
% %          NOTE: multiply w across trapezoidal face
repmat_area_floor_globe	= repmat(area_floor_globe, [1, 1, 1, num_t_ROMS]); % horizontal cell face area; (m2); (4D matrix: psi longitude-1 (49) X psi latitude-1 (79), num_w (43) X num_t_ROMS (366));
flux_UpDown             = omega(2:(end-1), 2:(end-1), :, :) .* repmat_area_floor_globe; % Up<-->Down volume flux; (m3/s); (4D matrix: psi longitude-1 (49) X psi latitude-1 (79), num_w (43) X num_t_ROMS (366));
% -------------------------------------------------------------------------


% step 4c: volume test -- calculate net flux into each ROMS cell ----------
%          NOTE: northward flux is "+"
%          NOTE: eastward flux is "+"
%          NOTE: upward flux is "+"
%          NOTE: face order matters; 
net_flux_NS_1           = flux_NS(:, 1:(end-1), :, :)     - flux_NS(:, 2:end, :, :);     % net N<-->S volume flux; (flux through SOUTH face minus flux through NORTH face); (m3/s); (4D matrix: psi longitude-1 (49) X psi latitude-1 (79) X num_z (42) X num_t_ROMS (366))
net_flux_EW_1           = flux_EW(1:(end-1), :, :, :)     - flux_EW(2:end, :, :, :);     % net E<-->W volume flux; (flux through WEST face minus flux through EAST face); (m3/s); (4D matrix: psi longitude-1 (49) X psi latitude-1 (79) X num_z (42) X num_t_ROMS (366))
net_flux_UpDown_1       = flux_UpDown(:, :, 1:(end-1), :) - flux_UpDown(:, :, 2:end, :); % net Up<-->Down volume flux; (flux through BOTTOM face minus flux through TOP face); (m3/s); (4D matrix: psi longitude-1 (49) X psi latitude-1 (79) X num_z (42) X num_t_ROMS (366))

net_flux_horizontal_1	= net_flux_NS_1 + net_flux_EW_1; % net HORIZONTAL volume flux into each ROMS cell; (m3/s); (4D matrix: psi longitude-1 (49) X psi latitude-1 (79) X num_z (42) X num_t_ROMS (366))
net_flux_horizontal_1	= sum(net_flux_horizontal_1, 3); % net HORIZONTAL volume flux into a column of water; (m3/s); (4D matrix: psi longitude-1 (49) X psi latitude-1 (79) X num_z (42) X num_t_ROMS (366))
            % NOTE: at this point, horizontal conservation is good. Most error is well within 50 m3/s. There are 6 instances greater than 200m3/s. Max error is 525m3/s.

net_flux_1              = net_flux_NS_1 + net_flux_EW_1 + net_flux_UpDown_1; % total net flux; should = 0; (m3/s); (4D matrix: psi longitude-1 (49) X psi latitude-1 (79) X num_z (42) X num_t_ROMS (366))
% hist(net_flux_1(:), 1000)
            % NOTE: at this point, 3D volume conservation is good. Most error is well within 20 m3/s. There are a few (476 out of 59,505,012) instances greater than 200m3/s. Max error is 506 m3/s.
% *************************************************************************





% *************************************************************************
% STEP 5: aggregate ROMS depth layers into ECOTRAN depth layers------------
%         NOTE: after this aggregation step, vertical flux matrices are re-ordered from surface to bottom

% % step 5a: calculate ROMS depth coordinates -------------------------------
% %          NOTE: these are depths to u, v, & w points, NOT cell edge depths (except for w points)
% % define ROMS parameters (from M. Jacox 3/18/2021)
% Vtransform      = 1; % 1 = original transformation
% Vstretching     = 1; % 1 = original (Song and Haidvogel, 1994)
% theta_s         = 5; % S-coordinate surface control parameter (scalar)
% theta_b         = 0.4; % S-coordinate bottom control parameter (scalar)
% report          = 0; % flag to report detailed information (OPTIONAL):
% N               = num_z; % number of layers; (N = 0 is bottom, N = 42 is the top)
% %	igrid: Staggered grid C-type (integer):
% %       igrid=1  => density points
% %       igrid=2  => streamfunction points
% %       igrid=3  => u-velocity points
% %       igrid=4  => v-velocity points
% %       igrid=5  => w-velocity points
% 
% z_w             = zeros(rows_rho, clms_rho, (num_z+1), num_t_ROMS);  % initialize; (4D matrix: rows_rho X clms_rho X (num_z+1) X num_t_ROMS);
% z_u             = zeros(rows_u, clms_u, num_z, num_t_ROMS);          % initialize; (4D matrix: rows_u   X clms_u   X num_z     X num_t_ROMS);
% z_v             = zeros(rows_v, clms_v, num_z, num_t_ROMS);          % initialize; (4D matrix: rows_v   X clms_v   X num_z     X num_t_ROMS);
% 
% for day_loop = 1:num_t_ROMS
%     % z at w-points
%     z_w(:, :, :, day_loop)	= set_depth(Vtransform, Vstretching, ...
%                               theta_s, theta_b, hc, N, ...
%                               5, h, zeta(:, :, day_loop), report); % depth at w-points; (m); (4D matrix: rho longitude (51) X rho latitude (81) X num_z+1 (43) X num_t_ROMS);
%                                                                    %    NOTE: z_w top layer = h
%                                                                    %    NOTE: z_w bottom layer = zeta 
%                                                                    
%     % z at u-points
%     z_u(:, :, :, day_loop)	= set_depth(Vtransform, Vstretching, ...
%                               theta_s, theta_b, hc, N, ...
%                               3, h, zeta(:, :, day_loop), report); % depth at u-points; (m); (4D matrix: u longitude (50) X u latitude (81) X num_z (42) X num_t_ROMS);
%                           
%     % z at v-points
%     z_v(:, :, :, day_loop)	= set_depth(Vtransform, Vstretching, ...
%                               theta_s, theta_b, hc, N, ...
%                               4, h, zeta(:, :, day_loop), report); % depth at v-points; (m); (4D matrix: v longitude (51) X v latitude (80) X num_z (42) X num_t_ROMS);
% end % day_loop
% % -------------------------------------------------------------------------


% step 5b: aggregate u, v, & w into ECOTRAN depth zones -------------------
z_wv                            = z_wv * (-1);  % set to negative; (m); (4D matrix: v longitude (51) X v latitude (80) X num_w (43) X num_t_ROMS (366))
z_wu                            = z_wu * (-1);  % set to negative; (m); (4D matrix: u longitude (50) X u latitude (81) X num_w (43) X num_t_ROMS (366))
z_w                             = z_w  * (-1);  % set to negative; (m); (4D matrix: rho longitude (51) X rho latitude (81) X num_w (43) X num_t_ROMS (366))

agg_flux_NS                     = zeros(rows_v-2, clms_v, num_agg_z, num_t_ROMS);         	% initialize; (m3/s); (4D matrix: v longitude-2 (49) X v latitude (80) X num_agg_z (7) X num_t_ROMS (366))
agg_flux_EW                     = zeros(rows_u, clms_u-2, num_agg_z, num_t_ROMS);        	% initialize; (m3/s); (4D matrix: u longitude (50) X u latitude-2 (79) X num_agg_z (7) X num_t_ROMS (366))
agg_flux_UpDown                 = zeros(rows_rho-2, clms_rho-2, (num_agg_z+1), num_t_ROMS);	% initialize; (m3/s); (4D matrix: rho longitude-2 (49) X rho latitude-2 (79) X num_agg_z+1 (8) X num_t_ROMS (366))

agg_v_height                    = zeros(rows_v, clms_v, num_agg_z);                  % initialize; height of cell @ v-point; (m); (3D matrix: v longitude (51) X v latitude (80) X num_agg_z (7))
agg_u_height                    = zeros(rows_u, clms_u, num_agg_z);                  % initialize; height of cell @ u-point; (m); (3D matrix: u longitude (50) X u latitude (81) X num_agg_z (7))

agg_temperature_vertical        = zeros(rows_rho, clms_rho, num_agg_z, num_t_ROMS); % initialize; (4D matrix: rho longitude (51) X rho latitude (81) X num_agg_z X num_t_ROMS)
agg_diatom_vertical             = zeros(rows_rho, clms_rho, num_agg_z, num_t_ROMS); % initialize; (4D matrix: rho longitude (51) X rho latitude (81) X num_agg_z X num_t_ROMS)
agg_nanophytoplankton_vertical	= zeros(rows_rho, clms_rho, num_agg_z, num_t_ROMS); % initialize; (4D matrix: rho longitude (51) X rho latitude (81) X num_agg_z X num_t_ROMS)

for depth_loop = 1:num_agg_z

	current_shallow                     = z_definition(depth_loop);     % (m); (top of ROMS cell)
    current_deep                        = z_definition(depth_loop + 1); % (m); (bottom of ROMS cell)

    % if current_shallow is the surface, set the value to the free surface (zeta)
    if current_shallow == 0
        current_shallow_v = z_wv(:, :, end, :);         % free surface (zeta); (m); (2D matrix: v longitude (51) X v latitude (80))
        current_shallow_u = z_wu(:, :, end, :);         % free surface (zeta); (m); (2D matrix: u longitude (50) X u latitude (81))
        current_shallow_w = z_w(:, :, end, :);          % free surface (zeta); (m); (2D matrix: rho longitude (51) X rho latitude (81))
    else
        current_shallow_v = z_definition(depth_loop);	% top of ROMS cell; (m); (scalar)
        current_shallow_u = z_definition(depth_loop);	% top of ROMS cell; (m); (scalar)
        current_shallow_w = z_definition(depth_loop);	% top of ROMS cell; (m); (scalar)
    end
    
    % calculate v aggregation ---------------------------------------------
    z_wv_distBELOWshallow               = current_shallow_v - z_wv;	% distance BELOW shallow boundary; (m); z_wv BELOW current_shallow are positive; (4D matrix: v longitude (51) X v latitude (80) X num_w (43) X num_t_ROMS (366))
    z_wv_distABOVEdeep                  = z_wv - current_deep;      % distance ABOVE deep boundary;    (m); z_wv ABOVE current_deep are positive;    (4D matrix: v longitude (51) X v latitude (80) X num_w (43) X num_t_ROMS (366))

    z_wv_distBELOWshallow(z_wv_distBELOWshallow < 0)	= 0;
    z_wv_distABOVEdeep(z_wv_distABOVEdeep < 0)          = 0;
    
    z_wv_fractionBELOWshallow           = z_wv_distBELOWshallow(:, :, 1:(end-1), :) ./ v_height; % fraction of cell BELOW shallow boundary; aligned with z_v; (4D matrix: v longitude (51) X v latitude (80) X num_z (42) X num_t_ROMS (366))
    z_wv_fractionABOVEdeep              = z_wv_distABOVEdeep(:, :, 2:end, :)        ./ v_height; % fraction of cell ABOVE deep boundary;    aligned with z_v; (4D matrix: v longitude (51) X v latitude (80) X num_z (42) X num_t_ROMS (366))

    z_wv_fractionBELOWshallow(z_wv_fractionBELOWshallow > 1)	= 1;
    z_wv_fractionABOVEdeep(z_wv_fractionABOVEdeep > 1)          = 1;

    z_wv_CellFraction                   = z_wv_fractionBELOWshallow .* z_wv_fractionABOVEdeep; % fraction of each ROMS cell in aggregate depth zone; (4D matrix: v longitude (51) X v latitude (80) X num_z (42) X num_t_ROMS (366))

    agg_flux_NS(:, :, depth_loop, :)	= sum((flux_NS .* z_wv_CellFraction(2:(end-1), :, :, :)), 3); % (m3/s); (4D matrix: v longitude-2 (49) X v latitude (80) X num_agg_z (7) X num_t_ROMS (366))
                                                                                                      %     NOTE: omitting far west & east border v terms
                                                                                                      %     NOTE: z_wu_CellFraction currently does not have a t-axis but the element-wise multiplication still works here

    current_v_height                    = v_height .* z_wv_CellFraction; % (m); (3D matrix: v longitude (51) X v latitude (80) X num_z (42))
    agg_v_height(:, :, depth_loop)      = sum(current_v_height, 3);      % height of cell @ v-point; (m); (3D matrix: v longitude-2 (49) X v latitude (80) X num_agg_z (7))
    % ---------------------------------------------------------------------


    % calculate u aggregation ---------------------------------------------
    z_wu_distBELOWshallow               = current_shallow_u - z_wu;	% distance BELOW shallow boundary; (m); (4D matrix: u longitude (50) X u latitude (81) X num_w (43) X num_t_ROMS (366)); NOTE: omitting west, east, south, & north boundary rho points  QQQ still need fix for zeta depth above 0m
    z_wu_distABOVEdeep                  = z_wu - current_deep;         % distance ABOVE deep boundary; (m); (4D matrix: u longitude (50) X u latitude (81) X num_w (43) X num_t_ROMS (366)); NOTE: omitting west, east, south, & north boundary rho points

    z_wu_distBELOWshallow(z_wu_distBELOWshallow < 0)	= 0;
    z_wu_distABOVEdeep(z_wu_distABOVEdeep < 0)          = 0;
    
    z_wu_fractionBELOWshallow           = z_wu_distBELOWshallow(:, :, 1:(end-1), :) ./ u_height; % fraction of cell BELOW shallow boundary; (4D matrix: u longitude (50) X u latitude (81) X num_z (42) X num_t_ROMS (366))
    z_wu_fractionABOVEdeep              = z_wu_distABOVEdeep(:, :, 2:end, :)        ./ u_height; % fraction of cell ABOVE deep boundary;    (4D matrix: u longitude (50) X u latitude (81) X num_z (42) X num_t_ROMS (366))
    
    z_wu_fractionBELOWshallow(z_wu_fractionBELOWshallow > 1)	= 1;
    z_wu_fractionABOVEdeep(z_wu_fractionABOVEdeep > 1)          = 1;
    
    z_wu_CellFraction                   = z_wu_fractionBELOWshallow .* z_wu_fractionABOVEdeep; % fraction of each ROMS cell in aggregate depth zone; (4D matrix: u longitude (50) X u latitude (81) X num_z (42) X num_t_ROMS (366))
            
    agg_flux_EW(:, :, depth_loop, :)	= sum((flux_EW .* z_wu_CellFraction(:, 2:(end-1), :, :)), 3); % (m3/s); (4D matrix: u longitude (50) X u latitude-2 (79) X num_agg_z (7) X num_t_ROMS (366))
                                                                                                      %     NOTE: omitting far south & north border u terms
                                                                                                      %     NOTE: z_wu_CellFraction currently does not have a t-axis but the element-wise multiplication still works here
    
	current_u_height                    = u_height .* z_wu_CellFraction; % (m); (3D matrix: u longitude (50) X u latitude (81) X num_z (42))
    agg_u_height(:, :, depth_loop)      = sum(current_u_height, 3);      % height of cell @ u-point; (m); (3D matrix: u longitude (50) X u latitude (81) X num_agg_z (7))
    % ---------------------------------------------------------------------

    
    % calculate BGC aggregation @ rho -------------------------------------
%     z_w_distBELOWdeep                   = current_deep - z_w(:, :, 1:(end-1), :);       % distance BELOW deep boundary;    (m); (4D matrix: rho longitude (51) X rho latitude (81) X num_z (42) X num_t_ROMS)
    z_w_distABOVEdeep                   = z_w(:, :, 2:end, :) - current_deep;           % distance ABOVE deep boundary;    (m); (4D matrix: rho longitude (51) X rho latitude (81) X num_z (42) X num_t_ROMS)
    z_w_distBELOWshallow                = current_shallow_w - z_w(:, :, 1:(end-1), :);  % distance BELOW shallow boundary; (m); (4D matrix: rho longitude (51) X rho latitude (81) X num_z (42) X num_t_ROMS)
%     z_w_distABOVEshallow                = z_w(:, :, 2:end, :) - current_shallow_w;      % distance ABOVE deep boundary;    (m); (4D matrix: rho longitude (51) X rho latitude (81) X num_z (42) X num_t_ROMS)
    
    ROMScell_height = z_w(:, :, 2:end, :) - z_w(:, :, 1:(end-1), :); % height of each ROMS cell @ rho (box center looking from above); (m); (4D matrix: rho longitude (51) X rho latitude (81) X num_z (42) X num_t_ROMS)
    
	z_w_distABOVEdeep(z_w_distABOVEdeep < 0)     	= 0;
	z_w_distBELOWshallow(z_w_distBELOWshallow < 0)	= 0;

    z_w_fractionABOVEdeep               = z_w_distABOVEdeep ./ ROMScell_height; % fraction of cell ABOVE deep boundary;    (4D matrix: rho longitude (51) X rho latitude (81) X num_z (42) X num_t_ROMS)
	z_w_fractionABOVEdeep(z_w_fractionABOVEdeep > 1)          = 1;
    
	z_w_fractionBELOWshallow           = z_w_distBELOWshallow ./ ROMScell_height; % fraction of cell BELOW shallow boundary; (4D matrix: rho longitude (51) X rho latitude (81) X num_z (42) X num_t_ROMS)
    z_w_fractionBELOWshallow(z_w_fractionBELOWshallow > 1)	= 1;

	z_w_CellFraction	= z_w_fractionBELOWshallow .* z_w_fractionABOVEdeep; % fraction of each ROMS cell in aggregate depth zone; (4D matrix: rho longitude (51) X rho latitude (81) X num_z (42) X num_t_ROMS)

    ROMScell_weighting      = (z_w_CellFraction .* ROMScell_height); % ROMS cell heights for weighted averages of BGC info; (4D matrix: rho longitude (51) X rho latitude (81) X num_z (42)); QQQ should there not be a time dimension here already?
    ROMScell_weighting      = repmat(ROMScell_weighting, [1, 1, 1, num_t_ROMS]); % (4D matrix: rho longitude (51) X rho latitude (81) X num_z (42)); QQQ should there not be a time dimension here already?
    
    % temperature
    agg_temperature                                     = sum((temperature .* ROMScell_weighting), 3) ./ sum(ROMScell_weighting, 3); % (4D matrix: rho longitude (51) X rho latitude (81) X 1 X num_z (42))
    agg_temperature(isnan(agg_temperature))             = 0; % correct for div/0 errors
    agg_temperature_vertical(:, :, depth_loop, :)       = agg_temperature; % (4D matrix: rho longitude (51) X rho latitude (81) X num_agg_z X num_t_ROMS)
    
    % diatom
    agg_diatom                                          = sum((diatom .* ROMScell_weighting), 3) ./ sum(ROMScell_weighting, 3); % (4D matrix: rho longitude (51) X rho latitude (81) X 1 X num_z (42))
    agg_diatom(isnan(agg_diatom))                       = 0; % correct for div/0 errors
    agg_diatom_vertical(:, :, depth_loop, :)            = agg_diatom; % (4D matrix: rho longitude (51) X rho latitude (81) X num_agg_z X num_t_ROMS)
    
    % nanophytoplankton
    agg_nanophytoplankton                           	= sum((nanophytoplankton .* ROMScell_weighting), 3) ./ sum(ROMScell_weighting, 3); % (4D matrix: rho longitude (51) X rho latitude (81) X 1 X num_z (42))
    agg_nanophytoplankton(isnan(agg_nanophytoplankton))	= 0; % correct for div/0 errors
    agg_nanophytoplankton_vertical(:, :, depth_loop, :)	= agg_nanophytoplankton; % (4D matrix: rho longitude (51) X rho latitude (81) X num_agg_z X num_t_ROMS)
    % ---------------------------------------------------------------------
    
    
    % calculate w aggregation ---------------------------------------------
    z_w_distBELOWdeep                   = current_deep - z_w(:, :, 1:(end-1), :);       % distance BELOW deep boundary;    (m); (4D matrix: rho longitude (51) X rho latitude (81) X num_z (42) X num_t_ROMS)
    z_w_distABOVEdeep                   = z_w(:, :, 2:end, :) - current_deep;           % distance ABOVE deep boundary;    (m); (4D matrix: rho longitude (51) X rho latitude (81) X num_z (42) X num_t_ROMS)
    z_w_distBELOWshallow                = current_shallow_w - z_w(:, :, 1:(end-1), :);  % distance BELOW shallow boundary; (m); (4D matrix: rho longitude (51) X rho latitude (81) X num_z (42) X num_t_ROMS)
    z_w_distABOVEshallow                = z_w(:, :, 2:end, :) - current_shallow_w;      % distance ABOVE deep boundary;    (m); (4D matrix: rho longitude (51) X rho latitude (81) X num_z (42) X num_t_ROMS)
    
    z_w_fracBELOWdeep                   = 1 - (z_w_distBELOWdeep    ./ w_height);       % fraction BELOW deep boundary;    (4D matrix: rho longitude (51) X rho latitude (81) X num_z (42) X num_t_ROMS)
    z_w_fracABOVEdeep                   = 1 - (z_w_distABOVEdeep    ./ w_height);       % fraction ABOVE deep boundary;    (4D matrix: rho longitude (51) X rho latitude (81) X num_z (42) X num_t_ROMS)
    z_w_fracBELOWshallow                = 1 - (z_w_distBELOWshallow ./ w_height);       % fraction BELOW shallow boundary; (4D matrix: rho longitude (51) X rho latitude (81) X num_z (42) X num_t_ROMS)
    z_w_fracABOVEshallow                = 1 - (z_w_distABOVEshallow ./ w_height);       % fraction ABOVE deep boundary;    (4D matrix: rho longitude (51) X rho latitude (81) X num_z (42) X num_t_ROMS)
    
    z_w_fracBELOWdeep(z_w_fracBELOWdeep < 0)        = 0;
    z_w_fracABOVEdeep(z_w_fracABOVEdeep < 0)        = 0;
    z_w_fracBELOWshallow(z_w_fracBELOWshallow < 0)  = 0;
    z_w_fracABOVEshallow(z_w_fracABOVEshallow < 0)  = 0;
    
    z_w_fracBELOWdeep(z_w_fracBELOWdeep > 1)        = 0;
    z_w_fracABOVEdeep(z_w_fracABOVEdeep > 1)        = 0;
    z_w_fracBELOWshallow(z_w_fracBELOWshallow > 1)  = 0;
    z_w_fracABOVEshallow(z_w_fracABOVEshallow > 1)  = 0;
    
    w_BELOWdeep                         = z_w_fracBELOWdeep(2:(end-1), 2:(end-1), :, :)    .* flux_UpDown(:, :, 1:(end-1), :); % contribution to w from BELOW deep boundary;    (m3/s); (4D matrix: rho longitude-2 (49) X rho latitude-2 (79) X num_z (42) X num_t_ROMS (366)); NOTE: omitting far west, east, south, & north border rho points
    w_ABOVEdeep                         = z_w_fracABOVEdeep(2:(end-1), 2:(end-1), :, :)    .* flux_UpDown(:, :, 2:end, :);     % contribution to w from ABOVE deep boundary;    (m3/s); (4D matrix: rho longitude-2 (49) X rho latitude-2 (79) X num_z (42) X num_t_ROMS (366)); NOTE: omitting far west, east, south, & north border rho points
    w_BELOWshallow                      = z_w_fracBELOWshallow(2:(end-1), 2:(end-1), :, :) .* flux_UpDown(:, :, 1:(end-1), :); % contribution to w from BELOW shallow boundary; (m3/s); (4D matrix: rho longitude-2 (49) X rho latitude-2 (79) X num_z (42) X num_t_ROMS (366)); NOTE: omitting far west, east, south, & north border rho points
    w_ABOVEshallow                      = z_w_fracABOVEshallow(2:(end-1), 2:(end-1), :, :) .* flux_UpDown(:, :, 2:end, :);	   % contribution to w from ABOVE deep boundary;    (m3/s); (4D matrix: rho longitude-2 (49) X rho latitude-2 (79) X num_z (42) X num_t_ROMS (366)); NOTE: omitting far west, east, south, & north border rho points
    
    w_deep                              = w_BELOWdeep + w_ABOVEdeep;      	% w across deep boundary;    (m3/s); (4D matrix: rho longitude-2 (49) X rho latitude-2 (79) X num_z (42) X num_t_ROMS)
    w_shallow                           = w_BELOWshallow + w_ABOVEshallow;	% w across shallow boundary; (m3/s); (4D matrix: rho longitude-2 (49) X rho latitude-2 (79) X num_z (42) X num_t_ROMS)
    
    agg_flux_UpDown(:, :, depth_loop, :)        = sum(w_shallow, 3); % vertical flux; ; (m3/s); (4D matrix: rho longitude-2 (49) X rho latitude-2 (79) X num_agg_z+1 (8) X num_t_ROMS (366))
    agg_flux_UpDown(:, :, (depth_loop+1), :)	= sum(w_deep, 3);    % vertical flux; ; (m3/s); (4D matrix: rho longitude-2 (49) X rho latitude-2 (79) X num_agg_z+1 (8) X num_t_ROMS (366))
                                                                     %   NOTE: agg_flux_UpDown is arranged shallowest layer on top and deepest on bottom (opposite to order of ROMS matrices)
	% ---------------------------------------------------------------------

end % depth_loop
% -------------------------------------------------------------------------


% step 5c: interpolate agg heights at corner (psi) points -----------------
agg_v_height_psi    = (agg_v_height(1:(end-1), :, :) + agg_v_height(2:end, :, :)) / 2; % height of cell @ psi point based on mean height @ neighbor v points; (m); (3D matrix: v longitude-1 (50) X v latitude (80) X num_agg_z (7))
agg_u_height_psi    = (agg_u_height(:, 1:(end-1), :) + agg_u_height(:, 2:end, :)) / 2; % height of cell @ psi point based on mean height @ neighbor u points; (m); (3D matrix: u longitude (50) X u latitude-1 (80) X num_agg_z (7))
agg_psi_height      = (agg_v_height_psi + agg_u_height_psi) / 2;                       % height of cell @ psi point based on mean height @ neighbor v & u points; (m); (3D matrix: rho longitude-1 (50) X rho latitude-1 (80) X num_agg_z (7))
% -------------------------------------------------------------------------


% step 5d: recalculate vertical fluxes to conserve volume -----------------
%          Required to account for the change in number of pre-aggregated ROMS cell depth zones now passing through each of the 4 horizontal faces
%          NOTE: northward flux is "+"
%          NOTE: eastward flux is "+"
%          NOTE: upward flux is "+"
%          NOTE: face order matters; 
net_flux_NS_2           = agg_flux_NS(:, 1:(end-1), :, :) - agg_flux_NS(:, 2:end, :, :); % net N<-->S volume flux; (flux through SOUTH face minus flux through NORTH face); (m3/s); (3D matrix: 49 X 79 X num_agg_z (7) X num_t_ROMS)
net_flux_EW_2           = agg_flux_EW(1:(end-1), :, :, :) - agg_flux_EW(2:end, :, :, :); % net E<-->W volume flux; (flux through WEST face minus flux through EAST face);   (m3/s); (3D matrix: 49 X 79 X num_agg_z (7) X num_t_ROMS)
net_flux_horizontal_2	= net_flux_NS_2 + net_flux_EW_2; % (m3/s); (3D matrix: psi longitude-1 (49) X psi latitude-1 (79) X num_agg_z (7) X num_t_ROMS (366))

required_w_top          = zeros((rows_rho-2), (clms_rho-2), (num_agg_z+1), num_t_ROMS); % (4D matrix: rho longitude-2 (49) X rho latitude-2 (79) X num_agg_z+1 X num_t_ROMS)
required_w_bottom       = zeros((rows_rho-2), (clms_rho-2), (num_agg_z+1), num_t_ROMS); % (4D matrix: rho longitude-2 (49) X rho latitude-2 (79) X num_agg_z+1 X num_t_ROMS)
agg_flux_UpDown_recalc	= zeros((rows_rho-2), (clms_rho-2), (num_agg_z+1), num_t_ROMS); % (4D matrix: rho longitude-2 (49) X rho latitude-2 (79) X num_agg_z+1 X num_t_ROMS)

for depth_loop = 1:num_agg_z % work from top down
    
    % work from top down
    required_w_top(:, :, (depth_loop+1), :)                 = (-1) * (net_flux_horizontal_2(:, :, depth_loop, :) - required_w_top(:, :, depth_loop, :)); % (m3/s); (4D matrix: rho longitude-2 (49) X rho latitude-2 (79) X num_agg_z+1 X num_t_ROMS)
    
    % work from bottom up
    required_w_bottom(:, :, (num_agg_z+1-depth_loop), :)	= (net_flux_horizontal_2(:, :, (num_agg_z+1-depth_loop), :) + required_w_bottom(:, :, (num_agg_z+2-depth_loop), :)); % (m3/s); (4D matrix: rho longitude-2 (49) X rho latitude-2 (79) X num_agg_z+1 X num_t_ROMS)

end % (depth_loop)

% use top-->down calcs for upper fluxes & bottom-->up calcs for lower fluxes
%       NOTE: this puts non-conservation error in the middle of the water-column where the depth zone heights are greatest and the relative error of the vertical flux (compared to horizontal fluxes) would likely be the smallest
looky_UpperZones	= 1:round((num_agg_z + 1) / 2);
looky_LowerZones	= (max(looky_UpperZones)+1):(num_agg_z + 1);
agg_flux_UpDown_recalc(:, :, looky_UpperZones, :) = required_w_top(:, :, looky_UpperZones, :);
agg_flux_UpDown_recalc(:, :, looky_LowerZones, :) = required_w_bottom(:, :, looky_LowerZones, :); % (m3/s); (4D matrix: rho longitude-2 (49) X rho latitude-2 (79) X num_agg_z+1 X num_t_ROMS)
% -------------------------------------------------------------------------


% step 5e: test volume conservation ---------------------------------------
%          NOTE: northward flux is "+"
%          NOTE: eastward flux is "+"
%          NOTE: upward flux is "+"
%          NOTE: face order matters; 
% net_flux_NS_2           = agg_flux_NS(:, 1:(end-1), :, :)    - agg_flux_NS(:, 2:end, :, :); % net N<-->S volume flux; (flux through SOUTH face minus flux through NORTH face); (m3/s); (4D matrix: 49 X 79 X num_agg_z (7) X num_t_ROMS)
% net_flux_EW_2           = agg_flux_EW(1:(end-1), :, :, :)    - agg_flux_EW(2:end, :, :, :); % net E<-->W volume flux; (flux through WEST face minus flux through EAST face);   (m3/s); (4D matrix: 49 X 79 X num_agg_z (7) X num_t_ROMS)
% net_flux_UpDown_2       = agg_flux_UpDown(:, :, 2:end, :)        - agg_flux_UpDown(:, :, 1:(end-1), :); % net TOP<-->BOTTOM flux; (flux through BOTTOM face minus flux through TOP face);  (m3/s); (4D matrix: 49 X 79 X num_agg_z (7) X num_t_ROMS)
net_flux_UpDown_2       = agg_flux_UpDown_recalc(:, :, 2:end, :) - agg_flux_UpDown_recalc(:, :, 1:(end-1), :); % net TOP<-->BOTTOM flux; (flux through BOTTOM face minus flux through TOP face); (m3/s); (4D matrix: rho longitude-2 (49) X rho latitude-2 (79) X num_agg_z (7) X num_t_ROMS (366))
                                        %   NOTE: agg_flux_UpDown is arranged shallowest layer on top and deepest on bottom (opposite to order of ROMS matrices)

net_flux_horizontal_2	= net_flux_NS_2 + net_flux_EW_2; % (m3/s); (3D matrix: psi longitude-1 (49) X psi latitude-1 (79) X num_agg_z (7) X num_t_ROMS (366))
net_flux_horizontal_2	= sum(net_flux_horizontal_2, 3); % net HORIZONTAL volume flux into a column of water; (m3/s); (4D matrix: psi longitude-1 (49) X psi latitude-1 (79) X 1 X num_t_ROMS (366))
% hist(net_flux_horizontal_2(:), 1000)
                        	% NOTE: at this point, horizontal conservation is good. Most error is well within 20 m3/s. There are 22 instances greater than 200m3/s. Max error is 525 m3/s.
                  
net_flux_2              = net_flux_NS_2 + net_flux_EW_2 + net_flux_UpDown_2; % (m3/s); (3D matrix: psi longitude-1 (49) X psi latitude-1 (79) X num_agg_z (7) X num_t_ROMS)
% hist(net_flux_2(:), 1000)
                            % NOTE: at this point, 3D volume conservation is good except for the middle depth level. There are 22 instances greater than 200m3/s. Max error is 525 m3/s
% -------------------------------------------------------------------------


% step 5f: calculate euclidean box volumes --------------------------------
%           NOTE: box volumes cannot be calculated from the dimensions of the 6 box faces through which ROMS fluxes are expressed (the edges of the NS and EW faces do not line up, faces are rectangular and NOT trapezoidal)
%           NOTE: the real-world box dimensions must be based on trapezoidal faces
BoxSize	= zeros(num_domains, 13, num_agg_z); % initialize; NOTE: BoxSize will still be a 3D matrix of stacked heights that will need to be flattened to a 2D matrix later)
            % 1) length_SouthFace
            % 2) length_NorthFace 
            % 3) length_WestFace
            % 4) length_EastFace
            % 5) height_SW 
            % 6) height_SE 
            % 7) height_NW 
            % 8) height_NE
            % 9) length_NS_globe
            % 10) area_NorthFace
            % 11) area_SouthFace
            % 12) area_floor
            % 13) BoxVolume

for domain_loop = 1:num_domains
    
    % SOUTH<-->NORTH face lengths
    current_NS_cells_u      = grid_addresses_u(domain_loop, south):grid_addresses_u(domain_loop, north); % (horizontal vector: 1 X num ROMS cells N<-->S)
    current_EW_cells_u      = grid_addresses_u(domain_loop, west):grid_addresses_u(domain_loop, east); % (horizontal vector: 1 X num ROMS cells E<-->W + 1)
    
    current_dist_NS_face	= dist_NS_face(current_EW_cells_u, current_NS_cells_u); % lengths of faces running N<-->S; (m); (2D matrix: num ROMS cells E<-->W + 1 X num ROMS cells N<-->S)
    current_dist_NS         = sum(current_dist_NS_face, 2); % (vertical vector: num ROMS cells E<-->W + 1 X 1)
    
    current_dist_NS_globe	= dist_NS_globe(current_EW_cells_u, current_NS_cells_u); % lengths of faces running N<-->S; (m); (2D matrix: num ROMS cells E<-->W + 1 X num ROMS cells N<-->S)
    current_dist_NS_globe	= sum(current_dist_NS_globe, 2); % (vertical vector: num ROMS cells E<-->W + 1 X 1)
    
    BoxSize(domain_loop, west, :)	= repmat(current_dist_NS(1), [1, 1, num_agg_z]);
    BoxSize(domain_loop, east, :)	= repmat(current_dist_NS(end), [1, 1, num_agg_z]);
    BoxSize(domain_loop, 9, :)      = repmat(current_dist_NS_globe(1), [1, 1, num_agg_z]); % distance between parallel north and south faces (identical for all longitudes)
    
    
    % WEST<-->EAST face lengths
    current_EW_cells_v      = grid_addresses_v(domain_loop, west):grid_addresses_v(domain_loop, east); % (horizontal vector: 1 X num ROMS cells E<-->W)
    current_NS_cells_v      = grid_addresses_v(domain_loop, south):grid_addresses_v(domain_loop, north); % (horizontal vector: 1 X num ROMS cells N<-->S + 1)
    
    current_dist_EW_face	= dist_EW_face(current_EW_cells_v, current_NS_cells_v); % (2D matrix)
    current_dist_EW         = sum(current_dist_EW_face, 1); % (horizontal vector)
    
    BoxSize(domain_loop, south, :) = repmat(current_dist_EW(1), [1, 1, num_agg_z]);
    BoxSize(domain_loop, north, :) = repmat(current_dist_EW(end), [1, 1, num_agg_z]);

    
    % box corner heights @ psi points
    current_agg_psi_height	= agg_psi_height(current_EW_cells_u, current_NS_cells_v, :); % height of cell @ psi point based on mean height @ neighbor v & u points; (m); (3D matrix: num ROMS cells E<-->W + 1 X num ROMS cells N<-->S + 1 X num_agg_z (7))
    psi_height_SW           = current_agg_psi_height(1,   1,   :); % SW corner height; (m); (3D matrix: 1 X 1 X num_agg_z (7))
    psi_height_SE           = current_agg_psi_height(end, 1,   :); % SE corner height; (m); (3D matrix: 1 X 1 X num_agg_z (7))
    psi_height_NW           = current_agg_psi_height(1,   end, :); % NW corner height; (m); (3D matrix: 1 X 1 X num_agg_z (7))
    psi_height_NE           = current_agg_psi_height(end, end, :); % NE corner height; (m); (3D matrix: 1 X 1 X num_agg_z (7))
    
    BoxSize(domain_loop, southwest, :) = psi_height_SW;
    BoxSize(domain_loop, southeast, :) = psi_height_SE;
    BoxSize(domain_loop, northwest, :) = psi_height_NW;
    BoxSize(domain_loop, northeast, :) = psi_height_NE;
    
    
    % calculate North & South trapezoidal face areas
    area_NorthFace          = ((psi_height_NW + psi_height_NE) / 2) .* BoxSize(domain_loop, north, :); % area of NORTH face; (m2); (3D matrix: 1 X 1 X num_agg_z (7))
    area_SouthFace          = ((psi_height_SW + psi_height_SE) / 2) .* BoxSize(domain_loop, south, :); % area of SOUTH face; (m2); (3D matrix: 1 X 1 X num_agg_z (7))
    BoxSize(domain_loop, 10, :) = area_NorthFace;
    BoxSize(domain_loop, 11, :)	= area_SouthFace;
    
    
    % domain floor area
	current_NS_cells_rho	= grid_addresses_rho(domain_loop, south):grid_addresses_rho(domain_loop, north); % (horizontal vector: 1 X num ROMS cells N<-->S)
    current_EW_cells_rho	= grid_addresses_rho(domain_loop, west):grid_addresses_rho(domain_loop, east); % (horizontal vector: 1 X num ROMS cells E<-->W)

    current_area_floor_face	= area_floor_face(current_EW_cells_rho, current_NS_cells_rho); % areas of cell floors within domian; (m2); (2D matrix: num ROMS cells E<-->W X num ROMS cells N<-->S)
    current_area_floor      = sum(sum(current_area_floor_face)); % floor area calculated from cell face lengths; (m2); (scalar)
    BoxSize(domain_loop, 12, :)	= repmat(current_area_floor, [1, 1, num_agg_z]); % distance between parallel north and south faces (identical for all longitudes)
    
    
    % calculate Box Volume as a truncated pyramid (frustrum) with bases being the parallel north & south faces
    BoxSize(domain_loop, 13, :) = (BoxSize(domain_loop, 9, :) / 3) .* (area_NorthFace + area_SouthFace + sqrt(area_NorthFace .* area_SouthFace)); % (m3); (3D matrix: num_domains X 13 X num_agg_z (7))
    
end % (domain_loop)
% -------------------------------------------------------------------------


% step 5g: clear temporary variables --------------------------------------
clear current_*
clear z_wv_distBELOWshallow z_wv_distABOVEdeep z_wv_fractionBELOWshallow z_wu_fractionABOVEdeep
clear repmat_dist_EW_v repmat_dist_NS_u repmat_area_*
clear required_w_top required_w_bottom
clear looky_*
clear agg_v_height_psi agg_u_height_psi psi_height_*
clear agg_temperature agg_NH4 agg_NO3 agg_DON agg_PON agg_diatom agg_nanophytoplankton
clear agg_microzooplankton agg_mesozooplankton agg_Pzooplankton
% *************************************************************************





% *************************************************************************
% STEP 6: calculate fluxes across ECOTRAN domain boundaries----------------
%         Re-express v & u flux in terms of DESTINY gains & losses rather than compass direction

% step 6a: latitudinal fluxes (v) -----------------------------------------
%          NOTE: north is positive (+)
[looky_neighbor, looky_DOI]	= find(LatitudinalConnectivity == 1); % find neighbors (looky_neighbor) of each domain (looky_DOI) within the ECOTRAN domain connectivity matrix (num_domainsWboundary X num_domainsWboundary)
num_connections             = length(looky_neighbor);
DS_flux_NS               	= zeros([num_domainsWboundary num_domainsWboundary num_agg_z num_t_ROMS]); % initialize; (m3/s); (4D matrix: num_domainsWboundary X num_domainsWboundary X num_agg_z X num_t_ROMS); NOTE use of num_agg_z

for v_loop = 1:num_connections
    
    current_DOI         = looky_DOI(v_loop); % current ECOTRAN domain of interest (DOI)
    current_neighbor	= looky_neighbor(v_loop);
    
    current_W_v         = LatitudeFluxConnectivity_W_v(current_neighbor, current_DOI); % row address of western ROMS cell in v; (scaler); NOTE: lon_v is in longitude center of cell, use lon_u for longitude boundaries
    current_E_v         = LatitudeFluxConnectivity_E_v(current_neighbor, current_DOI); % row address of eastern ROMS cell in v; (scaler); NOTE: lon_v is in longitude center of cell, use lon_u for longitude boundaries
    
    % find clm address of neighboring (shared) latitude 
	currentDOI_N_v      = grid_addresses_v(current_DOI, north);     % north v-point ROMS clm address of current ECOTRAN DOI; (scalar)
	currentDOI_S_v      = grid_addresses_v(current_DOI, south);     % south v-point ROMS clm address of current ECOTRAN DOI; (scalar)
    
	currentNeighbor_N_v	= grid_addresses_v(current_neighbor, north);     % north v-point ROMS clm address of current ECOTRAN neighbor; (scalar)
	currentNeighbor_S_v	= grid_addresses_v(current_neighbor, south);     % south v-point ROMS clm address of current ECOTRAN neighbor; (scalar)
    
    evaluate_matrix     = [currentDOI_N_v currentNeighbor_N_v; currentDOI_S_v currentNeighbor_S_v]; % (2D matrix: 2 X 2)
    
    diag_N2S            = diag(evaluate_matrix);
    diag_S2N            = fliplr(diag(flipud(evaluate_matrix)));
    
    evaluate_relation	= diff([diag_N2S diag_S2N]);
    
    DOIposition      = find(evaluate_relation == 0); % 1 = DOI is to SOUTH of neighbor; 2 = DOI is to NORTH of neighbor
    
    if DOIposition == 1
        use_DOIlatitude     = currentDOI_N_v; % northern latitude is shared by DOI & neighbor
        flux_sign           = -1; % sign to have boxes indicated as source (-) or destiny (+) of flux; northward ROMS v is given as "+"
    elseif DOIposition == 2
        use_DOIlatitude     = currentDOI_S_v; % southern latitude is shared by DOI & neighbor
        flux_sign           = +1; % sign to have boxes indicated as source (-) or destiny (+) of flux; northward ROMS v is given as "+"
    else
        error('No shared latitude where boxes were determined to be neighbors')
    end
    
    % calculate net DS_flux_NS across shared border
    current_v               = agg_flux_NS(current_W_v:current_E_v, use_DOIlatitude, :, :); % N<-->S volume flux; (m3/s); (4D matrix: num_ROMS cells E<-->W X 1 X num_agg_z X num_t_ROMS)
    looky_filler            = find(abs(current_v) > 10^20); % filter out "filler" values; QQQ filter out junk values earlier in code
    current_v(looky_filler)	= NaN; % v; (m2/s); (4D matrix: num_ROMS cells E<-->W X 1 X num_agg_z X num_t_ROMS)
    current_v               = nansum(current_v, 1); % (m3/s); (4D matrix: 1 X 1 X num_agg_z X num_t_ROMS)
    
    source_tag              = current_v * flux_sign; % POSITIVE when DOI is DESTINY; NEGATIVE when DOI is SOURCE; (4D matrix: 1 X 1 X num_agg_z X num_t_ROMS)
    source_tag(source_tag > 0) = +1; % source_tag is either +1 or -1; (4D matrix: 1 X 1 X num_agg_z X num_t_ROMS)
    source_tag(source_tag < 0) = -1; % source_tag is either +1 or -1; (4D matrix: 1 X 1 X num_agg_z X num_t_ROMS)
    
	current_v               = abs(current_v) .* source_tag; % v; NEGATIVE when DOI is SOURCE; POSITIVE when DOI is DESTINY; (m3/s); (4D matrix: 1 X 1 X num_agg_z X num_t_ROMS)
    
    % separate out source and destiny fluxes
    source_v                    = current_v; % (m3/s); (4D matrix: 1 X 1 X num_agg_z X num_t_ROMS)
	destiny_v                   = current_v; % (m3/s); (4D matrix: 1 X 1 X num_agg_z X num_t_ROMS)
    source_v(source_v > 0)      = 0; % keep only negative values for DOI as SOURCE; (m3/s); (4D matrix: 1 X 1 X num_agg_z X num_t_ROMS)
    destiny_v(destiny_v < 0)	= 0; % keep only positive values for DOI as destiny; (m3/s); (4D matrix: 1 X 1 X num_agg_z X num_t_ROMS)

	% put fluxes in proper location of DS_flux_NS matrix, v (N<-->S fluxes)
    %     NOTE: consecutive ADDITION to DS_flux_NS matrix to allow for cases where an ECOTRAN cell can receive or give to more than 1 other ECOTRAN cell
    DS_flux_NS(current_DOI, current_neighbor, :, :) = DS_flux_NS(current_DOI, current_neighbor, :, :) + destiny_v; % flux from neighbor INTO DOI; v (N<-->S fluxes); (m3/s); (4D matrix: num_domainsWboundary X num_domainsWboundary X num_agg_z X num_t_ROMS)
    DS_flux_NS(current_DOI, current_DOI, :, :)      = DS_flux_NS(current_DOI, current_DOI, :, :)      + source_v;  % flux OUT of DOI; v (N<-->S fluxes); (m3/s); (4D matrix: num_domainsWboundary X num_domainsWboundary X num_agg_z X num_t_ROMS)

end % (v_loop)
% -------------------------------------------------------------------------


% step 6b: longitudinal fluxes (u) ----------------------------------------
%          NOTE: east is positive (+)
[looky_neighbor, looky_DOI]	= find(LongitudinalConnectivity == 1); % find neighbors (looky_neighbor) of each domain (looky_DOI) within the ECOTRAN domain connectivity matrix (num_domainsWboundary X num_domainsWboundary)
num_connections             = length(looky_neighbor);
DS_flux_EW                	= zeros([num_domainsWboundary num_domainsWboundary num_agg_z num_t_ROMS]); % initialize; (4D matrix: num_domainsWboundary X num_domainsWboundary X num_agg_z X num_t_ROMS); NOTE use of num_agg_z

for u_loop = 1:num_connections
    
    current_DOI         = looky_DOI(u_loop); % current ECOTRAN domain of interest (DOI)
    current_neighbor	= looky_neighbor(u_loop);
    
    current_S_u         = LongitudeFluxConnectivity_S_u(current_neighbor, current_DOI); % clm address of southern ROMS cell in u; (scaler); NOTE: lat_u is in latitude center of cell, use lat_v for latitude boundaries
    current_N_u         = LongitudeFluxConnectivity_N_u(current_neighbor, current_DOI); % clm address of northern ROMS cell in u; (scaler); NOTE: lat_u is in latitude center of cell, use lat_v for latitude boundaries
    
    % find row address of neighboring (shared) longitude 
	currentDOI_W_u      = grid_addresses_u(current_DOI, west);     % west u-point ROMS row address of current ECOTRAN box (DOI)
	currentDOI_E_u      = grid_addresses_u(current_DOI, east);     % east u-point ROMS row address of current ECOTRAN box (DOI)
    
	currentNeighbor_W_u	= grid_addresses_u(current_neighbor, west);     % west u-point ROMS row address of current ECOTRAN box (neighbor)
	currentNeighbor_E_u	= grid_addresses_u(current_neighbor, east);     % east u-point ROMS row address of current ECOTRAN box (neighbor)
    
    evaluate_matrix     = [currentDOI_W_u currentNeighbor_W_u; currentDOI_E_u currentNeighbor_E_u]; % (2D matrix: 2 X 2)
    
    diag_W2E            = diag(evaluate_matrix); % (vertical vector: 2 X 1)
    diag_E2W            = fliplr(diag(flipud(evaluate_matrix))); % (vertical vector: 2 X 1)
    
    evaluate_relation	= diff([diag_W2E diag_E2W]); % (horizontal vector: 1 X 2)
    
    DOIposition         = find(evaluate_relation == 0); % 1 = DOI is to EAST of neighbor; 2 = DOI is to WEST of neighbor
    
    if DOIposition == 1
        use_DOIlongitude	= currentDOI_W_u;
        flux_sign           = 1; % sign to have boxes indicated as source (-) or destiny (+) of flux; eastward ROMS u is given as "+"
    elseif DOIposition == 2
        use_DOIlongitude	= currentDOI_E_u;
        flux_sign           = -1; % sign to have boxes indicated as source (-) or destiny (+) of flux; eastward ROMS u is given as "+"
    else
        error('No shared longitude where boxes were determined to be neighbors')
    end

    % calculate net DS_flux_EW across shared border
    current_u               = agg_flux_EW(use_DOIlongitude, current_S_u:current_N_u, :, :); % W<-->E volume flux; (m3/s); (4D matrix: 1 X num_ROMS cells S<-->N X num_agg_z X num_t_ROMS)
    looky_filler            = find(abs(current_u) > 10^20); % filter out "filler" values; QQQ filter out junk values earlier in code
    current_u(looky_filler)	= NaN; % u; (m2/s); (4D matrix: 1 X num_ROMS cells S<-->N X num_agg_z X num_t_ROMS)
    current_u               = nansum(current_u, 2); % (m3/s); (4D matrix: 1 X 1 X num_agg_z X num_t_ROMS)
    
    source_tag              = current_u * flux_sign; % POSITIVE when DOI is DESTINY; NEGATIVE when DOI is SOURCE; (4D matrix: 1 X 1 X num_agg_z X num_t_ROMS)
    source_tag(source_tag > 0) = +1; % source_tag is either +1 or -1
    source_tag(source_tag < 0) = -1; % source_tag is either +1 or -1
    
	current_u               = abs(current_u) .* source_tag; % u; NEGATIVE when DOI is SOURCE; POSITIVE when DOI is DESTINY; (m3/s); (4D matrix: 1 X 1 X num_agg_z X num_t_ROMS)
     
    % separate out source and destiny fluxes
    source_u                    = current_u; % (m3/s); (4D matrix: 1 X 1 X num_agg_z X num_t_ROMS)
	destiny_u                   = current_u; % (m3/s); (4D matrix: 1 X 1 X num_agg_z X num_t_ROMS)
    source_u(source_u > 0)      = 0; % keep only negative values for DOI as SOURCE; (m3/s); (4D matrix: 1 X 1 X num_agg_z X num_t_ROMS)
    destiny_u(destiny_u < 0)	= 0; % keep only positive values for DOI as destiny; (m3/s); (4D matrix: 1 X 1 X num_agg_z X num_t_ROMS)
    
    % put fluxes in proper location of DS_flux_EW matrix, u (E<-->W fluxes)
    %     NOTE: consecutive ADDITION to DS_flux_EW matrix to allow for cases where an ECOTRAN cell can receive or give to more than 1 other ECOTRAN cell
    DS_flux_EW(current_DOI, current_neighbor, :, :) = DS_flux_EW(current_DOI, current_neighbor, :, :) + destiny_u; % flux from neighbor INTO DOI; u (E<-->W fluxes); (m3/s); (4D matrix: num_domainsWboundary X num_domainsWboundary X num_agg_z X num_t_ROMS)
    DS_flux_EW(current_DOI, current_DOI, :, :)      = DS_flux_EW(current_DOI, current_DOI, :, :)      + source_u;  % flux OUT of DOI; u (E<-->W fluxes); (m3/s); (4D matrix: num_domainsWboundary X num_domainsWboundary X num_agg_z X num_t_ROMS)
    
end % (u_loop)
% -------------------------------------------------------------------------


% step 6c: recalculate vertical fluxes to conserve volume -----------------
%          Required to account for the change in number of pre-aggregated ROMS cell depth zones now passing through each of the 4 horizontal faces
%          NOTE: DS_flux_NS & DS_flux_EW are arranged in DESTINY<--SOURCE format
%          NOTE: flux_UpDown_recalc is arranged in vector format
%          NOTE: upward flux is "+"
%          NOTE: face order matters; 
DS_net_flux_horizontal	= DS_flux_NS + DS_flux_EW; % (m3/s); (4D matrix: num_domainsWboundary (41) X num_domainsWboundary (41) X num_agg_z (7) X num_t_ROMS (366))
DS_net_flux_horizontal	= sum(DS_net_flux_horizontal, 2); % net flux INTO each DESTINY box; (m3/s); (4D matrix: num_domainsWboundary (41) X 1 X num_agg_z (7) X num_t_ROMS (366))

required_w_top          = zeros(num_domainsWboundary, 1, (num_agg_z+1), num_t_ROMS); % (4D matrix: num_domainsWboundary (41) X 1 X num_agg_z+1 (8) X num_t_ROMS (366))
required_w_bottom       = zeros(num_domainsWboundary, 1, (num_agg_z+1), num_t_ROMS); % (4D matrix: num_domainsWboundary (41) X 1 X num_agg_z+1 (8) X num_t_ROMS (366))
flux_UpDown_recalc      = zeros(num_domainsWboundary, 1, (num_agg_z+1), num_t_ROMS); % (4D matrix: num_domainsWboundary (41) X 1 X num_agg_z+1 (8) X num_t_ROMS (366))

for depth_loop = 1:num_agg_z % work from top down
    
    % work from top down
    required_w_top(:, :, (depth_loop+1), :)                 = (-1) * (DS_net_flux_horizontal(:, :, depth_loop, :)        - required_w_top(:, :, depth_loop, :)); % (m3/s); (4D matrix: num_domainsWboundary (41) X 1 X num_agg_z+1 (8) X num_t_ROMS (366))
    
    % work from bottom up
    required_w_bottom(:, :, (num_agg_z+1-depth_loop), :)	= (DS_net_flux_horizontal(:, :, (num_agg_z+1-depth_loop), :) + required_w_bottom(:, :, (num_agg_z+2-depth_loop), :)); % (m3/s); (4D matrix: num_domainsWboundary (41) X 1 X num_agg_z+1 (8) X num_t_ROMS (366))

end % (depth_loop)

% use top-->down calcs for upper fluxes & bottom-->up calcs for lower fluxes
%       NOTE: this puts non-conservation error in the middle of the water-column where the depth zone heights are greatest and the relative error of the vertical flux (compared to horizontal fluxes) would likely be the smallest
looky_UpperZones        = 1:round((num_agg_z + 1) / 2);
looky_LowerZones        = (max(looky_UpperZones)+1):(num_agg_z + 1);
flux_UpDown_recalc(:, :, looky_UpperZones, :) = required_w_top(:, :, looky_UpperZones, :);
flux_UpDown_recalc(:, :, looky_LowerZones, :) = required_w_bottom(:, :, looky_LowerZones, :); % (m3/s); (4D matrix: num_domainsWboundary (41) X 1 X num_agg_z+1 (8) X num_t_ROMS (366))

% calculate net vertical fluxes
net_flux_UpDown_recalc	= flux_UpDown_recalc(:, :, 2:end, :) - flux_UpDown_recalc(:, :, 1:(end-1), :); % net TOP<-->BOTTOM flux; (flux through BOTTOM face minus flux through TOP face); (m3/s); (4D matrix: num_domainsWboundary (41) X 1 X num_agg_z (7) X num_t_ROMS (366))
                                        %   NOTE: net_flux_UpDown_recalc is arranged shallowest layer on top and deepest on bottom (opposite to order of ROMS matrices)
                                        %   NOTE: net_flux_UpDown_recalc is expressed as net flux into each ECOTRAN destiny domain
% -------------------------------------------------------------------------


% step 6d: evaluate conservation of volume --------------------------------
net_flux_3              = DS_net_flux_horizontal + net_flux_UpDown_recalc; % (m3/s); (3D matrix: num_domainsWboundary (41) X 1 X num_agg_z (7) X num_t_ROMS (366))
% hist(net_flux_3(:), 1000)
                            % NOTE: at this point, 3D volume conservation is good except for the middle depth level; there 15006 non-zero net flux instances, and max absolute net flux value is 6,913,207 (m3/s)
% -------------------------------------------------------------------------


% step 6e: clear temporary variables --------------------------------------
clear current* looky_* evaluate_* num_connections
clear use_DOI* flux_sign diag_* DOIposition
clear source_tag source_v destiny_v source_u destiny_u
clear required_w_top required_w_bottom
% *************************************************************************





% *************************************************************************
% STEP 7: force conservation of volume -----------------------------------
%          NOTE: conservation errors should be small relative to total fluxes

% step 7a: define conectivities between boxes -----------------------------
%           how much of DOI error to distribute to each destiny box

% Use HorizontalConnectivity for maximum connectivity to external boundaries QQQ move into grid code
HorizontalConnectivity                              = ROMSgrid.LongitudinalConnectivity + ROMSgrid.LatitudinalConnectivity; % 0 or 1; (2D matrix: num_domainsWboundary (41) X num_domainsWboundary (41))
HorizontalConnectivity(HorizontalConnectivity > 1)	= 1; % just in case any boxes border each other latitudinally AND longitudinally

% compress boundary cell rows & columns into 1 row & 1 column ----
HorizontalConnectivity((num_domains+1), :)        	= sum(HorizontalConnectivity((num_domains+1):end, :), 1);
HorizontalConnectivity((num_domains+2):end, :)    	= []; % (m3/s); (2D matrix: (num_domains+1) (22) X num_domainsWboundary (41))
HorizontalConnectivity(:, (num_domains+1))        	= sum(HorizontalConnectivity(:, (num_domains+1):end), 2);
HorizontalConnectivity(:, (num_domains+2):end)    	= []; % (m3/s); (2D matrix: num_domains+1 (22) X num_domains+1 (22))
HorizontalConnectivity(HorizontalConnectivity > 0)	= 1;

ExternalConnectivity_H                             	= HorizontalConnectivity(end, :); % (horizontal vector: 1 X SOURCE-->num_domains+1 (22))
looky_ExternalConnectivity                       	= find(ExternalConnectivity_H > 0); % SOURCE boxes that have connectivity to external boundary
looky_NonExternalConnectivity                       = find(ExternalConnectivity_H == 0); % SOURCE boxes that have NO connectivity to external boundary; use these addresses with FractionalDestinies at time-step
HorizontalConnectivity(1:num_domains, looky_ExternalConnectivity) = 0; % zero-out internal donations by SOURCE boxes that have connection to external boundary; leave external boundary row intact; (0 or 1); NOTE: this is done to speed distribution of error to the external environment 

repmat_HorizontalConnectivity                       = repmat(sum(HorizontalConnectivity, 1), [(num_domains+1), 1]); % (2D matrix: DESTINY-->num_domains+1 (22) X SOURCE-->num_domains+1 (22))
HorizontalConnectivity                              = HorizontalConnectivity ./ repmat_HorizontalConnectivity; % (fractions between 0 and 1); (2D matrix: DESTINY-->num_domains+1 (22) X SOURCE-->num_domains+1 (22))
HorizontalConnectivity(end, end)                    = 0.1; % terminal connectivity needed to allow solution
HorizontalConnectivity(isnan(HorizontalConnectivity)) = 0; % correct div/0 error (will probably not be necessary)


% connectivity between boxes as fractions, based on flow rates at each time-step
DS_flux_horizontal                                  = DS_flux_NS + DS_flux_EW; % horizontal flux; (m3/s); (4D matrix: DESTINY-->num_domainsWboundary (41) X SOURCE-->num_domainsWboundary (41) X num_agg_z (7) X num_t_ROMS (366));
DS_flux_horizontal(DS_flux_horizontal < 0)          = 0; % set loss from DOIs (diagonal) to 0

flux_horizontal_total                               = sum(DS_flux_horizontal, 1); % total flux out of each SOURCE box; (m3/s); (4D matrix: 1 X SOURCE-->num_domainsWboundary (41) X num_agg_z (7) X num_t_ROMS (366))
repmat_flux_horizontal_total                        = repmat(flux_horizontal_total, [num_domainsWboundary, 1, 1, 1]); % (m3/s); (4D matrix: num_domainsWboundary (41) X num_domainsWboundary (41) X num_agg_z (7) X num_t_ROMS (366))
FractionalDestinies_temp                            = DS_flux_horizontal ./ repmat_flux_horizontal_total; % fate of flux OUT of each SOURCE box; (unitless); (4D matrix: num_domainsWboundary (41) X num_domainsWboundary (41) X num_agg_z (7) X num_t_ROMS (366))
FractionalDestinies_temp(isnan(FractionalDestinies_temp)) = 0; % correct div/0 errors; (4D matrix: num_domainsWboundary (41) X num_domainsWboundary (41) X num_agg_z (7) X num_t_ROMS (366))
ExternalConnectivity_FD                             = sum(FractionalDestinies_temp((num_domains+1):end, :, :, :), 1); % fraction of flux flowing to external environment; (4D matrix: 1 X num_domainsWboundary (41) X num_agg_z (7) X num_t_ROMS (366))

FractionalDestinies                                 = zeros((num_domains+1), (num_domains+1), num_agg_z, num_t_ROMS); % (4D matrix: num_domains+1 (22) X num_domains+1 (22) X num_agg_z (7) X num_t_ROMS (366))
FractionalDestinies(1:num_domains, 1:num_domains, :, :) = FractionalDestinies_temp(1:num_domains, 1:num_domains, :, :);
FractionalDestinies(end, 1:num_domains, :, :)       = ExternalConnectivity_FD(1, 1:num_domains, :, :);
FractionalDestinies(end, end, :, :)                 = 0.1; % terminal connectivity needed to allow solution
% -------------------------------------------------------------------------


% step 7b: calculate net flux errors -------------------------------------
%           NOTE: include vertical fluxes
flux_error_1                                        = DS_net_flux_horizontal + net_flux_UpDown_recalc;  % net flux error; (m3/s); (4D matrix: num_domainsWboundary (41) X 1 X num_agg_z (7) X num_t_ROMS (366))
                                                                                                        % NOTE: all error should be restricted to the middle depth zone and 0 in all other zones
flux_error_1((num_domains+1):end, 1, :, :)          = 0; % zero out error in external boundary boxes (this error is does not need correction)
% hist(flux_error_1(:), 1000)
                            % NOTE: at this point, 3D volume conservation is good except for the middle depth level;
% -------------------------------------------------------------------------


% step 7c: redistribute the error according to FractionalDestinies --------
%          NOTE: Final changes to negative errors are made to inputs from boundary SOURCE boxes
NEW_DS_flux_horizontal      = DS_flux_NS + DS_flux_EW; % net horizontal flux; (m3/s); (4D matrix: DESTINY-->num_domainsWboundary (41) X SOURCE-->num_domainsWboundary (41) X num_agg_z (7) X num_t_ROMS (366))

% compress boundary cell rows & columns into 1 row & 1 column ----
NEW_DS_flux_horizontal((num_domains+1), :, :, :)                = sum(NEW_DS_flux_horizontal((num_domains+1):end, :, :, :), 1);
NEW_DS_flux_horizontal((num_domains+2):end, :, :, :)            = []; % (m3/s); (4D matrix: (num_domains+1) X num_domainsWboundary X num_agg_z X num_t_ROMS)
NEW_DS_flux_horizontal(:, (num_domains+1), :, :)                = sum(NEW_DS_flux_horizontal(:, (num_domains+1):end, :, :), 2);
NEW_DS_flux_horizontal(:, (num_domains+2):end, :, :)            = []; % (m3/s); (4D matrix: num_domains+1 (22) X num_domains+1 (22) X num_agg_z (7) X num_t_ROMS (366))
NEW_DS_flux_horizontal((num_domains+1), (num_domains+1), :, :)	= sum(NEW_DS_flux_horizontal(1:num_domains, (num_domains+1), :, :), 1) + NEW_DS_flux_horizontal((num_domains+1), (num_domains+1), :, :); % boundary imbalance (this value is not used) QQQ set to 0?

K                           = diag(ones(1, (num_domains+1)), 0); % 0 or 1; (2D matrix: DESTINY-->num_domains+1 (22) X SOURCE-->num_domains+1 (22))
K2                          = ones((num_domains+1)) - K;         % 0 or 1; (2D matrix: DESTINY-->num_domains+1 (22) X SOURCE-->num_domains+1 (22))

% find depth level of flux error(s); NOTE: error should restricted to 1 depth layer only
locate_error_depth          = flux_error_1(1:num_domains, 1, :, :); % (4D matrix: num_domains X 1 X num_agg_z X num_t_ROMS)
locate_error_depth          = squeeze(max(abs(locate_error_depth)));
locate_error_depth          = max(locate_error_depth, [], 2);
looky_error_depth           = find(locate_error_depth > 0);

for depth_loop = 1:length(looky_error_depth)
    
    current_depth = looky_error_depth(depth_loop);
    
    % save record of flux error that gets re-distributed
    total_flux_horizontal = sum(DS_flux_horizontal(1:num_domains, :, current_depth, :), 2); % (m3/s); (4D matrix: DESTINY-->num_domains (41) X SOURCE-->num_domainsWboundary (41) X num_agg_z (7) X num_t_ROMS (366))
    total_flux_horizontal = reshape(total_flux_horizontal, [num_domains, 1, num_t_ROMS]); % (m3/s); (3D matrix: DESTINY-->num_domains X 1 X num_t_ROMS (366))
    error_record(:, 1, :, depth_loop) = total_flux_horizontal; % flux INTO each DESTINY box (excludes transport OUT of box); (m3/s); (4D matrix: DESTINY-->num_domains X 2 X num_t_ROMS (366) X length(looky_error_depth))
    
    total_flux_error = flux_error_1(1:num_domains, 1, current_depth, :); % net flux error; (m3/s); (4D matrix: num_domains (41) X 1 X 1 X num_t_ROMS (366))
    total_flux_error = reshape(total_flux_error, [num_domains, 1, num_t_ROMS]); % (m3/s); (3D matrix: DESTINY-->num_domains X 1 X num_t_ROMS (366))
    error_record(:, 2, :, depth_loop) = total_flux_error; % total flux error; (m3/s); (4D matrix: DESTINY-->num_domains+1 X 2 X num_t_ROMS (366) X length(looky_error_depth))
    % -----------
    
	for time_loop = 1:num_t_ROMS

        current_DS_flux_horizontal      = NEW_DS_flux_horizontal(:, :, current_depth, time_loop); % (m3/s); (2D matrix: DESTINY-->num_domains+1 (22) X SOURCE-->num_domains+1 (22))
        
        current_FractionalDestinies     = FractionalDestinies(:, :, current_depth, time_loop); % fraction (unitless); (2D matrix: DESTINY-->num_domains+1 (22) X SOURCE-->num_domains+1 (22))
        current_HorizontalConnectivity	= HorizontalConnectivity; % 0 to 1; (2D matrix: DESTINY-->num_domains+1 (22) X SOURCE-->num_domains+1 (22))
        current_HorizontalConnectivity(:, looky_NonExternalConnectivity) = current_FractionalDestinies(:, looky_NonExternalConnectivity); % for boxes with NO external connectivity, use connectivities defined by fluxes @ t
        
        current_connectivity            = current_HorizontalConnectivity; % fraction (unitless); (2D matrix: DESTINY-->num_domains+1 (22) X SOURCE-->num_domains+1 (22))
        current_connectivity            = current_connectivity + (K * (-1)); % assign values of -1 along DOI-to-DOI diagonal
        current_connectivity(end, end)	= 0; % eliminate possibility of inclusion of boundary box as an isolated box

        % identify destiny boxes with no connectivity to other boxes, set rows in current_connectivity to 0, and renormalize columns to 1
        sum_current_connectivity        = sum(current_connectivity); % (horizontal vector: 1 X SOURCE-->num_domains+1)
        looky_IsolatedBox               = find(sum_current_connectivity == -1);
        
        if ~isempty(looky_IsolatedBox)
           
            current_connectivity(looky_IsolatedBox, :)  = 0; % allow no boxes to feed into the isolated box (NOTE: source & destiny boxes are not symmetrical, a destiny box may be isolated as a source box)
            current_connectivity(:, looky_IsolatedBox) = 0;
            current_connectivity(:, looky_IsolatedBox) = K(:, looky_IsolatedBox) * (-1);

            % renormalize current_connectivity to 1
            current_connectivity                        = current_connectivity + K; % assign values of 0 along DOI-to-DOI diagonal
            sum_current_connectivity                    = sum(current_connectivity(1:num_domains, :)); % (horizontal vector: 1 X SOURCE-->num_domains+1); exclude DESTINY boundary
            sum_current_connectivity                    = repmat(sum_current_connectivity, [num_domains, 1]); % (2D matrix: DESTINY-->num_domains (21) X SOURCE-->num_domains+1 (22))
            current_connectivity(1:num_domains, :)      = current_connectivity(1:num_domains, :) ./ sum_current_connectivity;
            current_connectivity(isnan(current_connectivity)) = 0; % correct for div/0 error
            current_HorizontalConnectivity              = current_connectivity;
            current_HorizontalConnectivity(end, end)	= 0.1; % necessary to be < 1 to solve
            current_connectivity                        = current_connectivity + (K * (-1)); % assign values of -1 along DOI-to-DOI diagonal; this also makes sure isolated box value on diagonal is also -1
        end % (~isempty(looky_IsolatedBox))
        % -------------------------
        
        
        current_connectivity(end, end)	= -1;

        flux_error_t                    = flux_error_1(1:num_domains, 1, current_depth, time_loop); % (m3/s); (vertical vector: num_domains (21) X 1)
        looky_positive_error            = find(flux_error_t > 0); % boxes with SURPLUS in-flow
        looky_negative_error            = find(flux_error_t < 0); % boxes with DEFICIT in-flow
        
        % process positive errors ----- (there is too much flow into DOI and not enough out of DOI)
        for PositiveDomain_loop = 1:length(looky_positive_error)
            
            current_domain                  = looky_positive_error(PositiveDomain_loop);
            current_flux_error              = flux_error_t(current_domain); % (m3/s); (scalar)
            
            error_vector                    = zeros((num_domains+1), 1); % (vertical vector: num_domains+1 (22) X 1)
            error_vector(current_domain)	= current_flux_error; % all zero except error in DOI column; (vertical vector: num_domains+1 (22) X 1)

            PositiveError_redistribute   	= inv(K - current_HorizontalConnectivity) * error_vector; % error to flow to/from each destiny box; (m3/s); (vertical vector: num_domains+1 (22) X 1); NOTE: use current_HorizontalConnectivity here and NOT current_connectivity
            PositiveError_redistribute(end)	= 0; % don't redistribute error to boundary as SOURCE columns
            PositiveError_redistribute    	= repmat(PositiveError_redistribute', [(num_domains+1), 1]); % (m3/s); (2D matrix: num_domains+1 (22) X num_domains+1 (22)); NOTE transpose

            PositiveFlux_correction         = PositiveError_redistribute .* current_connectivity;           % (m3/s); (2D matrix: num_domains+1 (22) X num_domains+1 (22)); NOTE: all columns should sum to 0; NOTE: use current_connectivity here and NOT current_HorizontalConnectivity
            current_DS_flux_horizontal      = current_DS_flux_horizontal + PositiveFlux_correction; % (m3/s); (2D matrix: num_domains+1 (22) X num_domains+1 (22))
            
        end % (PositiveDomain_loop)
        
        
        % process negative errors ----- (there is too little flow into DOI and too much out of DOI)
        for NegativeDomain_loop = 1:length(looky_negative_error)
            
            current_domain                  = looky_negative_error(NegativeDomain_loop);
            current_flux_error              = flux_error_t(current_domain); % (m3/s); (scalar)
                       
            error_vector                    = zeros((num_domains+1), 1); % (vertical vector: num_domains+1 (22) X 1)
            error_vector(current_domain)	= current_flux_error; % all zero except error in DOI column; (vertical vector: num_domains+1 (22) X 1)

            NegativeError_redistribute    	= inv(K - current_HorizontalConnectivity) * error_vector; % error to flow to/from each destiny box; (m3/s); (vertical vector: num_domains+1 (22) X 1);
            NegativeError_redistribute(end)	= 0; % don't redistribute error to boundary as SOURCE columns
            NegativeError_redistribute   	= repmat(NegativeError_redistribute', [(num_domains+1), 1]); % (m3/s); (2D matrix: num_domains+1 (22) X num_domains+1 (22)); NOTE transpose

            flux_correction                 = NegativeError_redistribute .* current_connectivity;   % (m3/s); (2D matrix: num_domains+1 (22) X num_domains+1 (22)); NOTE: all columns should sum to 0
            current_DS_flux_horizontal      = current_DS_flux_horizontal + flux_correction; % (m3/s); (2D matrix: num_domains+1 (22) X num_domains+1 (22))
            
            % As necessary, adjust flow OUT of SOURCE to external boundary to prevent negative flow to external destiny
            %       Move negative flow to external boundary DESTINY to flow out of SOURCE DOI
            %       Flow from external boundary column will be adjusted to account for negative net flux out of DOI a few steps further down
            looky_NegativeExport = find(current_DS_flux_horizontal(end, :) < 0);
            for error_loop = 1:length(looky_NegativeExport)
                current_DS_flux_horizontal(looky_NegativeExport(error_loop), looky_NegativeExport(error_loop)) = current_DS_flux_horizontal(looky_NegativeExport(error_loop), looky_NegativeExport(error_loop)) + current_DS_flux_horizontal(end, looky_NegativeExport(error_loop));
                current_DS_flux_horizontal(end, looky_NegativeExport(error_loop)) = 0;
            end
            

        end % (NegativeDomain_loop)
        
        % evaluate remaining error
        %       any error remaining after this point (due to isolated source box error) just get shunted to the boundary (assuming any remaining error is the need for exporting volume from a source box)
        remainder_error                                 = sum(current_DS_flux_horizontal(:, 1:num_domains), 1);
        current_DS_flux_horizontal(end, 1:num_domains)	= current_DS_flux_horizontal(end, 1:num_domains) - remainder_error;
        
        NEW_DS_flux_horizontal(:, 1:num_domains, current_depth, time_loop)	= current_DS_flux_horizontal(:, 1:num_domains); % net horizontal flux; (m3/s); (4D matrix: num_domainsWboundary X num_domainsWboundary X num_agg_z X num_t_ROMS)
        
    end % (time_loop)
    
end % (depth_loop)
% -------------------------------------------------------------------------

% QQQ make sure SOURCE external bondary box columns retain their original flux values


% step 7d: Re-evaluate error ----------------------------------------------
%          if error is negative error, adjust input from external boundary SOURCE column)
DS_net_flux_horizontal	= sum(NEW_DS_flux_horizontal, 2); % net flux INTO each DESTINY box; (m3/s); (4D matrix: num_domains+1 (22) X 1 X num_agg_z (7) X num_t_ROMS (366))
flux_error_2            = DS_net_flux_horizontal + net_flux_UpDown_recalc(1:(num_domains+1), 1, :, :);  % net flux error; (m3/s); (3D matrix: num_domainsWboundary (41) X 1 X num_agg_z (7) X num_t_ROMS (366))
                                                                                %   NOTE: all error should be restricted to the middle depth zone and 0 in all other zones
flux_error_2(end, 1, :, :) = 0; % zero out error in external boundary boxes (this error is does not need correction)
flux_error_2 = round(flux_error_2, 7);

% as necessary, adjust external boundary input for SOURCE boxes. This necessary for cases where SOURCE boxes with connectivity to external DESTINY boxes had negative net flux
for depth_loop = 1:length(looky_error_depth)
    
    current_depth	= looky_error_depth(depth_loop);
    
    for time_loop = 1:num_t_ROMS

        flux_error_t	= flux_error_2(1:num_domains, 1, current_depth, time_loop); % (m3/s); (vertical vector: num_domains (21) X 1)

        looky_positive_error = find(flux_error_t > 0); % boxes with SURPLUS in-flow
        if ~isempty(looky_positive_error)
            error('Error in flux correction @ final boundary input fix: positive net flux into SOURCE box')
        end

        looky_negative_error = find(flux_error_t < 0); % boxes with DEFICIT in-flow

        for domain_loop = 1:length(looky_negative_error)

            current_domain                  = looky_negative_error(domain_loop);
            current_flux_error              = flux_error_t(current_domain); % (m3/s); (scalar)

            NEW_DS_flux_horizontal(current_domain, end, current_depth, time_loop) = NEW_DS_flux_horizontal(current_domain, end, current_depth, time_loop) - current_flux_error;
            
        end % (domain_loop)
    end % (time_loop)
end % (depth_loop)
% -------------------------------------------------------------------------


% step 7e: Final re-evaluation of error -----------------------------------
% if it is negative error, try changing input from external boundary column
DS_net_flux_horizontal	= sum(NEW_DS_flux_horizontal, 2); % net flux INTO each DESTINY box; (m3/s); (4D matrix: num_domains+1 (22) X 1 X num_agg_z (7) X num_t_ROMS (366))
flux_error_3            = DS_net_flux_horizontal + net_flux_UpDown_recalc(1:(num_domains+1), 1, :, :);  % net flux error; (m3/s); (3D matrix: num_domainsWboundary (41) X 1 X num_agg_z (7) X num_t_ROMS (366))
                                                                                %   NOTE: all error should be restricted to the middle depth zone and 0 in all other zones
flux_error_3(end, 1, :, :) = 0; % zero out error in external boundary boxes (this error is does not need correction)
flux_error_3 = round(flux_error_3, 7);

if max(abs(flux_error_3)) > 0
    error('Error in flux correction: error persists after correction steps')
end
% -------------------------------------------------------------------------


% step 7f: clear temporary variables --------------------------------------

clear looky_* current_*
clear flux_horizontal_total repmat_flux_horizontal_total
clear ExternalConnectivity_* flux_error_t DS_net_flux_horizontal
clear PositiveDomain_loop PositiveError_redistribute 
clear NegativeDomain_loop NegativeError_redistribute
% *************************************************************************





% *************************************************************************
% STEP 8: reshape final flux matrix --------------------------------------
%          -- collapse boundary fluxes to 1 row & 1 column
%          -- express vertical fluxes in SOURCE-->>DESTINY format
%          -- combine horizontal & vertical fluxes in a single matrix
%          -- compact matrix

% step 8a: compress boundary cell rows & columns into 1 row & 1 column ----
NEW_DS_flux_horizontal((num_domains+1), :, :, :)                = sum(NEW_DS_flux_horizontal((num_domains+1):end, :, :, :), 1);
NEW_DS_flux_horizontal((num_domains+2):end, :, :, :)            = []; % (m3/s); (4D matrix: (num_domains+1) X num_domainsWboundary X num_agg_z X num_t_ROMS)
NEW_DS_flux_horizontal(:, (num_domains+1), :, :)                = sum(NEW_DS_flux_horizontal(:, (num_domains+1):end, :, :), 2);
NEW_DS_flux_horizontal(:, (num_domains+2):end, :, :)            = []; % (m3/s); (4D matrix: num_domains+1 (22) X num_domains+1 (22) X num_agg_z (7) X num_t_ROMS (366))
NEW_DS_flux_horizontal((num_domains+1), (num_domains+1), :, :)	= sum(NEW_DS_flux_horizontal(1:num_domains, (num_domains+1), :, :), 1) + NEW_DS_flux_horizontal((num_domains+1), (num_domains+1), :, :); % boundary imbalance (this value is not used)

net_flux_UpDown_recalc((num_domains+1), 1, :, :)                = sum(net_flux_UpDown_recalc((num_domains+1):end, 1, :, :), 1); % (m3/s); (4D matrix: num_domainsWboundary (41) X 1 X num_agg_z (7) X num_t_ROMS (366))
% NEW_FluxVertical                                                = net_flux_UpDown_recalc(1:(num_domains+1), 1, :, :); % (m3/s); (4D matrix: num_domains+1 (22) X 1 X num_agg_z (7) X num_t_ROMS (366))
% -------------------------------------------------------------------------


% step 8b: express vertical fluxes in SOURCE-->>DESTINY format ------------
DS_flux_UpDown = zeros((num_domains*num_agg_z+1), (num_domains*num_agg_z+1), num_t_ROMS); % initialize vertical flux in SOURCE-->>DESTINY format; (m3/s); (3D matrix: (num_domains*num_agg_z+1) (148) X (num_domains*num_agg_z+1) (148) X num_t_ROMS (366))

for time_loop = 1:num_t_ROMS
    for domain_loop = 1:num_domains
        
        current_agg_w = squeeze(flux_UpDown_recalc(domain_loop, 1, :, time_loop)); % (m3/s); (vertical vector: num_agg_z+1 (8) X 1); QQQ can get rid of border cells in flux_UpDown_recalc
        
        for flux_loop = 2:num_agg_z
            
            current_flux_vertical = current_agg_w(flux_loop); % (m3/s); (scalar)
            
            if current_flux_vertical > 0
                DS_flux_UpDown(((flux_loop-2)*num_domains+domain_loop), ((flux_loop-1)*num_domains+domain_loop), time_loop) = DS_flux_UpDown(((flux_loop-2)*num_domains+domain_loop), ((flux_loop-1)*num_domains+domain_loop), time_loop) + current_flux_vertical;
                DS_flux_UpDown(((flux_loop-1)*num_domains+domain_loop), ((flux_loop-1)*num_domains+domain_loop), time_loop) = DS_flux_UpDown(((flux_loop-1)*num_domains+domain_loop), ((flux_loop-1)*num_domains+domain_loop), time_loop) - current_flux_vertical;
            elseif current_flux_vertical < 0
                DS_flux_UpDown(((flux_loop-1)*num_domains+domain_loop), ((flux_loop-2)*num_domains+domain_loop), time_loop) = DS_flux_UpDown(((flux_loop-1)*num_domains+domain_loop), ((flux_loop-2)*num_domains+domain_loop), time_loop) - current_flux_vertical;
                DS_flux_UpDown(((flux_loop-2)*num_domains+domain_loop), ((flux_loop-2)*num_domains+domain_loop), time_loop) = DS_flux_UpDown(((flux_loop-2)*num_domains+domain_loop), ((flux_loop-2)*num_domains+domain_loop), time_loop) + current_flux_vertical;
            end
            
        end % (flux_loop)
    end % (domain_loop)
end % (time_loop)
% -------------------------------------------------------------------------


% step 8c: unstack depth layers in horizontal flux matrix ----------------
%           NOTE: in SOURCE-->>DESTINY format
reshaped_DS_flux_horizontal     = zeros((num_domains*num_agg_z+1), (num_domains*num_agg_z+1), num_t_ROMS); % initialize vertical flux in SOURCE-->>DESTINY format; (m3/s); (3D matrix: (num_domains*num_agg_z+1) (148) X (num_domains*num_agg_z+1) (148) X num_t_ROMS (366))

for time_loop = 1:num_t_ROMS
    for depth_loop = 1:num_agg_z
        
        current_flux_horizontal             = NEW_DS_flux_horizontal(:, :, depth_loop, time_loop); % (m3/s); (2D matrix: num_domains+1 (22) X num_domains+1 (22))
        current_border_destiny              = current_flux_horizontal(end, 1:(end-1)); % (m3/s); (horizontal vector: 1 X num_domains (21))
        current_border_source               = current_flux_horizontal(1:(end-1), end); % (m3/s); (vertical vector: num_domains (21) X 1)
        current_flux_horizontal(end, :)     = []; % delet border row; (m3/s); (2D matrix: num_domains (21) X num_domains+1 (22))
        current_flux_horizontal(:, end)     = []; % delet border column; (m3/s); (2D matrix: num_domains (21) X num_domains (21))
        
        current_RC                          = ((depth_loop-1)*num_domains+1):((depth_loop-1)*num_domains+num_domains);
        reshaped_DS_flux_horizontal(current_RC, current_RC, time_loop) = current_flux_horizontal;
        
        reshaped_DS_flux_horizontal(end, current_RC, time_loop) = current_border_destiny; % paste in cross-border fluxes
        reshaped_DS_flux_horizontal(current_RC, end, time_loop) = current_border_source;  % paste in cross-border fluxes; (m3/s); (3D matrix: (num_domains*num_agg_z+1) (148) X (num_domains*num_agg_z+1) (148) X num_t_ROMS (366))
        
    end % (depth_loop)
end % (time_loop)
% -------------------------------------------------------------------------


% step 8d: combine vertical & horizontal fluxes & reshape matrix ----------
DS_flux_total       = reshaped_DS_flux_horizontal + DS_flux_UpDown; % total flux in SOURCE-->>DESTINY format; (m3/s); (3D matrix: (num_domains*num_agg_z+1) (148) X (num_domains*num_agg_z+1) (148) X num_t_ROMS (366))
DS_flux_total       = DS_flux_total * (60 * 60 * 24); % convert time units; total flux in SOURCE-->>DESTINY format; (m3/d); (3D matrix: num_t_ROMS (366) X SOURCE (num_domains*num_agg_z+1) (148) X DESTINY (num_domains*num_agg_z+1) (148))
DS_flux_total(DS_flux_total < 0) = 0; % remove values on diagonals (diagonal values are all negative or 0)

% reshape matrix to (3D matrix: time X num_BoxesAndBoundary SOURCE X num_BoxesAndBoundary DESTINY)
ADVECTION           = permute(DS_flux_total, [3, 2, 1]); % total flux in SOURCE-->>DESTINY format; (m3/d); (3D matrix: num_t_ROMS (366) X SOURCE (num_boxes+1) (148) X DESTINY (num_boxes+1) (148))
% -------------------------------------------------------------------------


% step 8e: define HORIZONTALMIXING & VERTICALMIXING fluxes ----------------
%          NOTE: these are all 0 in this ROMS model
HORIZONTALMIXING	= zeros(num_t_ROMS, (num_boxes+1), (num_boxes+1)); % horizontal mixing rate; (m3/d); (3D matrix: num_t_ROMS X SOURCE (num_boxes+1) X DESTINY (num_boxes+1))
VERTICALMIXING      = zeros(num_t_ROMS, (num_boxes+1), (num_boxes+1)); % vertical mixing rate; (m3/d); (3D matrix: num_t_ROMS X SOURCE (num_boxes+1) X DESTINY (num_boxes+1))
% *************************************************************************





% *************************************************************************
% STEP 9: prepare sinking flux info----------------------------------------
% step 9a: unstack depths of BoxVolume & BoxFloorArea info -------------------
%          put boxes in same order as fluxes
%          PUT IN THIS ORDER-->> (2D matrix: time X domain);
BoxVolume               = BoxSize(:, 13, :); % extract box volume info; (m3)
BoxVolume               = BoxVolume(:); % stack in a single vertical vector
BoxVolume               = repmat(BoxVolume', [num_t_ROMS, 1]); % replicate over time vector; (2D matrix: num_t_ROMS (366) X num_boxes (147)); NOTE transpose; FFF allow for changing dynamic height

BoxFloorArea            = BoxSize(:, 12, :); % extract floor area info; (m2)
BoxFloorArea            = BoxFloorArea(:); % stack in a single vertical vector
BoxFloorArea            = BoxFloorArea'; % transpose; (m2); (horizontal vector: 1 X num_boxes (147)); NOTE transpose

BoxFloorArea(end+1)     = 0; % add dimensions (0m2) of boundary import box; (m2); (m2); (horizontal vector: 1 X SOURCE-->(num_boxes+1) (148))
repmat_BoxFloorArea     = repmat(BoxFloorArea, [1, 1, (num_boxes+1)]);	% (m2); (3D matrix: 1 X SOURCE-->(num_boxes+1) (148) X DESTINY-->(num_boxes+1) (148))
% -------------------------------------------------------------------------


% step 9b: reshape matrix -----------------------------------------------
looky_BottomBoxes       = ROMSgrid.looky_BottomBoxes; % addresses of bottom (benthic) boxes
VerticalConnectivity	= permute(VerticalConnectivity, [3, 2, 1]); % connectivity between source & destiny boxes; 0 or 1; (3D matrix: 1 X SOURCE-->(num_boxes+1) (148) X DESTINY-->(num_boxes+1) (148))
VerticalConnectivity(1, looky_BottomBoxes, end)	= 1; % sinking from box 1 to benthos; NOTE: not necessary if there is a specifically defined benthic box in the physical model (e.g., GoMexOcn); QQQ need to add benthic box to NCC and remove export sinking out of boundary
% -------------------------------------------------------------------------


% step 9c: multiply floor area by vertical connectivity -------------------
SINKING                 = VerticalConnectivity .* repmat_BoxFloorArea; % sinking source box floor area; (m2); (3D matrix: 1 X SOURCE-->(num_boxes+1) (148) X DESTINY-->(num_boxes+1) (148))
SINKING                 = repmat(SINKING, [num_t_ROMS, 1, 1]); % connectivity between source & destiny boxes; 0 or 1; (3D matrix: num_t_ROMS (366) X SOURCE-->(num_boxes+1) (148) X DESTINY-->(num_boxes+1) (148))
% *************************************************************************





% *************************************************************************
% STEP 10: compact fluxes-------------------------------------------------
%          compact flux time series when arranged as 3D matrix (num_t_ROMS X SOURCE (num_boxes+1) X DESTINY (num_boxes+1))
% returns:
%       CompactFlux
%           compact_flux            (2D matrix: num_t_ROMS X num_fluxes)
%                                       NOTE: fluxes include all linked boxes +1 for external links
%           looky_flux              (2D matrix: num_fluxes X 3)
%                                       clm 1: (destiny box) = list of boxes importing water volume (+ boundary)
%                                       clm 2: (source box) = list of boxes exporting water volume (+ boundary)
%                                       clm 3: flux address in non-compacted flux 2D matrix: destiny X source
%                                       NOTE: fluxes include all linked boxes +1 for external links
%                                       NOTE: values constant independent of t
%           looky_boundary_import	(2D matrix: num_fluxes_BoundaryImport X 3)
%                                       clm 1: (destiny box) = identity of boxes importing water volume
%                                       clm 2: (source box) = identity of boxes exporting water volume (always the boundary flux number)
%                                       clm 3: (import flux address) = addresses of import flux clm in compact_flux)
%           looky_boundary_export   (2D matrix: num_fluxes_BoundaryExport X 3)
%                                       clm 1: (destiny box) = identity of boxes importing water volume (always the boundary flux number)
%                                       clm 2: (source box) = identity of boxes exporting water volume
%                                       clm 3: (export flux address) = addresses of export fluxes clm in compact_flux
%           unique_source           (vertical vector: list of source boxes (+ boundary))
%           unique_destiny          (vertical vector: list of destiny boxes (+ boundary))
%           num_fluxes              number of realized fluxes between boxes (or boundary) over full time-series
%           fname_CompactFlux       name of this FuncName_CompactFlux function

% step 10a: compact fluxes ------------------------------------------------
CompactFlux_ADVECTION           = f_CompactFluxTimeSeries_11182019(ADVECTION);
CompactFlux_HORIZONTALMIXING	= f_CompactFluxTimeSeries_11182019(HORIZONTALMIXING);
CompactFlux_VERTICALMIXING      = f_CompactFluxTimeSeries_11182019(VERTICALMIXING);
CompactFlux_SINKING             = f_CompactFluxTimeSeries_11182019(SINKING); % compact SINKING as box floor areas and connectivity information; apply functional group sinking speeds in ECOTRANdynamic_ code
% FFF add migration conncectivity here?
% -------------------------------------------------------------------------


% step 10b: test volume conservation --------------------------------------
UnCompactFlux_ADVECTION         = f_UnCompactFluxTimeSeries_12112019(CompactFlux_ADVECTION);
display('Testing ADVECTION for flux imbalances...')
FluxImbalance_ADVECTION         = f_EvaluateFluxBalance_11262021(UnCompactFlux_ADVECTION); % evaluate for any flux imbalance
% *************************************************************************





% *************************************************************************
% STEP 11: prepare temperature & BGC info----------------------------------
% lat_rho_trim	= lat_rho(2:(end-1), 2:(end-1)); % prune EW & NW borders to align with agg_w; (2D matrix: rho longitude-2 (49) X rho latitude-2 (79));
% lon_rho_trim	= lon_rho(2:(end-1), 2:(end-1)); % prune EW & NW borders to align with agg_w; (2D matrix: rho longitude-2 (49) X rho latitude-2 (79));

build_mean_temperature          = zeros(num_t_ROMS, num_domains, num_agg_z); % temperature time-series (deg C); (3D matrix: num_t_ROMS X num_domains X num_agg_z)
build_mean_diatom               = zeros(num_t_ROMS, num_domains, num_agg_z); % diatom time-series (mmole N/m3); (3D matrix: num_t_ROMS X num_domains X num_agg_z)
build_mean_nanophytoplankton	= zeros(num_t_ROMS, num_domains, num_agg_z); % nanophytoplankton time-series (mmole N/m3); (3D matrix: num_t_ROMS X num_domains X num_agg_z)


agg_temperature_vertical_trim       = agg_temperature_vertical(2:(end-1), 2:(end-1), :, :); % prune EW & NW borders to align with agg_w; (4D matrix: rho longitude-2 (49) X rho latitude-2 (79) X num_agg_z X num_t_ROMS);
agg_temperature_vertical_trim(agg_temperature_vertical_trim == 0) = NaN; % make zeros (from mask) into NaN values (so that zeros aren't included in averaging)

agg_diatom_vertical_trim            = agg_diatom_vertical(2:(end-1), 2:(end-1), :, :); % prune EW & NW borders to align with agg_w; (4D matrix: rho longitude-2 (49) X rho latitude-2 (79) X num_agg_z X num_t_ROMS);
agg_diatom_vertical_trim(agg_diatom_vertical_trim == 0) = NaN; % make zeros (from mask) into NaN values (so that zeros aren't included in averaging)

agg_nanophytoplankton_vertical_trim	= agg_nanophytoplankton_vertical(2:(end-1), 2:(end-1), :, :); % prune EW & NW borders to align with agg_w; (4D matrix: rho longitude-2 (49) X rho latitude-2 (79) X num_agg_z X num_t_ROMS);
agg_nanophytoplankton_vertical_trim(agg_nanophytoplankton_vertical_trim == 0) = NaN; % make zeros (from mask) into NaN values (so that zeros aren't included in averaging)

for domain_loop = 1:num_domains
    
    current_grid_addresses_rho          = ROMSgrid.grid_addresses_rho(domain_loop, 1:4); % addresses of current ECOTRAN domain relative to lat_rho_trim & lon_rho_trim matrices; (horizontal matrix: 1 X 4)
    
    current_domain_temperature          = agg_temperature_vertical_trim(current_grid_addresses_rho(west):current_grid_addresses_rho(east), ...
                                                                current_grid_addresses_rho(south):current_grid_addresses_rho(north), ...
                                                                :, :); % potential temperature in ECOTRAN boxes; (Celsius); (4D matrix: num cells (longitude) X num cells (latitude) X num_agg_z X num_ROMS_t)
    current_domain_diatom               = agg_diatom_vertical_trim(current_grid_addresses_rho(west):current_grid_addresses_rho(east), ...
                                                                current_grid_addresses_rho(south):current_grid_addresses_rho(north), ...
                                                                :, :); % potential temperature in ECOTRAN boxes; (Celsius); (4D matrix: num cells (longitude) X num cells (latitude) X num_agg_z X num_ROMS_t)
    current_domain_nanophytoplankton	= agg_nanophytoplankton_vertical_trim(current_grid_addresses_rho(west):current_grid_addresses_rho(east), ...
                                                                current_grid_addresses_rho(south):current_grid_addresses_rho(north), ...
                                                                :, :); % potential temperature in ECOTRAN boxes; (Celsius); (4D matrix: num cells (longitude) X num cells (latitude) X num_agg_z X num_ROMS_t)

    current_area_floor_face             = area_floor_face(current_grid_addresses_rho(west):current_grid_addresses_rho(east), current_grid_addresses_rho(south):current_grid_addresses_rho(north), 1:num_agg_z); % (m2); (3D matrix: num cells (longitude) X num cells (latitude) X num_agg_z); NOTE: each depth layer has the same area
    current_area_floor_face             = repmat(current_area_floor_face, [1, 1, 1, num_t_ROMS]); % (m2); (4D matrix: num cells (longitude) X num cells (latitude) X num_agg_z X num_ROMS_t)

    % QQQ
    current_looky_NaN                   = find(isnan(current_domain_temperature)); % NOTE: same as find(isnan(current_domain_diatom)) and find(isnan(current_domain_nanophytoplankton))
    current_area_floor_face(current_looky_NaN) = NaN; % set NaN BGC cells to NaN in current_area_floor_face
    
    sum_weighted_temperature        = sum((current_domain_temperature .* current_area_floor_face), 1, 'omitnan'); % (deg C * m2); (4D matrix: 1 X num cells (latitude) X num_agg_z X num_ROMS_t)
    sum_weighted_temperature        = sum(sum_weighted_temperature, 2, 'omitnan'); % (4D matrix: 1 X 1 X num_agg_z X num_ROMS_t)
    
	sum_weighted_diatom             = sum((current_domain_diatom .* current_area_floor_face), 1, 'omitnan');
    sum_weighted_diatom             = sum(sum_weighted_diatom, 2, 'omitnan'); % (4D matrix: 1 X 1 X num_agg_z X num_ROMS_t)

	sum_weighted_nanophytoplankton	= sum((current_domain_nanophytoplankton .* current_area_floor_face), 1, 'omitnan');
    sum_weighted_nanophytoplankton	= sum(sum_weighted_nanophytoplankton, 2, 'omitnan'); % (4D matrix: 1 X 1 X num_agg_z X num_ROMS_t)

    sum_area_floor_face         = sum(current_area_floor_face, 1, 'omitnan'); % (m2); (4D matrix: 1 X num cells (latitude) X num_agg_z X num_ROMS_t)
    sum_area_floor_face         = sum(sum_area_floor_face,     2, 'omitnan'); % (4D matrix: 1 X 1 X num_agg_z X num_ROMS_t)
    
    current_mean_temperature        = sum_weighted_temperature ./ sum_area_floor_face; % (4D matrix: 1 X 1 X num_agg_z X num_t_ROMS)
    current_mean_temperature        = squeeze(current_mean_temperature)'; % (2D matrix: num_t_ROMS X num_agg_z)
    current_mean_temperature        = reshape(current_mean_temperature, [num_t_ROMS, 1, num_agg_z]); % (3D matrix: num_t_ROMS X 1 X num_agg_z)

	current_mean_diatom             = sum_weighted_diatom ./ sum_area_floor_face; % (4D matrix: 1 X 1 X num_agg_z X num_t_ROMS)
    current_mean_diatom             = squeeze(current_mean_diatom)'; % (2D matrix: num_t_ROMS X num_agg_z)
    current_mean_diatom             = reshape(current_mean_diatom, [num_t_ROMS, 1, num_agg_z]); % (3D matrix: num_t_ROMS X 1 X num_agg_z)

	current_mean_nanophytoplankton	= sum_weighted_nanophytoplankton ./ sum_area_floor_face; % (4D matrix: 1 X 1 X num_agg_z X num_t_ROMS)
    current_mean_nanophytoplankton	= squeeze(current_mean_nanophytoplankton)'; % (2D matrix: num_t_ROMS X num_agg_z)
    current_mean_nanophytoplankton	= reshape(current_mean_nanophytoplankton, [num_t_ROMS, 1, num_agg_z]); % (3D matrix: num_t_ROMS X 1 X num_agg_z)

    build_mean_temperature(:, domain_loop, :)       = current_mean_temperature; % temperature time-series (deg C); (2D matrix: num_t_ROMS X num_domains X num_agg_z)
    build_mean_diatom(:, domain_loop, :)            = current_mean_diatom; % diatom time-series (mmole N/m3); (2D matrix: num_t_ROMS X num_domains X num_agg_z)
    build_mean_nanophytoplankton(:, domain_loop, :)	= current_mean_nanophytoplankton; % nanophytoplankton time-series (mmole N/m3); (2D matrix: num_t_ROMS X num_domains X num_agg_z)
    
    
%     current_grid_addresses_rho_test(domain_loop, 1) = lat_rho_trim(1, current_grid_addresses_rho(1));
% 	current_grid_addresses_rho_test(domain_loop, 2) = lat_rho_trim(1, current_grid_addresses_rho(2));
%     current_grid_addresses_rho_test(domain_loop, 3) = lon_rho_trim(current_grid_addresses_rho(3), 1);
%     current_grid_addresses_rho_test(domain_loop, 4) = lon_rho_trim(current_grid_addresses_rho(4), 1);
% 
%     current_coordinates_rho(domain_loop, :)	= [ROMSgrid.ECOTRAN_grid_lat_rho(domain_loop, :) ROMSgrid.ECOTRAN_grid_lon_rho(domain_loop, :)];
    
end % (domain_loop)


% step 11b: unstack depths of build_mean_temperature & BGC ----------------
ROMS_temperature        = zeros(num_t_ROMS, num_boxes); % temperature time-series (deg C); (2D matrix: num_t_ROMS X num_boxes)
ROMS_diatom             = zeros(num_t_ROMS, num_boxes); % diatom time-series (mmole N/m3); (2D matrix: num_t_ROMS X num_boxes)
ROMS_nanophytoplankton	= zeros(num_t_ROMS, num_boxes); % nanophytoplankton time-series (mmole N/m3; (2D matrix: num_t_ROMS X num_boxes)


end_pointer = 0;
for depth_loop = 1:num_agg_z
    
    start_pointer	= end_pointer + 1;
    end_pointer     = num_domains * depth_loop; % 15 30 45 60

    ROMS_temperature(:, start_pointer:end_pointer)          = build_mean_temperature(:, :, depth_loop);
    ROMS_diatom(:, start_pointer:end_pointer)               = build_mean_diatom(:, :, depth_loop);
    ROMS_nanophytoplankton(:, start_pointer:end_pointer)	= build_mean_nanophytoplankton(:, :, depth_loop);

end % (depth_loop)

% *************************************************************************





% *************************************************************************
% STEP 12: pack results----------------------------------------------------
% step 12a: save names of this m-file and sub-functions -------------------
ROMSflux.fname_ROMS_FluxPrep            = fname_ROMS_FluxPrep;
ROMSflux.fname_CompactFlux           	= CompactFlux_ADVECTION.fname_CompactFlux;            % file name of f_CompactFluxTimeSeries function
ROMSflux.fname_UnCompactFluxTimeSeries	= UnCompactFlux_ADVECTION.fname_UnCompactFluxTimeSeries; % name of this f_UnCompactFluxTimeSeries sub-function
ROMSflux.fname_CalcNetFlux           	= UnCompactFlux_ADVECTION.fname_CalcNetFlux;          % name of this f_CalcNetFlux function
ROMSflux.fname_EvaluateFluxBalance      = FluxImbalance_ADVECTION.fname_EvaluateFluxBalance;	% name of this f_EvaluateFluxBalance sub-function
% -------------------------------------------------------------------------

% step 11b: time dimensions -----------------------------------------------
ROMSflux.ocean_time                     = ocean_time; % seconds since 1900-01-01 00:00:00; (vertical vector: num_t_ROMS (366) X 1)
ROMSflux.num_t_ROMS                     = num_t_ROMS;
% -------------------------------------------------------------------------

% step 11c: space dimensions ----------------------------------------------
ROMSflux.BoxVolume                      = BoxVolume; % box volume; (m3); (2D matrix: num_t_ROMS (366) X num_boxes (147)); FFF allow for changing dynamic height
ROMSflux.BoxFloorArea                   = BoxFloorArea(:, 1:(end-1)); % box floor area; (m2); (horizontal vector: 1 X num_boxes (147)); NOTE: excludes external boundary that was added earlier for SINKING algorithm
% QQQ mixed layer depth (MLD)?
% -------------------------------------------------------------------------

% step 11d: spatial relationships -----------------------------------------
ROMSflux.VerticalConnectivity           = VerticalConnectivity; % connectivity between source & destiny boxes; 0 or 1; (3D matrix: 1 X SOURCE-->(num_boxes+1) (148) X DESTINY-->(num_boxes+1) (148))
% -------------------------------------------------------------------------

% step 11d: physical driver time-series -----------------------------------
%           NOTE: fluxes are in SOURCE-->>DESTINY format
ROMSflux.ADVECTION                      = ADVECTION;        % advection rate;                (m3/d); (3D matrix: num_t_ROMS (366) X SOURCE (num_boxes+1) (148) X DESTINY (num_boxes+1) (148))
ROMSflux.HORIZONTALMIXING               = HORIZONTALMIXING; % horizontal mixing rate;        (m3/d); (3D matrix: num_t_ROMS (366) X SOURCE (num_boxes+1) X DESTINY (num_boxes+1))
ROMSflux.VERTICALMIXING                 = VERTICALMIXING;   % vertical mixing rate;          (m3/d); (3D matrix: num_t_ROMS (366) X SOURCE (num_boxes+1) X DESTINY (num_boxes+1))
ROMSflux.SINKING                        = SINKING;          % sinking source box floor area; (m2);   (3D matrix: num_t_ROMS (366) X SOURCE (num_boxes+1) X DESTINY (num_boxes+1))

ROMSflux.CompactFlux_ADVECTION          = CompactFlux_ADVECTION;
ROMSflux.CompactFlux_HORIZONTALMIXING	= CompactFlux_HORIZONTALMIXING;
ROMSflux.CompactFlux_VERTICALMIXING     = CompactFlux_VERTICALMIXING;
ROMSflux.CompactFlux_SINKING            = CompactFlux_SINKING; % compact SINKING as box floor areas and connectivity information; apply functional group sinking speeds in ECOTRANdynamic_ code
% -------------------------------------------------------------------------

% step 11e: flux error (before correction) & total horizontal flux --------
ROMSflux.error_record                   = error_record; % total flux error; (m3/s); (4D matrix: DESTINY-->num_domains X 2 [total horizontal flux; flux error] X num_t_ROMS (366) X length(looky_error_depth))
% -------------------------------------------------------------------------

% step 11f: temperature & BGC info ----------------------------------------
ROMSflux.ROMS_temperature               = ROMS_temperature; % temperature time-series (deg C); (2D matrix: num_t_ROMS X num_boxes)
ROMSflux.ROMS_diatom                    = ROMS_diatom; % diatom time-series (mmole N/m3; (2D matrix: num_t_ROMS X num_boxes)
ROMSflux.ROMS_nanophytoplankton         = ROMS_nanophytoplankton; % nanophytoplankton time-series (mmole N/m3); (2D matrix: num_t_ROMS X num_boxes)
% *************************************************************************


% end m-file***************************************************************