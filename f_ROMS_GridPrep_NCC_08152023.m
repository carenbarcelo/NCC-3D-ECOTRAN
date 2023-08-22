function ROMSgrid = f_ROMS_GridPrep_NCC_08152023(PHYSICSinput)
% Map ROMS grid to ECOTRAN grid
% takes:
%       PHYSICSinput
%           readFile_DepthLevels	file defining ROMS grid depth levels
%           readFile_grid           file defining ROMS grid
%           readFile_ExampleYear	example of 1 flux year, used to support grid re-mapping
%
% returns:
%       ROMSgrid structure variable defining how ROMS grid cells map into ECOTRAN grid and their connectivity with each other
% 
% calls:
%       none
%
% revision date: 8/15/2023
%       9/13/2022 processing temperature & BGD info
%       11/28/2022 added VerticalConnectivity_stack telling which boxes lay directly under each surface box; clm 1 is suface, clm 2 is depth layer 2, final column is the bottom
%       8/18/2023 cleaned up notation


% *************************************************************************
% STEP 1: set operating conditions-----------------------------------------
fname_ROMS_GridPrep     = mfilename; % save name of this m-file to keep in saved model results
display(['Running: ' fname_ROMS_GridPrep])

% step 1a: define ROMS NetCDF files to read -------------------------------
readFile_DepthLevels	= PHYSICSinput.readFile_DepthLevels;
% ncdisp(readFile_DepthLevels) % Display contents of NetCDF data source in Command Window

readFile_grid           = PHYSICSinput.readFile_grid;
% ncdisp(readFile_grid) % Display contents of NetCDF data source in Command Window

readFile_ExampleYear    = PHYSICSinput.readFile_ExampleYear; % example of 1 flux year, used to define grid re-mapping
% ncdisp(readFile_ExampleYear) % Display contents of NetCDF data source in Command Window
% -------------------------------------------------------------------------


% step 1b: define ECOTRAN grid --------------------------------------------
%          Boundaries: North    = 48 (47.8)
%                      South    = 40 (40.2)
%                      interval = 0.1

domain_definition = [

% GRID B
        % WA
        47.55 47.15 -15 -150 
        47.55 47.15 -151 -250
        47.55 47.15 -251 -1500

        % CR
        47.15 46.05 -15 -150
        47.15 46.05 -151 -250
        47.15 46.05 -251 -1500

        % NOR
        46.05 44.45 -15 -150
        46.05 44.45 -151 -250
        46.05 44.45 -251 -1500

        % SOR
        44.45 42.05 -15 -150
        44.45 42.05 -151 -250
        44.45 42.05 -251 -1500

        % NCal
        42.05 40.75 -15 -150
        42.05 40.75 -151 -250
        42.05 40.75 -251 -1500
                ]; % [north(lat) south(lat) east(inshore depth) west(offshore depth)]
                
z_definition =	[
    
        % surface
                    0
        % epipelgic bottom
                  -50
        % mesopelagic bottom
                 -100
        % bathypelagic bottom
                 -200
        % benthic bottom
                -5000 % <<--just a catch-all for really deep areas (that probably aren't in the geographic domain)
                ]; % (m)
            

z_definition        = sort(unique([0; z_definition]), 'descend'); % make sure depth zones start at 0 and are put into descending order  

[num_domains, ~]	= size(domain_definition);
num_agg_z           = length(z_definition) - 1; % number of ECOTRAN depth domains
num_boxes           = num_domains * num_agg_z;
% *************************************************************************





% *************************************************************************
% STEP 2: load ROMS variables from NetCDF files----------------------------
% step 2a: load ROMS grid variables ---------------------------------------
ncid        = netcdf.open(readFile_grid, 'NC_NOWRITE'); % Open readFile

% domain-wide terms
% varid       = netcdf.inqVarID(ncid, 'xl'); % Get variable ID of the first variable, given its name
% xl          = netcdf.getVar(ncid, varid); % 'domain length in the XI-direction'; (m); (scalar)

% varid       = netcdf.inqVarID(ncid, 'el');
% el          = netcdf.getVar(ncid, varid); % 'domain length in the ETA-direction'; (m); (scalar)

% varid       = netcdf.inqVarID(ncid, 'depthmin');
% depthmin 	= netcdf.getVar(ncid, varid); % 'Shallow bathymetry clipping depth'; (m); (scalar)

% varid       = netcdf.inqVarID(ncid, 'depthmax');
% depthmax 	= netcdf.getVar(ncid, varid); % 'Deep bathymetry clipping depth'; (m); (scalar)

% varid       = netcdf.inqVarID(ncid, 'spherical');
% spherical 	= netcdf.getVar(ncid, varid); % 'Grid type logical switch'; (scalar)


% values at rho-points (cell centers) (2D matrix: 186 X 181)
% varid       = netcdf.inqVarID(ncid, 'angle');
% angle       = netcdf.getVar(ncid, varid); % 'angle between xi axis and east'; (degree); (2D matrix: 186 X 181)

% varid       = netcdf.inqVarID(ncid, 'hraw');
% h_raw       = netcdf.getVar(ncid, varid); % 'Working bathymetry at RHO-points'; (m); (2D matrix: 186 X 181)

% varid       = netcdf.inqVarID(ncid, 'alpha');
% alpha       = netcdf.getVar(ncid, varid); % 'Weights between coarse and fine grids at RHO-points'; (no units); (2D matrix: 186 X 181)

% varid       = netcdf.inqVarID(ncid, 'f');
% Coriolis    = netcdf.getVar(ncid, varid); % 'Coriolis parameter at RHO-points'; (1/s); (2D matrix: 186 X 181)

% varid       = netcdf.inqVarID(ncid, 'pm');
% pm          = netcdf.getVar(ncid, varid); % 'curvilinear coordinate metric in XI'; (1/m); (2D matrix: 186 X 181)

% varid       = netcdf.inqVarID(ncid, 'pn');
% pn          = netcdf.getVar(ncid, varid); % 'curvilinear coordinate metric in eta'; (1/m); (2D matrix: 186 X 181)

% varid       = netcdf.inqVarID(ncid, 'dndx');
% dndx       	= netcdf.getVar(ncid, varid); % 'xi derivative of inverse metric factor pn'; (m); (2D matrix: 186 X 181)

% varid       = netcdf.inqVarID(ncid, 'dmde');
% dmde       	= netcdf.getVar(ncid, varid); % 'eta derivative of inverse metric factor pm'; (m); (2D matrix: 186 X 181)

varid       = netcdf.inqVarID(ncid, 'lat_rho');
lat_rho     = netcdf.getVar(ncid, varid); % 'latitude of RHO-points'; (degrees east); (2D matrix: 186 X 181)

varid       = netcdf.inqVarID(ncid, 'lon_rho');
lon_rho     = netcdf.getVar(ncid, varid); % 'longitude of RHO-points'; (degrees east); (2D matrix: 186 X 181)

varid       = netcdf.inqVarID(ncid, 'x_rho');
x_rho       = netcdf.getVar(ncid, varid); % 'x location of RHO-points'; (m); (2D matrix: 186 X 181)

varid       = netcdf.inqVarID(ncid, 'y_rho');
y_rho       = netcdf.getVar(ncid, varid); % 'y location of RHO-points'; (m); (2D matrix: 186 X 181)

varid       = netcdf.inqVarID(ncid, 'mask_rho');
mask_rho     = netcdf.getVar(ncid, varid); % 'mask on RHO-points'; (degrees east); (2D matrix: 186 X 181)

varid       = netcdf.inqVarID(ncid, 'h');
h           = netcdf.getVar(ncid, varid); % 'Final bathymetry at RHO-points'; (m); (2D matrix: 186 X 181)


% values at v-points (cell N/S boundaries) (2D matrix: 186 X 180)
varid       = netcdf.inqVarID(ncid, 'lat_v');
lat_v     	= netcdf.getVar(ncid, varid); % 'latitude of V-points'; (degrees east); (2D matrix: 186 X 180)

varid       = netcdf.inqVarID(ncid, 'lon_v');
lon_v     	= netcdf.getVar(ncid, varid); % 'longitude of V-points'; (degrees east); (2D matrix: 186 X 180)

varid       = netcdf.inqVarID(ncid, 'x_v');
x_v         = netcdf.getVar(ncid, varid); % 'x location of V-points'; (m); (2D matrix: 186 X 180)

varid       = netcdf.inqVarID(ncid, 'y_v');
y_v         = netcdf.getVar(ncid, varid); % 'y location of V-points'; (m); (2D matrix: 186 X 180)

varid       = netcdf.inqVarID(ncid, 'mask_v');
mask_v     	= netcdf.getVar(ncid, varid); % 'mask on V-points'; (degrees east); (2D matrix: 186 X 180)


% values at u-points (cell E/W boundaries) (2D matrix: 185 X 181)
varid       = netcdf.inqVarID(ncid, 'lat_u');
lat_u     	= netcdf.getVar(ncid, varid); % 'latitude of U-points'; (degrees east); (2D matrix: 185 X 181)

varid       = netcdf.inqVarID(ncid, 'lon_u');
lon_u     	= netcdf.getVar(ncid, varid); % 'longitude of U-points'; (degrees east); (2D matrix: 185 X 181)

varid       = netcdf.inqVarID(ncid, 'x_u');
x_u         = netcdf.getVar(ncid, varid); % 'x location of U-points'; (m); (2D matrix: 185 X 181)

varid       = netcdf.inqVarID(ncid, 'y_u');
y_u         = netcdf.getVar(ncid, varid); % 'y location of U-points'; (m); (2D matrix: 185 X 181)

varid       = netcdf.inqVarID(ncid, 'mask_u');
mask_u     	= netcdf.getVar(ncid, varid); % 'mask on U-points'; (degrees east); (2D matrix: 185 X 181)


% values at psi-points (cell E/W boundaries) (2D matrix: 185 X 180)
varid       = netcdf.inqVarID(ncid, 'lat_psi');
lat_psi     = netcdf.getVar(ncid, varid); % 'latitude of PSI-points'; (degrees east); (2D matrix: 185 X 180)

varid       = netcdf.inqVarID(ncid, 'lon_psi');
lon_psi     = netcdf.getVar(ncid, varid); % 'longitude of PSI-points'; (degrees east); (2D matrix: 185 X 180)

varid       = netcdf.inqVarID(ncid, 'x_psi');
x_psi     	= netcdf.getVar(ncid, varid); % 'x location of PSI-points'; (m); (2D matrix: 185 X 180)

varid       = netcdf.inqVarID(ncid, 'y_psi');
y_psi     	= netcdf.getVar(ncid, varid); % 'y location of PSI-points'; (m); (2D matrix: 185 X 180)

varid       = netcdf.inqVarID(ncid, 'mask_psi');
mask_psi	= netcdf.getVar(ncid, varid); % 'mask on PSI-points'; (degrees east); (2D matrix: 185 X 180)


% Close the NetCDF file
netcdf.close(ncid)
% -------------------------------------------------------------------------


% step 2b: open NetCDF depth file and load variables ----------------------
ncid        = netcdf.open(readFile_DepthLevels, 'NC_NOWRITE'); % Open readFile

varid       = netcdf.inqVarID(ncid, 'z_rho');
z_rho     	= netcdf.getVar(ncid, varid); % 'depth of rho points'; (m); (2D matrix: 186 X 181 X 42)

varid    	= netcdf.inqVarID(ncid, 'z_w');
z_w       	= netcdf.getVar(ncid, varid); % 'depth of w points'; (m); (2D matrix: 186 X 181 X 43)

varid     	= netcdf.inqVarID(ncid, 'z_v');
z_v       	= netcdf.getVar(ncid, varid); % 'depth of v points'; (m); (2D matrix: 186 X 180 X 42)

varid   	= netcdf.inqVarID(ncid, 'z_u');
z_u        	= netcdf.getVar(ncid, varid); % 'depth of u points'; (m); (2D matrix: 185 X 181 X 42)

varid      	= netcdf.inqVarID(ncid, 'z_psi');
z_psi     	= netcdf.getVar(ncid, varid); % 'depth of PSI-points'; (m); (2D matrix: 185 X 180 X 42)

% Close the NetCDF file
netcdf.close(ncid)
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% step 2c: load ROMS flow variable to get NCC grid extent  ----------------
ncid        = netcdf.open(readFile_ExampleYear, 'NC_NOWRITE'); % Open readFile

% values at rho-points (cell centers)
varid           = netcdf.inqVarID(ncid, 'lat_rho'); % Get variable ID of the first variable, given its name
lat_rho_flow  	= netcdf.getVar(ncid, varid);       % 'latitude of RHO-points'; (2D matrix: rho longitude (51 west:east) X rho latitude (81 south:north))
varid           = netcdf.inqVarID(ncid, 'lon_rho'); 
lon_rho_flow	= netcdf.getVar(ncid, varid);       % 'longitude of RHO-points'; (2D matrix: rho longitude (51 west:east) X rho latitude (81 south:north))

% Close the NetCDF file
netcdf.close(ncid)
% -------------------------------------------------------------------------


% step 2d: clear temporary variables --------------------------------------
clear ncid varid
% *************************************************************************





% *************************************************************************
% STEP 3: trim full grid ROMS terms to only those WITHIN the NCC ----------
% step 3a: round lats & lons (needed because of TINY rounding errors) -----
lat_rho   	= round(lat_rho, 2);
lon_rho  	= round(lon_rho, 2);
lat_v   	= round(lat_v, 2);
lon_v    	= round(lon_v, 2);
lat_u    	= round(lat_u, 2);
lon_u   	= round(lon_u, 2);
lat_psi     = round(lat_psi, 2);
lon_psi     = round(lon_psi, 2);
% -------------------------------------------------------------------------


% step 3b: make sure values align to Arakawa-C grid -----------------------
%          NOTE: using 1 year of the ROMS flow files to define NCC Arakawa-C grid boundaries

% latitude range (outer N/S boundaries defined by rho & u)
NorthBoundary	= max(lat_rho_flow(1, :)); % NOTE: assumes latitudes run across columns; NOTE: assumes location entirely in northern hemisphere
SouthBoundary	= min(lat_rho_flow(1, :)); % NOTE: assumes latitudes run across columns; NOTE: assumes location entirely in northern hemisphere

% longitude range (outer E/W boundaries defined by rho & v)
WestBoundary	= min(lon_rho_flow(:, 1)); % NOTE: assumes longitudes run down rows; NOTE: assumes location entirely in western hemisphere (negative longtitudes)
EastBoundary	= max(lon_rho_flow(:, 1)); % NOTE: assumes longitudes run down rows; NOTE: assumes location entirely in western hemisphere (negative longtitudes)

% rho terms
looky_TooNorth = find(lat_rho(1, :) > NorthBoundary); % NOTE: assumes latitudes run across columns; 
if ~isempty(looky_TooNorth)
    lat_rho(:, looky_TooNorth)      = []; % delete columns north of Arakawa-C longitude domain
    lon_rho(:, looky_TooNorth)      = []; % delete columns north of Arakawa-C longitude domain
    x_rho(:, looky_TooNorth)        = []; % delete columns north of Arakawa-C longitude domain
    y_rho(:, looky_TooNorth)        = []; % delete columns north of Arakawa-C longitude domain
    mask_rho(:, looky_TooNorth)     = []; % delete columns north of Arakawa-C longitude domain
    h(:, looky_TooNorth)            = []; % delete columns north of Arakawa-C longitude domain
    z_rho(:, looky_TooNorth, :)     = []; % delete columns north of Arakawa-C longitude domain
    z_w(:, looky_TooNorth, :)       = []; % delete columns north of Arakawa-C longitude domain
end % (~isempty(looky_TooNorth))
    
looky_TooSouth = find(lat_rho(1, :) < SouthBoundary); % NOTE: assumes latitudes run across columns; 
if ~isempty(looky_TooSouth)
    lat_rho(:, looky_TooSouth)      = []; % delete columns south of Arakawa-C longitude domain
    lon_rho(:, looky_TooSouth)      = []; % delete columns south of Arakawa-C longitude domain
    x_rho(:, looky_TooSouth)        = []; % delete columns south of Arakawa-C longitude domain
    y_rho(:, looky_TooSouth)        = []; % delete columns south of Arakawa-C longitude domain
    mask_rho(:, looky_TooSouth)     = []; % delete columns south of Arakawa-C longitude domain
    h(:, looky_TooSouth)            = []; % delete columns south of Arakawa-C longitude domain
	z_rho(:, looky_TooSouth, :)     = []; % delete columns south of Arakawa-C longitude domain
    z_w(:, looky_TooSouth, :)       = []; % delete columns south of Arakawa-C longitude domain
end % (~isempty(looky_TooSouth))

looky_TooWest	= find(lon_rho(:, 1) < WestBoundary); % NOTE: assumes longitudes run down rows; 
if ~isempty(looky_TooWest)
	lat_rho(looky_TooWest, :)       = []; % delete rows west of Arakawa-C longitude domain
    lon_rho(looky_TooWest, :)       = []; % delete rows west of Arakawa-C longitude domain
    x_rho(looky_TooWest, :)         = []; % delete rows west of Arakawa-C longitude domain
    y_rho(looky_TooWest, :)         = []; % delete rows west of Arakawa-C longitude domain
    mask_rho(looky_TooWest, :)      = []; % delete rows west of Arakawa-C longitude domain
    h(looky_TooWest, :)            	= []; % delete rows west of Arakawa-C longitude domain
	z_rho(looky_TooWest, :, :)      = []; % delete rows west of Arakawa-C longitude domain
    z_w(looky_TooWest, :, :)        = []; % delete rows west of Arakawa-C longitude domain
end % (~isempty(looky_TooWest))

looky_TooEast	= find(lon_rho(:, 1) > EastBoundary); % NOTE: assumes longitudes run down rows; 
if ~isempty(looky_TooEast)
    lat_rho(looky_TooEast, :)       = []; % delete rows east of Arakawa-C longitude domain
    lon_rho(looky_TooEast, :)       = []; % delete rows east of Arakawa-C longitude domain
    x_rho(looky_TooEast, :)         = []; % delete rows east of Arakawa-C longitude domain
    y_rho(looky_TooEast, :)         = []; % delete rows east of Arakawa-C longitude domain
    mask_rho(looky_TooEast, :)      = []; % delete rows east of Arakawa-C longitude domain
    h(looky_TooEast, :)             = []; % delete rows east of Arakawa-C longitude domain
	z_rho(looky_TooEast, :, :)      = []; % delete rows east of Arakawa-C longitude domain
    z_w(looky_TooEast, :, :)        = []; % delete rows east of Arakawa-C longitude domain
end % (~isempty(looky_TooEast))
% -----


% v terms
looky_TooNorth = find(lat_v(1, :) > NorthBoundary); % NOTE: assumes latitudes run across columns; 
if ~isempty(looky_TooNorth)
    lat_v(:, looky_TooNorth)        = []; % delete columns north of Arakawa-C longitude domain
    lon_v(:, looky_TooNorth)        = []; % delete columns north of Arakawa-C longitude domain
    x_v(:, looky_TooNorth)          = []; % delete columns north of Arakawa-C longitude domain
    y_v(:, looky_TooNorth)          = []; % delete columns north of Arakawa-C longitude domain
    mask_v(:, looky_TooNorth)       = []; % delete columns north of Arakawa-C longitude domain
    z_v(:, looky_TooNorth, :)       = []; % delete columns north of Arakawa-C longitude domain
end % (~isempty(looky_TooNorth))
    
looky_TooSouth = find(lat_v(1, :) < SouthBoundary); % NOTE: assumes latitudes run across columns; 
if ~isempty(looky_TooSouth)
    lat_v(:, looky_TooSouth)        = []; % delete columns south of Arakawa-C longitude domain
    lon_v(:, looky_TooSouth)        = []; % delete columns south of Arakawa-C longitude domain
    x_v(:, looky_TooSouth)          = []; % delete columns south of Arakawa-C longitude domain
    y_v(:, looky_TooSouth)          = []; % delete columns south of Arakawa-C longitude domain
    mask_v(:, looky_TooSouth)       = []; % delete columns south of Arakawa-C longitude domain
    z_v(:, looky_TooSouth, :)       = []; % delete columns south of Arakawa-C longitude domain
end % (~isempty(looky_TooSouth))

looky_TooWest	= find(lon_v(:, 1) < WestBoundary); % NOTE: assumes longitudes run down rows; 
if ~isempty(looky_TooWest)
    lat_v(looky_TooWest, :)         = []; % delete rows west of Arakawa-C longitude domain
    lon_v(looky_TooWest, :)         = []; % delete rows west of Arakawa-C longitude domain
    x_v(looky_TooWest, :)           = []; % delete rows west of Arakawa-C longitude domain
    y_v(looky_TooWest, :)           = []; % delete rows west of Arakawa-C longitude domain
    mask_v(looky_TooWest, :)        = []; % delete rows west of Arakawa-C longitude domain
    z_v(looky_TooWest, :, :)        = []; % delete rows west of Arakawa-C longitude domain
end % (~isempty(looky_TooWest))

looky_TooEast	= find(lon_v(:, 1) > EastBoundary); % NOTE: assumes longitudes run down rows; 
if ~isempty(looky_TooEast)
    lat_v(looky_TooEast, :)         = []; % delete rows east of Arakawa-C longitude domain
    lon_v(looky_TooEast, :)         = []; % delete rows east of Arakawa-C longitude domain
    x_v(looky_TooEast, :)           = []; % delete rows east of Arakawa-C longitude domain
    y_v(looky_TooEast, :)           = []; % delete rows east of Arakawa-C longitude domain
    mask_v(looky_TooEast, :)        = []; % delete rows east of Arakawa-C longitude domain
	z_v(looky_TooEast, :, :)        = []; % delete rows east of Arakawa-C longitude domain
end % (~isempty(looky_TooEast))
% -----


% u terms
looky_TooNorth = find(lat_u(1, :) > NorthBoundary); % NOTE: assumes latitudes run across columns; 
if ~isempty(looky_TooNorth)
    lat_u(:, looky_TooNorth)        = []; % delete columns north of Arakawa-C longitude domain
    lon_u(:, looky_TooNorth)        = []; % delete columns north of Arakawa-C longitude domain
    x_u(:, looky_TooNorth)          = []; % delete columns north of Arakawa-C longitude domain
    y_u(:, looky_TooNorth)          = []; % delete columns north of Arakawa-C longitude domain
    mask_u(:, looky_TooNorth)       = []; % delete columns north of Arakawa-C longitude domain
    z_u(:, looky_TooNorth, :)       = []; % delete columns north of Arakawa-C longitude domain
end % (~isempty(looky_TooNorth))
    
looky_TooSouth = find(lat_u(1, :) < SouthBoundary); % NOTE: assumes latitudes run across columns; 
if ~isempty(looky_TooSouth)
    lat_u(:, looky_TooSouth)        = []; % delete columns south of Arakawa-C longitude domain
    lon_u(:, looky_TooSouth)        = []; % delete columns south of Arakawa-C longitude domain
    x_u(:, looky_TooSouth)          = []; % delete columns south of Arakawa-C longitude domain
    y_u(:, looky_TooSouth)          = []; % delete columns south of Arakawa-C longitude domain
    mask_u(:, looky_TooSouth)       = []; % delete columns south of Arakawa-C longitude domain
	z_u(:, looky_TooSouth, :)       = []; % delete columns south of Arakawa-C longitude domain
end % (~isempty(looky_TooSouth))

looky_TooWest	= find(lon_u(:, 1) < WestBoundary); % NOTE: assumes longitudes run down rows; 
if ~isempty(looky_TooWest)
    lat_u(looky_TooWest, :)         = []; % delete rows west of Arakawa-C longitude domain
    lon_u(looky_TooWest, :)         = []; % delete rows west of Arakawa-C longitude domain
    x_u(looky_TooWest, :)           = []; % delete rows west of Arakawa-C longitude domain
    y_u(looky_TooWest, :)           = []; % delete rows west of Arakawa-C longitude domain
    mask_u(looky_TooWest, :)        = []; % delete rows west of Arakawa-C longitude domain
    z_u(looky_TooWest, :, :)        = []; % delete rows west of Arakawa-C longitude domain
end % (~isempty(looky_TooWest))

looky_TooEast	= find(lon_u(:, 1) > EastBoundary); % NOTE: assumes longitudes run down rows; 
if ~isempty(looky_TooEast)
    lat_u(looky_TooEast, :)         = []; % delete rows east of Arakawa-C longitude domain
    lon_u(looky_TooEast, :)         = []; % delete rows east of Arakawa-C longitude domain
    x_u(looky_TooEast, :)           = []; % delete rows east of Arakawa-C longitude domain
    y_u(looky_TooEast, :)           = []; % delete rows east of Arakawa-C longitude domain
    mask_u(looky_TooEast, :)        = []; % delete rows east of Arakawa-C longitude domain
	z_u(looky_TooEast, :, :)        = []; % delete rows east of Arakawa-C longitude domain
end % (~isempty(looky_TooEast))
% -----


% psi terms
looky_TooNorth = find(lat_psi(1, :) > NorthBoundary); % NOTE: assumes latitudes run across columns; 
if ~isempty(looky_TooNorth)
    lat_psi(:, looky_TooNorth)   	= []; % delete columns north of Arakawa-C longitude domain
    lon_psi(:, looky_TooNorth)   	= []; % delete columns north of Arakawa-C longitude domain
    x_psi(:, looky_TooNorth)     	= []; % delete columns north of Arakawa-C longitude domain
    y_psi(:, looky_TooNorth)     	= []; % delete columns north of Arakawa-C longitude domain
    mask_psi(:, looky_TooNorth)   	= []; % delete columns north of Arakawa-C longitude domain
    z_psi(:, looky_TooNorth, :)   	= []; % delete columns north of Arakawa-C longitude domain
end % (~isempty(looky_TooNorth))
    
looky_TooSouth = find(lat_psi(1, :) < SouthBoundary); % NOTE: assumes latitudes run across columns; 
if ~isempty(looky_TooSouth)
    lat_psi(:, looky_TooSouth)    	= []; % delete columns south of Arakawa-C longitude domain
    lon_psi(:, looky_TooSouth)    	= []; % delete columns south of Arakawa-C longitude domain
    x_psi(:, looky_TooSouth)       	= []; % delete columns south of Arakawa-C longitude domain
    y_psi(:, looky_TooSouth)      	= []; % delete columns south of Arakawa-C longitude domain
    mask_psi(:, looky_TooSouth)    	= []; % delete columns south of Arakawa-C longitude domain
	z_psi(:, looky_TooSouth, :)   	= []; % delete columns south of Arakawa-C longitude domain
end % (~isempty(looky_TooSouth))

looky_TooWest	= find(lon_psi(:, 1) < WestBoundary); % NOTE: assumes longitudes run down rows; 
if ~isempty(looky_TooWest)
    lat_psi(looky_TooWest, :)     	= []; % delete rows west of Arakawa-C longitude domain
    lon_psi(looky_TooWest, :)      	= []; % delete rows west of Arakawa-C longitude domain
    x_psi(looky_TooWest, :)       	= []; % delete rows west of Arakawa-C longitude domain
    y_psi(looky_TooWest, :)       	= []; % delete rows west of Arakawa-C longitude domain
    mask_psi(looky_TooWest, :)    	= []; % delete rows west of Arakawa-C longitude domain
	z_psi(looky_TooWest, :, :)    	= []; % delete rows west of Arakawa-C longitude domain
end % (~isempty(looky_TooWest))

looky_TooEast	= find(lon_psi(:, 1) > EastBoundary); % NOTE: assumes longitudes run down rows; 
if ~isempty(looky_TooEast)
    lat_psi(looky_TooEast, :)   	= []; % delete rows east of Arakawa-C longitude domain
    lon_psi(looky_TooEast, :)    	= []; % delete rows east of Arakawa-C longitude domain
    x_psi(looky_TooEast, :)       	= []; % delete rows east of Arakawa-C longitude domain
    y_psi(looky_TooEast, :)       	= []; % delete rows east of Arakawa-C longitude domain
    mask_psi(looky_TooEast, :)     	= []; % delete rows east of Arakawa-C longitude domain
	z_psi(looky_TooEast, :, :)    	= []; % delete rows east of Arakawa-C longitude domain
end % (~isempty(looky_TooEast))
% -------------------------------------------------------------------------


% step 3c: get final matrix sizes -----------------------------------------
[rows_rho, clms_rho]	= size(lat_rho);	% rows = w longitudes (rho=v, 51); clms = w latitudes (rho=u, 81)
[~, ~, num_z]           = size(z_rho);
[rows_h, clms_h]        = size(h);          % rows = longitudes (h=rho); clms = latitudes (h=rho)
[~, ~, num_w]           = size(z_w);        % num_w = num_z + 1;
[rows_v, clms_v]        = size(lat_v);      % rows = v longitudes (v=rho, 51); clms = v latitudes (v=rho-1, u-1, 80)
[rows_u, clms_u]        = size(lat_u);      % rows = u longitudes (u=rho+1, v+1, 50); clms = latitudes (u=rho, 81)

if (rows_rho ~= rows_h) || (clms_rho ~= clms_h)
    error('rho & h matrices are different sizes')
end
% -------------------------------------------------------------------------


% step 3d: clear temporary variables --------------------------------------
clear looky_*
% *************************************************************************





% *************************************************************************
% STEP 4: calculate ROMS cell dimensions & volume fluxes-------------------
% step 4a: calculate ROMS cell dimensions (m) from psi points -------------
% distances relative to the globe, due N<->S & E<->W (to be used for volume flux rate calculations)
dist_NS_globe	= y_psi(:, 2:end) - y_psi(:, 1:(end-1)); % N<-->S lengths; (m); abs(NORTH minus SOUTH); (2D matrix: u longitude (50) X v latitude-1 (79));
dist_EW_globe	= x_psi(2:end, :) - x_psi(1:(end-1), :); % W<-->E lengths; (m); abs(EAST minus WEST); (2D matrix: u longitude-1 (49) X v latitude (80));

% distances that are face lengths (to be used for box volume calculations)
dist_NS_face	= sqrt(dist_NS_globe.^2 + (x_psi(:, 2:end) - x_psi(:, 1:(end-1))).^2); % N<-->S face lengths; (m); (2D matrix: u longitude (50) X v latitude-1 (79)); NOTE: x_psi(south) - x_psi(north)
dist_EW_face	= sqrt(dist_EW_globe.^2 + (y_psi(2:end, :) - y_psi(1:(end-1), :)).^2); % E<-->W face lengths; (m); (2D matrix: u longitude-1 (49) X v latitude (80)); NOTE: y_psi(east) - y_psi(west)
% -------------------------------------------------------------------------


% step 4b: calculate depths of top & bottom bounds of each cell face ------
%          NOTE: z_w are the depths @ each horizontal depth-layer face at the w point (top & bottom face depths at the cell center)
%          NOTE: z_wv is the depth of each horizontal depth-layer face at the v point
%          NOTE: z_wu is the depth of each horizontal depth-layer face at the u point
h_v             = (h(:, 2:end) + h(:, 1:(end-1))) / 2;	 % bathymetry (h) @ v-point; (m); (2D matrix: v longitude (51) X v latitude (80))
h_u             = (h(2:end, :) + h(1:(end-1), :)) / 2;	 % bathymetry (h) @ u-point; (m); (2D matrix: u longitude (50) X u latitude (81))
h_psi           = (h_v(2:end, :) + h_v(1:(end-1), :)) / 2; % bathymetry (h) @ psi-point; (m); (2D matrix: psi longitude (50) X psi latitude (80)); NOTE: calculated from v points (confirmed to be same as h when calculated from u points)
z_wv            = (z_w(:, 2:end, :) + z_w(:, 1:(end-1), :)) / 2;	% depth of w @ v-point; (m); (3D matrix: v longitude (51) X v latitude (80) X num_w (43))
                                                                    % NOTE: z_wv top layer = h_v (bathymetry @ v) (JR confirmed this)
z_wu            = (z_w(2:end, :, :) + z_w(1:(end-1), :, :)) / 2;	% depth of w @ u-point; (m); (3D matrix: u longitude (50) X u latitude (81) X num_w (43))
                                                                    % NOTE: z_wu top layer = h_u (bathymetry @ v) (JR confirmed this at t=1)
z_wpsi          = (z_wv(2:end, :, :) + z_wv(1:(end-1), :, :)) / 2;  % depth of w @ psi-point; (m); (3D matrix: psi longitude (50) X psi latitude (80) X num_w (43)); NOTE: calculated from v points (confirmed to be same as h when calculated from u points)
% -------------------------------------------------------------------------


% step 4c: calculate ROMS cell-face heights -------------------------------
v_height        = z_wv(:, :, 1:(end-1))   - z_wv(:, :, 2:end);   % height of cell @ v-point; (bottom minus top);   (m); (3D matrix: v longitude (51) X v latitude (80) X num_z (42)); NOTE: JR confirms sum of these heights = bathymetry @ v
u_height        = z_wu(:, :, 1:(end-1))   - z_wu(:, :, 2:end);   % height of cell @ u-point; (bottom minus top);   (m); (3D matrix: u longitude (50) X u latitude (81) X num_z (42))
w_height        = z_w(:, :, 1:(end-1))    - z_w(:, :, 2:end);    % height of cell @ center; (bottom minus top);    (m); (3D matrix: rho longitude (51) X rho latitude (81) X num_z (42))
psi_height      = z_wpsi(:, :, 1:(end-1)) - z_wpsi(:, :, 2:end); % height of cell @ psi-point; (bottom minus top); (m); (3D matrix: psi longitude (50) X psi latitude (80) X num_z (42))
% -------------------------------------------------------------------------


% step 4d: calculate lateral cell face area -------------------------------
%          NOTE: rectangular face based on EW & NS distances relative to the globe
%          NOTE: trapezoidal face doesn't work
repmat_dist_EW_globe	= repmat(dist_EW_globe, [1, 1, num_z]); % (m); (3D matrix: v longitude-2 (49) X psi & v latitude (80), num_z (42));
area_EW_globe           = repmat_dist_EW_globe .* v_height(2:(end-1), :, :); % E<-->W ROMS cell area; (m2); (3D matrix: v longitude-2 (49) X v latitude (80) X num_z (42)); NOTE: omitting west & east border v_height outside of ROMS cells

repmat_dist_NS_globe	= repmat(dist_NS_globe, [1, 1, num_z]); % (m); (3D matrix: psi & u longitude (50) X u latitude-2 (79), num_z (42));
area_NS_globe         	= repmat_dist_NS_globe .* u_height(:, 2:(end-1), :); % N<-->S ROMS cell area; (m2); (3D matrix: u longitude (50) X u latitude-2 (79) X num_z (42)); NOTE: omitting south & north border u_height outside of ROMS cells
% -------------------------------------------------------------------------


% step 4e: calculate horizontal (top & bottom) cell face area -------------
%          NOTE: trapezoidal face

% floor area based on distances relative to the globe, due N<->S & E<->W (to be used for volume flux rate calculations)
dist_EW_mean            = (dist_EW_globe(:, 2:end) + dist_EW_globe(:, 1:(end-1))) / 2; % mean E<-->W distance; mean(north EW dist + south EW dist); (m); (2D matrix: psi longitude-1 (49) X psi latitude-1 (79))
repmat_dist_EW_globe	= repmat(dist_EW_mean, [1, 1, num_w]);          % E<-->W distance; (m); (3D matrix: psi longitude-1 (49) X psi latitude-1 (79), num_w (43));
dist_NS_mean            = (dist_NS_globe(2:end, :) + dist_NS_globe(1:(end-1), :)) / 2; % mean N<-->S distance; (m); mean(east NS dist + west NS dist); (2D matrix: psi longitude-1 (49) X psi latitude-1 (79))
repmat_dist_NS_globe	= repmat(dist_NS_mean, [1, 1, num_w]);          % N<-->S distance; (m); (3D matrix: psi longitude-1 (49) X psi latitude-1 (79), num_w (43));
area_floor_globe        = repmat_dist_NS_globe .* repmat_dist_EW_globe;	% cell floor area based on EW & NS distances relative to the globe; (m2); (3D matrix: rho longitude-2 (49) X rho latitude-2 (79), num_w (43)); QQQ try using the area_floor_face instead & compare error level

% floor area based on distances that are face lengths (to be used for box volume calculations)
dist_EW_mean            = (dist_EW_face(:, 2:end) + dist_EW_face(:, 1:(end-1))) / 2; % mean E<-->W distance; mean(north EW dist + south EW dist); (m); (2D matrix: psi longitude-1 (49) X psi latitude-1 (79))
repmat_dist_EW_face     = repmat(dist_EW_mean, [1, 1, num_w]);          % E<-->W distance; (m); (3D matrix: psi longitude-1 (49) X psi latitude-1 (79), num_w (43));
dist_NS_mean            = (dist_NS_face(2:end, :) + dist_NS_face(1:(end-1), :)) / 2; % mean N<-->S distance; (m); mean(east NS dist + west NS dist); (2D matrix: psi longitude-1 (49) X psi latitude-1 (79))
repmat_dist_NS_face     = repmat(dist_NS_mean, [1, 1, num_w]);          % N<-->S distance; (m); (3D matrix: psi longitude-1 (49) X psi latitude-1 (79), num_w (43));
area_floor_face         = repmat_dist_NS_face .* repmat_dist_EW_face;	% cell floor area based on EW & NS distances relative to face lengths; (m2); (3D matrix: rho longitude-2 (49) X rho latitude-2 (79), num_w (43));
% -------------------------------------------------------------------------


% step 4f: clear temporary variables --------------------------------------
clear repmat_dist_EW_globe repmat_dist_NS_globe
clear dist_EW_mean dist_NS_mean
% *************************************************************************





% *************************************************************************
% STEP 5: generate ECOTRAN grid geometry in terms of ROMS cells------------
%       NOTES: grid_addresses define the south, north, west, & east addresses of each ECOTRAN domain within each of the ROMS matrices
%              grid_addresses_rho & grid_addresses_v are identical for WEST & EAST & NORTH (rho SOUTH is +1 higher) (true for main grid and longitudinal boundary cells)
%              grid_addresses_rho & grid_addresses_u are identical for SOUTH & NORTH & EAST (rho WEST is +1 higher for main grid) (rho WEST is +1 higher for offshore boundary; rho EAST is +1 higher for offshore boundary) (rho & u are identical for SOUTH & NORTH & WEST & EAST for inshore boundary)

% step 5a: initialize and shape variables ---------------------------------
grid_addresses_v	= zeros(num_domains, 7); % initialize; (2D matrix: num_domains X 7); [south north west east (num cells lat) (num cells lon) (num cells total)]
grid_addresses_u	= zeros(num_domains, 7); % initialize; (2D matrix: num_domains X 7); [south north west east (num cells lat) (num cells lon) (num cells total)]
grid_addresses_rho	= zeros(num_domains, 7); % initialize; (2D matrix: num_domains X 7); [south north west east (num cells lat) (num cells lon) (num cells total)]
bathymetry_range    = zeros(num_domains, 3); % initialize; (2D matrix: num_domains X 3); [max min (num cells)]

south	= 1; % grid_addresses column definitions
north	= 2;
west	= 3;
east	= 4;

% trim bounding cells from lat & lon terms
lat_v_trim      = lat_v(2:(end-1), :); % prune EW borders to align with agg_v; (2D matrix: v longitude-2 (49) X v latitude (80));
lon_v_trim      = lon_v(2:(end-1), :); % prune EW borders to align with agg_v; (2D matrix: v longitude-2 (49) X v latitude (80));
lat_u_trim      = lat_u(:, 2:(end-1)); % prune NW borders to align with agg_u; (2D matrix: u longitude (50) X u latitude-2 (79));
lon_u_trim      = lon_u(:, 2:(end-1)); % prune NW borders to align with agg_u; (2D matrix: u longitude (50) X u latitude-2 (79));
lat_rho_trim	= lat_rho(2:(end-1), 2:(end-1)); % prune EW & NW borders to align with agg_w; (2D matrix: rho longitude-2 (49) X rho latitude-2 (79));
lon_rho_trim	= lon_rho(2:(end-1), 2:(end-1)); % prune EW & NW borders to align with agg_w; (2D matrix: rho longitude-2 (49) X rho latitude-2 (79));
h_trim          = h(2:(end-1), 2:(end-1)); % prune EW & NW borders to align with agg_w; (2D matrix: rho longitude-2 (49) X rho latitude-2 (79));
% -------------------------------------------------------------------------

for domain_loop = 1:num_domains
    
    current_north       = domain_definition(domain_loop, 1);
	current_south       = domain_definition(domain_loop, 2);
	current_shallow     = domain_definition(domain_loop, 3);
	current_deep        = domain_definition(domain_loop, 4);

    % step 5b: find North/South boundaries of current domian --------------
    % latitude domain in v
    looky_north         = find(lat_v_trim(1, :) <= current_north); % column address in lat_v at or south of northern bound; (horizontal vector: 1 X ???)
    looky_south     	= find(lat_v_trim(1, :) >= current_south); % column address in lat_v at or north of southern bound; (horizontal vector: 1 X ???)
    grid_addresses_v(domain_loop, south:north)	= [min(looky_south), max(looky_north)]; % [south north west east];
    
    % latitude domain in u (rho & u latitudes should be the same)
    looky_north         = find(lat_u_trim(1, :) <= current_north); % column address at or south of northern bound; (horizontal vector: 1 X ???)
    looky_south         = find(lat_u_trim(1, :) >= current_south); % column address at or north of southern bound; (horizontal vector: 1 X ???)
    grid_addresses_u(domain_loop, south:north) = [min(looky_south) max(looky_north)]; % [south north west east];
    
    % latitude domain in rho (rho & u latitudes should be the same)
    looky_north         = find(lat_rho_trim(1, :) <= current_north); % column address at or south of northern bound; (horizontal vector: 1 X ???)
    looky_south         = find(lat_rho_trim(1, :) >= current_south); % column address at or north of southern bound; (horizontal vector: 1 X ???)
    grid_addresses_rho(domain_loop, south:north) = [min(looky_south) max(looky_north)]; % [south north west east ? ?];
    % ---------------------------------------------------------------------
    
    
    % step 5c: find East/West boundaries based on bathymetry (h) ----------
    current_addresses_rho       = grid_addresses_rho(domain_loop, :); % (horizontal vector: 1 X 6)
    current_latBand_h           = h_trim(:, current_addresses_rho(south):current_addresses_rho(north)); % (2D matrix: rows_rho X (num ROMS cells in lat zone)); top row = west, bottom row = east;
    
    median_h                    = median(current_latBand_h, 2) * (-1); % median bathymetry along each longitude band within the current latitude zone; (vertical vector: rows_h X 1); NOTE: median seems to work better for NCC than mean or mode (reduces on land cells)
    looky_east                  = find(median_h <= current_shallow); % (vertical vector: number_of_qualifying_deep_depths X 1);
    looky_west                  = find(median_h >= current_deep); % (vertical vector: number_of_qualifying_shallow_depths X 1);
    
    grid_addresses_rho(domain_loop, west:east) = [min(looky_west), max(looky_east)]; % (rho & v longitudes should be the same)
    
    current_west_rho            = lon_rho_trim(grid_addresses_rho(domain_loop, west)); % longitude; (scalar)
    current_east_rho            = lon_rho_trim(grid_addresses_rho(domain_loop, east)); % longitudes; (scalar)

	% longitude domain in v (rho & v longitudes should be the same)
    looky_east                  = find(lon_v_trim(:, 1) <= current_east_rho); % column address at or south of northern bound; (horizontal vector: 1 X ???); NOTE: DON'T need to add +1 to this since rho and v lie along the same longitudes
    looky_west                  = find(lon_v_trim(:, 1) >= current_west_rho); % column address at or south of northern bound; (horizontal vector: 1 X ???); NOTE: DON'T need to subtract -1 to this since rho and v lie along the same longitudes
    grid_addresses_v(domain_loop, west:east)	= [min(looky_west), max(looky_east)]; % [south north west east]

	% longitude domain in u
    looky_east                  = find(lon_u_trim(:, 1) <= current_east_rho); % column address at or south of northern bound; (horizontal vector: 1 X ???) NOTE: need to add +1 to this since rho are internal points of boxes
    looky_west                  = find(lon_u_trim(:, 1) >= current_west_rho); % column address at or south of northern bound; (horizontal vector: 1 X ???) NOTE: need to subtract -1 to this since rho are internal points of boxes
    grid_addresses_u(domain_loop, west:east) = [min(looky_west)-1, max(looky_east)+1]; % [south north west east]
	% ---------------------------------------------------------------------
    
    
	% step 5d: record bathymetry (h) range for grid evaluation ------------
    bathymetry_range(domain_loop, 1) = min(min(h_trim(grid_addresses_rho(domain_loop, west):grid_addresses_rho(domain_loop, east), grid_addresses_rho(domain_loop, south):grid_addresses_rho(domain_loop, north)))) * (-1);
	bathymetry_range(domain_loop, 2) = max(max(h_trim(grid_addresses_rho(domain_loop, west):grid_addresses_rho(domain_loop, east), grid_addresses_rho(domain_loop, south):grid_addresses_rho(domain_loop, north)))) * (-1);
	bathymetry_range(domain_loop, 3) = numel(h_trim(grid_addresses_rho(domain_loop,   west):grid_addresses_rho(domain_loop, east), grid_addresses_rho(domain_loop, south):grid_addresses_rho(domain_loop, north)));
    % ---------------------------------------------------------------------
    
    
	% step 5e: record number of cells in the ECOTRAN grid -----------------
    grid_addresses_rho(domain_loop, 7)	= numel(h_trim(grid_addresses_rho(domain_loop, west):grid_addresses_rho(domain_loop, east), grid_addresses_rho(domain_loop, south):grid_addresses_rho(domain_loop, north)));
    % ---------------------------------------------------------------------
    
end % end domain_loop
% -------------------------------------------------------------------------


% step 5f: find inshore & offshore bounding ROMS cells for each latitude band
%          NOTE: rho & v longitude grid_addresses should have EAST = WEST and clm 6 = 1 (rho & v longitudes are in the MIDDLE of the boundary cell so EAST = WEST longitude address)
unique_latitude_bands_rho	= unique(grid_addresses_rho(:, south:north), 'rows'); % (2D matrix: num_latitude_bands X 2)
unique_latitude_bands_v     = unique(grid_addresses_v(:,   south:north), 'rows'); % (2D matrix: num_latitude_bands X 2)
unique_latitude_bands_u     = unique(grid_addresses_u(:,   south:north), 'rows'); % (u & rho latitude should be the same); (2D matrix: num_latitude_bands X 2)

[num_latitude_bands, ~]     = size(unique_latitude_bands_rho);
num_EWbounds                = num_latitude_bands * 2; % number of ECOTRAN East & West external boundary cells

grid_addresses_EWbound_rho	= zeros((num_latitude_bands*2), 7); % initialize; (2D matrix: (num_latitude_bands*2) X 7); [south north west east (num cells lat) (num cells lon) (num cells total)]
grid_addresses_EWbound_u	= zeros((num_latitude_bands*2), 7); % initialize; (2D matrix: (num_latitude_bands*2) X 7); [south north west east (num cells lat) (num cells lon) (num cells total)]
grid_addresses_EWbound_v	= zeros((num_latitude_bands*2), 7); % initialize; (2D matrix: (num_latitude_bands*2) X 7); [south north west east (num cells lat) (num cells lon) (num cells total)]

% find westmost and east most ROMS cells that will bound the ECOTRAN grid
global_westmost_rho     	= min(grid_addresses_rho(:, west)) - 1; % NOTE: rho & v & u are the same; QQQ add test that these addresses remain within the provided ROMS grid
global_westmost_u       	= min(grid_addresses_u(:,   west)) - 1;
global_westmost_v         	= min(grid_addresses_v(:,   west)) - 1;

global_eastmost_rho      	= max(grid_addresses_rho(:, east)) + 1; % NOTE: rho & v & u are the same
global_eastmost_u       	= max(grid_addresses_u(:,   east)) + 1;
global_eastmost_v       	= max(grid_addresses_v(:,   east)) + 1;

for latitude_loop = 1:num_latitude_bands
    
    current_latitude_band_rho	= unique_latitude_bands_rho(latitude_loop, :); % (u & rho latitude should be the same)
    current_latitude_band_u     = unique_latitude_bands_u(latitude_loop, :);   % (u & rho latitude should be the same)
    current_latitude_band_v     = unique_latitude_bands_v(latitude_loop, :);
    
	looky_latitude_band_rho     = find(ismember(grid_addresses_rho(:, south:north), current_latitude_band_rho, 'rows')); % (vertical vector: number of ECOTRAN domains in current_latitude_band_rho)
    looky_latitude_band_u       = find(ismember(grid_addresses_u(:, south:north),   current_latitude_band_u, 'rows'));   % (vertical vector: number of ECOTRAN domains in current_latitude_band_u); (should be same as looky_latitude_band_rho)
    looky_latitude_band_v       = find(ismember(grid_addresses_v(:, south:north),   current_latitude_band_v, 'rows'));   % (vertical vector: number of ECOTRAN domains in current_latitude_band_v); (should be same as looky_latitude_band_rho)
    
    westmost_rho                = min(grid_addresses_rho(looky_latitude_band_rho, west)); % (QQQ rho & v are the same???)
    westmost_u                  = min(grid_addresses_u(looky_latitude_band_u, west));
    westmost_v                  = min(grid_addresses_v(looky_latitude_band_v, west));
    
    eastmost_rho                = max(grid_addresses_rho(looky_latitude_band_rho, east)); % (QQQ rho & v & u are the same???)
    eastmost_u                  = max(grid_addresses_u(looky_latitude_band_u, east));
    eastmost_v                  = max(grid_addresses_v(looky_latitude_band_v, east));
    
    BoundingCell_west_rho       = [global_westmost_rho (westmost_rho-1)]; % (QQQ rho & v are the same???)
    BoundingCell_west_u         = [global_westmost_u    westmost_u];
	BoundingCell_west_v         = [global_westmost_v   (westmost_v-1)]; % (QQQ rho & v are the same???)

    BoundingCell_east_rho       = [(eastmost_rho+1)     global_eastmost_rho]; % (QQQ rho & v & u are the same???)
	BoundingCell_east_u         = [ eastmost_u          global_eastmost_u];
    BoundingCell_east_v         = [(eastmost_v+1)       global_eastmost_v]; % (QQQ rho & v are the same???)

    % offshore boundary
    grid_addresses_EWbound_rho(latitude_loop, 1:4)	= [current_latitude_band_rho BoundingCell_west_rho];
    grid_addresses_EWbound_u(latitude_loop, 1:4)	= [current_latitude_band_u   BoundingCell_west_u];
    grid_addresses_EWbound_v(latitude_loop, 1:4)	= [current_latitude_band_v   BoundingCell_west_v];
    
    % inshore boundary
	grid_addresses_EWbound_rho((num_latitude_bands+latitude_loop), 1:4)	= [current_latitude_band_rho BoundingCell_east_rho];
	grid_addresses_EWbound_u((num_latitude_bands+latitude_loop), 1:4)	= [current_latitude_band_u BoundingCell_east_u];
	grid_addresses_EWbound_v((num_latitude_bands+latitude_loop), 1:4)	= [current_latitude_band_v BoundingCell_east_v];

end % (latitude_loop)

% QQQ add some check to make sure boundary cells are still within the provided ROMS grid
% -------------------------------------------------------------------------


% step 5g: find north & south ROMS boundary cells -------------------------
northmost_rho           = max(grid_addresses_rho(:, north));
northmost_v             = max(grid_addresses_v(:, north));
northmost_u             = max(grid_addresses_u(:, north));

southmost_rho           = min(grid_addresses_rho(:, south));
southmost_v             = min(grid_addresses_v(:, south));
southmost_u             = min(grid_addresses_u(:, south));

looky_northmost_rho     = find(grid_addresses_rho(:, north) == northmost_rho);
looky_northmost_v       = find(grid_addresses_v(:, north) == northmost_v);
looky_northmost_u       = find(grid_addresses_u(:, north) == northmost_u);

looky_southmost_rho     = find(grid_addresses_rho(:, south) == southmost_rho);
looky_southmost_v       = find(grid_addresses_v(:, south) == southmost_v);
looky_southmost_u       = find(grid_addresses_u(:, south) == southmost_u);

num_NSbounds            = length(looky_northmost_rho) + length(looky_southmost_rho); % number of ECOTRAN North & South external boundary cells

% NOTE: rho & u are in middle of NS boundary (their North grid address will have the same as their South grid address)
grid_addresses_Nbound_rho	= grid_addresses_rho(looky_northmost_rho, :);
grid_addresses_Nbound_rho(:, south:north)	= northmost_rho + 1;
grid_addresses_Sbound_rho	= grid_addresses_rho(looky_southmost_rho, :);
grid_addresses_Sbound_rho(:, south:north)	= southmost_rho - 1;
grid_addresses_NSbound_rho	= [grid_addresses_Nbound_rho; grid_addresses_Sbound_rho];

grid_addresses_Nbound_u     = grid_addresses_u(looky_northmost_u, :);
grid_addresses_Nbound_u(:, south:north)     = northmost_u + 1;
grid_addresses_Sbound_u     = grid_addresses_u(looky_southmost_u, :);
grid_addresses_Sbound_u(:, south:north)     = southmost_u - 1;
grid_addresses_NSbound_u	= [grid_addresses_Nbound_u; grid_addresses_Sbound_u];

grid_addresses_Nbound_v     = grid_addresses_v(looky_northmost_v, :);
grid_addresses_Nbound_v(:, north)           = northmost_v + 1;
grid_addresses_Nbound_v(:, south)           = northmost_v;
grid_addresses_Sbound_v     = grid_addresses_v(looky_southmost_v, :);
grid_addresses_Sbound_v(:, north)           = southmost_v;
grid_addresses_Sbound_v(:, south)           = southmost_v - 1;
grid_addresses_NSbound_v	= [grid_addresses_Nbound_v; grid_addresses_Sbound_v];
% -------------------------------------------------------------------------


% step 5h: append grid_addresses_boundary to grid_addresses ---------------
grid_addresses_rho          = [grid_addresses_rho; grid_addresses_EWbound_rho; grid_addresses_NSbound_rho];
grid_addresses_u            = [grid_addresses_u;   grid_addresses_EWbound_u;   grid_addresses_NSbound_u];
grid_addresses_v            = [grid_addresses_v;   grid_addresses_EWbound_v;   grid_addresses_NSbound_v];
% -------------------------------------------------------------------------


% step 5i: get ROMS lat, lon, & cell count for ECOTRAN grid ---------------
grid_addresses_rho(:, 5)	= grid_addresses_rho(:, north) - grid_addresses_rho(:, south) + 1;  % number of ROMS cells latitudinally
grid_addresses_u(:, 5)      = grid_addresses_u(:, north)   - grid_addresses_u(:, south)   + 1;  % number of ROMS cells latitudinally
grid_addresses_v(:, 5)      = grid_addresses_v(:, north)   - grid_addresses_v(:, south);        % number of ROMS cells latitudinally

grid_addresses_rho(:, 6)	= grid_addresses_rho(:, east) - grid_addresses_rho(:, west) + 1;	% number of ROMS cells longitudinally
grid_addresses_u(:, 6)      = grid_addresses_u(:, east)   - grid_addresses_u(:, west);          % number of ROMS cells longitudinally
grid_addresses_v(:, 6)      = grid_addresses_v(:, east)   - grid_addresses_v(:, west)   + 1;	% number of ROMS cells longitudinally

grid_addresses_rho((num_domains+1):end, 7)	= grid_addresses_rho((num_domains+1):end, 5) .* grid_addresses_rho((num_domains+1):end, 6); % number of ROMS cells; NOTE: count for the internal ECOTRAN domain was already calculated at step 4d)
grid_addresses_u(:, 7)      = grid_addresses_u(:, 5)   .* grid_addresses_u(:, 6);               % number of ROMS cells
grid_addresses_v(:, 7)      = grid_addresses_v(:, 5)   .* grid_addresses_v(:, 6);            	% number of ROMS cells

ECOTRAN_grid_lat_rho        = [lat_rho_trim(1, grid_addresses_rho(:, south))' lat_rho_trim(1, grid_addresses_rho(:, north))'];	% use lat_v or lat_rho
ECOTRAN_grid_lat_u          = [lat_u_trim(1, grid_addresses_u(:, south))'     lat_u_trim(1, grid_addresses_u(:, north))'];
ECOTRAN_grid_lat_v          = [lat_v_trim(1, grid_addresses_v(:, south))'     lat_v_trim(1, grid_addresses_v(:, north))'];        % use lat_v or lat_rho

ECOTRAN_grid_lon_rho        = [lon_rho_trim(grid_addresses_rho(:, west), 1) lon_rho_trim(grid_addresses_rho(:, east), 1)];
ECOTRAN_grid_lon_u          = [lon_u_trim(grid_addresses_u(:, west), 1)     lon_u_trim(grid_addresses_u(:, east), 1)];            % use lon_u
ECOTRAN_grid_lon_v          = [lon_v_trim(grid_addresses_v(:, west), 1)     lon_v_trim(grid_addresses_v(:, east), 1)];

ECOTRAN_gridCoordinates     = [lat_v_trim(1, grid_addresses_v(:, north))' lat_v_trim(1, grid_addresses_v(:, south))' lon_u_trim(grid_addresses_u(:, west), 1) lon_u_trim(grid_addresses_u(:, east), 1)]; % [north south west east]
% -------------------------------------------------------------------------


% step 5j: clear temporary variables --------------------------------------
clear current_* looky_*
clear mean_h median_h mode_h
clear grid_addresses_Nbound_* grid_addresses_Sbound_* southmost* northmost*
clear global_* unique_* BoundingCell*
% *************************************************************************





% *************************************************************************
% STEP 6: aggregate & organize ROMS cells into ECOTRAN grid----------------
%         ROMS cells are defined by rho-points
% step 6a: initialize variables -------------------------------------------
num_domainsWboundary            = num_domains + num_EWbounds + num_NSbounds;

overlapAddress_S_rho            = zeros(num_domainsWboundary); % initialize; (2D matrix: num_domainsWboundary X num_domainsWboundary)
overlapAddress_N_rho            = zeros(num_domainsWboundary); % initialize; (2D matrix: num_domainsWboundary X num_domainsWboundary)
overlapAddress_W_rho            = zeros(num_domainsWboundary); % initialize; (2D matrix: num_domainsWboundary X num_domainsWboundary)
overlapAddress_E_rho            = zeros(num_domainsWboundary); % initialize; (2D matrix: num_domainsWboundary X num_domainsWboundary)

overlapAddress_S_v              = zeros(num_domainsWboundary); % initialize; (2D matrix: num_domainsWboundary X num_domainsWboundary)
overlapAddress_N_v              = zeros(num_domainsWboundary); % initialize; (2D matrix: num_domainsWboundary X num_domainsWboundary)
overlapAddress_W_v              = zeros(num_domainsWboundary); % initialize; (2D matrix: num_domainsWboundary X num_domainsWboundary)
overlapAddress_E_v              = zeros(num_domainsWboundary); % initialize; (2D matrix: num_domainsWboundary X num_domainsWboundary)

overlapAddress_S_u              = zeros(num_domainsWboundary); % initialize; (2D matrix: num_domainsWboundary X num_domainsWboundary)
overlapAddress_N_u              = zeros(num_domainsWboundary); % initialize; (2D matrix: num_domainsWboundary X num_domainsWboundary)
overlapAddress_W_u              = zeros(num_domainsWboundary); % initialize; (2D matrix: num_domainsWboundary X num_domainsWboundary)
overlapAddress_E_u              = zeros(num_domainsWboundary); % initialize; (2D matrix: num_domainsWboundary X num_domainsWboundary)

evaluate_latitude_matrix        = zeros(num_domainsWboundary); % initialize; (2D matrix: num_domainsWboundary X num_domainsWboundary)
evaluate_longitude_matrix       = zeros(num_domainsWboundary); % initialize; (2D matrix: num_domainsWboundary X num_domainsWboundary)

evaluate_SouthNeighbor_matrix	= zeros(num_domainsWboundary); % initialize; (2D matrix: num_domainsWboundary X num_domainsWboundary)
evaluate_NorthNeighbor_matrix	= zeros(num_domainsWboundary); % initialize; (2D matrix: num_domainsWboundary X num_domainsWboundary)
evaluate_WestNeighbor_matrix	= zeros(num_domainsWboundary); % initialize; (2D matrix: num_domainsWboundary X num_domainsWboundary)
evaluate_EastNeighbor_matrix	= zeros(num_domainsWboundary); % initialize; (2D matrix: num_domainsWboundary X num_domainsWboundary)
% -------------------------------------------------------------------------


for domain_loop = 1:num_domainsWboundary
    
    % step 6b: ROMS cell addresses of current ECOTRAN domain --------------
    currentBox_S_rho	= grid_addresses_rho(domain_loop, south);	% south rho-point address of current ECOTRAN box; (scalar)
    currentBox_N_rho	= grid_addresses_rho(domain_loop, north);	% north rho-point address of current ECOTRAN box; (scalar)
    currentBox_W_rho	= grid_addresses_rho(domain_loop, west);	% west rho-point address of current ECOTRAN box; (scalar)
    currentBox_E_rho	= grid_addresses_rho(domain_loop, east);	% east rho-point address of current ECOTRAN box; (scalar)
    
    currentBox_S_v      = grid_addresses_v(domain_loop, south);     % south v-point address of current ECOTRAN box; (scalar)
    currentBox_N_v      = grid_addresses_v(domain_loop, north);     % north v-point address of current ECOTRAN box; (scalar)
    currentBox_W_v      = grid_addresses_v(domain_loop, west);      % west v-point address of current ECOTRAN box; (scalar)
    currentBox_E_v      = grid_addresses_v(domain_loop, east);      % east v-point address of current ECOTRAN box; (scalar)
    
    currentBox_S_u      = grid_addresses_u(domain_loop, south);     % south u-point address of current ECOTRAN box; (scalar)
    currentBox_N_u      = grid_addresses_u(domain_loop, north);     % north u-point address of current ECOTRAN box; (scalar)
    currentBox_W_u      = grid_addresses_u(domain_loop, west);      % west u-point address of current ECOTRAN box; (scalar)
    currentBox_E_u      = grid_addresses_u(domain_loop, east);      % east u-point address of current ECOTRAN box; (scalar)
    % ---------------------------------------------------------------------

    
    % step 6c: find overlapping latitude & longitude rho-points -----------
    %          (will be same for v & u)
    evaluate_southward	= grid_addresses_rho(:, north) - currentBox_S_rho; % all potential overlaps are 0 or '+'; (vertical vector: num_domainsWboundary X 1)
    evaluate_southward(evaluate_southward >= 0) = 1; % (vertical vector: num_domainsWboundary X 1)
    evaluate_southward(evaluate_southward < 0)	= 0;
    
	evaluate_northward	= currentBox_N_rho - grid_addresses_rho(:, south); % all potential overlaps are 0 or '+'; (vertical vector: num_domainsWboundary X 1)
    evaluate_northward(evaluate_northward >= 0) = 1; % (vertical vector: num_domainsWboundary X 1)
    evaluate_northward(evaluate_northward < 0)	= 0;
    
    evaluate_westward	= grid_addresses_rho(:, east) - currentBox_W_rho; % all potential overlaps are 0 or '+'; (vertical vector: num_domainsWboundary X 1)
    evaluate_westward(evaluate_westward >= 0)	= 1; % (vertical vector: num_domainsWboundary X 1)
    evaluate_westward(evaluate_westward < 0)	= 0;

	evaluate_eastward	= currentBox_E_rho - grid_addresses_rho(:, west); % all potential overlaps are 0 or '+'; (vertical vector: num_domainsWboundary X 1)
    evaluate_eastward(evaluate_eastward >= 0)	= 1; % (vertical vector: num_domainsWboundary X 1)
    evaluate_eastward(evaluate_eastward < 0)	= 0;

  	evaluate_latitude	= evaluate_southward .* evaluate_northward; % BOTH conditions must be true for overlap; (vertical vector: num_domainsWboundary X 1)
	evaluate_longitude	= evaluate_westward  .* evaluate_eastward;  % BOTH conditions must be true for overlap; (vertical vector: num_domainsWboundary X 1)
 
    evaluate_latitude_matrix(:, domain_loop)	= evaluate_latitude;  % values = 0 or 1; overlap = 1; (2D matrix: num_domainsWboundary X num_domainsWboundary)
    evaluate_longitude_matrix(:, domain_loop)	= evaluate_longitude; % values = 0 or 1; overlap = 1; (2D matrix: num_domainsWboundary X num_domainsWboundary)

    
    % find addresses of southern overlapping rho, v, & u-points
    overlapAddress_S_rho(:, domain_loop)	= repmat(currentBox_S_rho, [num_domainsWboundary, 1]); % (2D matrix: num_domainsWboundary X num_domainsWboundary)
    overlapAddress_S_rho(:, domain_loop)	= max([overlapAddress_S_rho(:, domain_loop) grid_addresses_rho(:, south)], [], 2); % pick north-most of overlaping southern rho-points
    overlapAddress_S_rho(:, domain_loop)	= overlapAddress_S_rho(:, domain_loop) .* evaluate_latitude; % southern rho-point for ALL potentially overlapping ECOTRAN boxes
    
    overlapAddress_S_v(:, domain_loop)      = repmat(currentBox_S_v, [num_domainsWboundary, 1]); % (2D matrix: num_domainsWboundary X num_domainsWboundary)
    overlapAddress_S_v(:, domain_loop)      = max([overlapAddress_S_v(:, domain_loop) grid_addresses_v(:, south)], [], 2); % pick north-most of overlaping southern v-points
    overlapAddress_S_v(:, domain_loop)      = overlapAddress_S_v(:, domain_loop) .* evaluate_latitude; % southern v-point for ALL potentially overlapping ECOTRAN boxes
    
    overlapAddress_S_u(:, domain_loop)      = repmat(currentBox_S_u, [num_domainsWboundary, 1]); % (2D matrix: num_domainsWboundary X num_domainsWboundary)
    overlapAddress_S_u(:, domain_loop)      = max([overlapAddress_S_u(:, domain_loop) grid_addresses_u(:, south)], [], 2); % pick north-most of overlaping southern u-points
    overlapAddress_S_u(:, domain_loop)      = overlapAddress_S_u(:, domain_loop) .* evaluate_latitude; % southern u-point for ALL potentially overlapping ECOTRAN boxes
    
    
    % find addresses of northern overlapping rho, v, & u-points
    overlapAddress_N_rho(:, domain_loop)	= repmat(currentBox_N_rho, [num_domainsWboundary, 1]); % (2D matrix: num_domainsWboundary X num_domainsWboundary)
    overlapAddress_N_rho(:, domain_loop)	= min([overlapAddress_N_rho(:, domain_loop) grid_addresses_rho(:, north)], [], 2); % pick south-most of overlaping northern rho-points
    overlapAddress_N_rho(:, domain_loop)	= overlapAddress_N_rho(:, domain_loop) .* evaluate_latitude; % northern rho-point for ALL potentially overlapping ECOTRAN boxes
    
    overlapAddress_N_v(:, domain_loop)      = repmat(currentBox_N_v, [num_domainsWboundary, 1]); % (2D matrix: num_domainsWboundary X num_domainsWboundary)
    overlapAddress_N_v(:, domain_loop)      = min([overlapAddress_N_v(:, domain_loop) grid_addresses_v(:, north)], [], 2); % pick south-most of overlaping northern v-points
    overlapAddress_N_v(:, domain_loop)      = overlapAddress_N_v(:, domain_loop) .* evaluate_latitude; % northern v-point for ALL potentially overlapping ECOTRAN boxes
    
    overlapAddress_N_u(:, domain_loop)      = repmat(currentBox_N_u, [num_domainsWboundary, 1]); % (2D matrix: num_domainsWboundary X num_domainsWboundary)
    overlapAddress_N_u(:, domain_loop)      = min([overlapAddress_N_u(:, domain_loop) grid_addresses_u(:, north)], [], 2); % pick south-most of overlaping northern u-points
    overlapAddress_N_u(:, domain_loop)      = overlapAddress_N_u(:, domain_loop) .* evaluate_latitude; % northern u-point for ALL potentially overlapping ECOTRAN boxes
    
    
    % find addresses of western overlapping rho, v, & u-points
    overlapAddress_W_rho(:, domain_loop)    = repmat(currentBox_W_rho, [num_domainsWboundary, 1]); % (2D matrix: num_domainsWboundary X num_domainsWboundary)
    overlapAddress_W_rho(:, domain_loop)    = max([overlapAddress_W_rho(:, domain_loop) grid_addresses_rho(:, west)], [], 2); % pick east-most of overlaping western rho-points
    overlapAddress_W_rho(:, domain_loop)    = overlapAddress_W_rho(:, domain_loop) .* evaluate_longitude; % western rho-point for ALL potentially overlapping ECOTRAN boxes

    overlapAddress_W_v(:, domain_loop)      = repmat(currentBox_W_v, [num_domainsWboundary, 1]); % (2D matrix: num_domainsWboundary X num_domainsWboundary)
    overlapAddress_W_v(:, domain_loop)      = max([overlapAddress_W_v(:, domain_loop) grid_addresses_v(:, west)], [], 2); % pick east-most of overlaping western v-points
    overlapAddress_W_v(:, domain_loop)      = overlapAddress_W_v(:, domain_loop) .* evaluate_longitude; % western v-point for ALL potentially overlapping ECOTRAN boxes
    
    overlapAddress_W_u(:, domain_loop)      = repmat(currentBox_W_u, [num_domainsWboundary, 1]); % (2D matrix: num_domainsWboundary X num_domainsWboundary)
    overlapAddress_W_u(:, domain_loop)      = max([overlapAddress_W_u(:, domain_loop) grid_addresses_u(:, west)], [], 2); % pick east-most of overlaping western u-points
    overlapAddress_W_u(:, domain_loop)      = overlapAddress_W_u(:, domain_loop) .* evaluate_longitude; % western u-point for ALL potentially overlapping ECOTRAN boxes
    
    
    % find addresses of eastern overlapping rho, v, & u-points
    overlapAddress_E_rho(:, domain_loop)    = repmat(currentBox_E_rho, [num_domainsWboundary, 1]);
    overlapAddress_E_rho(:, domain_loop)    = min([overlapAddress_E_rho(:, domain_loop) grid_addresses_rho(:, east)], [], 2); % pick west-most of overlaping eastern rho-points
    overlapAddress_E_rho(:, domain_loop)    = overlapAddress_E_rho(:, domain_loop) .* evaluate_longitude; % eastern rho-point for ALL potentially overlapping ECOTRAN boxes
    
    overlapAddress_E_v(:, domain_loop)      = repmat(currentBox_E_v, [num_domainsWboundary, 1]);
    overlapAddress_E_v(:, domain_loop)      = min([overlapAddress_E_v(:, domain_loop) grid_addresses_v(:, east)], [], 2); % pick west-most of overlaping eastern v-points
    overlapAddress_E_v(:, domain_loop)      = overlapAddress_E_v(:, domain_loop) .* evaluate_longitude; % eastern v-point for ALL potentially overlapping ECOTRAN boxes
    
    overlapAddress_E_u(:, domain_loop)      = repmat(currentBox_E_u, [num_domainsWboundary, 1]);
    overlapAddress_E_u(:, domain_loop)      = min([overlapAddress_E_u(:, domain_loop) grid_addresses_u(:, east)], [], 2); % pick west-most of overlaping eastern u-points
    overlapAddress_E_u(:, domain_loop)      = overlapAddress_E_u(:, domain_loop) .* evaluate_longitude; % eastern u-point for ALL potentially overlapping ECOTRAN boxes
    % ---------------------------------------------------------------------
    
    
    % step 6d: find neighboring boxes -------------------------------------
    evaluate_SouthNeighbor = currentBox_S_rho - grid_addresses_rho(:, north); % will be "+1" if a potential neighbor; (vertical vector: num_domainsWboundary X 1)
    evaluate_SouthNeighbor(evaluate_SouthNeighbor ~= 1) = 0;
    
    evaluate_NorthNeighbor = grid_addresses_rho(:, south) - currentBox_N_rho; % will be "+1" if a potential neighbor; (vertical vector: num_domainsWboundary X 1)
    evaluate_NorthNeighbor(evaluate_NorthNeighbor ~= 1) = 0;
    
    evaluate_WestNeighbor = currentBox_W_rho - grid_addresses_rho(:, east); % will be "+1" if a potential neighbor; (vertical vector: num_domainsWboundary X 1)
    evaluate_WestNeighbor(evaluate_WestNeighbor ~= 1) = 0;
    
    evaluate_EastNeighbor = grid_addresses_rho(:, west) - currentBox_E_rho; % will be "+1" if a potential neighbor; (vertical vector: num_domainsWboundary X 1)
    evaluate_EastNeighbor(evaluate_EastNeighbor ~= 1) = 0;
    
    evaluate_SouthNeighbor_matrix(:, domain_loop)	= evaluate_SouthNeighbor; % values = 0 ro 1; potential neighbor = 1; (2D matrix: num_domainsWboundary X num_domainsWboundary)
    evaluate_NorthNeighbor_matrix(:, domain_loop)	= evaluate_NorthNeighbor; % values = 0 ro 1; potential neighbor = 1; (2D matrix: num_domainsWboundary X num_domainsWboundary)
    evaluate_WestNeighbor_matrix(:, domain_loop)	= evaluate_WestNeighbor;  % values = 0 ro 1; potential neighbor = 1; (2D matrix: num_domainsWboundary X num_domainsWboundary)
    evaluate_EastNeighbor_matrix(:, domain_loop)	= evaluate_EastNeighbor;  % values = 0 ro 1; potential neighbor = 1; (2D matrix: num_domainsWboundary X num_domainsWboundary)
    % ---------------------------------------------------------------------

end % domain_loop


% step 6e: evaluate latitudinal & longitudinal connectivity ---------------
evaluate_LatitudeNeighbor_matrix	= evaluate_SouthNeighbor_matrix + evaluate_NorthNeighbor_matrix; % EITHER neighbor to north or to south may be true
evaluate_LatitudeNeighbor_matrix(evaluate_LatitudeNeighbor_matrix ~= 0) = 1; % values = 0 ro 1; potential neighbor = 1; (2D matrix: num_domainsWboundary X num_domainsWboundary)

evaluate_LongitudeNeighbor_matrix	= evaluate_WestNeighbor_matrix + evaluate_EastNeighbor_matrix; % EITHER neighbor to north or to south may be true
evaluate_LongitudeNeighbor_matrix(evaluate_LongitudeNeighbor_matrix ~= 0) = 1; % values = 0 ro 1; potential neighbor = 1; (2D matrix: num_domainsWboundary X num_domainsWboundary)

% conditions for longitudinal connectivity:
%   boxes have N/S overlap
%   boxes have a neighboring E/W boundary
LongitudinalConnectivity            = evaluate_latitude_matrix .* evaluate_LongitudeNeighbor_matrix; % BOTH conditions must be true; values = 0 or 1; neighbor = 1; (2D matrix: num_domainsWboundary X num_domainsWboundary)

% conditions for latitudinal connectivity:
%   boxes have E/W overlap
%   boxes have a neighboring N/S boundary
LatitudinalConnectivity             = evaluate_longitude_matrix .* evaluate_LatitudeNeighbor_matrix; % BOTH conditions must be true; values = 0 or 1; neighbor = 1; (2D matrix: num_domainsWboundary X num_domainsWboundary)
% -------------------------------------------------------------------------


% step 6f: apply connectivity matrices to overlap boundaries --------------
% boundary addresses for N/S fluxes -----
LatitudeFluxConnectivity_E_rho = overlapAddress_E_rho .* LatitudinalConnectivity; % eastern rho-point for ALL overlapping ECOTRAN boxes; (2D matrix: num_domainsWboundary X num_domainsWboundary)
LatitudeFluxConnectivity_W_rho = overlapAddress_W_rho .* LatitudinalConnectivity; % western rho-point for ALL overlapping ECOTRAN boxes; (2D matrix: num_domainsWboundary X num_domainsWboundary)

LatitudeFluxConnectivity_E_v = overlapAddress_E_v .* LatitudinalConnectivity; % eastern v-point for ALL overlapping ECOTRAN boxes; (2D matrix: num_domainsWboundary X num_domainsWboundary)
LatitudeFluxConnectivity_W_v = overlapAddress_W_v .* LatitudinalConnectivity; % western v-point for ALL overlapping ECOTRAN boxes; (2D matrix: num_domainsWboundary X num_domainsWboundary)

LatitudeFluxConnectivity_E_u = overlapAddress_E_u .* LatitudinalConnectivity; % eastern u-point for ALL overlapping ECOTRAN boxes; (2D matrix: num_domainsWboundary X num_domainsWboundary)
LatitudeFluxConnectivity_W_u = overlapAddress_W_u .* LatitudinalConnectivity; % western u-point for ALL overlapping ECOTRAN boxes; (2D matrix: num_domainsWboundary X num_domainsWboundary)


% boundary addresses for E/W fluxes -----
LongitudeFluxConnectivity_N_rho = overlapAddress_N_rho .* LongitudinalConnectivity; % northern rho-point for ALL overlapping ECOTRAN boxes; (2D matrix: num_domainsWboundary X num_domainsWboundary)
LongitudeFluxConnectivity_S_rho = overlapAddress_S_rho .* LongitudinalConnectivity; % southern rho-point for ALL overlapping ECOTRAN boxes; (2D matrix: num_domainsWboundary X num_domainsWboundary)

LongitudeFluxConnectivity_N_v = overlapAddress_N_v .* LongitudinalConnectivity; % northern v-point for ALL overlapping ECOTRAN boxes; (2D matrix: num_domainsWboundary X num_domainsWboundary)
LongitudeFluxConnectivity_S_v = overlapAddress_S_v .* LongitudinalConnectivity; % southern v-point for ALL overlapping ECOTRAN boxes; (2D matrix: num_domainsWboundary X num_domainsWboundary)

LongitudeFluxConnectivity_N_u = overlapAddress_N_u .* LongitudinalConnectivity; % northern u-point for ALL overlapping ECOTRAN boxes; (2D matrix: num_domainsWboundary X num_domainsWboundary)
LongitudeFluxConnectivity_S_u = overlapAddress_S_u .* LongitudinalConnectivity; % southern u-point for ALL overlapping ECOTRAN boxes; (2D matrix: num_domainsWboundary X num_domainsWboundary)
% -------------------------------------------------------------------------


% step 6g: connectivity within the aggregated grid ------------------------
% QQQ WHY NO EXTERNAL BOUNDARY BOXES???
LatitudinalConnectivity_ECOTRAN     = LatitudinalConnectivity(1:num_domains, 1:num_domains);
LongitudinalConnectivity_ECOTRAN    = LongitudinalConnectivity(1:num_domains, 1:num_domains);
HorizontalConnectivity_ECOTRAN      = LatitudinalConnectivity_ECOTRAN + LongitudinalConnectivity_ECOTRAN;
HorizontalConnectivity_ECOTRAN(HorizontalConnectivity_ECOTRAN > 1) = 1; % for weird cases where cells can border both horizontally and vertically (won't happen with a regular rectangular grid)
% -------------------------------------------------------------------------


% step 6h: evaluate directional east, west, north, or south connectivity --
%          (for volume conservation forcing)
WesternConnectivity     = evaluate_WestNeighbor_matrix  .* evaluate_latitude_matrix;  % neighbors directly to the WEST of DOI; (2D matrix: num_domainsWboundary X num_domainsWboundary)
EasternConnectivity     = evaluate_EastNeighbor_matrix  .* evaluate_latitude_matrix;  % neighbors directly to the EAST of DOI; (2D matrix: num_domainsWboundary X num_domainsWboundary)
NorthernConnectivity	= evaluate_NorthNeighbor_matrix .* evaluate_longitude_matrix; % neighbors directly to the NORTH of DOI; (2D matrix: num_domainsWboundary X num_domainsWboundary)
SouthernConnectivity	= evaluate_SouthNeighbor_matrix .* evaluate_longitude_matrix; % neighbors directly to the SOUTH of DOI; (2D matrix: num_domainsWboundary X num_domainsWboundary)

% merge external boundary cells into 1 row & column
WesternConnectivity((num_domains+1), :)         = sum(WesternConnectivity((num_domains+1):end, :), 1);
WesternConnectivity((num_domains+2):end, :)     = []; % (2D matrix: (num_domains+1) X num_domainsWboundary)
WesternConnectivity(:, (num_domains+1))         = sum(WesternConnectivity(:, (num_domains+1):end), 2);
WesternConnectivity(:, (num_domains+2):end)     = []; % (2D matrix: (num_domains+1) X (num_domains+1))
WesternConnectivity(:, end)                     = 0; % (2D matrix: (num_domains+1) X (num_domains+1))

EasternConnectivity((num_domains+1), :)         = sum(EasternConnectivity((num_domains+1):end, :), 1);
EasternConnectivity((num_domains+2):end, :)     = []; % (2D matrix: (num_domains+1) X num_domainsWboundary)
EasternConnectivity(:, (num_domains+1))         = sum(EasternConnectivity(:, (num_domains+1):end), 2);
EasternConnectivity(:, (num_domains+2):end)     = []; % (2D matrix: (num_domains+1) X (num_domains+1))
EasternConnectivity(:, end)                     = 0; % (2D matrix: (num_domains+1) X (num_domains+1))

NorthernConnectivity((num_domains+1), :)        = sum(NorthernConnectivity((num_domains+1):end, :), 1);
NorthernConnectivity((num_domains+2):end, :)	= []; % (2D matrix: (num_domains+1) X num_domainsWboundary)
NorthernConnectivity(:, (num_domains+1))        = sum(NorthernConnectivity(:, (num_domains+1):end), 2);
NorthernConnectivity(:, (num_domains+2):end)	= []; % (2D matrix: (num_domains+1) X (num_domains+1))
NorthernConnectivity(:, end)                    = 0; % (2D matrix: (num_domains+1) X (num_domains+1))

SouthernConnectivity((num_domains+1), :)        = sum(SouthernConnectivity((num_domains+1):end, :), 1);
SouthernConnectivity((num_domains+2):end, :)	= []; % (2D matrix: (num_domains+1) X num_domainsWboundary)
SouthernConnectivity(:, (num_domains+1))        = sum(SouthernConnectivity(:, (num_domains+1):end), 2);
SouthernConnectivity(:, (num_domains+2):end)	= []; % (2D matrix: (num_domains+1) X (num_domains+1))
SouthernConnectivity(:, end)                    = 0; % (2D matrix: (num_domains+1) X (num_domains+1))
% -------------------------------------------------------------------------


% step 6i: define vertical connectivity -----------------------------------
%          for particle sinking & DVM
%          NOTE: this is an unstacked matrix (all depth layers are present in the 2D matrix)
%          NOTE: no -1 values on diagonal in ECOTRAN2
VerticalConnectivity	= zeros(num_domains*num_agg_z); % (2D matrix: (num_domains*num_agg_z) (147) X (num_domains*num_agg_z) (147))
temp_diag               = diag(ones(1, (num_domains*(num_agg_z-1))), 0);
VerticalConnectivity((num_domains+1):end, 1:(num_domains*(num_agg_z-1))) = temp_diag; % 0 or 1; (2D matrix: DESTINY-->(num_domains*num_agg_z) (147) X SOURCE-->(num_domains*num_agg_z) (147))
VerticalConnectivity((end+1), (end+1))	= 0; % add cells for inport/export cross-boundary fluxes (usually 0 for vertical mixing and sinking); (2D matrix: DESTINY-->(num_boxes+1) (148) X SOURCE-->(num_boxes+1) (148))

VerticalConnectivity_stack	= repmat(1:num_domains:num_boxes, [num_domains, 1]) + repmat((0:num_domains-1)', [1, 4]); % (2D matrix: num_domains X num_agg_z); tells which boxes lay directly under each surface box; clm 1 is suface, clm 2 is depth layer 2, final column is the bottom
looky_BottomBoxes           = (num_boxes - (num_domains + 1)):num_boxes; % addresses of bottom boxes
% -------------------------------------------------------------------------


% step 6j: clean up temporary variables -----------------------------------
clear currentBox_* overlapAddress_* evaluate_*
clear LatitudinalConnectivity_ECOTRAN LongitudinalConnectivity_ECOTRAN
clear temp_diag
% *************************************************************************





% *************************************************************************
% STEP 7: pack up results-------------------------------------------------
ROMSgrid.fname_ROMS_GridPrep            = fname_ROMS_GridPrep;
ROMSgrid.readFile_grid                  = readFile_grid;
ROMSgrid.readFile_DepthLevels           = readFile_DepthLevels;
ROMSgrid.readFile_ExampleYear           = readFile_ExampleYear;

ROMSgrid.z_definition                   = z_definition;
ROMSgrid.domain_definition              = domain_definition;

ROMSgrid.mask_v                         = mask_v;
ROMSgrid.mask_u                         = mask_u;
ROMSgrid.mask_rho                       = mask_rho;
ROMSgrid.mask_psi                       = mask_psi;

ROMSgrid.rows_v                         = rows_v;
ROMSgrid.clms_v                         = clms_v;
ROMSgrid.rows_u                         = rows_u;
ROMSgrid.clms_u                         = clms_u;
ROMSgrid.rows_rho                       = rows_rho;
ROMSgrid.clms_rho                       = clms_rho;
ROMSgrid.rows_h                         = rows_h;
ROMSgrid.clms_h                         = clms_h;
ROMSgrid.num_z                          = num_z;
ROMSgrid.num_w                          = num_w;
ROMSgrid.num_domains                    = num_domains;
ROMSgrid.num_domainsWboundary           = num_domainsWboundary;
ROMSgrid.num_agg_z                      = num_agg_z;
ROMSgrid.num_boxes                      = num_boxes;

ROMSgrid.z_w                            = z_w;  % z_w are the depths @ each horizontal depth-layer face at the w point (top & bottom face depths at the cell center)
ROMSgrid.z_wv                           = z_wv; % z_wv is the depth of each horizontal depth-layer face at the v point
ROMSgrid.z_wu                           = z_wu; % z_wu is the depth of each horizontal depth-layer face at the u point

ROMSgrid.v_height                       = v_height; % height of cell @ v-point; (bottom minus top);   (m); (3D matrix: v longitude (51) X v latitude (80) X num_z (42)); NOTE: JR confirms sum of these heights = bathymetry @ v
ROMSgrid.u_height                       = u_height; % height of cell @ u-point; (bottom minus top);   (m); (3D matrix: u longitude (50) X u latitude (81) X num_z (42))
ROMSgrid.w_height                       = w_height; % height of cell @ center; (bottom minus top);    (m); (3D matrix: rho longitude (51) X rho latitude (81) X num_z (42))
ROMSgrid.h                              = h; % bathymetry at RHO-points; (2D matrix: rho longitude (51 west:east) X rho latitude (81 south:north))

ROMSgrid.dist_NS_globe              	= dist_NS_globe; % N<-->S lengths relative to globe; (m); abs(NORTH minus SOUTH); (2D matrix: u longitude (50) X v latitude-1 (79));
ROMSgrid.dist_EW_globe              	= dist_EW_globe; % W<-->E lengths relative to globe; (m); abs(EAST minus WEST); (2D matrix: u longitude-1 (49) X v latitude (80));
ROMSgrid.dist_NS_face                   = dist_NS_face;  % N<-->S face lengths; (m); (2D matrix: u longitude (50) X v latitude-1 (79)); NOTE: x_psi(south) - x_psi(north)
ROMSgrid.dist_EW_face                   = dist_EW_face;  % E<-->W face lengths; (m); (2D matrix: u longitude-1 (49) X v latitude (80)); NOTE: y_psi(east) - y_psi(west)

ROMSgrid.area_NS_globe                  = area_NS_globe;    % N<-->S ROMS lateral cell area based on EW & NS distances relative to the globe; (m2); (3D matrix: v longitude-2 (49) X v latitude (80) X num_z (42))
ROMSgrid.area_EW_globe                  = area_EW_globe;    % E<-->W ROMS lateral cell area based on EW & NS distances relative to the globe; (m2); (3D matrix: v longitude-2 (49) X v latitude (80) X num_z (42))
ROMSgrid.area_floor_globe               = area_floor_globe; % cell floor area based on EW & NS distances relative to the globe; (m2); (3D matrix: rho longitude-2 (49) X rho latitude-2 (79) X num_w (43)); (includes only values at rho points the inside domain)
ROMSgrid.area_floor_face                = area_floor_face;	% cell floor area based on EW & NS distances relative to face lengths; (m2); (3D matrix: rho longitude-2 (49) X rho latitude-2 (79) X num_w (43)); (includes only values at rho points the inside domain)

ROMSgrid.LatitudinalConnectivity      	= LatitudinalConnectivity; % values = 0 or 1; neighbor = 1; (2D matrix: num_domainsWboundary X num_domainsWboundary)
ROMSgrid.LatitudeFluxConnectivity_W_v	= LatitudeFluxConnectivity_W_v; % western v-point for ALL overlapping ECOTRAN boxes; (2D matrix: num_domainsWboundary X num_domainsWboundary)
ROMSgrid.LatitudeFluxConnectivity_E_v	= LatitudeFluxConnectivity_E_v; % eastern v-point for ALL overlapping ECOTRAN boxes; (2D matrix: num_domainsWboundary X num_domainsWboundary)
ROMSgrid.grid_addresses_v              	= grid_addresses_v; % (2D matrix: num_domains X 7); [south north west east (num cells lat) (num cells lon) (num cells total)]

ROMSgrid.LongitudinalConnectivity     	= LongitudinalConnectivity; % values = 0 or 1; neighbor = 1; (2D matrix: num_domainsWboundary X num_domainsWboundary)
ROMSgrid.LongitudeFluxConnectivity_S_u	= LongitudeFluxConnectivity_S_u; % southern u-point for ALL overlapping ECOTRAN boxes; (2D matrix: num_domainsWboundary X num_domainsWboundary)
ROMSgrid.LongitudeFluxConnectivity_N_u	= LongitudeFluxConnectivity_N_u; % northern u-point for ALL overlapping ECOTRAN boxes; (2D matrix: num_domainsWboundary X num_domainsWboundary)
ROMSgrid.grid_addresses_u            	= grid_addresses_u; % (2D matrix: num_domains X 7); [south north west east (num cells lat) (num cells lon) (num cells total)]

ROMSgrid.grid_addresses_rho            	= grid_addresses_rho; % (2D matrix: num_domains X 7); [south north west east (num cells lat) (num cells lon) (num cells total)]

% these variables might be useful but are not currently used
ROMSgrid.bathymetry_range               = bathymetry_range; % bathymetry range (h) for each ECOTRAN domain; (2D matrix: num_domains X 3); [max min (num cells)]

ROMSgrid.ECOTRAN_grid_lat_rho           = ECOTRAN_grid_lat_rho; % (2D matrix: num_domains X 2); [South North]
ROMSgrid.ECOTRAN_grid_lon_rho           = ECOTRAN_grid_lon_rho; % (2D matrix: num_domains X 2); [West East]
ROMSgrid.ECOTRAN_grid_lat_v             = ECOTRAN_grid_lat_v;   % (2D matrix: num_domains X 2); [South North]
ROMSgrid.ECOTRAN_grid_lon_v             = ECOTRAN_grid_lon_v;   % (2D matrix: num_domains X 2); [West East]
ROMSgrid.ECOTRAN_grid_lat_u             = ECOTRAN_grid_lat_u;   % (2D matrix: num_domains X 2); [South North]
ROMSgrid.ECOTRAN_grid_lon_u             = ECOTRAN_grid_lon_u;   % (2D matrix: num_domains X 2); [West East]

ROMSgrid.HorizontalConnectivity_ECOTRAN	= HorizontalConnectivity_ECOTRAN; % (2D matrix: num_domains X num_domains); [DESTINY SOURCE]; NOTE: If wanted, go into code to save LatitudinalConnectivity_ECOTRAN & LongitudinalConnectivity_ECOTRAN.
ROMSgrid.VerticalConnectivity           = VerticalConnectivity; % 0 or 1; (2D matrix: DESTINY-->(num_boxes+1) (148) X SOURCE-->(num_boxes+1) (148))
ROMSgrid.VerticalConnectivity_stack     = VerticalConnectivity_stack; % (2D matrix: num_domains X num_agg_z); tells which boxes lay directly under each surface box; clm 1 is suface, clm 2 is depth layer 2, final column is the bottom

ROMSgrid.looky_BottomBoxes              = looky_BottomBoxes; % addresses of bottom (benthic) boxes
% *************************************************************************


% end m-file***************************************************************