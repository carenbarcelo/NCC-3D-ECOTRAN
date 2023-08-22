function ECOTRANmigration = f_ECOTRANmigration_NCC_02222022(ECOTRANphysics)
% returns area of overlap (m2) between sub-domain boxes to which are later
% applied swimming speeds between boxes to get migration volumetric flux
% rates (m3/d).
% Overlaps are defined manually for each model.
% This code works for NCC ROMS and is applied to DVM only
% by Jim Ruzicka
%
% calls:
%       none
%
% takes:
%       ECOTRANphysics
%           num_t                   number of time steps	(vertical vector: num_t X 1)
%           num_boxes               1-MLD, 2-epipelagic, 3-mesopelagic, 4-bathypelagic, 5-benthic (QQQ delete MLD)
%           BoxArea_y               side or face area; North-South (y) orientation; (m2)	(2D matrix: time X num_boxes)
%           BoxArea_x               side or face area; East-West (x) orientation; (m2)      (2D matrix: time X num_boxes)
%           BoxFloorArea            floor area; (m2)        (2D matrix: time X num_boxes SOURCE)
%
% returns:
%       ECOTRANmigration
%           num_boxes
%           num_fluxes_migration	migration fluxes (potential number per group)
%           MIGRATION               migration as boundary area between boxes; (m2); (3D matrix: num_t X source (num_boxes+1) X destiny (num_boxes+1))
%           fname_ECOTRANmigration	name of this f_ECOTRANmigration function
%
% revision date: 2-22-2022


% *************************************************************************
% STEP 1: basic domain parameters (NCC)-------------------------------
% step 1a: function code name ---------------------------------------------
fname_ECOTRANmigration	= mfilename; % name of this f_ECOTRANmigration function
display(['Running: ' fname_ECOTRANmigration])

% step 1b: unpack input terms ---------------------------------------------
t_grid                  = ECOTRANphysics.t_grid; % time points (days)
ROMS_time               = ECOTRANphysics.ROMS_time; % seconds since 1900-01-01 00:00:00; (vertical vector: 1 X days)
num_t                   = ECOTRANphysics.num_t;         % number of time steps; (vertical vector: num_t X 1)
num_boxes               = ECOTRANphysics.num_boxes;     % (scalar: num_boxes); 1-MLD, 2-epipelagic, 3-mesopelagic, 4-bathypelagic, 5-benthic (QQQ delete MLD)
% BoxArea_y               = ECOTRANphysics.BoxArea_y;	  % side or face area; North-South (y) orientation; (m2); (2D matrix: time X num_boxes)
% BoxArea_x               = ECOTRANphysics.BoxArea_x;	  % side or face area; East-West (x) orientation; (m2); (2D matrix: time X num_boxes)
BoxFloorArea            = ECOTRANphysics.BoxFloorArea;	% floor area; (m2); (horizontal vector: 1 X SOURCE (num_boxes))

VerticalConnectivity	= ECOTRANphysics.VerticalConnectivity; % connectivity between source & destiny boxes; 0 or 1; (3D matrix: 1 X SOURCE-->(num_boxes+1) (148) X DESTINY-->(num_boxes+1) (148))
SINKING                 = ECOTRANphysics.SINKING;       % sinking source box floor area (m2);   (3D matrix: num_t_ROMS X SOURCE-->(num_boxes+1) X DESTINY-->(num_boxes+1))
% *************************************************************************



% *************************************************************************
% STEP 2: Physical Structure (GoMexOCN)------------------------------------
%         define spatial relationships between boxes

% step 2a: vertical migration ---------------------------------------------
looky_connectivity      = find(VerticalConnectivity == 1);
num_fluxes_migration	= length(looky_connectivity); % number of potential migration fluxes per group); number of 1s in the "links" matrices; so far, this is just vertical migration fluxes

% vertical_links       	= repmat(VerticalConnectivity, [num_t, 1, 1]); % (3D matrix: num_t X SOURCE (num_boxes + boundary)  X DESTINY (num_boxes + boundary))
% BoxFloorArea(end+1)     = 0; % add dimensions (0m2) of boundary import box; (m2); (m2); (horizontal vector: 1 X SOURCE-->(num_boxes+1) (148))
% repmat_BoxFloorArea     = repmat(BoxFloorArea, [1, 1, (num_boxes+1)]);	% (m2); (3D matrix: 1 X SOURCE-->(num_boxes+1) (148) X DESTINY-->(num_boxes+1) (148))
% repmat_BoxFloorArea     = repmat(repmat_BoxFloorArea, [num_t, 1]);	% floor area; (m2); (3D matrix: num_t (366) X SOURCE-->(num_boxes+1) (148) X DESTINY-->(num_boxes+1) (148))
% vertical_links          = vertical_links .* repmat_BoxFloorArea; % QQQ still needs boundary box floor area, or zeros?

% use SINKING connectivity directly for DVM
vertical_links          = SINKING; % sinking source box floor area (m2);   (3D matrix: num_t_ROMS X SOURCE-->(num_boxes+1) X DESTINY-->(num_boxes+1))

% FFF I'll need a vector to tell which depth level a box is on in order to
% later define migration between specific depth layers in cases where I
% have many geographic domains
% -------------------------------------------------------------------------

% step 2b: east-west migration --------------------------------------------

% step 2c: north-south migration ------------------------------------------

% step 2d: combine all migration links ------------------------------------
migration_links         = vertical_links; % add other links here (east-west face overlap & north-south face overlap); (3D matrix: num_t X SOURCE (num_boxes + boundary)  X DESTINY (num_boxes + boundary))
% *************************************************************************



% *************************************************************************
% STEP 3: define flux time-series------------------------------------------
% step 3a: define flux time-series ----------------------------------------
MIGRATION               = migration_links;	% migration as boundary area between boxes; (m2); (3D matrix: num_t X SOURCE (num_boxes + boundary)  X DESTINY (num_boxes + boundary))

% step 3b: compact fluxes -------------------------------------------------
CompactFlux_MIGRATION	= f_CompactFluxTimeSeries_11182019(MIGRATION);	% compact MIGRATION 3D matrix


% step 3c: interpolate ROMS_time to t_grid --------------------------------
temp_ROMS_time                           	= t_grid(1):1:(ROMS_time(end) - ROMS_time(1))/86400 + 1;

old_MIGRATION_flux                          = CompactFlux_MIGRATION.compact_flux;
CompactFlux_MIGRATION.compact_flux      	= interp1(temp_ROMS_time, old_MIGRATION_flux, t_grid, 'pchip');
% *************************************************************************



% *************************************************************************
% STEP 4: pack results for export-----------------------------------------
ECOTRANmigration.fname_ECOTRANmigration     = fname_ECOTRANmigration;	% name of this f_ECOTRANmigration function
ECOTRANphysics.fname_CompactFlux            = CompactFlux_MIGRATION.fname_CompactFlux;           % file name of f_CompactFluxTimeSeries function

ECOTRANmigration.num_boxes              	= num_boxes;
ECOTRANmigration.num_fluxes_migration       = num_fluxes_migration;     % migration fluxes (potential number per group)

ECOTRANmigration.MIGRATION                  = MIGRATION;                % migration as boundary area between boxes; (m2); (3D matrix: num_t X SOURCE (num_boxes+boundary) X DESTINY (num_boxes+boundary))
ECOTRANmigration.CompactFlux_MIGRATION      = CompactFlux_MIGRATION; % (structure)
% *************************************************************************


% end m-file***************************************************************