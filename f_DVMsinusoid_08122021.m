function DVM = f_DVMsinusoid_08122021(ECOTRANphysics, ODEinput, DVMinput, t_grid)
% Calculate DVM flux rates between all model domain boxes for each functional group at each time point
% Uses a daily sinusoidal migration pattern
% by Jim Ruzicka
%
% calls:
%   none
%
% takes:
%	ECOTRANphysics:
%       BoxVolume                       (m3); (2D matrix: num_t X num_boxes)
%       num_boxes
%   ODEinput:
%       num_grps                        number of functional groups
%       num_fluxes_migration            migration fluxes (potential number of fluxes for each group); includes external fluxes (external flux count defaults to 2 if all external fluxes equal 0)
%       looky_MigrationFlux             definitions of linkages between boxes; (2D matrix: num_fluxes_migration X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)])
%                                           NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
%                                           FFF: will eventually need distinct variables for DVM and for other migration types
%	DVMinput:
%       MigrationArea_compact           migration as boundary area between boxes; (m2); (3D matrix: num_t X num_fluxes_migration)
%       MigrationSpeed                  maximum migration speed; (m/d); (3D matrix: num_grps, source (num_boxes+1) X destiny (num_boxes+1))
%                                           NOTE: dusk speeds are negative (in "standard" DVM)
%   t_grid
%   dt                  time-step;  (d); (scalar)
%
% returns:
%   DVM:
%       MIGRATION_compact           DVM migration rate of each group between each box; (m3/d); (3D matrix: num_t X num_grps X num_fluxes_migration)
%       dawn_phase                  migration speed scaler; dawn phase scalers are all positive; (vertical vector: num_t X 1)
%       dusk_phase                  migration speed scaler; dusk phase scalers are all negative; (vertical vector: num_t X 1)
%       num_fluxes_migration        migration fluxes (potential number of fluxes for each group); includes external fluxes (external flux count defaults to 2 if all external fluxes equal 0)
%       looky_MigrationFlux     	(2D matrix: num_fluxes_migration X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
%
% revision date: 8/12/2021


% *************************************************************************
% STEP 1: prepare operating terms------------------------------------------
fname_DVMmodel          = mfilename; % save name of this m-file to keep in saved model results
display(['Running: ' fname_DVMmodel])
DVM.fname_DVMmodel      = fname_DVMmodel;


% step 1a: unpack & prepare inputs ----------------------------------------
julian_day                      = floor(t_grid);
num_t                           = length(t_grid);
num_days                        = max(julian_day);

MigrationArea_compact           = DVMinput.MigrationArea_compact;	% migration as boundary area between boxes; (m2); (3D matrix: num_t X num_fluxes_migration)
MigrationSpeed                  = DVMinput.MigrationSpeed;   % (m/d); (3D matrix: num_grps, source (num_boxes+1) X destiny (num_boxes+1))

BoxVolume                       = ECOTRANphysics.BoxVolume; % (m3); (2D matrix: num_t X num_boxes)
num_boxes                       = ECOTRANphysics.num_boxes;

num_grps                        = ODEinput.num_grps;
num_fluxes_migration            = ODEinput.num_fluxes_migration; % migration fluxes (potential number of fluxes for each group); includes external fluxes (external flux count defaults to 2 if all external fluxes equal 0)
looky_MigrationFlux             = ODEinput.looky_MigrationFlux; % (2D matrix: num_fluxes_migration X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t

MIGRATION_compact               = zeros(num_t, num_grps, num_fluxes_migration); % initialize; (m3/d); (3D matrix: num_t X num_grps X num_fluxes_migration)
% -------------------------------------------------------------------------


% step 1b: prepare terms for individual functional groups -----------------
% separate sunrise & sunset migration speeds -----
MigrationSpeed_dawn     = MigrationSpeed;   % maximum migration speed at sunrise; (m/d); (3D matrix: num_grps, source (num_boxes+1) X destiny (num_boxes+1))
MigrationSpeed_dusk     = MigrationSpeed;   % maximum migration speed at sunset; (m/d); (3D matrix: num_grps, source (num_boxes+1) X destiny (num_boxes+1))
MigrationSpeed_dawn(MigrationSpeed_dawn < 0) = 0;       % (m/d); (3D matrix: num_grps, SOURCE (num_boxes+1) X DESTINY (num_boxes+1))
MigrationSpeed_dusk(MigrationSpeed_dusk > 0) = 0;       % (m/d); (3D matrix: num_grps, SOURCE (num_boxes+1) X DESTINY (num_boxes+1))
% -------------------------------------------------------------------------


% step 1c: get number of migrators and their addresses --------------------
temp_MigrationSpeed     = sum(MigrationSpeed, 3);
[temp_rows, ~]          = find(MigrationSpeed ~= 0);
looky_DVMgrps           = unique(temp_rows);     % row addresses of groups with DVM
num_DVMgrps             = length(looky_DVMgrps); % number of DVM groups
% *************************************************************************





% *************************************************************************
% STEP 2: prepare sinusoidal function--------------------------------------
amplitude                   = 1;
period_length               = 1; % period length in days
period                      = (1/period_length)*2*pi(); % NOTE: 1 day period = 2*pi
phase                       = period * (0/4);
diel_cycle                  = amplitude * sin((period * t_grid) + phase);

dawn_phase                  = diel_cycle; % migration speed scaler for DAWN; (vertical vector: num_t X 1)
dusk_phase                  = diel_cycle; % migration speed scaler for DUSK; (vertical vector: num_t X 1)
dawn_phase(dawn_phase < 0)	= 0; % dawn phase scalers are all positive
dusk_phase(dusk_phase > 0)	= 0; % dusk phase scalers are all negative
% *************************************************************************





% *************************************************************************
% STEP 3: calculate migration speed for each group at each time point------
% step 3a: repeat steps for each migrating group --------------------------
if num_DVMgrps > 0

    for DVMgrp_loop = 1:num_DVMgrps

        current_DVMgrp_RC           = looky_DVMgrps(DVMgrp_loop); % address (row number) of current DVM migrator group
        
        % step 3b: calculate migration speed for each group at each time point------
        %         generate MIGRATION_compact matrix defining DVM migration rate (m3/d) of each group between each box
        %         NOTE: dusk phase speeds are negative (in "standard" DVM)
        for flux_loop = 1:num_fluxes_migration

            current_SOURCE          = looky_MigrationFlux(flux_loop, 2); % current SOURCE box address
            current_DESTINY         = looky_MigrationFlux(flux_loop, 1); % current DESTINY box address
            
            if (current_SOURCE <= num_boxes) && (current_DESTINY <= num_boxes) % if steps over DVM from outside of domain (leaving these rates as 0)

                current_MigrationArea   = MigrationArea_compact(:, flux_loop); % migration as boundary area between boxes; (m2); (vertical vector: num_t X 1)
                
                current_SOURCE_volume	= BoxVolume(:, current_SOURCE);  % (m3); (vertical vector: num_t X 1)
                current_DESTINY_volume	= BoxVolume(:, current_DESTINY); % (m3); (vertical vector: num_t X 1)
                
                SOURCE_DESTINY_VolumeRatio = (current_SOURCE_volume ./ current_DESTINY_volume); % (proportion); (vertical vector: num_t X 1)
                SOURCE_DESTINY_VolumeRatio(isnan(SOURCE_DESTINY_VolumeRatio)) = 0; % correct for div/0 error; assume ratio is 0 (div/0 means DESINY volum = 0 and we might as well scale DVM to 0)
                
                current_MigrationSpeed_dawn = MigrationSpeed_dawn(current_DVMgrp_RC, current_SOURCE, current_DESTINY) * dawn_phase; % (m/d); (vertical vector: num_t X 1)
                current_MigrationSpeed_dusk = MigrationSpeed_dusk(current_DVMgrp_RC, current_SOURCE, current_DESTINY) * dusk_phase; % (m/d); (vertical vector: num_t X 1)
                    
                current_DVMgrp_dawn	= current_MigrationSpeed_dawn .* current_MigrationArea; % (m3/d); (vertical vector: num_t X 1)
                current_DVMgrp_dusk	= current_MigrationSpeed_dusk .* current_MigrationArea; % (m3/d); (vertical vector: num_t X 1)
                    
                MIGRATION_compact(:, current_DVMgrp_RC, flux_loop)     = current_DVMgrp_dawn + current_DVMgrp_dusk; % (m3/d); (3D matrix: num_t X num_grps X num_fluxes_migration)
                
            end % end if (current_SOURCE <= num_boxes) && (current_DESTINY <= num_boxes)
        end % flux_loop
        % -----------------------------------------------------------------

    end % DVMgrp_loop
end % if num_DVMgrps > 0
% *************************************************************************





% *************************************************************************
% STEP 4: pack up results--------------------------------------------------
DVM.MIGRATION_compact       = MIGRATION_compact; % DVM migration rate of each group between each box; (m3/d); (3D matrix: num_t X num_grps X num_fluxes_migration)
DVM.dawn_phase              = dawn_phase;        % dawn phase migration speed scaler; dawn phase scalers are all positive; (vertical vector: num_t X 1)
DVM.dusk_phase              = dusk_phase;        % dusk phase migration speed scaler; dusk phase scalers are all negative; (vertical vector: num_t X 1)
DVM.num_fluxes_migration	= ODEinput.num_fluxes_migration;	% migration fluxes (potential number of fluxes for each group); includes external fluxes (external flux count defaults to 2 if all external fluxes equal 0)
DVM.looky_MigrationFlux     = ODEinput.looky_MigrationFlux;     % (2D matrix: num_fluxes_migration X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
% *************************************************************************


% end m-file***************************************************************