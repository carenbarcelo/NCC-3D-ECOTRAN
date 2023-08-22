function [AddressStruct, TrophicStruct, FuncRespStruct, PhysicsStruct, fname_PrepMexODE] = f_PrepMexODE_09192022(ODEinput)
% by Jim Ruzicka
% prepare ECOTRAN variables for using C++ ODE solver.
% pack parameters & drivers for C++ mex function
%
% calls:
%       f_unspoolMATRIX_04282020        linearize multidimenional matrices up to 4-D for use in C++
%
% takes:
%       ODEinput
%
% returns:
%       AddressStruct
%       TrophicStruct
%       FuncRespStruct
%       PhysicsStruct
%
% revision date: 9-19-2022
%       8/21/2022 - added pb parameter to PhysicsStruct
%       8/21/2022 - added Thornton Lessem temperature response variable
%       9/19/2022 - added external forcing rate variables


% *************************************************************************
% STEP 1: pack parameters & drivers for C++ mex function-------------------
fname_PrepMexODE	= mfilename; % save name of this m-file to keep in saved model results
display(['Running: ' fname_PrepMexODE])


% step 1a: group addresses & counts ---------------------------------------
%          NOTE: data type uint8 limits number of groups to 255 (change to uint16 to allow for 65,535 groups)
%          NOTE: uint32 data type for num_t (allows for 4,294,967,295 time steps)
AddressStruct.num_grps                       	= uint8(ODEinput.num_grps);
AddressStruct.num_nutrients                    	= uint8(ODEinput.num_nutrients);
AddressStruct.num_NH4                           = uint8(ODEinput.num_NH4);
AddressStruct.num_plgcNH4                    	= uint8(ODEinput.num_plgcNH4);
AddressStruct.num_bnthNH4                       = uint8(ODEinput.num_bnthNH4);
AddressStruct.num_ANYPrimaryProd               	= uint8(ODEinput.num_ANYPrimaryProd);
AddressStruct.num_phytoplankton               	= uint8(ODEinput.num_phytoplankton);
AddressStruct.num_macroalgae                 	= uint8(ODEinput.num_macroalgae);
AddressStruct.num_eggs                         	= uint8(ODEinput.num_eggs);
AddressStruct.num_ANYdetritus                 	= uint8(ODEinput.num_ANYdetritus);
AddressStruct.num_livingANDfleets               = uint8(ODEinput.num_livingANDfleets);
AddressStruct.num_drivers                     	= uint8(ODEinput.num_drivers);
AddressStruct.num_boxes                        	= uint8(ODEinput.num_boxes);        % number of spatial boxes in model
AddressStruct.num_t                            	= uint32(ODEinput.num_t);           % length of t_grid (note uint32 data type)

AddressStruct.looky_driver                   	= uint8(ODEinput.looky_driver);     % row address(es) of driver group(s) (e.g., NO3); vector
AddressStruct.looky_nutrients                  	= uint8(ODEinput.looky_nutrients);
AddressStruct.looky_NO3                     	= uint8(ODEinput.looky_NO3);        % NOTE: use with Michaelis-Menten option
AddressStruct.looky_NH4                     	= uint8(ODEinput.looky_NH4);
AddressStruct.looky_plgcNH4                   	= uint8(ODEinput.looky_plgcNH4);
AddressStruct.looky_bnthNH4                    	= uint8(ODEinput.looky_bnthNH4);
% AddressStruct.looky_ANYPrimaryProducer       	  = uint8(ODEinput.looky_ANYPrimaryProducer); % NOTE: use with Michaelis-Menten option
AddressStruct.looky_eggs                      	= uint8(ODEinput.looky_eggs);
% AddressStruct.looky_fleets                 	  = uint8(ODEinput.looky_fleets;
AddressStruct.looky_livingANDfleets            	= uint8(ODEinput.looky_livingANDfleets);
AddressStruct.looky_ANYdetritus                 = uint8(ODEinput.looky_ANYdetritus);

AddressStruct.looky_externalForcing           	= uint8(ODEinput.looky_externalForcing); % NEW!!
AddressStruct.num_externalForcing_grps       	= uint8(ODEinput.num_externalForcing_grps); % NEW!!
% -------------------------------------------------------------------------



% step 1b: EnergyBudget ---------------------------------------------------
[TrophicStruct.EnergyBudget, ~]                         = f_unspoolMATRIX_04282020(ODEinput.EnergyBudget); % (3D matrix: num_grps X num_grps X num_boxes)
% -------------------------------------------------------------------------

% step 1c: ConsumptionBudget terms ----------------------------------------
[TrophicStruct.ConsumptionBudget_feces, ~]              = f_unspoolMATRIX_04282020(ODEinput.ConsumptionBudget_feces);        % (3D matrix: num_t X num_grps X num_boxes)
[TrophicStruct.ConsumptionBudget_metabolism, ~]         = f_unspoolMATRIX_04282020(ODEinput.ConsumptionBudget_metabolism);	% (3D matrix: num_t X num_grps X num_boxes)
[TrophicStruct.ConsumptionBudget_eggs, ~]           	= f_unspoolMATRIX_04282020(ODEinput.ConsumptionBudget_eggs);         % (3D matrix: num_t X num_grps X num_boxes)
[TrophicStruct.ConsumptionBudget_predation, ~]      	= f_unspoolMATRIX_04282020(ODEinput.ConsumptionBudget_predation);    % (3D matrix: num_t X num_grps X num_boxes)
[TrophicStruct.ConsumptionBudget_senescence, ~]     	= f_unspoolMATRIX_04282020(ODEinput.ConsumptionBudget_senescence);	% (3D matrix: num_t X num_grps X num_boxes)
[TrophicStruct.ConsumptionBudget_ba, ~]             	= f_unspoolMATRIX_04282020(ODEinput.ConsumptionBudget_ba);          	% (3D matrix: num_t X num_grps X num_boxes)
[TrophicStruct.ConsumptionBudget_em, ~]                 = f_unspoolMATRIX_04282020(ODEinput.ConsumptionBudget_em);          	% (3D matrix: num_t X num_grps X num_boxes)
% -------------------------------------------------------------------------

% step 1d: TransferEfficiency ---------------------------------------------
[TrophicStruct.TransferEfficiency, ~]                   = f_unspoolMATRIX_04282020(ODEinput.TransferEfficiency);              % (3D matrix: 1 X num_grps X num_boxes)
% -------------------------------------------------------------------------

% step 1e: functional group fate parameters -------------------------------
[TrophicStruct.fate_feces, ~]                         = f_unspoolMATRIX_04282020(ODEinput.fate_feces);   	% (3D matrix: num_ANYdetritus X num_grps X num_boxes)
[TrophicStruct.fate_metabolism, ~]                    = f_unspoolMATRIX_04282020(ODEinput.fate_metabolism);	% (3D matrix: num_nutrients X num_grps X num_boxes)
[TrophicStruct.fate_eggs, ~]                          = f_unspoolMATRIX_04282020(ODEinput.fate_eggs);    	% (3D matrix: num_eggs X num_grps X num_boxes)
[TrophicStruct.fate_predation, ~]                     = f_unspoolMATRIX_04282020(ODEinput.fate_predation);  % (3D matrix: num_livingANDfleets X num_grps X num_boxes)
[TrophicStruct.fate_senescence, ~]                    = f_unspoolMATRIX_04282020(ODEinput.fate_senescence);	% (3D matrix: num_ANYdetritus X num_grps X num_boxes)
% -------------------------------------------------------------------------


% step 1f: initial & boundary conditions, external driver -----------------
[TrophicStruct.production_initial, ~]               = f_unspoolMATRIX_04282020(ODEinput.production_initial);	% production rates to use as initial conditions; (mmole N/m3/d); (3D matrix: num_grps X 1 X num_boxes)

% QQQ don't need to replicate in C++ code loops
[TrophicStruct.productionC_initial_repmat, ~]      	= f_unspoolMATRIX_04282020(ODEinput.productionC_initial_repmat);	% consumption inflow rate; (mmole N/m3/d); (3D matrix: num_grps X num_grps X num_boxes); NOTE: replicated vertical vectors across columns
% [TrophicStruct.biomass_boundary, ~]                   = f_unspoolMATRIX_04282020(ODEinput.biomass_boundary);            % NOTE: use this line only for NON-REFLECTIVE boundary option; (mmoles N/m3); (3D matrix: num_t X num_grps X num_boxes+1)
[TrophicStruct.external_driver, ~]                 	= f_unspoolMATRIX_04282020(ODEinput.external_driver);             % external boundary biomass conditions for each box; external NO3 driver; (mmole N/m3); (3D matrix: num_t X num_drivers X num_boxes+1)
%   FFF external_driver defined for all boxes now, but will try to compact to just the boxes with fluxes (3D matrix: num_t X num_drivers X num_fluxes_BoundaryImport)
%   FFF could allow for allow multiple external_driver grps?

[TrophicStruct.externalForcing, ~]                 	= f_unspoolMATRIX_04282020(ODEinput.externalForcing);             % NEW!!! external input rate(s) of select group(s) to each box; (mmole N/m3/d); (3D matrix: num_t X num_externalForcing_grps X num_boxes)

[TrophicStruct.biomass_plus1, ~]                 	= f_unspoolMATRIX_04282020(ODEinput.biomass_plus1);               % initialized as all zeros; add 1 clm for boundary fluxes; (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1)
% -------------------------------------------------------------------------


% % step 1g: scenario rules -------------------------------------------------
% %          NOTE: this could be anything, using a vertical vector for scaling NO3 levels for now
% % NNN come back to this - find dynmic scenario code in other file
% [TrophicStruct.Scenario_scaler, dims_Scenario_scaler]         = f_unspoolMATRIX_04282020(ODEinput.Scenario_scaler);     % scenario time-series (scaler); (3D matrix: time X num_grps X num_boxes)
% -------------------------------------------------------------------------

% step 1h: functional response parameters -------------------------------
% QQQ no need to replicate for C++ loops?
[FuncRespStruct.FunctionalResponseParams, ~]     	= f_unspoolMATRIX_04282020(ODEinput.FunctionalResponseParams); % producer vulnerability (m_p); (3D matrix: CONSUMERS X prey group X num_boxes) replicated across clms (= producers)
% -------------------------------------------------------------------------

% step 1i: Thornton Lessem consumption rate scaler for temperature response
% FFF keeping name of variable q_TemperatureScaler for now, but it is a scaler driven by anything that affects q rate (temperature, mining plumes, ...)
[FuncRespStruct.q_TemperatureScaler, ~]             = f_unspoolMATRIX_04282020(ODEinput.q_TemperatureScaler);	% Thornton-Lessem temperature adjustment to consumption rate; value between 0 and 1; (3D matrix: num_t X num_grps X num_boxes)
% -------------------------------------------------------------------------

% % step 1j: terms used for Michaelis-Menten option -------------------------
% [FuncRespStruct.MichaelisMenten_Vmax, ~]        	= f_unspoolMATRIX_04282020(ODEinput.MichaelisMenten_Vmax);      % Vmax  = maximum nutrient uptake rate; (1/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes); NOTE: use for Michaelis-Menten function
% [FuncRespStruct.MichaelisMenten_KNO3, ~]       	= f_unspoolMATRIX_04282020(ODEinput.MichaelisMenten_KNO3);      % KNO3  = NO3 half-saturation constant; (mmol N/m3); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes); NOTE: use for Michaelis-Menten function
% [FuncRespStruct.MichaelisMenten_KNH4, ~]         	= f_unspoolMATRIX_04282020(ODEinput.MichaelisMenten_KNH4);      % KNH4  = NH4 half-saturation constant; (mmol N/m3); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes); NOTE: use for Michaelis-Menten function
% [FuncRespStruct.MichaelisMenten_alpha, ~]          = f_unspoolMATRIX_04282020(ODEinput.MichaelisMenten_alpha);     % alpha = initial slope of light response curve; (m2/W/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes); NOTE: use for Michaelis-Menten function
% [FuncRespStruct.MichaelisMenten_psi, ~]            = f_unspoolMATRIX_04282020(ODEinput.MichaelisMenten_psi);       % psi   = NO3 uptake inhibition by NH4; (m3/mmole N); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes); NOTE: use for Michaelis-Menten function
% [FuncRespStruct.MichaelisMenten_w, ~]              = f_unspoolMATRIX_04282020(ODEinput.MichaelisMenten_w);         % w	= sinking rate; (m/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes); NOTE: use for Michaelis-Menten function
% [FuncRespStruct.MichaelisMenten_eta, ~]            = f_unspoolMATRIX_04282020(ODEinput.MichaelisMenten_eta);       % eta = non-grazing mortality; (1/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes); NOTE: use for Michaelis-Menten function
% 
% FuncRespStruct.MLD                                 = ODEinput.MLD;                       % time-series of Mixed Layer Depth (for light intensity calculations, NOT advection calcs); (m); (vertical vector: length = num_t); NOTE: use for Michaelis-Menten option
% FuncRespStruct.Io                                  = ODEinput.Io;                        % time-series of surface PAR light intensity; daily integrated solar raditation at ocean surface averaged across 24 hours; (W m^-2 h^-1); NOTE: the daily average is also what was used by Spitz in her NPZD models; NOTE: use for Michaelis-Menten function
% FuncRespStruct.Kw                                  = ODEinput.Kw;                        % Light attenuation of seawater; (scaler); NOTE: use for Michaelis-Menten function
% FuncRespStruct.Kp                                  = ODEinput.Kp;                        % Light attenuation of phytoplankton (Newberger et al., 2003); (m2/mmol N); NOTE: use for Michaelis-Menten function
% 
% [FuncRespStruct.ProductionFraction_macroalgae, ~]	= ODEinput.ProductionFraction_macroalgae; % tie any macroalgae production to a fixed multiple of total phytoplankton production in each box; (3D matrix: num_macroalgae X 1 X num_boxes)
% -------------------------------------------------------------------------



% step 1k: physics terms --------------------------------------------------
% QQQ delete t_grid altogether because it may not be used in C++ interpolation?
PhysicsStruct.t_grid                            = ODEinput.t_grid;                                                  % (vertical vector: num_t x 1)
PhysicsStruct.dt                                = ODEinput.dt;                                                      % time-step
[PhysicsStruct.BoxVolume, ~]                    = f_unspoolMATRIX_04282020(ODEinput.BoxVolume);                  % (m3); (2D matrix: num_t X num_boxes)
[PhysicsStruct.qb, ~]                        	= f_unspoolMATRIX_04282020(ODEinput.qb);                         % intrinsic (weight-specific) consumption rates; (1/d); (3D matrix: num_t X num_grps X num_boxes)
[PhysicsStruct.pb, ~]                        	= f_unspoolMATRIX_04282020(ODEinput.pb);                         % intrinsic (weight-specific) production rates;  (1/d); (3D matrix: num_t X num_grps X num_boxes)
% [PhysicsStruct.biomass_migrator, ~]       	  = f_unspoolMATRIX_04282020(ODEinput.biomass_migrator);  % special definition of boundary biomasses for migrators; (mmoles N/m3); (3D matrix: num_t X num_grps X num_boxes); NOTE: deactived by default
[PhysicsStruct.RetentionScaler, ~]              = f_unspoolMATRIX_04282020(ODEinput.RetentionScaler); % (scaler 0-1); (2D matrix: num_grps X num_boxes DESTINY)

[PhysicsStruct.ADVECTION_compact, ~]                        = f_unspoolMATRIX_04282020(ODEinput.ADVECTION_compact);          % (m3/d); (2D matrix: num_t X num_fluxes_advection)
[num_fluxes_AdvectionBoundary_import, ~]                    = size(ODEinput.looky_AdvectionBoundary_import);
[num_fluxes_AdvectionBoundary_export, ~]                    = size(ODEinput.looky_AdvectionBoundary_export);
PhysicsStruct.num_fluxes_advection                          = uint32(ODEinput.num_fluxes_advection);
PhysicsStruct.num_fluxes_AdvectionBoundary_import           = uint32(num_fluxes_AdvectionBoundary_import);
PhysicsStruct.num_fluxes_AdvectionBoundary_export           = uint32(num_fluxes_AdvectionBoundary_export);
[PhysicsStruct.looky_AdvectionFlux, ~]                      = f_unspoolMATRIX_04282020(ODEinput.looky_AdvectionFlux);            % (2D matrix: num_fluxes_advection X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
[PhysicsStruct.looky_AdvectionBoundary_import, ~]           = f_unspoolMATRIX_04282020(ODEinput.looky_AdvectionBoundary_import); % (2D matrix: num_fluxes_BoundaryImport X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
[PhysicsStruct.looky_AdvectionBoundary_export, ~]           = f_unspoolMATRIX_04282020(ODEinput.looky_AdvectionBoundary_export); % (2D matrix: num_fluxes_BoundaryExport X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)])

[PhysicsStruct.HORIZONTALMIXING_compact, ~]                 = f_unspoolMATRIX_04282020(ODEinput.HORIZONTALMIXING_compact);	% (m3/d); (2D matrix: num_t X num_fluxes_HorizontalMixing)
[num_fluxes_HorizontalMixingBoundary_import, ~]             = size(ODEinput.looky_HorizontalMixingBoundary_import);
[num_fluxes_HorizontalMixingBoundary_export, ~]             = size(ODEinput.looky_HorizontalMixingBoundary_export);
PhysicsStruct.num_fluxes_HorizontalMixing                   = uint32(ODEinput.num_fluxes_HorizontalMixing);
PhysicsStruct.num_fluxes_HorizontalMixingBoundary_import	= uint32(num_fluxes_HorizontalMixingBoundary_import);
PhysicsStruct.num_fluxes_HorizontalMixingBoundary_export	= uint32(num_fluxes_HorizontalMixingBoundary_export);
[PhysicsStruct.looky_HorizontalMixingFlux, ~]               = f_unspoolMATRIX_04282020(ODEinput.looky_HorizontalMixingFlux);            % (2D matrix: num_fluxes_HorizontalMixing X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
[PhysicsStruct.looky_HorizontalMixingBoundary_import, ~]	= f_unspoolMATRIX_04282020(ODEinput.looky_HorizontalMixingBoundary_import); % (2D matrix: num_fluxes_HorizontalMixingBoundary_import X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
[PhysicsStruct.looky_HorizontalMixingBoundary_export, ~]	= f_unspoolMATRIX_04282020(ODEinput.looky_HorizontalMixingBoundary_export); % (2D matrix: num_fluxes_HorizontalMixingBoundary_export X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)])

[PhysicsStruct.VERTICALMIXING_compact, ~]                   = f_unspoolMATRIX_04282020(ODEinput.VERTICALMIXING_compact);     % (m3/d); (2D matrix: num_t X num_fluxes_VerticalMixing)
[num_fluxes_VerticalMixingBoundary_import, ~]               = size(ODEinput.looky_VerticalMixingBoundary_import);
[num_fluxes_VerticalMixingBoundary_export, ~]               = size(ODEinput.looky_VerticalMixingBoundary_export);
PhysicsStruct.num_fluxes_VerticalMixing                     = uint32(ODEinput.num_fluxes_VerticalMixing);
PhysicsStruct.num_fluxes_VerticalMixingBoundary_import      = uint32(num_fluxes_VerticalMixingBoundary_import);
PhysicsStruct.num_fluxes_VerticalMixingBoundary_export      = uint32(num_fluxes_VerticalMixingBoundary_export);
[PhysicsStruct.looky_VerticalMixingFlux, ~]                 = f_unspoolMATRIX_04282020(ODEinput.looky_VerticalMixingFlux);            % (2D matrix: num_fluxes_VerticalMixing X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
[PhysicsStruct.looky_VerticalMixingBoundary_import, ~]      = f_unspoolMATRIX_04282020(ODEinput.looky_VerticalMixingBoundary_import); % (2D matrix: num_fluxes_VerticalMixingBoundary_import X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
[PhysicsStruct.looky_VerticalMixingBoundary_export, ~]      = f_unspoolMATRIX_04282020(ODEinput.looky_VerticalMixingBoundary_export); % (2D matrix: num_fluxes_VerticalMixingBoundary_export X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)])

[PhysicsStruct.SINKING_compact, ~]                          = f_unspoolMATRIX_04282020(ODEinput.SINKING_compact);            % (m3/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)
[num_fluxes_SinkingBoundary_import, ~]                      = size(ODEinput.looky_SinkingBoundary_import);
[num_fluxes_SinkingBoundary_export, ~]                      = size(ODEinput.looky_SinkingBoundary_export);
PhysicsStruct.num_fluxes_sinking                            = uint32(ODEinput.num_fluxes_sinking);
PhysicsStruct.num_fluxes_SinkingBoundary_import             = uint32(num_fluxes_SinkingBoundary_import);
PhysicsStruct.num_fluxes_SinkingBoundary_export             = uint32(num_fluxes_SinkingBoundary_export);
[PhysicsStruct.looky_SinkingFlux, ~]                        = f_unspoolMATRIX_04282020(ODEinput.looky_SinkingFlux);            % (2D matrix: num_fluxes_sinking X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
[PhysicsStruct.looky_SinkingBoundary_import, ~]             = f_unspoolMATRIX_04282020(ODEinput.looky_SinkingBoundary_import); % (2D matrix: num_fluxes_SinkingBoundary_import X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
[PhysicsStruct.looky_SinkingBoundary_export, ~]             = f_unspoolMATRIX_04282020(ODEinput.looky_SinkingBoundary_export); % (2D matrix: num_fluxes_SinkingBoundary_export X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)])

[PhysicsStruct.MIGRATION_compact, ~]                        = f_unspoolMATRIX_04282020(ODEinput.MIGRATION_compact);          % (m3/d); (3D matrix: num_t X num_grps X num_fluxes_migration)
[num_fluxes_MigrationBoundary_import, ~]                    = size(ODEinput.looky_MigrationBoundary_import);
[num_fluxes_MigrationBoundary_export, ~]                    = size(ODEinput.looky_MigrationBoundary_export);
PhysicsStruct.num_fluxes_migration                          = uint32(ODEinput.num_fluxes_migration);
PhysicsStruct.num_fluxes_MigrationBoundary_import           = uint32(num_fluxes_MigrationBoundary_import);
PhysicsStruct.num_fluxes_MigrationBoundary_export           = uint32(num_fluxes_MigrationBoundary_export);
[PhysicsStruct.looky_MigrationFlux, ~]                      = f_unspoolMATRIX_04282020(ODEinput.looky_MigrationFlux);            % (2D matrix: num_fluxes_migration X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
[PhysicsStruct.looky_MigrationBoundary_import, ~]           = f_unspoolMATRIX_04282020(ODEinput.looky_MigrationBoundary_import); % (2D matrix: num_fluxes_BoundaryImport X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
[PhysicsStruct.looky_MigrationBoundary_export, ~]           = f_unspoolMATRIX_04282020(ODEinput.looky_MigrationBoundary_export); % (2D matrix: num_fluxes_BoundaryExport X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)])
% *************************************************************************


% end m-file***************************************************************