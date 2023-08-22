function NetFlux = f_calcNetFlux_12112019(compact_flux_t, CompactFlux)
% calculate net flux into and net flux out of each model box and across outer domain boundaries
% by Jim Ruzicka
%
% calls:
%       none
%
% takes:
%       compact_flux_t              (horizontal vector: 1 X num_fluxes)
%       CompactFlux
%           looky_flux              (2D matrix: num_fluxes X 5)
%                                       clm 1: (destiny box) = list of boxes importing water volume (+ boundary)
%                                       clm 2: (source box) = list of boxes exporting water volume (+ boundary)
%                                       clm 3: (destiny box address) = addresses of destiny boxes in compact_flux, tells you to which destiny box each flux is to be added, (order number of each unique destiny box within looky_flux)
%                                       clm 4: (source box address) = addresses of source boxes in compact_flux, tells you to which source box each flux is to be added, (order number of each unique source box within looky_flux)
%                                       clm 5: flux address in non-compacted flux 2D matrix: destiny X source
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
%
% returns:
%       NetFlux
%           flux_import_t           net flux INTO each box                 (horizontal vector: 1 X num_boxes+1); NOTE: +1 to account for fluxes across external domain boundary; NOTE: incudes (adds in) all fluxes across domain boundary
%           flux_export_t           net flux OUT OF each box               (horizontal vector: 1 X num_boxes+1); NOTE: +1 to account for fluxes across external domain boundary; NOTE: incudes (adds in) all fluxes across domain boundary
%           flux_domain_import_t	flux INTO each box from outside domain (horizontal vector: 1 X num_boxes+1); NOTE: +1 to account for fluxes across external domain boundary
%           flux_domain_export_t	flux OUT OF each box to outside domain (horizontal vector: 1 X num_boxes+1); NOTE: +1 to account for fluxes across external domain boundary
%           fname_CalcNetFlux       name of this f_CalcNetFlux function
%
% revision: 12-11-2019


% *************************************************************************
% STEP 1: unpack variables-------------------------------------------------
fname_calcNetFlux       = mfilename; % name of this f_calcNetFlux function
% display(['   Running: ' fname_CalcNetFlux])

looky_flux              = CompactFlux.looky_flux;            % (2D matrix: num_fluxes_advection X 3-->[(destiny box) (source box) (flux address in non-compacted flux 2D matrix: destiny X source)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
looky_boundary_import	= CompactFlux.looky_boundary_import;
looky_boundary_export	= CompactFlux.looky_boundary_export;
% *************************************************************************



% *************************************************************************
% STEP 2: initialize output vectors----------------------------------------
num_boxes             	= max(looky_flux(:, 1)) - 1; % the largest box number is the boundary, subtract 1 to get num_boxes
flux_import_t           = zeros(1, (num_boxes+1));   % initialize; (horizontal vector: 1 X num_boxes+1); NOTE: +1 to account for fluxes across external domain boundary
flux_export_t           = zeros(1, (num_boxes+1));   % initialize; (horizontal vector: 1 X num_boxes+1); NOTE: +1 to account for fluxes across external domain boundary
flux_domain_import_t	= zeros(1, (num_boxes+1));   % initialize; (horizontal vector: 1 X num_boxes+1); NOTE: +1 to account for fluxes across external domain boundary
flux_domain_export_t	= zeros(1, (num_boxes+1));   % initialize; (horizontal vector: 1 X num_boxes+1); NOTE: +1 to account for fluxes across external domain boundary
% *************************************************************************



% *************************************************************************
% STEP 3: initialize output vectors----------------------------------------
[looky_destiny, ~, c1]	= unique(looky_flux(:, 1)); % looky_destiny = identify each destiny box with a flux IN; c1 = order number of each unique destiny box within looky_flux, same as clm 3 in looky_flux
[looky_source,  ~, c2]	= unique(looky_flux(:, 2)); % looky_source = identify each source box with a flux out; c2 = order number of each unique source box within looky_flux, same as clm 4 in looky_flux
% *************************************************************************



% *************************************************************************
% STEP 4: calculate net fluxes for each source & destiny box---------------
sum_import              = [looky_destiny, accumarray(c1, compact_flux_t(1, :))]; % total flux INTO each destiny box; (m3/d); (2D matrix: num_fluxes X [(destiny box) (flux rate)])
sum_export              = [looky_source,  accumarray(c2, compact_flux_t(1, :))]; % total flux OUT OF each source box; (m3/d); (2D matrix: num_fluxes X [(source box) (flux rate)])
% *************************************************************************



% *************************************************************************
% STEP 5: build up import & export flux vectors----------------------------
%         put fluxes of sum_import & sum_export into the appropriate box order and set non-flux boxes to 0
flux_import_t(sum_import(:, 1))	= sum_import(:, 2); % net flux INTO each box; (m3/d); (horizontal vector: 1 X num_boxes+1); NOTE: +1 to account for fluxes across external domain boundary
flux_export_t(sum_export(:, 1))	= sum_export(:, 2); % net flux OUT OF each box;	(m3/d); (horizontal vector: 1 X num_boxes+1); NOTE: +1 to account for fluxes across external domain boundary
% *************************************************************************



% *************************************************************************
% STEP 6: external boundary import & export flux vectors-------------------
if ~isempty(looky_boundary_import)
    flux_domain_import_t(looky_boundary_import(:, 1))	= compact_flux_t(1, looky_boundary_import(:, 3)); % net flux across domain boundary INTO each box; (m3/d); (horizontal vector: 1 X num_boxes+1); NOTE: +1 to account for fluxes across external domain boundary
end % (~isempty(looky_boundary_import))

if ~isempty(looky_boundary_export)
    flux_domain_export_t(looky_boundary_export(:, 2))	= compact_flux_t(1, looky_boundary_export(:, 3)); % net flux across domain boundary OUT OF each box; (m3/d); (horizontal vector: 1 X num_boxes+1); NOTE: +1 to account for fluxes across external domain boundary
end % (~isempty(looky_boundary_export))
% *************************************************************************



% *************************************************************************
% STEP 7: pack results-----------------------------------------------------
NetFlux.flux_import_t           = flux_import_t;      	% net flux INTO each box                 (horizontal vector: 1 X num_boxes+1); NOTE: +1 to account for fluxes across external domain boundary; NOTE: incudes (adds in) all fluxes across domain boundary
NetFlux.flux_export_t           = flux_export_t;      	% net flux OUT OF each box               (horizontal vector: 1 X num_boxes+1); NOTE: +1 to account for fluxes across external domain boundary; NOTE: incudes (adds in) all fluxes across domain boundary
NetFlux.flux_domain_import_t	= flux_domain_import_t;	% flux INTO each box from outside domain (horizontal vector: 1 X num_boxes+1); NOTE: +1 to account for fluxes across external domain boundary
NetFlux.flux_domain_export_t	= flux_domain_export_t;	% flux OUT OF each box to outside domain (horizontal vector: 1 X num_boxes+1); NOTE: +1 to account for fluxes across external domain boundary
NetFlux.fname_calcNetFlux       = fname_calcNetFlux;	% name of this f_CalcNetFlux function
% *************************************************************************


% end m-file***************************************************************