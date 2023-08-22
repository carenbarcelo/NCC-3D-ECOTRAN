function UnCompactFlux = f_UnCompactFluxTimeSeries_12112019(CompactFlux)
% UnCompact a flux time series to provide IMPORT & EXPORT fluxes for each box and the domain as a whole
% by Jim Ruzicka
%
% calls:
%       f_calcNetFlux_12112019      calculate net flux into and net flux out of each model box and across outer domain boundaries
%
% takes:
%       CompactFlux
%           compact_flux            (2D matrix: num_t X num_fluxes)
%                                       NOTE: fluxes include all linked boxes +1 for external links
%           looky_flux              (2D matrix: num_fluxes X 5)
%                                       clm 1: (destiny box) = list of boxes importing water volume (+ boundary)
%                                       clm 2: (source box) = list of boxes exporting water volume (+ boundary)
%                                       clm 3: (destiny box address) = addresses of destiny boxes in compact_flux
%                                       clm 4: (source box address) = addresses of source boxes in compact_flux
%                                       clm 5: flux address in non-compacted flux 2D matrix: destiny X source
%                                       NOTE: fluxes include all linked boxes +1 for external links
%                                       NOTE: values constant independent of t
%           looky_boundary_import	(2D matrix: num_fluxes_BoundaryImport X 3)
%                                       clm 1: (destiny box) = identity of boxes importing water volume
%                                       clm 2: (source box) = identity of boxes exporting water volume (always the boundary flux number)
%                                       clm 3: (destiny box address) = addresses of destiny boxes in compact_flux
%           looky_boundary_export   (2D matrix: num_fluxes_BoundaryExport X 3)
%                                       clm 1: (destiny box) = identity of boxes importing water volume (always the boundary flux number)
%                                       clm 2: (source box) = identity of boxes exporting water volume
%                                       clm 3: (source box address) = addresses of source boxes in compact_flux
%           num_fluxes              number of realized fluxes between boxes (or boundary) over full time-series
%
% returns:
%       UnCompactFlux
%           flux_import                     flux INTO each box (includes external domain fluxes) (2D matrix: num_t X num_boxes)
%           flux_export                     flux OUT OF each box (includes external domain fluxes) (2D matrix: num_t X num_boxes)
%           flux_imbalance                  imbalance of fluxes for each box (includes and external domain fluxes added into each box-specific flux) (2D matrix: num_t X num_boxes)
%           flux_domain_import              total flux INTO domain from outside (vertical vector: num_t X 1)
%           flux_domain_export              total flux OUT OF domain to outside (vertical vector: num_t X 1)
%           flux_domain_imbalance           total flux imbalance across all domain boundaries (vertical vector: num_t X 1)
%           num_t                           number of time steps (scaler)
%           num_boxes                       number of boxes (scaler)
%           fname_UnCompactFluxTimeSeries	name of this f_UnCompactFluxTimeSeries sub-function
%           fname_CalcNetFlux               name of this f_CalcNetFlux function
%
% revision date: 12-11-2019


% *************************************************************************
% STEP 1: unpack variables-------------------------------------------------
fname_UnCompactFluxTimeSeries	= mfilename; % save name of this m-file to keep in saved model results
display(['   Running: ' fname_UnCompactFluxTimeSeries])

compact_flux                        = CompactFlux.compact_flux; % (2D matrix: num_t X num_fluxes); NOTE: fluxes include all linked boxes +1 for external links
looky_flux                          = CompactFlux.looky_flux;   % (2D matrix: num_fluxes_advection X 3-->[(destiny box) (source box) (flux address in non-compacted flux 2D matrix: destiny X source)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
% *************************************************************************



% *************************************************************************
% STEP 2: intialize variables----------------------------------------------
num_boxes                   = max(looky_flux(:, 1)) - 1; % the largest box number is the boundary, subtract 1 to get num_boxes
[num_t, ~]                  = size(compact_flux);

flux_import                 = zeros(num_t, (num_boxes+1)); % flux INTO each box (includes external domain fluxes)
flux_export                 = zeros(num_t, (num_boxes+1)); % flux OUT OF each box (includes external domain fluxes)
flux_imbalance              = zeros(num_t, (num_boxes+1)); % imbalance of fluxes for each box (includes external domain fluxes)

flux_domain_import          = zeros(num_t, 1); % total flux INTO domain from outside
flux_domain_export          = zeros(num_t, 1); % total flux OUT OF domain to outside
flux_domain_imbalance       = zeros(num_t, 1); % total flux imbalance across all domain boundaries
% *************************************************************************



% *************************************************************************
% STEP 3: build flux timeseries--------------------------------------------

for time_loop = 1:num_t
    
    compact_flux_t                      = compact_flux(time_loop, :); % flux @ t; (horizontal vector: 1 X num_fluxes); NOTE: fluxes include all linked boxes +1 for external links
    
    [NetFlux]                           = f_calcNetFlux_12112019(compact_flux_t, CompactFlux); % calculate net flux into and net flux out of each model box and across outer domain boundaries
    
    flux_import(time_loop, :)           = NetFlux.flux_import_t;
    flux_export(time_loop, :)           = NetFlux.flux_export_t;
    flux_imbalance(time_loop, :)        = NetFlux.flux_import_t - NetFlux.flux_export_t;
    
    flux_domain_import(time_loop)       = sum(NetFlux.flux_domain_import_t);
    flux_domain_export(time_loop)       = sum(NetFlux.flux_domain_export_t);
    flux_domain_imbalance(time_loop)	= flux_domain_import(time_loop) - flux_domain_export(time_loop);
    
end
% *************************************************************************



% *************************************************************************
% STEP 4: pack UnCompact timeseries----------------------------------------
%         NOTE: trimming off the boundary "box" column
UnCompactFlux.flux_import           = flux_import(:, 1:(end-1));  % flux INTO each box (includes external domain fluxes) (2D matrix: num_t X num_boxes)
UnCompactFlux.flux_export           = flux_export(:, 1:(end-1));  % flux OUT OF each box (includes external domain fluxes) (2D matrix: num_t X num_boxes)
UnCompactFlux.flux_imbalance        = flux_imbalance(:, 1:(end-1)); % imbalance of fluxes for each box (includes and external domain fluxes added into each box-specific flux) (2D matrix: num_t X num_boxes)

UnCompactFlux.flux_domain_import    = flux_domain_import;       % total flux INTO domain from outside (vertical vector: num_t X 1)
UnCompactFlux.flux_domain_export	= flux_domain_export;       % total flux OUT OF domain to outside (vertical vector: num_t X 1)
UnCompactFlux.flux_domain_imbalance	= flux_domain_imbalance;    % total flux imbalance across all domain boundaries (vertical vector: num_t X 1)

UnCompactFlux.num_t                 = num_t;                    % number of time steps (scaler)
UnCompactFlux.num_boxes             = num_boxes;                % number of boxes (scaler)

UnCompactFlux.fname_UnCompactFluxTimeSeries	= fname_UnCompactFluxTimeSeries; % name of this f_UnCompactFluxTimeSeries sub-function
UnCompactFlux.fname_CalcNetFlux     = NetFlux.fname_calcNetFlux;	% name of this fname_calcNetFlux function
% *************************************************************************


% end m-file***************************************************************