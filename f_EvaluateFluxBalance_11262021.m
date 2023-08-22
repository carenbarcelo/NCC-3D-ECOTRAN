function FluxImbalance = f_EvaluateFluxBalance_11262021(UnCompactFlux)
% examine for flux time-series for imbalances IN & OUT of invidual boxes and IN & OUT of the overall domain
% by Jim Ruzicka
%
% calls:
%       none
%
% takes:
%   UnCompactFlux
%       flux_imbalance          time-series of flux imbalances for individual boxes     (2D matrix: num_t X num_boxes)
%       flux_domain_imbalance   time-series of flux imbalances for whole domain         (vertical vector: num_t X 1)
%
% returns:
%   FluxImbalance
%       import_imbalance_time           list of time points of IMPORT imbalances
%       import_imbalance_box            list of boxes points with IMPORT imbalances
%       export_imbalance_time           list of time points of EXPORT imbalances
%       export_imbalance_box            list of boxes points with EXPPORT imbalances
%       domain_import_imbalance_time	list of time points of whole domain IMPORT imbalances
%       domain_export_imbalance_time	list of time points of whole domain EXPORT imbalances
%       fname_EvaluateFluxBalance       name of this f_EvaluateFluxBalance sub-function
%
% revision date: 11-26-2021


% *************************************************************************
% STEP 1: unpack variables-------------------------------------------------
fname_EvaluateFluxBalance	= mfilename; % save name of this m-file to keep in saved model results
display(['   Running: ' fname_EvaluateFluxBalance])

flux_imbalance              = UnCompactFlux.flux_imbalance;
flux_domain_imbalance       = UnCompactFlux.flux_domain_imbalance;
% *************************************************************************



% *************************************************************************
% STEP 2: find any flux imbalances-----------------------------------------
[import_imbalance_time, import_imbalance_box]	= find(round(flux_imbalance, 2) > 0);
[export_imbalance_time, export_imbalance_box]	= find(round(flux_imbalance, 2) < 0);

import_imbalance_time                           = unique(import_imbalance_time);
import_imbalance_box                            = unique(import_imbalance_box);
export_imbalance_time                           = unique(export_imbalance_time);
export_imbalance_box                            = unique(export_imbalance_box);


[domain_import_imbalance_time]                  = find(round(flux_domain_imbalance, 1) > 0);
[domain_export_imbalance_time]                  = find(round(flux_domain_imbalance, 1) < 0);

domain_import_imbalance_time                    = unique(domain_import_imbalance_time);
domain_export_imbalance_time                    = unique(domain_export_imbalance_time);
% *************************************************************************



% *************************************************************************
% STEP 3: display imbalance warnings---------------------------------------

if (isempty(import_imbalance_time)  && isempty(export_imbalance_time) && isempty(domain_import_imbalance_time) && isempty(domain_export_imbalance_time))
    display('   -->Flux time series is OK. No imbalances calculated')
end


if ~isempty(import_imbalance_time)
    display(['   -->WARNING: flux import imbalance at t = ' num2str(import_imbalance_time')])
    display(['                                  and box = ' num2str(import_imbalance_box')])

end


if ~isempty(export_imbalance_time)
    display(['   -->WARNING: flux export imbalance at t = ' num2str(export_imbalance_time')])
    display(['                                  and box = ' num2str(export_imbalance_box')])
end


if ~isempty(domain_import_imbalance_time)
    display(['   -->WARNING: total import imbalance INTO DOMAIN at t = ' num2str(domain_import_imbalance_time')])
end


if ~isempty(domain_export_imbalance_time)
    display(['   -->WARNING: total export imbalance OUT OF DOMAIN at t = ' num2str(domain_export_imbalance_time')])
end
% *************************************************************************



% *************************************************************************
% STEP 4: pack up error information----------------------------------------
FluxImbalance.import_imbalance_time         = import_imbalance_time;
FluxImbalance.import_imbalance_box          = import_imbalance_box;
FluxImbalance.export_imbalance_time         = export_imbalance_time;
FluxImbalance.export_imbalance_box          = export_imbalance_box;

FluxImbalance.domain_import_imbalance_time	= domain_import_imbalance_time;
FluxImbalance.domain_export_imbalance_time	= domain_export_imbalance_time;
FluxImbalance.fname_EvaluateFluxBalance     = fname_EvaluateFluxBalance;	% name of this f_EvaluateFluxBalance sub-function
% *************************************************************************


% end m-file***************************************************************