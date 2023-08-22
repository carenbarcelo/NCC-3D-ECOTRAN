function CompactFlux = f_CompactFluxTimeSeries_11182019(FLUX)
% compact flux time series when arranged as 3D matrix (time X source box X destiny box)
% by Jim Ruzicka
%
% calls:
%       none
%
% takes:
%       FLUX                (3D matrix: time X num_BoxesAndBoundary SOURCE X num_BoxesAndBoundary DESTINY)
%
% returns:
%       CompactFlux
%           compact_flux            (2D matrix: num_t X num_fluxes)
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
%
% revision: 11-18-2019


% *************************************************************************
% STEP 1: function code name ----------------------------------------------
fname_CompactFlux        = mfilename; % save name of this m-file to keep in saved model results
display(['   Running: ' fname_CompactFlux])
% *************************************************************************



% *************************************************************************
% STEP 2: identify linked boxes across full time-series--------------------
[num_t, num_BoxesAndBoundary, ~]	= size(FLUX);
test_sum                        	= sum(FLUX, 1);         % (3D matrix: 1 X num_BoxesAndBoundary X num_BoxesAndBoundary); NOTE: num_BoxesAndBoundary = num_boxes + 1 for boundary
test_sum                           	= squeeze(test_sum);	% (2D matrix: num_BoxesAndBoundary SOURCE X num_BoxesAndBoundary DESTINY)
[source_boxes, destiny_boxes]       = find(test_sum > 0);	%  NOTE: transposed order for source & destiny boxes; (vertical vectors: num_fluxes X 1)
[looky_links_linear]               	= find(test_sum > 0);	% (vertical vector: num_fluxes X 1)
num_fluxes                         	= length(looky_links_linear);

external_flux                       = num_BoxesAndBoundary; % boundary fluxes are in the outer column & layer
looky_domain_import                 = find(source_boxes == external_flux);
looky_domain_export                 = find(destiny_boxes == external_flux);
% *************************************************************************



% *************************************************************************
% test if fluxes exist
if num_fluxes > 0 % if fluxes exist for this physical process, do this...
    
    % *********************************************************************
    % STEP 3: compact the FLUX matrix--------------------------------------
    compact_flux	= zeros(num_t, num_fluxes); % initialize; (2D matrix: num_t X num_fluxes)
    
    for time_loop = 1:num_t
        current_flux                = FLUX(time_loop, :, :);
        compact_flux(time_loop, :)	= current_flux(looky_links_linear); % use looky_links_linear
    end
    % *********************************************************************
    
    
    
    % *********************************************************************
    % STEP 4: identify unique destiny & source boxes-----------------------
    [unique_source, ~, looky_source_boxes]      = unique(source_boxes);   % unique_source = identify each source box with a flux out; looky_source_boxes = order number of each source box within looky_flux
    [unique_destiny, ~, looky_destiny_boxes]	= unique(destiny_boxes);  % unique_destiny = identify each destiny box with a flux in; looky_destiny_boxes = order number of each destiny box within looky_flux 
	% *********************************************************************

    
    
    % *********************************************************************
    % STEP 5: external boundary import & export flux vectors---------------
    if ~isempty(looky_domain_import) % if there are defined import fluxes, do this...
        
        looky_boundary_import       = [destiny_boxes(looky_domain_import) source_boxes(looky_domain_import) looky_domain_import];
        
    else  % if there are NO defined import fluxes, do this...
        
        % add a null import flux
        num_fluxes                                  = num_fluxes + 1;
        compact_flux(:, num_fluxes)                 = 0; % add import flux = 0; (2D matrix: num_t X num_fluxes)
        looky_boundary_import                       = [1 external_flux num_fluxes]; % NOTE: identity of destination box doesn't matter becasue flux = 0
        source_boxes(num_fluxes)                    = external_flux;
        destiny_boxes(num_fluxes)                   = 1;
        [unique_source, ~, looky_source_boxes]      = unique(source_boxes);   % unique_source = identify each source box with a flux out; looky_source_boxes = order number of each source box within looky_flux
        [unique_destiny, ~, looky_destiny_boxes]	= unique(destiny_boxes);  % unique_destiny = identify each destiny box with a flux in; looky_destiny_boxes = order number of each destiny box within looky_flux 
        looky_links_linear(num_fluxes)              = num_BoxesAndBoundary*(num_BoxesAndBoundary-1) + 1;
        
    end % (~isempty(looky_domain_import))
    
    
    if ~isempty(looky_domain_export) % if there are defined export fluxes, do this...
        
        looky_boundary_export       = [destiny_boxes(looky_domain_export) source_boxes(looky_domain_export) looky_domain_export];
        
    else  % if there are NO defined export fluxes, do this...
        
        % add a null export flux
        num_fluxes                                  = num_fluxes + 1;
        compact_flux(:, num_fluxes)                 = 0; % add export flux = 0; (2D matrix: num_t X num_fluxes)
        looky_boundary_export                       = [external_flux 1 num_fluxes];
        source_boxes(num_fluxes)                    = 1;
        destiny_boxes(num_fluxes)                   = external_flux;
        [unique_source, ~, looky_source_boxes]      = unique(source_boxes);   % unique_source = identify each source box with a flux out; looky_source_boxes = order number of each source box within looky_flux
        [unique_destiny, ~, looky_destiny_boxes]	= unique(destiny_boxes);  % unique_destiny = identify each destiny box with a flux in; looky_destiny_boxes = order number of each destiny box within looky_flux 
        looky_links_linear(num_fluxes)              = num_BoxesAndBoundary;

    end % (~isempty(looky_domain_export))
    % *********************************************************************
    
    
else % if fluxes DO NOT exist for this physical process, do this...
    
    % *********************************************************************
    % step 3: make compact_flux a time-series of 0; (2D matrix: num_t X 2)
    num_fluxes = 2;
    compact_flux        = zeros(num_t, 2);          % if no fluxes,     
    source_boxes        = [external_flux; 1];    % (vertical vector: num_fluxes X 1) [import flux, export flux]
    destiny_boxes       = [1; external_flux];    % (vertical vector: num_fluxes X 1) [import flux, export flux]
    looky_links_linear	= [(num_BoxesAndBoundary*(num_BoxesAndBoundary-1) + 1); num_BoxesAndBoundary]; % doesn't matter??
    % *********************************************************************

    
    % *********************************************************************
    % STEP 4: identify unique destiny & source boxes-----------------------
    unique_destiny      = [1; external_flux];
    unique_source       = [1; external_flux];
    looky_destiny_boxes	= [1; 2];
    looky_source_boxes	= [2; 1];
    % *********************************************************************
    
    
    % *********************************************************************
    % STEP 5: external boundary import & export flux vectors---------------
    looky_boundary_import = [1 external_flux 1];
    looky_boundary_export = [external_flux 1 2];
    % *********************************************************************
    
end % (if num_fluxes > 0)



% *************************************************************************
% STEP 6: finalize compact_flux & looky_flux-------------------------------
% filter out negative values 
%        NOTE: these are very small rounding errors in 11th decimal place
looky_negative                  = find(compact_flux < 0);
compact_flux(looky_negative)	= 0;


% build compact_flux address information 
looky_flux	= [destiny_boxes source_boxes looky_links_linear]; % addresses in compact_flux
%             (2D matrix: num_fluxes X 3-->[(destiny box) (source box) (flux address in original destiny X source matrix)])
% *************************************************************************



% *************************************************************************
% STEP 7: pack compact_flux address information----------------------------
CompactFlux.compact_flux            = compact_flux;
CompactFlux.looky_flux              = looky_flux;
CompactFlux.looky_boundary_import	= looky_boundary_import;
CompactFlux.looky_boundary_export	= looky_boundary_export;
CompactFlux.unique_source           = unique_source;
CompactFlux.unique_destiny          = unique_destiny;
CompactFlux.num_fluxes              = num_fluxes;
CompactFlux.fname_CompactFlux       = fname_CompactFlux; % name of this f_CompactFluxTimeSeries function
% *************************************************************************


% end m-file***************************************************************