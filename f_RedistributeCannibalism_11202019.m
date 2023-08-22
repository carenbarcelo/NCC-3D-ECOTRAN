function [DIET_NoCannibalism, cannibalism_vector] = f_RedistributeCannibalism_11202019(DIET)
% Get rid of the non-zero diag terms; redistribute cannibalism terms among
% all other non-zero terms for each predator (columns)
%
% calls:
%       none
%
% takes:
%       DIET
%
% returns:
%       DIET_NoCannibalism          (matrix)
%       cannibalism_vector          (horizontal vector)
%
% revision date: 11-20-2019


% *************************************************************************
% STEP 1: measure size of DIET matrix--------------------------------------
fname_RedistributeCannibalism     = mfilename; % save name of this m-file to keep in saved model results
display(['   Running: ' fname_RedistributeCannibalism])

[rows, clms] = size(DIET);
if rows ~= clms
    error('ERROR: diet matrix must be square')
end
% *************************************************************************



% *************************************************************************
% STEP 2: retain cannibalism vector----------------------------------------
cannibalism_vector = diag(DIET)'; % (horizontal vector: 1 X clms)
% *************************************************************************



% *************************************************************************
% STEP 3: set cannibalism diagonal to 0 in DIET_NoCannibalism--------------
DIET_NoCannibalism = DIET;
for row_loop = 1:rows
    DIET_NoCannibalism(row_loop, row_loop) = 0;
end
% *************************************************************************



% *************************************************************************
% STEP 4: rescale each clm in DIET_NoCannibalism to 1----------------------
DIET_NoCannibalism                              = DIET_NoCannibalism ./ repmat(sum(DIET_NoCannibalism, 1), [rows 1]); % (2D matrix: rows X clms)
DIET_NoCannibalism(isnan(DIET_NoCannibalism))   = 0; % fix div/0 errors (phyto & detritus diets)
% *************************************************************************


% end m-file***************************************************************