function ThorntonLessem_parameters = f_TLparameterization_08182023(temp_ref, temp_ref_std, ...
                                                                   looky_grpType1, looky_grpType2, looky_grpType3)
% prepare Thornton-Lessem parameters for living groups
%
% Written by Caren Barcel— & Jim Ruzicka
%
% Variables needed for physiology~temperature relationship with TL.
% The definitions of the eight Thornton-Lessem parameters:
%     (1)	te1     Temperature for xk1 (in ºC)
%     (2)	te2     Temperature for xk2 (in ºC)
%     (3)	te3     Temperature for xk3 (in ºC)
%     (4)	te4     Temperature for xk4 (in ºC)
%     (5)	xk1     Proportion of CMAX at te1
%     (6)	xk2     Proportion of CMAX at te2
%     (7)	xk3     Proportion of CMAX at te3
%     (8)	xk4     Proportion of CMAX at te4
%
% calls: none
%
% takes:
%       temp_ref         	mean temperature across desired reference time period                	(3D matrix: 1 X num_grps X num_domains) NOTE: we only need to consider geographic regions 1-15
%       temp_ref_std        standard deviation of temperatures across desired referencetime period	(3D matrix: 1 X num_grps X num_domains) NOTE: we only need to consider geographic regions
%       looky_grpType1      addresses of groups classified as surface dwellers                      (horizontal vector)
%       looky_grpType2      addresses of groups classified as mid-water column dwellers             (horizontal vector)
%       looky_grpType3      addresses of groups classified as benthic dwellers                      (horizontal vector)
%
% returns:
%       ThorntonLessem_parameters    The eight Thornton-Lessem parameters                           (2D matrix: 8 X num_grps)
%
% revision date: 5/18/2023 
%       5/18/2023 cleaned up notation


%% *************************************************************************
% STEP 1: T-L set up based on reference temperatures during first part of time series (1980-2005)

% Mean depth specific domain wide temperature and STD
temp_ref        = mean(temp_ref, 3, 'omitnan');     % mean across domains; (horizontal vector: 1 X num_grps)   
temp_ref_std	= mean(temp_ref_std, 3, 'omitnan');	% mean across domains; (horizontal vector: 1 X num_grps)

% Calculate range of SD +/- around mean temperature
for ii = 1:1:10
    temp_ref_std_test       = temp_ref_std * ii;
    temp_ref_plus(ii, :)	= temp_ref + temp_ref_std_test;
    temp_ref_minus(ii, :)	= temp_ref - temp_ref_std_test;
end
% *************************************************************************



%% *************************************************************************
% STEP 2: Set up T-L parameter set based on above specified conditions
%         First: Define targeted T-L set up for each group occupying each depth

% Default set up used in Barcelo et al.
TLsizeBase_Type1 = 4;
TLsizeTop_Type1 = 1;

TLsizeBase_Type2 = 4;
TLsizeTop_Type2 = 1;

TLsizeBase_Type3 = 6;
TLsizeTop_Type3 = 3;


% Alternate parameter sets for sensitivty analyses

% % 1sd less at base means
% TLsizeBase_Type1 = 3;
% TLsizeTop_Type1 = 1;
% 
% TLsizeBase_Type2 = 3;
% TLsizeTop_Type2 = 1;
% 
% TLsizeBase_Type3 = 5;
% TLsizeTop_Type3 = 3;


% % 2sd less at base means
% TLsizeBase_Type1 = 2;
% TLsizeTop_Type1 = 1;
% 
% TLsizeBase_Type2 = 2;
% TLsizeTop_Type2 = 1;
% 
% TLsizeBase_Type3 = 4;
% TLsizeTop_Type3 = 3;


% % 3sd less at base means
% TLsizeBase_Type1 = 1;
% TLsizeTop_Type1 = 1;
% 
% TLsizeBase_Type2 = 1;
% TLsizeTop_Type2 = 1;
% 
% TLsizeBase_Type3 = 3;
% TLsizeTop_Type3 = 3;


% % 1std more at base means
% TLsizeBase_Type1 = 5;
% TLsizeTop_Type1 = 1;
% 
% TLsizeBase_Type2 = 5;
% TLsizeTop_Type2 = 1;
% 
% TLsizeBase_Type3 = 7;
% TLsizeTop_Type3 = 3;


% % 2sd more at base means
% TLsizeBase_Type1 = 6;
% TLsizeTop_Type1 = 1;
% 
% TLsizeBase_Type2 = 6;
% TLsizeTop_Type2 = 1;
% 
% TLsizeBase_Type3 = 8;
% TLsizeTop_Type3 = 3;


% % 3sd more at base means
% TLsizeBase_Type1 = 7;
% TLsizeTop_Type1 = 1;
% 
% TLsizeBase_Type2 = 7;
% TLsizeTop_Type2 = 1;
% 
% TLsizeBase_Type3 = 9;
% TLsizeTop_Type3 = 3;
% *************************************************************************



%% *************************************************************************
% Step 3: 
% Set all groups to 1SD at top and 4SD at bottom
ThorntonLessem_parameters                       = cat(1, ...
                                                        temp_ref_minus(TLsizeBase_Type1, :), ...
                                                        temp_ref_minus(TLsizeTop_Type1, :), ...
                                                        temp_ref_plus(TLsizeTop_Type1, :), ...
                                                        temp_ref_plus(TLsizeBase_Type1, :), ...
                                                        repelem(0.01, 106), ...
                                                        repelem(0.98,106), ...
                                                        repelem(0.98, 106), ...
                                                        repelem(0.01, 106));

% **If needed** assign group specific T-L curves, specify with these calls:
% SURFACE GROUPS: Broaden size of T-L curve to XXSD at top and XXSD at bottom for benthic groups (indexed by looky_grpType3)
ThorntonLessem_parameters(:, looky_grpType1)    = cat(1, ...
                                                        temp_ref_minus(TLsizeBase_Type1, looky_grpType1), ...
                                                        temp_ref_minus(TLsizeTop_Type1, looky_grpType1), ...
                                                        temp_ref_plus(TLsizeTop_Type1, looky_grpType1), ...
                                                        temp_ref_plus(TLsizeBase_Type1, looky_grpType1), ...
                                                        repelem(0.01, size(looky_grpType1, 2)), ...
                                                        repelem(0.98, size(looky_grpType1, 2)), ...
                                                        repelem(0.98, size(looky_grpType1, 2)), ...
                                                        repelem(0.01, size(looky_grpType1, 2)));

% MIDWATER GROUPS: Broaden size of T-L curve to XXSD at top and XXSD at bottom for benthic groups (indexed by looky_grpType3)
ThorntonLessem_parameters(:, looky_grpType2)    = cat(1, temp_ref_minus(TLsizeBase_Type2, looky_grpType2), ...
                                                        temp_ref_minus(TLsizeTop_Type2, looky_grpType2), ...
                                                        temp_ref_plus(TLsizeTop_Type2, looky_grpType2), ...
                                                        temp_ref_plus(TLsizeBase_Type2, looky_grpType2), ...
                                                        repelem(0.01, size(looky_grpType2, 2)), ...
                                                        repelem(0.98, size(looky_grpType2, 2)), ...
                                                        repelem(0.98, size(looky_grpType2, 2)), ...
                                                        repelem(0.01, size(looky_grpType2, 2)));

% BENTHIC GROUPS: Broaden size of T-L curve to XxSD at top and XXSD at bottom for benthic groups (indexed by looky_grpType3)
ThorntonLessem_parameters(:, looky_grpType3)    = cat(1, temp_ref_minus(TLsizeBase_Type3, looky_grpType3), ...
                                                        temp_ref_minus(TLsizeTop_Type3, looky_grpType3), ...
                                                        temp_ref_plus(TLsizeTop_Type3, looky_grpType3), ...
                                                        temp_ref_plus(TLsizeBase_Type3, looky_grpType3), ...
                                                        repelem(0.01, size(looky_grpType3, 2)), ...
                                                        repelem(0.98, size(looky_grpType3, 2)), ...
                                                        repelem(0.98, size(looky_grpType3, 2)), ...
                                                        repelem(0.01, size(looky_grpType3, 2)));

% ThLe *JUST* FOR FISH
ThorntonLessem_parameters(:, [1:20, 59:106])    = 1; % set non-fish groups to 1

end
