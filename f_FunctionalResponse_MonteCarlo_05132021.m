function [FunctionalResponseParams, fname_FunctionalResponse_MonteCarlo] = f_FunctionalResponse_MonteCarlo_05132021(ECOTRAN, ODEinput, FunctionalResponse_CV)
% prepare vulnerability array
% by Jim Ruzicka
% NOTE: producer "vulnerabilities", m_p, are a characteristic of the CONSUMER
% NOTE: MC layer 1 is the mean "type" set of parameters
%
% calls:
%       none
%
% takes:
% 	ECOTRAN
%       FunctionalResponseParams	(2D matrix: num_grps->producers X 4)
%       num_MC
%       num_grps
%   ODEinput
%       looky_nutrients
%       looky_ANYdetritus
%	FunctionalResponse_CV           (2D matrix: num_grps->CONSUMERS X 4)
%
% returns:
%	FunctionalResponseParams        (producer "vulnerabilities", m_p) (3D matrix: CONSUMERS X prey group X num_MC) replicated across clms (= producers)
%
% revision date: 5-13-2021
%               (5-13-2021 JR added function file name to returned variables)


% *************************************************************************
% STEP 1: set operating conditions-----------------------------------------

% step 1a: log function file name -----------------------------------------
fname_FunctionalResponse_MonteCarlo	= mfilename; % name of this f_FunctionalResponse_MonteCarlo function
display(['Running: ' fname_FunctionalResponse_MonteCarlo])
% -------------------------------------------------------------------------

% step 1b: set distribution for random variables --------------------------
distribution_type           = 'normal';	% 'uniform'
m_p_min                     = 0;    	% minimum vulnerability parameter (m_p); arbitrary
m_p_max                     = 100;    	% maximum vulnerability parameter (m_p); arbitrary
% *************************************************************************



% *************************************************************************
% STEP 2: unpack variables-------------------------------------------------
FunctionalResponseParams    = ECOTRAN.FunctionalResponseParams; % vulnerability parameter (m_p); (2D matrix: num_grps->producers X 4); NOTE: first column is used (other columns for future applications)
num_MC                      = ECOTRAN.num_MC;
num_grps                    = ECOTRAN.num_grps;
looky_nutrients             = ODEinput.looky_nutrients;
looky_ANYdetritus          	= ODEinput.looky_ANYdetritus;
% *************************************************************************



% *************************************************************************
% STEP 2: define functional predator-prey relations------------------------
% step 2a: select a "linear" or "ECOSIM" prey vulnerability ---------------
%          NOTE: prey vulnerability (m_p) is the first column of the four potential FunctionalResponseParams;
%          NOTE: to change one half-sat constant for individual groups, change m_p([Grp_rowclm])
%                m_p = 0 is "linear"
%                m_p = 1 is "non-linear" ECOSIM default
%          FFF: other three columns are for future options
m_p                         = FunctionalResponseParams(:, 1);	% (vertical vector: num_grps X 1)
m_p_CV                      = FunctionalResponse_CV(:, 1);      % (CV); (vertical vector: num_grps X 1)
% -------------------------------------------------------------------------

% step 2b: replicate for each MC model ------------------------------------
m_p                         = repmat(m_p, [1, 1, num_MC]);      % (3D matrix: num_grps X 1 X num_MC)
m_p_CV                      = repmat(m_p_CV, [1, 1, num_MC]);	% (CV); (3D matrix: num_grps X 1 X num_MC)
% *************************************************************************



% *************************************************************************
% STEP 3: generate random functional response vector-----------------------

% OPTION 1: random UNIFORM distribution -----------------------------------
if strcmp(distribution_type, 'uniform')
    
    % step 3a: interval over which to draw random vulnerabilities ---------
    interval_low_m_p                            = m_p - (m_p .* m_p_CV);            % (3D matrix: num_grps X 1 X num_MC)
    interval_high_m_p                           = m_p + (m_p .* m_p_CV);            % (3D matrix: num_grps X 1 X num_MC)
    % ---------------------------------------------------------------------

    % step 3b: filter out extreme values (m_p) ----------------------------
    looky_small                                 = find(interval_low_m_p < m_p_min);
    interval_low_m_p(looky_small)               = m_p_min;                          % (3D matrix: num_grps X 1 X num_MC)
    looky_big                                   = find(interval_high_m_p > m_p_max); 
    interval_high_m_p(looky_big)                = m_p_max;                          % (3D matrix: num_grps X 1 X num_MC)
    % ---------------------------------------------------------------------

    % step 3c: random vulnerabilities-- -----------------------------------
    random_vector                               = rand(num_grps, 1, num_MC);	% draw from UNIFORM distribution; (3D matrix: num_grps X 1 X num_MC)
    MonteCarloVector                            = interval_low_m_p + ((interval_high_m_p - interval_low_m_p) .* random_vector); % produce the new random m_p vector; (3D matrix: num_grps X 1 X num_MC)
    % ---------------------------------------------------------------------

    % step 3d: set nutrient & detritus vulnerabilites (m_p) to linear -----
    MonteCarloVector(looky_nutrients, 1, :) 	= 0;	% set so that flow to nutrients is linear (should already be done, but just in case); (3D matrix: num_grps X 1 X num_MC)
    MonteCarloVector(looky_ANYdetritus, 1, :)	= 0;	% set so that flow to detritus is linear (should already be done, but just in case); (3D matrix: num_grps X 1 X num_MC)
    % ---------------------------------------------------------------------

% OPTION 2: random NORMAL distribution ------------------------------------
elseif strcmp(distribution_type, 'normal')
    
    % step 3a: interval over which to draw random vulnerabilities ---------
    interval_low_m_p                            = zeros(num_grps, 1, num_MC);	% (3D matrix: num_grps X 1 X num_MC)
    interval_high_m_p                           = (m_p .* m_p_CV);              % (3D matrix: num_grps X 1 X num_MC); NOTE: this is the STDev when using NORMAL distribution
    % ---------------------------------------------------------------------

    % step 3b: random vulnerabilities (m_p) -------------------------------
    random_vector                            	= randn(num_grps, 1, num_MC);	% draw from NORMAL distribution; (3D matrix: num_grps X 1 X num_MC)
    MonteCarloVector                          	= m_p + (interval_high_m_p .* random_vector); % produce the new random m_p vector; (3D matrix: num_grps X 1 X num_MC)
    % ---------------------------------------------------------------------

    % step 3c: filter out extreme values (m_p) ----------------------------
    looky_small                               	= find(MonteCarloVector < m_p_min);
    MonteCarloVector(looky_small)             	= m_p(looky_small);	% set too small values to the defined mean
    looky_big                                	= find(MonteCarloVector > m_p_max); 
    MonteCarloVector(looky_big)               	= m_p(looky_big);	% set too large values to defined mean
    % ---------------------------------------------------------------------

    % step 3d: set nutrient & detritus vulnerabilites (m_p) to linear -----
    MonteCarloVector(looky_nutrients, 1, :) 	= 0;	% set so that flow to nutrients is linear (should already be done, but just in case); (3D matrix: num_grps X 1 X num_MC)
    MonteCarloVector(looky_ANYdetritus, 1, :)	= 0;	% set so that flow to detritus is linear (should already be done, but just in case); (3D matrix: num_grps X 1 X num_MC)
    % ---------------------------------------------------------------------

end % end (if strcmp(distribution_type, 'uniform'))
% *************************************************************************



% *************************************************************************
% STEP 4: finalize vulnerabilities (m_p)-----------------------------------
m_p_repmat                          = repmat(MonteCarloVector, [1, num_grps, 1]);	% replicate m_p across clms; (3D matrix: num_grps->CONSUMERS X num_grps->PRODUCERS X num_MC)
m_p_repmat(:, :, 1)                 = repmat(m_p(:, 1, 1), [1, num_grps, 1]);       % define MC layer 1 as the mean "type" set of parameters; (3D matrix: num_grps->CONSUMERS X num_grps->PRODUCERS X num_MC)
m_p_repmat(:, looky_nutrients, :)   = 0;            % set so that flow from nutrients is linear; (3D matrix: num_grps->CONSUMERS X num_grps->PRODUCERS X num_MC)
m_p_repmat(:, looky_ANYdetritus, :)	= 0;            % set so that flow from detritus pools is linear; (3D matrix: num_grps->CONSUMERS X num_grps->PRODUCERS X num_MC)
FunctionalResponseParams            = m_p_repmat;	% pack to send to ODE; (3D matrix: num_grps->CONSUMERS X num_grps->PRODUCERS X num_MC)
% *************************************************************************


% end m-file***************************************************************