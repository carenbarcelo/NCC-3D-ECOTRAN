function PredationBudget = f_calcPredationBudget_12102019(TotalConsumptionVector, DIET)
% Calculate the predation budget (a.k.a.: bottom-up matrix, inverse transposition, production matrix)
% A_cp = (D_pc * q_c) / sigma(D_pc * q_c) for all c; Eq. 9 in Steele & Ruzicka (2011)
%
% 	A_cp						= PredationMatrix; (2D matrix)
% 	D_pc						= DIET matrix; (2D matrix: num_grps X num_grps)
% 	q_c							= total consumption by c; (vector)
%   (D_pc * q_c)                = CONSUMPTION matrix: (t/km2/y); (2D matrix: num_grps X num_grps; producers X consumers)
% 	sigma(D_pc * q_c) for all c	= predation_total upon p by all consumers c
% 	p							= producer index
% 	c							= consumer index
%
% calls:
%       none
%
% takes:
%       DIET                    DIET matrix (D_pc)                          (2D matrix: num_grps X num_grps)
%       TotalConsumptionVector  total consumption by each consumer c (q_c)  (horizontal vector: 1 X num_grps)
%
% returns:
%       PredationBudget         (A_cp) fraction of p production flowing to each consumer c (2D matrix: num_grps X num_grps; consumers X prodcers)
%
% revision date: 11-20-2019


% *************************************************************************
% STEP 1: measure matrix & vector sizes------------------------------------
fname_CalcPredationBudget	= mfilename; % save name of this f_CalcPredationBudget sub-function
display(['   Running: ' fname_CalcPredationBudget])

[rows_diet, clms_diet]	= size(DIET);
[rows_cons, clms_cons]	= size(TotalConsumptionVector);
num_grps                = rows_diet;

if clms_diet ~= clms_cons
    error('ERROR: consumption vector and DIET matrix are of mis-matched size')
end
% *************************************************************************



% *************************************************************************
% STEP 2: calculate total predation on each producer p---------------------
%         NOTE: predation_total = sigma(D_pc * q_c) for all c.
%         NOTE: We do NOT exclude any flow to nutrients, eggs, or detritus
%               when calculating the PredationBudget. BUT, we DO exclude 
%               any flow nutrients, eggs, or detritus in other functions
%               when we are specifically calculating predation terms in
%               ProductionBudget & ConsumptionBudget
%         NOTE: This is the denominator in Eq. 9 of Steele & Ruzicka (2011)
%         NOTE: We could just supply the function with the CONSUMPTION
%               matrix directly. We don't do that here simply in order to show
%               the calculations step-by-step as described by Eq. 9
%               in Steele & Ruzicka (2011)
%         NOTE: predation_total includes flow to eggs & detritus pools, but
%               this is later accounted for in the calculation of the
%               EnergyBudget by scaling consumer & fleet rows by the predation
%               terms of the ConsumptionBudget
CONSUMPTION     = DIET .* repmat(TotalConsumptionVector, [num_grps, 1]);	% (t/km2/y); (2D matrix: num_grps X num_grps)
predation_total = sum(CONSUMPTION, 2);                                      % total predation on each producer p; (t/km2/y); (vertical vector: num_grps X 1)
% *************************************************************************



% *************************************************************************
% STEP 3: calculate the bottom-up, PredationBudget-------------------------
%         NOTE: PredationBudget = A_cp; consumers c run down rows
PredationBudget                         = CONSUMPTION' ./ repmat(predation_total', [num_grps 1]); % (proportions); (2D matrix: num_grps X num_grps); NOTE transposes
PredationBudget(isnan(PredationBudget)) = 0;        % fix div/0 errors for groups with no predation, e.g. mammals, (PredationBudget = 0)
% *************************************************************************


% end m-file***************************************************************