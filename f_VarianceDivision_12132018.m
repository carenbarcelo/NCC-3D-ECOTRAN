function quotient_VAR = f_VarianceDivision_12132018(mean_A, VAR_A, mean_B, VAR_B)
% calculate the variance of two products
%
%           RULE: VAR     = STD^2 = (mean * CV)^2
%
%           RULE: CV      = STD / mean
%
%           RULE: adding & subtracting uncertainties: 
%                       mean(A+B) = mean(A) + mean(B) OR mean(A-B) = mean(A) - mean(B) 
%                       var(A+B)  = var(A) + var(B) (variances are additive whether adding or subtracting terms)
%
%           RULE: multiplying uncertainties of INDEPENDENT variables (i.e., we can ignore covariance term): 
%                 (from: https://chem.libretexts.org/Core/Analytical_Chemistry/Quantifying_Nature/Significant_Digits/Propagation_of_Error)
%                       mean(A*B) = mean(A) * mean(B)
%                       var(A*B)  = [mean(A*B)^2] * [(var(A)/mean(A)^2) + (var(B)/mean(B)^2)]    NOTE: Goodman (1960) derives the same equation but also adds the term [... + var(A)*var(B)]
%                                                                                                      (Goodman, L.A. 1960. On the exact variance of products. Journal of the American Statistical Association. December 1960, 708?713. DOI: 10.2307/2281592)
%           RULE: dividing uncertainties of INDEPENDENT variables (i.e., we ignore covariance term):
%                 (from: Seltman, H. Approximations for Mean and Variance of a Ratio. Carnegie Mellon University)
%                       mean(A/B) = mean(A) / mean(B) 
%                       var(A/B)  = [mean(A)^2/mean(B)^2] * [(var(A)/mean(A)^2) + (var(B)/mean(B)^2)]
%
%
% revision date: 12-13-2018


% *************************************************************************
% STEP 1: calculate product of two variances-------------------------------
quotient_VAR                        = ((mean_A.^2) ./ (mean_B.^2)) .* ((VAR_A ./ (mean_A.^2)) + (VAR_B ./ (mean_B.^2)));
quotient_VAR(isnan(quotient_VAR))	= 0;                % fix div/0 errors
quotient_VAR                        = abs(quotient_VAR);	% make absolute values, variance terms are all positive


% end m-file---------------------------------------------------------------