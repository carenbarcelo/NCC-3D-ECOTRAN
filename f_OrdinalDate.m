function OrdinalDate = f_OrdinalDate(DatesOfInterest)
% calculate ordinal date of DatesOfInterest vector
% ordinal date is day of year with January 1 any year = 1
% by Jim Ruzicka
% calls:
%       none
% revision date: 2-8-2011

% *************************************************************************
% STEP 1: extract year info------------------------------------------------
temp_date_vector      = datevec(DatesOfInterest);
temp_year             = temp_date_vector(:, 1);
temp_month            = temp_date_vector(:, 2);
temp_day              = temp_date_vector(:, 3);





% *************************************************************************
% STEP 2: create vector of January 1sts------------------------------------
JanFirst_day_vector   = ones(length(temp_year), 1);
JanFirst_month_vector = ones(length(temp_year), 1);
JanFirst_vector       = datenum(temp_year, JanFirst_month_vector, JanFirst_day_vector);





% *************************************************************************
% STEP 3: calculate ordinal date-------------------------------------------
OrdinalDate           = datenum(DatesOfInterest) - JanFirst_vector + 1;


% end m-file***************************************************************