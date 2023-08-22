function [PSUEDOARRAY, PSUEDOARRAY_dims] = f_unspoolMATRIX_04282020(MATRIX)
% linearize multidimenional matrices up to 4-D for use in C++
% by Jim Ruzicka
% takes:
%           MATRIX              up to 4 dimensions
% returns:
%           PSUEDOARRAY         1-dimension vector
%           PSUEDOARRAY_dims    dimensions of original MATRIX (horizontal vector: 1 X 4 [num_rows, num_clms, num_layers, num_fourthdim])
% How to index the vectorized matrix:
%       3D matrix indexing is: PSUEDOARRAY[(layer_index*num_rows*num_clms) + (row_index*num_clms) + clm_index])
% revision date: 4-28-2020

% STEP 1: get matrix dimensions-------------------------------------------
[num_rows, num_clms, num_layers, num_fourthdim] = size(MATRIX);

% STEP 2: initilaize output-----------------------------------------------
PSUEDOARRAY         = zeros([1, (num_rows * num_clms * num_layers * num_fourthdim)]); % initialize the PSUDEOARRAY
PSUEDOARRAY_dims	= [num_rows, num_clms, num_layers, num_fourthdim];

% STEP 3: build the PSUDEOARRAY-------------------------------------------
for fourthdim_index = 1:num_fourthdim
    for layer_index = 1:num_layers
        for clm_index = 1:num_clms
            for row_index = 1:num_rows
                
                PSUEDOARRAY(1, (((fourthdim_index-1) * num_layers * num_rows * num_clms) + ((layer_index-1) * num_rows * num_clms) + ((row_index-1) * num_clms) + clm_index))= ...
                    MATRIX(row_index, clm_index, layer_index, fourthdim_index); % note subtraction of 1 from each index because C++ addressing is 0-based, clm_index on left is really = (clm_index-1) + 1 (+1 here because we are still in matlab and first index = 1;
                
%                 % for debugging...
%                 PSUEDOARRAY(1, (((fourthdim_index-1) * num_layers * num_rows * num_clms) + ((layer_index-1) * num_rows * num_clms) + ((row_index-1) * num_clms) + clm_index))= ...
%                     (((fourthdim_index-1) * num_layers * num_rows *
%                     num_clms) + ((layer_index-1) * num_rows * num_clms) +
%                     ((row_index-1) * num_clms) + clm_index); % using
%                     index numbers for debugging
                
            end % end for row_index
        end % end for clm_index
    end % end for layer_index
end % end for fourthdim_index


% end function************************************************************