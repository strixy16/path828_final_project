function [cell_data, mat_data] = RemoveSamples(cell_data,samples_to_remove, row_start_ind, col_start_ind)
%Name: RemoveSamples
%Description: Function to remove indicated samples (columns) from a
%             provided cell array dataset. Returns the dataset as a cell 
%             array and converted to a matrix.
%
%INPUT:  - cell_data: data to remove samples from, cell array, assumes 
%                     first column is feature labels
%        - samples_to_remove: index of samples to be removed from cell_data,
%                             double matrix
%        - row_start_ind: row index to start from for cell2mat conversion
%        - col_start_ind: column index to start from for cell2mat
%        conversion
%
%OUTPUT: - cell_data: input cell_data with indicated samples removed, cell
%                     array 
%        - mat_data: cell_data converted to a double matrix
%
%Environment: MATLAB R2020b
%
%Notes:
%
%Author: Katy Scott
%
%Last edited: 2 December 2020
%Revison: 2 Dec 2020 -> replaced hardcoded indices in mat_data conversion
%                       with row_start_ind and col_start_ind

    % Assuming input has feature IDs as first column, need to add 1 to the
    % indices of samples_to_remove to correct for this
    idx_samples_to_remove = samples_to_remove + 1;

    cell_data(:,idx_samples_to_remove) = [];

    mat_data = cell2mat(cell_data(row_start_ind:end, col_start_ind:end));
end

