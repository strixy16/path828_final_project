function normalized_data = SampleNormalizationRF(datain)
%Name: Sample Normalization RF
%Description: Computes sample normalization on a matrix. 
%             Computes the sum of values for each column in the matrix and 
%             divides each value in a column of the matrix by its sum.
%
%INPUT:  - datain: matrix to perform sample normalization on, double
%
%OUTPUT: - normalized_data: datain after sample normalization has been
%                            performed
%
%Environment: MATLAB R2020b
%
%Notes: This is a mini-assignment for PATH828 
%
%Author: Katy Scott
%
%Revisions: 2020-10-05 -> changed output variable from result to normalized_data
%           2020-10-07 -> added input checking 
    % Check that input data is correct type (a matrix)
    if ~ismatrix(datain)
        disp("Input must be a matrix");
        % Return invalid result if input is incorrect
        normalized_data = -1;
        return
    % If input is a cell matrix, convert to a double matrix and continue
    elseif iscell(datain)
        disp("Input is cell matrix. Converting to double matrix");
        datain = cell2mat(datain);
    end
    
    % Sum each column of the datain matrix to use for normalization
    sum_cols = sum(datain);

    % Perform normalization by dividing each value in column by the sum of that column
    normalized_data = datain ./ sum_cols;
end