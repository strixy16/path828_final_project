function avg_correlations = SampleCorrelation(datain, correlation)
%Name: Sample Correlation
%Description: For each sample (column) in an input matrix, compute an average
%             correlation to all other samples (columns). Correlation type 
%             is specified as second input.
%             Returns a vector of the average correlation for each column.
%
%INPUT:  - datain: matrix of samples to compute correlations from, samples
%                  are columns, type double
%        - correlation: string or char array indicating correlation type 
%                       to use
%
%OUTPUT: - avg_correlations: vector of the average correlations for each
%                            column, has same length as number of columns 
%                            in datain
%
%Environment: MATLAB R2020b
%
%Notes: This is a mini-assignment for PATH828 
%
%Author: Katy Scott
%
%Revisions: 2020-10-05 -> changed output variable name from result to
%                         avg_correlations
%                      -> changed nanmean to mean with omitnan flag
%                      -> changed dimension in mean from 2 to 1, row mean
%                         to column mean
%           2020-10-07 -> Added input checking 
%                      -> Made char string an acceptable input for
%                         correlation

    % Check that input data is correct type 
    if ~ismatrix(datain) 
        disp("First argument must be a matrix");
        % Return invalid result if input is incorrect
        avg_correlations = -1;
        return
    elseif ~ischar(correlation) && ~isstring(correlation)
        disp("Second argument must be a string or char array");
        % Return invalid result if input is incorrect
        avg_correlations = -1;
        return
    end
    
    % Checking which correlation input was passed
    switch correlation
        case {"Spearman", 'Spearman'}
            % Spearman correlation was chosen, perform correlation between
            % each sample(column) of input matrix
            rho = corr(datain, 'Type', 'Spearman');
            
        case {"Pearson", 'Pearson'}
            % Pearson correlation was chosen, perform correlation between 
            % each sample(column) of input matrix
            rho = corr(datain, 'Type', 'Pearson');
            
        otherwise
            % Incorrect correlation input was received.
            disp("Invalid correlation input.");
            % Return invalid result.
            avg_correlations = -1;
            return
    end
    
    % Remove correlation between column and itself
    
    % Make diagonal matrix of NaN for correlation matrix 
    nans = diag(NaN(1,size(rho,1)));
    
    % Set diagonal to 0 and add NaN so it can be ignored in mean
    % calculation
    rho = rho - diag(diag(rho)) + nans;
    
    % omitnan flag ignores all NaN values in mean calculation
    avg_correlations = mean(rho,1,'omitnan');
end