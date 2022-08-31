function marked_features = MarkLowCounts(datain, Q)
%Name: Mark Low Counts
%Description: Function that computes the threshold level based on Q.
%             For each miRNA, checks and marks if expression is above the 
%             threshold in at least one sample. Returns a boolean vector
%             indicating which features have no samples above the threshold 
%             calculated from the quantile Q.
%
%INPUT:  - input: matrix to filter, type double 
%        - Q: 0 <= Q <= 1, indicates quartile level to establish the
%         filtering threshold
%
%OUTPUT: - marked_features: boolean column vector, indicating which features 
%                           have no samples above the provided quantile 
%                           with a 1
%
%Environment: MATLAB R2020b
%
%Notes: This is a mini-assignment for PATH828 
%
%Author: Katy Scott
%
%Last edited: 30 October 2020

    % Input checking for datain 
    % Can't use datain if not a matrix
    if ~ismatrix(datain)
        disp("Input must be a 2D matrix");
        marked_features = -1;
        return
        
    % Can't use datain if values aren't double
    elseif ~isa(datain, 'double')
        disp("Input must be a double matrix");
        marked_features = -1;
        return
        
    % Can't use cell matrix, but can convert to double and continue on
    elseif iscell(datain)
        disp("Input is cell matrix. Converting to double matrix");
        datain = cell2mat(datain); 
    end 
    
    % Input checking for Q, since:
    % quantile function takes in a numeric argument for Q
    if ~isnumeric(Q)
        disp("Q must be a numeric value");
        marked_features = -1;
        return
        
    % quantile function needs a value between 0 and 1
    elseif Q < 0 || Q > 1
        disp("Q value must be between 0 and 1");
        marked_features = -1;
        return
    end

    % Need data as a single row vector so quantile is calculated over all values
    input_as_vector = reshape(datain, [1 size(datain,1)*size(datain,2)]);
    
    % Making a threshold value to compare to samples
    threshold_level = quantile(input_as_vector, Q);
    
    samples_above_threshold = (datain > threshold_level);
    
    % Checking which features (rows) have no samples above the threshold
    marked_features = all(samples_above_threshold == 0, 2);
end