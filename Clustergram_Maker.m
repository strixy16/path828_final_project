function Clustergram_Maker(data, col_labels, rowpdist, colpdist, linkage, title)
% Name: Clustergram_Maker
% Description: Function that makes clustergram based on input arguments
%
% Example input: Clustergram_Maker(data, col_labels, 'Euclidean',
%                                 'Spearman', 'average', 'My Clustergram');
%
%INPUT:  - data: Data to make clustergram with, 2D matrix
%        - col_labels: Column labels, cell_array, same length as data
%        - rowpdist: Distance metric for the data, string
%        - colpdist: Distance metric for the column labels, string
%        - linkage: Algorithm for computing the distance between clusters,
%                   string
%        - title: Title for plot, string
%
%OUTPUT: - None
%
%Environment: MATLAB R2020b
%
%Notes: Created for use in PATH828 Assignment 3
%
%Author: Katy Scott
%
%Last edited: 30 November 2020
%
%Revisions: 30 November 2020  -> rotated column labels (only works for small
%                            labels)
%Revisions: 5 November 2020   -> commented out rotating column labels
%TODO: - input checking
%      - set some defaults
%      - Make so tempTitle is the default, but takes title otherwise

%% INPUT CHECKING 
    % Can't use datain if not a matrix
    if ~ismatrix(data)
        disp("Data must be a 2D matrix");
        return
    % Can't use datain if values aren't double
    elseif ~isa(data, 'double')
        disp("Data must be a double matrix");
        return  
    end
    
%% CLUSTERGRAM GENERATION
cg = clustergram(data, ...
            'ColumnLabels', col_labels,...
            'RowPDist', rowpdist,...
            'ColumnPDist', colpdist,...
            'ColorMap', 'parula',...
            'Linkage', linkage);
%             'ColumnLabelsRotate', 0);

tempTitle =['Linkage: ' linkage ', Row: ' rowpdist ', Column: ' colpdist];
        
addTitle(cg, {title, tempTitle});

end