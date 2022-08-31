function [outliers, high_outliers, low_outliers] = scatterplot_mark_outliers(data, options)
% SCATTERPLOT_MARK_OUTLIERS(data) makes scatter plot of data, plots horizontal line at the upper and lower
% outlier boundries defined as :
%    25th/75th percentile +/- Alpha x IQR. 
%    see <a href="matlab:web('https://doi.org/10.1080/01621459.1987.10478551')">Hoaglin and Iglewicz (1986)</a>
% Marks outliers and returns logical vectors with outlier cols as 'true'
%
% Syntax:
%   [outliers, high_outliers, low_outliers] = SCATTERPLOT_MARK_OUTLIERS(data[, options])
%
% Inputs:
%   data: 1xN numerical matrix
% Options
%   Alpha: float, defaults to 1.5
%   Plot
%   ColumnLabels: cell array of strings; must be same dimenstions as data
%   PlotTitle: string
%   ylabel: string
%   
% output:
%   [outliers, high_outliers, low_outliers]: all logical arrays indicating
%       all outliers, outliers above the upper threshold, and outliers below 
%       the lower threshold, respectively. 

% Author: Kathrin Tyryshkin
% Last edited: 27 September 2020 by Justin Wong
%              1  December 2020 by Katy Scott - updated default alpha to
%              use computeAlphaOutliers function and added ShowPlot variable to
%              options

arguments
    data (1,:) {mustBeNumeric}
    options.Alpha {mustBeNumeric} = computeAlphaOutliers(size(data,2))
    options.ShowPlot = false
    options.ColumnLabels = {}
    options.PlotTitle string = ''
    options.ylabel string = ''
end


% Identify outler boundaries
qupper = quantile(data, 0.75) + options.Alpha * iqr(data);
qlower = quantile(data, 0.25) - options.Alpha * iqr(data);

% Identify outliers
high_outliers = data > qupper;
low_outliers = data < qlower;
outliers = high_outliers | low_outliers; 

if options.ShowPlot
    x = 1:length(data);

    figure, hold on

    % Plot data
    plot(x, data, 'b.', 'markersize', 15);
    % Plot outlier boundaries
    plot([0 length(data)+1], [qupper qupper], '-r', 'linewidth', 2);
    plot([0 length(data)+1], [qlower qlower], '-r', 'linewidth', 2);
    % Mark and label outliers
    plot(x(outliers), data(outliers), '*m', 'markersize', 15); 
    if ~isempty(options.ColumnLabels)
        text(x(outliers)+0.3, data(outliers), options.ColumnLabels(outliers));
    end


    % Plot aesthetics
    if ~isempty(options.PlotTitle)
        title(options.PlotTitle, 'FontSize',16, 'FontName', 'Helvetica');
    end
    if ~isempty(options.ColumnLabels)
        set(gca, 'XTick', 1:length(options.ColumnLabels));
        set(gca, 'XTickLabel', options.ColumnLabels);
        ax = gca;
        ax.XTickLabelRotation = 45;
    end
    if ~isempty(options.ylabel)
        ylabel(options.ylabel, 'FontSize',14, 'FontName', 'Helvetica');
    end

    hold off;
end
end 