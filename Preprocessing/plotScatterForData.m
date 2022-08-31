function plotScatterForData(data, plot_title, ylbl, datalbl)
%Name: plotScatterForData
%Description: Creating scatterplot for input data, with outlier
%             boundaries and outlier points identified. 
%
%INPUT:  - data: vector to make scatterplot of, double (i.e. mean
%                correlations)
%        - plot_title: title for scatter plot
%        - ylbl: y-axis label
%        - datalbl: x-axis labels
%
%OUTPUT: None
%
%Environment: MATLAB R2020b
%
%Notes: 
%
%Author: Kathrin Tyryshkin
%
%Revisions: Oct. 2020 -> changed alpha from 1.5 to using
%                        computeAlphaOutliers
%                     -> added comments

alpha = computeAlphaOutliers(size(data,2));
% Calculate outlier boundaries
qupper = quantile(data, 0.75)+alpha*iqr(data);
qlower = quantile(data, 0.25)-alpha*iqr(data);
% Marking outliers
flag1 = data<qlower;
flag2 = data>qupper;

x = 1:length(data);
figure, hold on
% Plot data
plot(x, data, 'b.', 'markersize', 15);
% Plotting the outlier boundaries
plot([0 length(data)+1], [qupper qupper], '-r', 'linewidth', 2);
plot([0 length(data)+1], [qlower qlower], '-r', 'linewidth', 2); hold on;
% Plot outliers
plot(x(flag1), data(flag1), '*m', 'markersize', 15);
plot(x(flag2), data(flag2), '*m', 'markersize', 15);
% Labelling outliers with gene name
text(x(flag1)+0.7, data(flag1), datalbl(flag1));
text(x(flag2)+0.7, data(flag2), datalbl(flag2));

% Centering the data in the plot
axis([0 length(data)+1 nanmin(data)-nanstd(data) nanmax(data)+nanstd(data)]);
% Plot aesthetics
title(plot_title, 'FontSize',16, 'FontName', 'Helvetica');
set(gca, 'XTick', 1:length(datalbl));
set(gca, 'XTickLabel', datalbl);
ax = gca;
ax.XTickLabelRotation=45;
ylabel(ylbl, 'FontSize',14, 'FontName', 'Helvetica');
hold off;
