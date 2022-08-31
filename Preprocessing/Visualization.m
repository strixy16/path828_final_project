function Visualization(data, data_title, sample_labels, y_label, plots, log_flag) 
%Name: Visualization
%Description: Function to create plots for data visualization. See inputs
%             for plot options.
%
%INPUT:  - data: double matrix, data to be plotted, samples treated as
%                columns
%        - data_title: string, title to appear on plots 
%        - sample_labels: cell array, labels the x-axis in the scatter
%                         plots with each of the samples. Must be the same 
%                         size as the number of columns in data.
%        - y_label: y-axis label for boxplot
%        - plots: boolean array, indicates which plots to make
%               - index 1: Boxplot
%               - index 2: Total Counts vs. Samples scatter plot
%               - index 3: IQR vs. Samples scatter plot
%               - index 4: Correlation vs. Samples scatter plot, uses
%                          Spearman correlation
%
%               - Example: [0 1 1 0] makes total count and IQR plots
%        - log_flag: boolean value, indicates whether to perform log2
%                    transformation on data for the box plots
%
%OUTPUT: None
%
%Environment: MATLAB R2020b
%
%Author: Katy Scott
%
%Last edited: 30 November 2020
%Revisions:  30 Nov. 2020 -> added plot_hist section

    % Creating understandable boolean flags
    plot_box = plots(1);
    plot_scatter_total = plots(2);
    plot_scatter_iqr = plots(3);
    plot_scatter_corr = plots(4);
    plot_hist = plots(5);
    
    if plot_box
        figure
        if log_flag
            % Log transforming data to handle skew in data
            log_data = log2(replaceZeros(data, 'lowval'));
            boxplot(log_data);
            title(append(data_title, ', log2 transformed'));
        else
            boxplot(data);
            title(data_title);
        end
        xlabel('Samples');
        ylabel(y_label);
    end
    
    if plot_scatter_total
        % Sum expression values for each sample
        total_counts = sum(data,1);
        plotScatterForData(total_counts, data_title,'Total Counts', sample_labels);
    end
        
    if plot_scatter_iqr
        % Calculate IQR for each sample
        features_IQR = iqr(data,1);
        plotScatterForData(features_IQR, data_title,'IQR', sample_labels);
    end
            
    if plot_scatter_corr
        % Calculate average correlation between samples
        avg_correlations = SampleCorrelation(data, 'Spearman');
        plotScatterForData(avg_correlations, data_title,'Correlation', sample_labels);
    end
    
    if plot_hist
        % Normal distribution check plot
        figure
        histogram(data);
        title(data_title);
        xlabel('Expression Level');
        ylabel('Counts');
    end
end