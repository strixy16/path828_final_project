% Script for PATH project

% Switch for preprocessing visualization
SHOW_PREPROCESS_VISUALS = false;
% Switch for quality control visualization
SHOW_QC_VISUALS = false;
% Switch for filtering visualization
SHOW_FILTER_VISUALS = false;
% Threshold for filtering
FILTER_THRESHOLD = 0.99;
% Years of survival
SURVIVAL_YEARS = 2;


addpath('Preprocessing');
% close all force;

%% Data Loading
load('../Data/AllCDEsfixed.mat');
% load('../Data/mRNA_raw_counts.mat');
load('../Data/mRNA_RSEM_counts.mat');

% Get just ICC patients
icc_cdes = AllCDEsfixed(AllCDEsfixed.icd_10=='c22.1',:);

% Convert patient IDs to strings, remove dashes, and make all uppercase
pat_ids = string(icc_cdes.bcr_patient_barcode);
pat_ids = upper(erase(pat_ids, "-"));

% Select mRNA counts for ICC patients ids only
icc_mRNA_RSEM_counts = mRNA_RSEM_counts(:,contains(mRNA_RSEM_counts.Properties.VariableNames, pat_ids));
icc_mRNA_RSEM_counts = [mRNA_RSEM_counts(:,1) icc_mRNA_RSEM_counts];

cell_icc_mRNA_RSEM = table2cell(icc_mRNA_RSEM_counts);
mat_icc_mRNA_RSEM = cell2mat(cell_icc_mRNA_RSEM(:,2:end));

if SHOW_PREPROCESS_VISUALS
    Visualization(mat_icc_mRNA_RSEM, 'ICC RSEM mRNA Data', pat_ids, 'Expression Levels', [1 1 0 0 1], 0);
end
%% Normalization and Transformation
% norm_log_icc_mRNA = NormalizeAndTransform(mat_icc_mRNA);
log_icc_mRNA_RSEM = log2(replaceZeros(mat_icc_mRNA_RSEM, 'lowval'));

if SHOW_PREPROCESS_VISUALS
    Visualization(log_icc_mRNA_RSEM, 'ICC RSEM mRNA Data, log2 transformed', pat_ids, 'Expression Levels', [1 0 1 1 1], 0);
end
%% QC for batch effects and outliers
% Total counts outliers
% From Katherine: don't have to remove these, can just mark them
tc_icc_mRNA_RSEM = sum(mat_icc_mRNA_RSEM,1);

[tc_outliers, tc_high_outliers, tc_low_outliers] = scatterplot_mark_outliers(tc_icc_mRNA_RSEM, 'ShowPlot', false, 'ColumnLabels', pat_ids, 'PlotTitle',...
                                                                                        'ICC RSEM mRNA Data, log transformed');

% Marking total counts outliers
idx_marked_outliers = find(tc_outliers==1);

% Removing total counts outliers
% idx_outliers_to_remove = find(high_outliers==1);
% [cell_icc_mRNA_RSEM, mat_remout_icc_mRNA_RSEM] = RemoveSamples(cell_icc_mRNA_RSEM,idx_outliers_to_remove, 1, 2);
% % Removing from patient ids
% pat_ids(idx_outliers_to_remove) = [];
% % Removing from clinical data
% icc_cdes(idx_outliers_to_remove, :) =[];
% 
% log_icc_mRNA_RSEM = log2(replaceZeros(mat_remout_icc_mRNA_RSEM, 'lowval'));
% 
% if SHOW_QC_VISUALS
%     Visualization(mat_remout_icc_mRNA_RSEM, 'ICC RSEM mRNA Data, outliers removed', pat_ids, 'Expression Levels', [0 1 0 0 0], 0);
%     Visualization(log_icc_mRNA_RSEM, 'ICC RSEM mRNA Data, outliers removed, log transformed', pat_ids, 'Expression Levels', [1 0 1 1 1], 0);
% end
%% Filtering

features_to_remove = MarkLowCounts(log_icc_mRNA_RSEM, FILTER_THRESHOLD);
idx_features_to_remove = features_to_remove == 1;

cell_filtered_icc_mRNA_RSEM = cell_icc_mRNA_RSEM;
filtered_log_icc_mRNA_RSEM = log_icc_mRNA_RSEM;

cell_filtered_icc_mRNA_RSEM(idx_features_to_remove,:) = []; 
filtered_log_icc_mRNA_RSEM(idx_features_to_remove,:) = [];

if SHOW_FILTER_VISUALS
    Visualization(filtered_log_icc_mRNA_RSEM, ...
                 {'ICC RSEM mRNA Data, log2 transformed', ['filter threshold = ' num2str(FILTER_THRESHOLD)]},...
                  pat_ids, 'Expression Levels', [1 0 0 0 1], 0);
    Visualization(cell2mat(cell_filtered_icc_mRNA_RSEM(:,2:end)), 'ICC RSEM mRNA Data, filtered',...
                  pat_ids, [], [0 1 0 0 0], 0);
end

%% UNSUPERVISED LEARNING

% Setting up labels

% Getting years to death/followup instead of day
years_to_death = icc_cdes.days_to_death / 365;
survive_past = double(years_to_death > SURVIVAL_YEARS);

years_to_followup = icc_cdes.days_to_last_followup /365;
followup_past = double(years_to_followup > SURVIVAL_YEARS) * 3;
followup_before = double(years_to_followup <= SURVIVAL_YEARS) * 2;

% Key: 0 = death before SURVIVAL YEARS
%      1 = death after SURVIVAL_YEARS
%      2 = followup before SURVIVAL_YEARS
%      3 = follwoup after SURVIVAL_YEARS
surv_follow_past = survive_past + followup_past + followup_before;

%% Hierarchical Clustering
% median centering the data
medcent_log_icc_mRNA = filtered_log_icc_mRNA_RSEM - median(median(filtered_log_icc_mRNA_RSEM));

% cluster_title = ['ICC RSEM mRNA ' num2str(FILTER_THRESHOLD) ' filtered vs. ' num2str(SURVIVAL_YEARS) ' year survival/followup'];
% Clustergram_Maker(medcent_log_icc_mRNA, surv_follow_past, 'Chebychev','Spearman', 'average', cluster_title) 
% Clustergram_Maker(medcent_log_icc_mRNA, surv_follow_past, 'Euclidean', 'Spearman', 'average', cluster_title)
% 
% outlier_cluster_title = ['ICC RSEM mRNA ' num2str(FILTER_THRESHOLD) ' filtered vs. total counts outlier'];
% % Clustering with outliers as labels
% Clustergram_Maker(medcent_log_icc_mRNA, double(tc_outliers), 'Chebychev','Spearman', 'average', outlier_cluster_title) 
% Clustergram_Maker(medcent_log_icc_mRNA, double(tc_outliers), 'Euclidean', 'Spearman', 'average', outlier_cluster_title)


%% Kmeans Clustering

% TODO: make this a function
% Running k-means with k-values 2-10
K_MIN = 2;
K_MAX = 3;
% Need to try different distance metrics

medcent_log_icc_mRNA = medcent_log_icc_mRNA';

cluster_indices = zeros(size(medcent_log_icc_mRNA,1),(K_MAX-K_MIN + 1));
avg_dist_to_centroids = zeros(1, (K_MAX-K_MIN + 1));
for k=K_MIN:K_MAX
    [idx, C, sumd] = kmeans(medcent_log_icc_mRNA,k);
    
    avg_dist_to_centroids(:,k-K_MIN+1) = mean(sumd);
    figure
    silhouette(medcent_log_icc_mRNA, idx);
    title(['Silhouette for ' num2str(k) ' clusters']);
    
%     figure
%     gscatter(X(:,1), X(:,2), idx)
    
end

figure
plot(K_MIN:K_MAX, avg_dist_to_centroids)
xlim([1 K_MAX]);
xlabel('Number of clusters (k)');
ylabel('Avg. distance to centroids');




%% SUPERVISED LEARNING
%% Feature Selection
%%% Selecting from literature %%%
% Jarnagin: IDH1, IDH2, FGFR2, TP53, KRAS, CDKN2A, CDKN2B
% Andersen additions: BRAF, EGFR, PIK3CA

lit_genes = {'IDH1'; 'IDH2'; 'FGFR2'; 'TP53'; 'KRAS'; 'CDKN2A'; 'CDKN2B'; ...
             'BRAF'; 'EGFR'; 'PIK3CA'};

cell_lit_selectgenes = cell_icc_mRNA_RSEM(contains(cell_icc_mRNA_RSEM(:,1), lit_genes), :);
mat_lit_selectgenes = cell2mat(cell_lit_selectgenes(:,2:end));
log_mat_lit_selectgenes = log2(replaceZeros(mat_lit_selectgenes, 'lowval'));
medcent_log_lit_selectgenes = log_mat_lit_selectgenes - median(median(log_mat_lit_selectgenes));


% Setting up 2 label sets for feature selection functions
% using 4 classes, death before and after YEARS and follow up before and
% after YEARS
labels_4 = surv_follow_past;

% using 2 classes, combining death and followup, still before and after YEARS
labels_2 = survive_past + double(years_to_followup > SURVIVAL_YEARS);

% Using sequentialfs

% Setting up display output for sequentialfs
% Will show output for each iteration
% opts = statset('Display','iter');
% 
% num_iter = 25;
% num_features = size(medcent_log_icc_mRNA, 1);
% fs_totals_2 = zeros(1, num_features);
% cell_fs_2 = {num_iter, 2}; % store a cell containing the features each iteration
% 
% fs_totals_4 = zeros(1, num_features);
% cell_fs_4 = {num_iter, 2}; 
% for i = 1:num_iter
%     % Creating k partitions with approximately the same number of observations
%     % and class proportions 
%     c_2 = cvpartition(labels_2,'k',5);
%     c_4 = cvpartition(labels_4,'k',5);
%     
%     opts = statset('Display','iter');
%     % using neighbourhood component analysis for classification
%     fun = @(XT,yT,Xt,yt)loss(fscnca(XT,yT),Xt,yt, 'LossFun', 'classiferror');
%     
%     [fs_2, history_2] = sequentialfs(fun,medcent_log_icc_mRNA',labels_2,'cv',c_2,'options',opts);
%     cell_fs_2{i,1} = fs_2;
%     cell_fs_2{i,2} = cell_filtered_icc_mRNA_RSEM(fs_2, 1);
%     fs_totals_2 = fs_totals_2 + fs_2;
%     
%     [fs_4, history_4] = sequentialfs(fun,medcent_log_icc_mRNA',labels_4,'cv',c_4,'options',opts);
%     cell_fs_4{i,1} = fs_4;
%     cell_fs_4{i,2} = cell_filtered_icc_mRNA_RSEM(fs_4, 1);
%     fs_totals_4 = fs_totals_4 + fs_4;
% end

load('fscnca_run_dec10_25iter.mat')

% Thresholds here chosen so that there are more than 2 genes selected 
fs_final_2 = fs_totals_2 >= 3;
fs_final_4 = fs_totals_4 >= 3;

cell_fscnca_selectgenes_class2 = cell_filtered_icc_mRNA_RSEM(fs_final_2, :);
mat_fscnca_selectgenes_class2 = cell2mat(cell_fscnca_selectgenes_class2(:,2:end));

cell_fscnca_selectgenes_class4 = cell_filtered_icc_mRNA_RSEM(fs_final_4, :);
mat_fscnca_selectgenes_class4 = cell2mat(cell_fscnca_selectgenes_class4(:,2:end));

%%% Using fscmrmr %%%
tbl_medcent_log_icc_mRNA = array2table(medcent_log_icc_mRNA');

% Using 4 classes
[idx4, scores4] = fscmrmr(tbl_medcent_log_icc_mRNA, labels_4);
figure
bar(scores4(idx4(1:25)))
xlabel('Predictor rank')
ylabel('Predictor importance score')
title("MRMR Feature Selection for 4 classes")
% Selecting out all genes with any predictive power to compare
idx4_mrmr_selectgenes = idx4(1:10);
cell_mrmr_selectgenes_class4 = cell_filtered_icc_mRNA_RSEM(idx4_mrmr_selectgenes, :);
mat_mrmr_selectgenes_class4 = cell2mat(cell_mrmr_selectgenes_class4(:,2:end));

% Using 2 classes
[idx2, scores2] = fscmrmr(tbl_medcent_log_icc_mRNA, labels_2);
figure
bar(scores2(idx2(1:25)))
xlabel('Predictor rank')
ylabel('Predictor importance score')
title("MRMR Feature Selection for 2 classes")
% Selecting out all genes with any predictive power to compare
idx2_mrmr_selectgenes = idx2(1:6);
cell_mrmr_selectgenes_class2 = cell_filtered_icc_mRNA_RSEM(idx2_mrmr_selectgenes, :);
mat_mrmr_selectgenes_class2 = cell2mat(cell_mrmr_selectgenes_class2(:,2:end));

% Clustergrams with selected features     
% litselect_title_4 = 'ICC RSEM mRNA - literature selected features - 4 classes';
% Clustergram_Maker(medcent_log_lit_selectgenes, labels_4, 'Euclidean', 'Spearman', 'average', litselect_title_4)
% 
% litselect_title_2 = 'ICC RSEM mRNA - literature selected features - 2 classes';
% Clustergram_Maker(medcent_log_lit_selectgenes, labels_2, 'Euclidean', 'Spearman', 'average', litselect_title_2)
% 
% fscncaselect_title_4 = 'ICC RSEM mRNA - fscnca selected features - 4 classes';
% Clustergram_Maker(mat_fscnca_selectgenes_class4, labels_4, 'Euclidean', 'Spearman', 'average', fscncaselect_title_4)
% 
% fscncaselect_title_2 = 'ICC RSEM mRNA - fscnca selected features - 2 classes';
% Clustergram_Maker(mat_fscnca_selectgenes_class2, labels_2, 'Euclidean', 'Spearman', 'average', fscncaselect_title_2)
% 
% mrmrselect_title_4 = 'ICC RSEM mRNA - MRMR selected features - 4 classes';
% Clustergram_Maker(mat_mrmr_selectgenes_class4, labels_4, 'Euclidean', 'Spearman', 'average', mrmrselect_title_4)
% 
% mrmrselect_title_2 = 'ICC RSEM mRNA - MRMR selected features - 2 classes';
% Clustergram_Maker(mat_mrmr_selectgenes_class2, labels_2, 'Euclidean', 'Spearman', 'average', mrmrselect_title_2)

%% Classification Learner App
% Setting up data for use in the classification learner app
% Literature selected genes
cl_litsel_ICC_mRNA_2 = [medcent_log_lit_selectgenes; labels_2']';
cl_litsel_ICC_mRNA_4 = [medcent_log_lit_selectgenes; labels_4']';

% fscnca selected genes
cl_fscnca_ICC_mRNA_2 = [mat_fscnca_selectgenes_class2; labels_2']';
cl_fscnca_ICC_mRNA_4 = [mat_fscnca_selectgenes_class4; labels_4']';

% mrmr selected genes
cl_mrmr_ICC_mRNA_2 = [mat_mrmr_selectgenes_class2; labels_2']';
cl_mrmr_ICC_mRNA_4 = [mat_mrmr_selectgenes_class4; labels_4']';

%% SPSS data setup
pat_id_col = [0; pat_ids];

spss_fscnca2_mRNA = cell_fscnca_selectgenes_class2';
spss_fscnca2_mRNA = cell2table(spss_fscnca2_mRNA(2:end,:), 'VariableNames', spss_fscnca2_mRNA(1,:));

spss_mrmr2_mRNA = cell_mrmr_selectgenes_class2';
spss_mrmr2_mRNA = cell2table(spss_mrmr2_mRNA(2:end,:), 'VariableNames', spss_mrmr2_mRNA(1,:));

spss_litsel_mRNA = cell_lit_selectgenes';
spss_litsel_mRNA = cell2table(spss_litsel_mRNA(2:end,:), 'VariableNames', spss_litsel_mRNA(1,:));

save spss_selectgenes.mat spss_fscnca2_mRNA spss_mrmr2_mRNA spss_litsel_mRNA