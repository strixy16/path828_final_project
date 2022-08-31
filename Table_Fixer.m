% Fixing tables script
% Code used to put the CDEs data from firebrowse into the desired setup in
% an excel spreadsheet
load('../Data/AllCDEs.mat');
% Features are in a column instead of row, need to switch this
temp = table2cell(AllCDEs_total)';
fixed_AllCDEs = cell2table(temp(2:end,:), 'VariableNames', temp(1,:));
writetable(fixed_AllCDEs, '../Data/AllCDEs_fixed.xlsx');

load('../Data/selectCDEs.mat');
temp = table2cell(selectCDEs)';
fixed_selectCDEs = cell2table(temp(2:end,:), 'VariableNames', temp(1,:));
writetable(fixed_selectCDEs, '../Data/selectCDEs_fixed.xlsx');
