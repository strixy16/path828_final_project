% Adding select genes to Clinical spreadsheet to use in SPSS

load('../Data/selectCDEsfixed.mat');
load('spss_selectgenes.mat')

% Get just ICC patients
icc_cdes = selectCDEsfixed(selectCDEsfixed.icd_10=='c22.1',:);

days_to_death_or_followup = sum([icc_cdes.days_to_death icc_cdes.days_to_last_followup],2, 'omitnan');
censored = isnan(icc_cdes.days_to_death);
age = floor(icc_cdes.days_to_birth/-365);

icc_cdes = addvars(icc_cdes, days_to_death_or_followup, 'Before', 'days_to_birth'); 
icc_cdes = addvars(icc_cdes, censored, 'Before', 'days_to_birth');
icc_cdes = addvars(icc_cdes, age, 'Before', 'days_to_birth');

% icc_cdes.days_to_death_or_followup = days_to_death_or_follow_up;
% icc_cdes.censored = censored;

spss_fscnca2_mRNA_table = [icc_cdes spss_fscnca2_mRNA];
spss_mrmr2_mRNA_table = [icc_cdes spss_mrmr2_mRNA];
spss_litsel_mRNA_table = [icc_cdes spss_litsel_mRNA];

writetable(spss_fscnca2_mRNA_table, '../Data/spss_fscnca2_mRNA.xlsx')
writetable(spss_mrmr2_mRNA_table, '../Data/spss_mrmr2_mRNA.xlsx')
writetable(spss_litsel_mRNA_table, '../Data/spss_litsel_mRNA.xlsx')