% script for determining how omics help discriminate cancerous/healthy samples according to the clinical data
clear, clc
results_path = "path";  % where to save results

addpath(genpath('path\SVM results'));

%% non-permuted labels
% load fluxes only
load('SVMl2_results_FVA.mat');
FVA_results = fullConversion(SVM_results);
% load genes only
load('SVMl2_results_genes.mat');
genes_results = fullConversion(SVM_results);
% load clinical data only
load('Int2SVMl2_results_clinical_only');
clinical = fullConversion(SVM_results);
% load fluxes with clinical data
load('Int3SVMl2_results_FVA');
FVA_clinical_results = fullConversion(SVM_results);
% load genes with clinical data
load('Int3SVMl2_results_genes');
genes_clinical_results = fullConversion(SVM_results);
% load fluxes and genes only
load('Int2SVMl2_results');
FVA_genes_results = fullConversion(SVM_results);
% load all data
load('Int4SVMl2_results');
all_results = fullConversion(SVM_results);
clear SVM_results;
load('path\combined_dataset.mat');
clear geneData_combat_p;
% select clinical data columns
tumour_data = sampleData(strcmp(sampleData.Status, 'Tumor'), {'Age', 'Sex', 'Status'});
healthy_data = sampleData(~strcmp(sampleData.Status, 'Tumor'), {'Age', 'Sex', 'Status'});

% data is the same for every run
healthy_test = FVA_results.testSamples(1:5, :);
tumour_test = FVA_results.testSamples(6:end, :);
healthy_test = reshape(healthy_test, [1, numel(healthy_test)]);
tumour_test = reshape(tumour_test, [1, numel(tumour_test)]); 

healthy_true = FVA_results.true(1:5, :);
tumour_true = FVA_results.true(6:end, :);
healthy_true = reshape(healthy_true, [1, numel(healthy_true)]);
tumour_true = reshape(tumour_true, [1, numel(tumour_true)]);

tumour_data = tumour_data(tumour_test, :);
healthy_data = healthy_data(healthy_test, :);
data = vertcat(healthy_data, tumour_data);
data.true = transpose(horzcat(healthy_true, tumour_true));

omic_combinations = {FVA_results, genes_results, clinical, FVA_clinical_results, genes_clinical_results, FVA_genes_results, all_results};
combination_names = {'FVA', 'genes', 'clinical', 'FVA_clinical', 'genes_clinical', 'FVA_genes', 'all'};

%compute and save comparisons for age
[~, ~, ~, ~] = stratifiedAnalyses(data, omic_combinations, ...
    combination_names, {'Age'}, strcat(results_path, "\age.csv"));
% compute and save comparisons for gender
[~, ~, ~, ~] = stratifiedAnalyses(data, omic_combinations, ...
    combination_names, {'Sex'}, strcat(results_path, "\sex.csv"));
% compute and save comparisons for health status
[~, ~, ~, ~] = stratifiedAnalyses(data, omic_combinations, ...
    combination_names, {'Status'}, strcat(results_path, "\status.csv"));
% compute and save comparisons for age and gender
[~, ~, ~, ~] = stratifiedAnalyses(data, omic_combinations, ...
    combination_names, {'Age', 'Sex'}, strcat(results_path, "\age_sex.csv"));
% compute and save comparisons for age and health status
[~, ~, ~, ~] = stratifiedAnalyses(data, omic_combinations, ...
    combination_names, {'Age', 'Status'}, strcat(results_path, "\age_status.csv"));
% compute and save comparisons for gender and health status
[~, ~, ~, ~] = stratifiedAnalyses(data, omic_combinations, ...
    combination_names, {'Sex', 'Status'}, strcat(results_path, "\sex_status.csv"));
% compute and save comparisons for all clinical data
[~, ~, ~, ~]= stratifiedAnalyses(data, omic_combinations, ...
    combination_names, {'Age', 'Sex', 'Status'}, strcat(results_path, "\all.csv"));

%% permuted labels
% load fluxes only
load('SVMl2_results_FVA_label_permutation.mat');
FVA_results_label_permutation = fullConversion(SVM_results);
% load genes only
load('SVMl2_results_genes_label_permutation.mat');
genes_results_label_permutation = fullConversion(SVM_results);
% load clinical data only
load('Int2SVMl2_results_clinical_only_label_permutation');
clinical_label_permutation = fullConversion(SVM_results);
% load fluxes with clinical data
load('Int3SVMl2_results_FVA_label_permutation');
FVA_clinical_results_label_permutation = fullConversion(SVM_results);
% load genes with clinical data
load('Int3SVMl2_results_genes_label_permutation');
genes_clinical_results_label_permutation = fullConversion(SVM_results);
% load fluxes and genes only
load('Int2SVMl2_results_label_permutation');
FVA_genes_results_label_permutation = fullConversion(SVM_results);
% load all data
load('Int4SVMl2_results_label_permutation');
all_results_label_permutation = fullConversion(SVM_results);
clear SVM_results;

omic_combinations_label_permutation = {FVA_results_label_permutation, genes_results_label_permutation, ... 
    clinical_label_permutation, FVA_clinical_results_label_permutation, ... 
    genes_clinical_results_label_permutation, FVA_genes_results_label_permutation, all_results_label_permutation};
combination_names_label_permutation = {'FVA_label_permutation', 'genes_label_permutation', ... 
    'clinical_label_permutation', 'FVA_clinical_label_permutation', 'genes_clinical_label_permutation', ... 
    'FVA_genes_label_permutation', 'all_label_permutation'};

%compute and save comparisons for age
[~, ~, ~, ~] = stratifiedAnalyses(data, omic_combinations_label_permutation, ...
    combination_names_label_permutation, {'Age'}, strcat(results_path, "\age_label_permutation.csv"));
% compute and save comparisons for gender
[~, ~, ~, ~] = stratifiedAnalyses(data, omic_combinations_label_permutation, ...
    combination_names_label_permutation, {'Sex'}, strcat(results_path, "\sex_label_permutation.csv"));
% compute and save comparisons for health status
[~, ~, ~, ~] = stratifiedAnalyses(data, omic_combinations_label_permutation, ...
    combination_names_label_permutation, {'Status'}, strcat(results_path, "\status_label_permutation.csv"));
% compute and save comparisons for age and gender
[~, ~, ~, ~] = stratifiedAnalyses(data, omic_combinations_label_permutation, ...
    combination_names_label_permutation, {'Age', 'Sex'}, strcat(results_path, "\age_sex_label_permutation.csv"));
% compute and save comparisons for age and health status
[~, ~, ~, ~] = stratifiedAnalyses(data, omic_combinations_label_permutation, ...
    combination_names_label_permutation, {'Age', 'Status'}, strcat(results_path, "\age_status_label_permutation.csv"));
% compute and save comparisons for gender and health status
[~, ~, ~, ~] = stratifiedAnalyses(data, omic_combinations_label_permutation, ...
    combination_names_label_permutation, {'Sex', 'Status'}, strcat(results_path, "\sex_status_label_permutation.csv"));
% compute and save comparisons for all clinical data
[~, ~, ~, ~]= stratifiedAnalyses(data, omic_combinations_label_permutation, ...
    combination_names_label_permutation, {'Age', 'Sex', 'Status'}, strcat(results_path, "\all_label_permutation.csv"));


function res = acc(pred_Y, Ytest)
    res = sum(pred_Y == Ytest) / length(Ytest);
end

function res = mcc(pred_Y, Ytest)
    if length(unique([pred_Y; Ytest])) > 1
        C = confusionmat(Ytest, pred_Y);
        TN = C(1, 1);
        TP = C(2, 2);
        FN = C(2, 1);
        FP = C(1, 2);
        res = (TP*TN-FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    else
        res = 0;
    end
end

function [unique_vals, accuracies, MCCs, sampleSizes] = stratifiedAnalyses(data, omics, names, columns, path)
    [unique_vals, ~, ~] = unique(data(:, columns));
    unique_vals = rmmissing(unique_vals);   % discard rows with Nans
    sampleSizes = zeros(1, size(unique_vals, 1));
    T = unique_vals;
    
    for i=1:length(omics)
        healthy_predicted = omics{i}.predicted(1:5, :);
        tumour_predicted = omics{i}.predicted(6:end, :);
        healthy_predicted = reshape(healthy_predicted, [1, numel(healthy_predicted)]);
        tumour_predicted = reshape(tumour_predicted, [1, numel(tumour_predicted)]);

        data_omic = data;
        data_omic.predicted = transpose(horzcat(healthy_predicted, tumour_predicted));
    
        accuracies = zeros(1, size(unique_vals, 1));
        MCCs = zeros(1, size(unique_vals, 1));
        for j=1:size(unique_vals, 1)
            filtered = data_omic(ismember(data_omic(:, columns), unique_vals(j, :)), {'predicted', 'true'});
            accuracies(j) = acc(filtered.predicted, filtered.true);
            MCCs(j) = mcc(filtered.predicted, filtered.true);
            if i==1   % only for the first omic, since the samples are the same for all of them
                sampleSizes(j) = size(filtered, 1);
            end
        end
        if i==1
            T.sampleSizes = transpose(sampleSizes);
        end
        T.(strcat('accuracy_', names{i})) = transpose(accuracies);
        T.(strcat('MCC_', names{i})) = transpose(MCCs);
    end
    writetable(T,path);
end