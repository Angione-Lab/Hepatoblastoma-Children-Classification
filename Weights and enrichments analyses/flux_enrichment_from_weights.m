% flux enrichment analysis based on flux weights from the experiment runs
clc, clear

data_path = 'path\';

load('path\recon2.2.mat');

metric = "median_"; % "", "median_"

FVA_weights = readtable(strcat(data_path, metric, 'FVA_weights.csv'));
FVA_clinical_weights = readtable(strcat(data_path, metric, 'FVA_clinical_weights.csv'));
FVA_genes_weights = readtable(strcat(data_path, metric, 'FVA_genes_weights.csv'));
FVA_all_weights = readtable(strcat(data_path, metric, 'FVA_all_weights.csv'));
omics = {FVA_weights, FVA_clinical_weights, FVA_genes_weights, FVA_all_weights};

C = {};
k =1;
max_el = 0;

for i=1:size(omics, 2)
    percentile = prctile(table2array(omics{i}(:, 3)), 97);
    mask = table2array(omics{i}(:, 3)) >= percentile;
    [~,indices] = ismember(table2cell(omics{i}(mask, 1)), model.rxns);
    flux_enrichment = FEA(model, indices, 'subSystems');
    enriched_indices = cell2mat(flux_enrichment(2:end, 2)) <= 0.05;
    enriched_indices = logical([0; enriched_indices]);
    if sum(enriched_indices) > max_el
        max_el = sum(enriched_indices);
    end
    C{k} = flux_enrichment(enriched_indices, 2);
    C{k+1} = flux_enrichment(enriched_indices, 3);
    k = k +2;
end
T = cell(max_el, 2 * size(omics, 2));
T(:) = {'N/A'};
for i=1:size(T, 2)
    temp = C(i);
    T(1:numel(temp{:}), i) = temp{:};
end
column_names = {'Adjusted P-values for FVA_weights' 'Pathways for FVA_weights' ... 
    'Adjusted P-values for FVA_clinical_weights' 'Pathways for FVA_clinical_weights' ... 
    'Adjusted P-values for FVA_genes_weights' 'Pathways for FVA_genes_weights' ... 
    'Adjusted P-values for FVA_all_weights' 'Pathways for FVA_all_weights'};
T = vertcat(column_names, T);
writetable(cell2table(T), strcat(data_path, metric, "FVA_enrichment.csv"));

