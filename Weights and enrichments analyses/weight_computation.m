% weights computation for reactions, genes and pathways
clc, clear

save_path = 'path\';

addpath(genpath('path\SVM results'));
addpath(genpath('path\countmember'));

%% load data
% fluxes and pathways
load('path\recon2.2.mat');
load('SVMl2_results_FVA.mat');
FVA_results = fullConversion(SVM_results);   
load('Int3SVMl2_results_FVA.mat');
FVA_clinical_results = fullConversion(SVM_results);
% genes
load('common_gene_ids.mat');
load('SVMl2_results_genes.mat');
genes_results = fullConversion(SVM_results);   
load('Int3SVMl2_results_genes.mat');
genes_clinical_results = fullConversion(SVM_results);
% fluxes and genes
load('Int2SVMl2_results.mat');
FVA_genes_results = fullConversion(SVM_results);
load('Int4SVMl2_results.mat');
all_results = fullConversion(SVM_results);
clear SVM_results;

rxn_weights = zeros(size(model.rxns));
int3rxn_weights = zeros(size(model.rxns));
int2rxn_weights = zeros(size(model.rxns));
int4rxn_weights = zeros(size(model.rxns));

pathways_weights = zeros(size(unique(model.subSystems)));
int3pathways_weights = zeros(size(unique(model.subSystems)));
int2pathways_weights = zeros(size(unique(model.subSystems)));
int4pathways_weights = zeros(size(unique(model.subSystems)));

genes_weights = zeros(size(common_gene_ids));
int3genes_weights = zeros(size(common_gene_ids));
int2genes_weights = zeros(size(common_gene_ids));
int4genes_weights = zeros(size(common_gene_ids));
%% compute weights
start_time = tic;
for i=1:size(FVA_results.mdl, 1) % loop across the repetitions
    for j=1:size(FVA_results.mdl{i}.Beta) % loop across the weights/variables
        % reactions
        [indices_r] = ismember(model.rxns, FVA_results.predictors{i}{j});  % same for all simulations
        rxn_weights(indices_r) = rxn_weights(indices_r) + abs(FVA_results.mdl{i}.Beta(j)) ./ numel(FVA_results.predictors{i}{j}); % division in equal parts
        % reactions from integration with clinical data
        int3FVA_betas = abs(FVA_clinical_results.plsda_score_f{i} * FVA_clinical_results.mdl{i}.Beta(1:2));
        int3rxn_weights(indices_r) = int3rxn_weights(indices_r) + int3FVA_betas(j) ./ numel(FVA_results.predictors{i}{j}); % apart from the betas, 
                                                                                                                                % the rest is the same for all simulations
        % pathways
        unique_pathways = unique(FVA_results.pathways{i}{j});
        [indices_p] = ismember(unique(model.subSystems), unique_pathways);   % same for all simulations
        pathways_occurrences = countmember(unique_pathways, FVA_results.pathways{i}{j});
        pathways_weights(indices_p) = pathways_weights(indices_p) + abs(FVA_results.mdl{i}.Beta(j)) ./ numel(FVA_results.pathways{i}{j}) .* pathways_occurrences;
        % pathways from integration with clinical data
        int3pathways_weights(indices_p) = int3pathways_weights(indices_p) + int3FVA_betas(j) ./ numel(FVA_results.pathways{i}{j}) .* pathways_occurrences;
        
        % genes
        [indices_g] = ismember(common_gene_ids, genes_results.genes{i}{j});  % same for all simulations
        genes_weights(indices_g) = genes_weights(indices_g) + abs(genes_results.mdl{i}.Beta(j)) ./ numel(genes_results.genes{i}{j});
        % genes from integration with clinical data
        int3genes_betas = abs(genes_clinical_results.plsda_score_g{i} * genes_clinical_results.mdl{i}.Beta(1:2));
        int3genes_weights(indices_g) = int3genes_weights(indices_g) + int3genes_betas(j) ./ numel(genes_results.genes{i}{j});        
        
        % fluxes and genes
        int2FVA_betas = abs(FVA_genes_results.plsda_score_f{i} * FVA_genes_results.mdl{i}.Beta(1:2));
        int2genes_betas = abs(FVA_genes_results.plsda_score_g{i} * FVA_genes_results.mdl{i}.Beta(3:4));
        int2rxn_weights(indices_r) = int2rxn_weights(indices_r) + int2FVA_betas(j) ./ numel(FVA_results.predictors{i}{j}); 
        int2pathways_weights(indices_p) = int2pathways_weights(indices_p) + int2FVA_betas(j) ./ numel(FVA_results.pathways{i}{j}) .* pathways_occurrences;
        int2genes_weights(indices_g) = int2genes_weights(indices_g) + int2genes_betas(j) ./ numel(genes_results.genes{i}{j});        

        int4FVA_betas = abs(all_results.plsda_score_f{i} * all_results.mdl{i}.Beta(1:2));
        int4genes_betas = abs(all_results.plsda_score_g{i} * all_results.mdl{i}.Beta(3:4));
        int4rxn_weights(indices_r) = int4rxn_weights(indices_r) + int4FVA_betas(j) ./ numel(FVA_results.predictors{i}{j}); 
        int4pathways_weights(indices_p) = int4pathways_weights(indices_p) + int4FVA_betas(j) ./ numel(FVA_results.pathways{i}{j}) .* pathways_occurrences;
        int4genes_weights(indices_g) = int4genes_weights(indices_g) + int4genes_betas(j) ./ numel(genes_results.genes{i}{j});
        
        % safety check
        if or(sum(cellfun(@iscell, genes_results.genes{i}{j})), or(sum(cellfun(@iscell, FVA_results.predictors{i}{j})), sum(cellfun(@iscell, FVA_results.pathways{i}{j}))))
            fprintf("\nNested cell array! i=%d, j=%d", i, j);
        end
    end
end
end_time = toc(start_time);
fprintf("\nTotal duration of weight computation: %f minutes\n", end_time./60);
%% save data
% reactions
idx_r = find(rxn_weights);
writetable(cell2table(horzcat(model.rxns(idx_r), model.rxnNames(idx_r), num2cell(rxn_weights(idx_r)), ... 
    repmat({'FVA'}, size(idx_r))), 'VariableNames', ... 
    {'Reaction' 'Reaction_Name' 'Weight' 'Omic_combination'}), strcat(save_path, 'FVA_weights.csv'));

writetable(cell2table(horzcat(model.rxns(idx_r), model.rxnNames(idx_r), num2cell(int3rxn_weights(idx_r)), ... 
    repmat({'FVA_clinical'}, size(idx_r))), 'VariableNames', ... 
    {'Reaction' 'Reaction_Name' 'Weight' 'Omic_combination'}), strcat(save_path, 'FVA_clinical_weights.csv'));

writetable(cell2table(horzcat(model.rxns(idx_r), model.rxnNames(idx_r), num2cell(int2rxn_weights(idx_r)), ... 
    repmat({'FVA_genes'}, size(idx_r))), 'VariableNames', ... 
    {'Reaction' 'Reaction_Name' 'Weight' 'Omic_combination'}), strcat(save_path, 'FVA_genes_weights.csv'));

writetable(cell2table(horzcat(model.rxns(idx_r), model.rxnNames(idx_r), num2cell(int4rxn_weights(idx_r)), ... 
    repmat({'all'}, size(idx_r))), 'VariableNames', ... 
    {'Reaction' 'Reaction_Name' 'Weight' 'Omic_combination'}), strcat(save_path, 'FVA_all_weights.csv'));
% pathways
unique_pathways = unique(model.subSystems);
idx_p = find(pathways_weights);
writetable(cell2table(horzcat(unique_pathways(idx_p), num2cell(pathways_weights(idx_p)), ... 
    repmat({'FVA'}, size(idx_p))), 'VariableNames', ... 
    {'Pathway' 'Weight' 'Omic_combination'}), strcat(save_path, 'pathways_weights.csv'));

writetable(cell2table(horzcat(unique_pathways(idx_p), num2cell(int3pathways_weights(idx_p)), ... 
    repmat({'FVA_clinical'}, size(idx_p))), 'VariableNames', ... 
    {'Pathway' 'Weight' 'Omic_combination'}), strcat(save_path, 'pathways_clinical_weights.csv'));

writetable(cell2table(horzcat(unique_pathways(idx_p), num2cell(int2pathways_weights(idx_p)), ... 
    repmat({'FVA_genes'}, size(idx_p))), 'VariableNames', ... 
    {'Pathway' 'Weight' 'Omic_combination'}), strcat(save_path, 'pathways_genes_weights.csv'));

writetable(cell2table(horzcat(unique_pathways(idx_p), num2cell(int4pathways_weights(idx_p)), ... 
    repmat({'all'}, size(idx_p))), 'VariableNames', ... 
    {'Pathway' 'Weight' 'Omic_combination'}), strcat(save_path, 'pathways_all_weights.csv'));
% genes
idx_g = find(genes_weights);
genes_array = cellstr(common_gene_ids);
writetable(cell2table(horzcat(genes_array(idx_g), num2cell(genes_weights(idx_g)), ... 
    repmat({'genes'}, size(idx_g))), 'VariableNames', ... 
    {'Gene' 'Weight' 'Omic_combination'}), strcat(save_path, 'genes_weights.csv'));

writetable(cell2table(horzcat(genes_array(idx_g), num2cell(int3genes_weights(idx_g)), ... 
    repmat({'genes_clinical'}, size(idx_g))), 'VariableNames', ... 
    {'Gene' 'Weight' 'Omic_combination'}), strcat(save_path, 'genes_clinical_weights.csv'));

writetable(cell2table(horzcat(genes_array(idx_g), num2cell(int2genes_weights(idx_g)), ... 
    repmat({'FVA_genes'}, size(idx_g))), 'VariableNames', ... 
    {'Gene' 'Weight' 'Omic_combination'}), strcat(save_path, 'genes_FVA_weights.csv'));

writetable(cell2table(horzcat(genes_array(idx_g), num2cell(int4genes_weights(idx_g)), ... 
    repmat({'all'}, size(idx_g))), 'VariableNames', ... 
    {'Gene' 'Weight' 'Omic_combination'}), strcat(save_path, 'genes_all_weights.csv'));


