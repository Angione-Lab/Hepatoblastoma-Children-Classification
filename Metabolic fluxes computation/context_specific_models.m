clear, clc

load('recon2.2.mat');

f = zeros(length(model.rxns), 1);
model.c = f;                                
f(strcmp(model.rxns, 'biomass_reaction')) = 1; % biomass index 
model.c = f;                                   % set biomass as objective "classical" objective

gamma = 2;

addpath(genpath('path\MATLAB\cobratoolbox'));
% initCobraToolbox

[reaction_expression, pos_genes_in_react_expr, ixs_geni_sorted_by_length] = compute_reaction_expression(model);

% load gene expression and gene ids
addpath(char("path\data"));
load('Combined data\gene_exp');  
load('common_gene_ids');  
load('sample_ids');  

genes = model.genes;
genes_in_dataset = common_gene_ids; 
gene_exp = geneexp; 
Cmin{size(gene_exp, 1), length(model.rxns)} = [];
Cmax{size(gene_exp, 1), length(model.rxns)} = [];

GeneExpressionArray = ones(numel(genes),1); 

k = 1;
changeCobraSolver('gurobi', 'all');
%% Generate and save metabolic fluxes

% compute fold change
data = gene_exp./mean(gene_exp);

% apply the bounds and compute the fluxes for each gene espression profile
for t=1:size(data,1)          % select a unique profile
    expr_profile = data(t,:);
    pos_genes_in_dataset = zeros(numel(genes),1);
    for i=1:length(pos_genes_in_dataset)
        position = find(strcmp(genes{i},genes_in_dataset),1); 
        if ~isempty(position)                                   
            pos_genes_in_dataset(i) = position(1);              
            GeneExpressionArray(i) = expr_profile(pos_genes_in_dataset(i));         
        end
    end
    if or(sum(isinf(GeneExpressionArray)) >= 1, sum(isnan(GeneExpressionArray)) >= 1)
        fprintf("\nError in the gene expression data!");
    end
    [minFlux, maxFlux] = compute_fluxes(gamma, GeneExpressionArray, model, genes, reaction_expression, pos_genes_in_react_expr, ixs_geni_sorted_by_length);
	Cmin(k,:) = [sample_ids(k), num2cell(transpose(minFlux))];
	Cmax(k,:) = [sample_ids(k), num2cell(transpose(maxFlux))];
    k = k +1;
end
T = cell2table(Cmin);
T.Properties.VariableNames{1} = 'sample_id';  
writetable(T,"path\minFluxes.csv");
T = cell2table(Cmax);
T.Properties.VariableNames{1} = 'sample_id'; 
writetable(T,"path\maxFluxes.csv");
