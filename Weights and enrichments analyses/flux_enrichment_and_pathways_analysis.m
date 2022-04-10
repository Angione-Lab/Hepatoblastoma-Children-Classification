function [flux_enrichment, pathways, significant_pathways] = flux_enrichment_and_pathways_analysis(rows1, rows2, fluxes, threshold, model, pvalue_f, pvalue_p, save_path)
% Function to compute flux enrichment and pathways (statistically significant and non-significant) on the fluxes which are deemed statistically significantly different among the two groups
% Input:
%       rows1: indices of the samples of interest for the first group
%       rows2: indices of the samples of interest for the second group
%       fluxes: flux data to use for the Wilcoxon rank-sum test and the flux enrichment analysis
%       threshold: fluxes whose value is lower then threshold for all samples will be discarded before the analysis
%       model: model to use
%       pvalue_f: p_value below which the fluxes are considered statistically significantly EQUAL among the two groups
%       pvalue_p: p_value below which the pathways are considered statistically significant
%       save_path: path for saving the flux_enrichment table
%
% N.B. pvalue_f defines the threshold ABOVE WHICH the fluxes are considered DIFFERENT
%
    
    indices = find(all(abs(fluxes)> threshold, 1)); % select only the reactions which are > threshold for all the samples
    group1 = fluxes(rows1, indices);
    group2 = fluxes(rows2, indices);
    pvalues_w = zeros(length(indices), 1);
    
    for i=1:length(indices)       % compute p-values for each reaction pair among the two groups
        [p,~] = ranksum(group1(:,i),group2(:,i));
        pvalues_w(i) = p;
    end
    adjusted_pvalues = mafdr(pvalues_w, 'BHFDR', true);
    final_indices = indices((adjusted_pvalues > pvalue_f)); % the null hypothesis is rejected, i.e. the two distributions are different 

    flux_enrichment = FEA(model, final_indices, 'subSystems');

    p_values = cell2mat(flux_enrichment(2:end, 2));  
    pathways = string(flux_enrichment(2:end, 3));
    significant_pathways = pathways(p_values <= pvalue_p);

    writetable(cell2table(flux_enrichment), save_path, 'WriteVariableNames', false);
end