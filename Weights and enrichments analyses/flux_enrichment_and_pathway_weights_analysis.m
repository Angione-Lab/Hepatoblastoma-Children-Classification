clear, clc

min_max = "min";     % select whether to run the analysis on the min or the max fluxes, either "min" or "max"
remove_irrev = ""; % valid only if min_max = "min", if "remove" than remove all the irreversible reactions from the model (only for the enrichment)

load("path\recon2.2.mat");
load("path\metadata.mat");

file_path = "path";

path_rev = ""; % changes only if min_max = "min" and remove_irrev = "remove"
if min_max == "max"
    reactions = readtable(strcat(file_path, '\FVA\maxFluxes.csv'));
    reactions = reactions(:, 2:end);
    reactions = table2array(reactions);
    enrichment_save_path = '\FVA\Flux_enrichment_max';
    pathway_weigths_save_path = '\FVA\Pathway weights_max\';
else
    reactions = readtable(strcat(file_path, '\FVA\minFluxes.csv'));
    reactions = reactions(:, 2:end);
    reactions = table2array(reactions);
    enrichment_save_path = '\FVA\Flux_enrichment_min';
    pathway_weigths_save_path = '\FVA\Pathway weights_min\';
        
    if remove_irrev == "remove"     % remove non-reversible reactions
        indices_non_rev = model.lb >= 0;
        reactions_rev = reactions(:, ~indices_non_rev);
        model_rev = removeRxns(model, model.rxns(indices_non_rev));
        path_rev = "_rev";
    else
        reactions_rev = reactions;
        model_rev = model;
    end
end

    %% Flux enrichment according to the health status
tumour = table2cell(metadata(:, 2)) == "Tumor";  % rows associated to tumour samples

[health_flux_enrichment, health_pathways, health_significant_pathways] = ...
    flux_enrichment_and_pathways_analysis(tumour, ~tumour, reactions_rev, 1e-7, model_rev, 0.05, 0.05, ...
    strcat(file_path, strcat(enrichment_save_path, '\Health status\flux_enrichment_analysis', path_rev, '.csv')));
%% Flux enrichment according to the gender
male = table2cell(metadata(:, 4)) == "Male";  % rows associated to male samples

[gender_flux_enrichment, gender_pathways, gender_significant_pathways] = ...
    flux_enrichment_and_pathways_analysis(male, ~male, reactions_rev, 1e-7, model_rev, 0.05, 0.05, ...
    strcat(file_path, strcat(enrichment_save_path, '\Gender\flux_enrichment_analysis', path_rev, '.csv')));
%% Flux enrichment according to the age
rows_to_exclude = find(cellfun(@isnan,table2cell(metadata(:, 3))));  % rows with no age information
age_metadata = metadata;
age_metadata(rows_to_exclude,:) = [];
age = table2array(age_metadata(:, 3)) > 4.5; % rows associated to individuals older than 4.5 years
age2 = table2array(age_metadata(:, 3)) > 3; % rows associated to individuals older than 3 years

% first threshold (4.5 years old)
[age_flux_enrichment, age_pathways, age_significant_pathways] = ...
    flux_enrichment_and_pathways_analysis(age, ~age, reactions_rev, 1e-7, model_rev, 0.05, 0.05, ...
    strcat(file_path, strcat(enrichment_save_path, '\Age\age_flux_enrichment_analysis', path_rev, '.csv')));

% second threshold (3 years old)
[age2_flux_enrichment, age2_pathways, age2_significant_pathways] = ...
    flux_enrichment_and_pathways_analysis(age2, ~age2, reactions_rev, 1e-7, model_rev, 0.05, 0.05, ...
    strcat(file_path, strcat(enrichment_save_path, '\Age\age2_flux_enrichment_analysis', path_rev, '.csv')));
%% Flux enrichment according to the clinical course
alive = table2cell(metadata(:, 8)) == "Alive";  % rows associated to alive samples
dead = table2cell(metadata(:, 8)) == "Dead";  % rows associated to dead samples

[clinical_flux_enrichment, clinical_pathways, clinical_significant_pathways] = ...
    flux_enrichment_and_pathways_analysis(alive, dead, reactions_rev, 1e-7, model_rev, 0.05, 0.05, ...
    strcat(file_path, strcat(enrichment_save_path, '\Clinical course\flux_enrichment_analysis', path_rev, '.csv')));
%% Shared enriched pathways
fprintf("\nShared enriched pathways among the different stratifications\n");
% total intersection 
intersection1 = intersect(health_significant_pathways, gender_significant_pathways); 
intersection2 = intersect(age_significant_pathways, age2_significant_pathways);
intersection3 = intersect(intersection1, clinical_significant_pathways);
intersection_final = intersect(intersection3, intersection2);
fprintf("\Total intersection:\n");
fprintf("\nNum dimensions of intersection: %d\n", length(intersection_final));
fprintf(1, '%s \n ', intersection_final{:})
%% Average pathway weights
stratifications = {tumour, ~tumour, male, ~male, age, ~age, age2, ~age2, alive, dead};
stratifications_p_sig = {health_significant_pathways, gender_significant_pathways, ...
    age_significant_pathways, age2_significant_pathways, clinical_significant_pathways};
stratifications_p_no_sig = {health_pathways, gender_pathways, age_pathways, age2_pathways, clinical_pathways};
stratifications_files = {"Health status\tumour", "Health status\control", "Gender\male", "Gender\female", ...
    "Age\old", "Age\young", "Age\old2", "Age\young2", "Clinical course\alive", "Clinical course\dead"};

for st=1:length(stratifications)
% compute the mean weight for each relevant pathway and save
    mean_reactions = mean(abs(reactions(stratifications{st}, :)));        % mean along each column
    if mod(st, 2)
        strat = stratifications_p_sig{fix(st/2)+1};
    else
        strat = stratifications_p_sig{st/2}; 
    end
    average_weigth = zeros(size(strat)); 
    for j=1:size(strat, 1)
        k = 0;                                                            % number of times we find the pathway
        for i=1:size(model.subSystems, 1) 
            flat = model.subSystems(i);
            flat = [flat{:}];                                             
            if ismember(strat(j), flat)               % if the flux belongs to the pathway
                average_weigth(j) = average_weigth(j) + mean_reactions(i);   
                k = k + 1;
            end
        end
        average_weigth(j) = average_weigth(j) / k;                                  
    end
    writetable(cell2table(cellstr(horzcat(strat, num2cell(average_weigth)))), strcat(file_path, pathway_weigths_save_path, stratifications_files{st}, '_average_pathway_weights.csv'));
end

% all the pathways
for st=1:length(stratifications)
    mean_reactions = mean(abs(reactions(stratifications{st}, :)));        % mean along each column
    if mod(st, 2)
        strat = stratifications_p_no_sig{fix(st/2)+1};
    else
        strat = stratifications_p_no_sig{st/2}; 
    end
    average_weigth = zeros(size(strat)); 
    for j=1:size(strat, 1)
        k = 0;                                                            % number of times we find the pathway
        for i=1:size(model.subSystems, 1) 
            flat = model.subSystems(i);
            flat = [flat{:}];                                             
            if ismember(strat(j), flat)               % if the flux belongs to the pathway
                average_weigth(j) = average_weigth(j) + mean_reactions(i);   
                k = k + 1;
            end
        end
        average_weigth(j) = average_weigth(j) / k;                                  
    end
    writetable(cell2table(cellstr(horzcat(strat, num2cell(average_weigth)))), strcat(file_path, pathway_weigths_save_path, stratifications_files{st}, '_no_sig_average_pathway_weights.csv'));
end
