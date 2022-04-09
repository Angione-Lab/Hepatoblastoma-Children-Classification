function [minFlux, maxFlux] = compute_fluxes(gamma, gene_expression, model, genes, reaction_expression, pos_genes_in_react_expr, ixs_geni_sorted_by_length)
% function modified from evaluate_objective.m, for code explanation refer to that file
yt=gene_expression';      
eval_reaction_expression = reaction_expression;

for i=ixs_geni_sorted_by_length 
    posizioni_gene = pos_genes_in_react_expr{i};
    for j=1:length(posizioni_gene)  
        eval_reaction_expression{posizioni_gene(j)} = strrep(eval_reaction_expression{posizioni_gene(j)}, genes{i}, num2str(yt(i),'%.15f'));  
    end
end

eval_reaction_expression( cellfun(@isempty, eval_reaction_expression) ) = {'1.0'};  
num_reaction_expression = zeros(1,length(eval_reaction_expression));

for i=1:length(num_reaction_expression)
    str = eval_reaction_expression{i};
    num_parenthesis = numel(strfind(str,')'));
    while (num_parenthesis > 32) 
        to_replace = 'min.\d*+\.+\d*,\d*+\.+\d*.|max.\d*+\.+\d*,\d*+\.+\d*.|min..\d*+\.+\d*.,\d*+\.+\d*.|max..\d*+\.+\d*.,\d*+\.+\d*.|min..\d*+\.+\d*.,.\d*+\.+\d*..|max..\d*+\.+\d*.,.\d*+\.+\d*..|min.\d*+\.+\d*,.\d*+\.+\d*..|max.\d*+\.+\d*,.\d*+\.+\d*..';  %searches for all the strings of kind min(NUM.NUM,NUM.NUM) or max(NUM.NUM,NUM.NUM) or  min((NUM.NUM),NUM.NUM) or max((NUM.NUM),NUM.NUM) or  min((NUM.NUM),(NUM.NUM)) or max(NUM.NUM,(NUM.NUM)) or  min(NUM.NUM,(NUM.NUM)) or max((NUM.NUM),(NUM.NUM))
        substrings_to_replace = regexp(str, to_replace, 'match');
        if isempty(substrings_to_replace)
            num_parenthesis = 0; 
        else
            for j = 1:numel(substrings_to_replace)
                ss_rep = substrings_to_replace{j};
                str = strrep(str,ss_rep,num2str(eval(ss_rep),'%.15f'));
            end
            num_parenthesis = numel(strfind(str,')'));
        end
    end
    str = regexprep(str,'/','');
    try
        num_reaction_expression(i) = eval(str);
    catch
        num_reaction_expression(i) = 1;
    end
end

reaction_expressionUpCorrected = prctile(num_reaction_expression,99);
num_reaction_expression(num_reaction_expression > reaction_expressionUpCorrected) = reaction_expressionUpCorrected;

if or(sum(isinf(num_reaction_expression)) >= 1, sum(isnan(num_reaction_expression)) >= 1)
    fprintf("\nError in the evaluated and corrected gene expression data!");
end

for i=1:length(num_reaction_expression)   %loop over the array of the geneset expressions
    if num_reaction_expression(i)>=1
        model.lb(i) = model.lb(i)*(1+gamma*log(num_reaction_expression(i)));
        model.ub(i) = model.ub(i)*(1+gamma*log(num_reaction_expression(i)));
    else
        model.lb(i) = model.lb(i)/(1+gamma*abs(log(num_reaction_expression(i))));
        model.ub(i) = model.ub(i)/(1+gamma*abs(log(num_reaction_expression(i))));        
    end
end

if or(sum(isinf(model.lb)) >= 1, sum(isnan(model.lb)) >= 1)
    fprintf("\nError in the new lower bounds!");
end
if or(sum(isinf(model.ub)) >= 1, sum(isnan(model.ub)) >= 1)
    fprintf("\nError in the new upper bounds!");
end

%%Setting of the manual bounds from the literature
% we set them here (except for the import reactions) because the mapping is done with an arbitrary function, which means that it cannot override the experimental values

model = changeRxnBounds(model, 'EX_glc(e)', 2.025, 'l'); 
model = changeRxnBounds(model, 'EX_his_L(e)', -0.04425, 'l');
model = changeRxnBounds(model, 'EX_ile_L(e)', -0.0585, 'l');
model = changeRxnBounds(model, 'EX_leu_L(e)', -0.0825, 'l');
model = changeRxnBounds(model, 'EX_lys_L(e)', -0.2325, 'l');
model = changeRxnBounds(model, 'EX_met_L(e)', -0.12, 'l');
model = changeRxnBounds(model, 'EX_phe_L(e)', -0.202774, 'l'); 
model = changeRxnBounds(model, 'EX_thr_L(e)', -0.12, 'l');
model = changeRxnBounds(model, 'EX_trp_L(e)', -0.0075, 'l');
model = changeRxnBounds(model, 'EX_val_L(e)', -0.04125, 'l'); 
model = changeRxnBounds(model, 'EX_h2o(e)', 25.3228, 'l');
model = changeRxnBounds(model, 'EX_o2(e)', -28.05, 'l');
model = changeRxnBounds(model, 'EX_co2(e)', 21.7219, 'l');
model = changeRxnBounds(model, 'EX_ala_L(e)', -0.02325, 'l');
model = changeRxnBounds(model, 'EX_asn_L(e)', -0.00135, 'l');
model = changeRxnBounds(model, 'EX_gln_L(e)', -2.325, 'l');
model = changeRxnBounds(model, 'EX_tyr_L(e)', -0.05775, 'l');
model = changeRxnBounds(model, 'EX_cys_L(e)', -0.0555, 'l');
model = changeRxnBounds(model, 'EX_arg_L(e)', -0.2175, 'l');
model = changeRxnBounds(model, 'EX_gly(e)', -0.2625, 'l');
model = changeRxnBounds(model, 'EX_pro_L(e)', 0.02925, 'l'); 
model = changeRxnBounds(model, 'EX_ser_L(e)', -0.1425, 'l');
model = changeRxnBounds(model, 'EX_asp_L(e)', 0.00825, 'l');
model = changeRxnBounds(model, 'EX_glu_L(e)', 0.15, 'l');
model = changeRxnBounds(model, 'EX_nh4(e)', -0.165, 'l');
model = changeRxnBounds(model, 'EX_so4(e)', 0.16121, 'l');
model = changeRxnBounds(model, 'EX_h(e)', -0.42825, 'l'); 
model = changeRxnBounds(model, 'EX_glyc(e)', -6.675, 'l');
model = changeRxnBounds(model, 'EX_orn(e)', 0.125, 'l');
model = changeRxnBounds(model, 'EX_acac(e)', 0.1275, 'l'); 
model = changeRxnBounds(model, 'EX_bhb(e)', 0.05775, 'l'); 
model = changeRxnBounds(model, 'EX_lac_L(e)', -0.063, 'l');
model = changeRxnBounds(model, 'EX_urea(e)', 3.375, 'l');

for i=1:length(model.lb)
    if model.lb(i) > model.ub(i)
        fprintf("Index: %d, lower bound: %f, upper bound: %f", i, model.lb(i), model.ub(i));
    end
end

%% Compute Flux Variability Analysis
minFlux, maxFlux] = fluxVariability(model, 100, 'max', model.rxns);

format longG; format compact;
