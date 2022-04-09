function [reaction_expression, pos_genes_in_react_expr, ixs_geni_sorted_by_length] = compute_reaction_expression(fbamodel)

genesets = fbamodel.grRules;

%the genesets are already written nicely and ordered by reaction, so there
%is no need to assign each geneset to a reaction, so there is no need to
%use the following code. We only need to convert "AND" into "and", "OR"
%into "or" because of how associate_genes_reactions.m is written

genesets = regexprep(genesets,' AND ',' and '); 
genesets = regexprep(genesets,' OR ',' or '); 

reaction_expression = cell(length(genesets),1);
reaction_expression(:) = {''};

for i = 1:length(genesets)
    str_geneset = genesets{i};
    aux = associate_genes_reactions(str_geneset);
    reaction_expression{i} = aux;  
end

reaction_expression=strrep(reaction_expression,' ',''); %removes white spaces
geni = fbamodel.genes;

for i=1:size(geni)
    lung(i)=length(geni{i});
end

[~, ixs_geni_sorted_by_length] = sort(lung,'descend');
reaction_expression_aux = reaction_expression;

for i = 1:numel(ixs_geni_sorted_by_length)
    j = ixs_geni_sorted_by_length(i);
    matches = strfind(reaction_expression_aux,geni{j});     %this and the following instruction find the locations of the gene 'bXXXX' in the array reaction_expression
    pos_genes_in_react_expr{j} = find(~cellfun('isempty', matches));
    reaction_expression_aux(pos_genes_in_react_expr{j}) = strrep(reaction_expression_aux(pos_genes_in_react_expr{j}),geni{j},'');   %replaced with empty char so it's not found again later (this avoid for instance looking for 'HGNC:111' and replacing also partially another gene named for instance 'HGNC:1113'
end

end