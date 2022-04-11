function converted = fullConversion(data)
    converted.predicted = transpose(data.predicted); 
    converted.true = transpose(data.true);
    converted.mdl = data.mdl;
    converted.accuracies = transpose(data.accuracies);
    converted.MCCs = transpose(data.MCCs);
    converted.testSamples = transpose(data.testSamples);
    if isfield(data, 'predictors')
        converted.predictors = transpose(data.predictors);
        converted.pathways = transpose(data.pathways);
    end
    if isfield(data, 'genes')
        converted.genes = transpose(data.genes);
    end
    if isfield(data, 'plsda_scores_f')
        converted.plsda_score_f = data.plsda_scores_f;
    end
    if isfield(data, 'plsda_scores_g')
        converted.plsda_score_g = data.plsda_scores_g;       
    end
end
