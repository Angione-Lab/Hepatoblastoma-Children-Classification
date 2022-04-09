
clear
packages = 'path';
addpath(genpath([packages '\cbrewer']));
addpath(genpath([packages '\fix_xticklabels']));
addpath(genpath([packages '\ihstevenson-beeswarm']));

fig = 'tumour_groups_PCA'; % 'sample_data', 'batch_correction_PCA', 'tumour_groups_PCA'

if strcmp(fig, 'sample_data')
    load('path\combined_dataset.mat')
    sampleData.Status = replace(sampleData.Status, 'Normal Liver', 'Control');
    sampleData.Status = replace(sampleData.Status, 'Tumor', 'Tumour');
    sampleData.Race = replace(sampleData.Race, 'Unknown', 'Unavailable');
    sampleData.Race = replace(sampleData.Race, 'Â NA', 'Unavailable');
    sampleData.Race(strcmp(sampleData.Race, '')) = {'Unavailable'};
    sampleData.Stage = replace(sampleData.Stage, 'N/A', 'Unavailable');
    sampleData.Stage(strcmp(sampleData.Stage, '')) = {'Unavailable'};
    sampleData.ClinicalCourse = replace(sampleData.ClinicalCourse, 'N/A', 'Unavailable');
    sampleData.ClinicalCourse(strcmp(sampleData.ClinicalCourse, '')) = {'Unavailable'};
    
    label_position_shrink_factor = 0.9;
    cb1 = cbrewer('qual', 'Set2', 3);
    cb2 = cbrewer('qual', 'Pastel2', 6);
    cb3 = cbrewer('qual', 'Set1', 3);
    
    ax1 = subplot(2, 3, 1);
    h = pie(categorical(sampleData.Status));
    h(2).Position = label_position_shrink_factor*h(2).Position;
    h(4).Position = label_position_shrink_factor*h(4).Position;
    title('Status')
    colormap(ax1, cb1(1:2, :))
    
    ax2 = subplot(2, 3, 2);
    h = pie(categorical(sampleData.Sex));
    h(2).Position = label_position_shrink_factor*h(2).Position;
    h(4).Position = [1 -1 0];
    title('Sex')
    colormap(ax2, cb2([1,3], :))
    
    ax3 = subplot(2, 3, 3);
    h = pie(categorical(sampleData.Race));
    h(2).Position = [0.65 1.1 0];
    h(4).Position = [-0.55 1.08 0];
    h(6).Position = 0.85*h(6).Position;
    h(8).Position = 0.85*h(8).Position;
    h(10).Position = 0.85*h(10).Position;
    h(12).Position = [1.2 -1 0];
    title('Race')
    colormap(ax3, cb2)
    
    ax4 = subplot(2, 3, 4);
    h = pie(categorical(sampleData.Stage));
    h(2).Position = 0.95*h(2).Position;
    h(4).Position = 0.95*h(4).Position;
    h(6).Position = 0.95*h(6).Position;
    h(8).Position = 0.95*h(8).Position;
    h(10).Position = [1.2 -1 0];
    title('Stage')
    colormap(ax4, cb2)
    
    ax5 = subplot(2, 3, 5);
    h = pie(categorical(sampleData.ClinicalCourse));
    h(2).Position = label_position_shrink_factor*h(2).Position;
    h(4).Position = label_position_shrink_factor*h(4).Position;
    h(6).Position = [1 1 0];
    title('Clinical course')
    colormap(ax5, cb2)
    
    ax6 = subplot(2, 3, 6);
    h = beeswarm(ones(sum(~isnan(sampleData.Age)), 1), sampleData.Age(~isnan(sampleData.Age)), 'corral_style', 'random', 'dot_size', .5);
    set(gca, 'xtick', [], 'xlim', [0.5 1.5])
    title('Age (yr)')
    box on
    
    hold on
    l = line([0.8 1.2], [1 1]*median(sampleData.Age(~isnan(sampleData.Age))), 'LineWidth', 1.5, 'Color', cb3(1, :));
    hold off
    legend(l, 'Median');
    
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperSize = [6.5 10.5];
    fig.PaperPosition = [-1.15 -0.5 12 7];
    fig.PaperOrientation = 'landscape';
    print('sample data', '-dpdf', '-r600');
    
elseif strcmp(fig, 'batch_correction_PCA')
    load('path\batch_correction.mat')
    batch_labels = [repmat({'GSE75271'}, size(sampleData1, 1), 1); repmat({'GSE131329'}, size(sampleData2, 1), 1); repmat({'E-MEXP-1851'}, size(sampleData3, 1), 1)];
    cb1 = cbrewer('qual', 'Set1', 3);
    
    [coeff,score,latent,tsquared,explained,mu] = pca(zscore(geneData'));
    subplot(1, 2, 1)
    gscatter(score(:,1), score(:,2), batch_labels, cb1)
    legend('Location', 'northwest')
    legend({}, 'FontSize', 20)
    xlabel(['PC1 (', num2str(round(explained(1), 1)), '%)'], 'FontSize', 20, 'FontWeight', 'bold')
    ylabel(['PC2 (', num2str(round(explained(2), 1)), '%)'], 'FontSize', 20, 'FontWeight', 'bold')
    title('Original data', 'FontSize', 20)
    
    [coeff,score,latent,tsquared,explained,mu] = pca(zscore(geneData_combat_p'));
    subplot(1, 2, 2)
    gscatter(score(:,1), score(:,2), batch_labels, cb1)
    legend('off')
    xlabel(['PC1 (', num2str(round(explained(1), 1)), '%)'], 'FontSize', 20, 'FontWeight', 'bold')
    ylabel(['PC2 (', num2str(round(explained(2), 1)), '%)'], 'FontSize', 20, 'FontWeight', 'bold')
    title('Batch-corrected data', 'FontSize', 20)
    
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperSize = [5 10.5];
    fig.PaperPosition = [-.85 .01 12 5];
    fig.PaperOrientation = 'landscape';
    print('batch correction PCA', '-dpdf', '-r600');
    
elseif strcmp(fig, 'tumour_groups_PCA')
    fluxes = readtable('path\maxFluxes.csv', 'ReadRowNames', true);
    load('path\combined_dataset.mat')
    fluxes = fluxes{:, :};
    sampleData.Status(strcmp(sampleData.Status, 'Tumor')) = {'Tumour'};
    cb2 = cbrewer('qual', 'Set2', 3);
    idx = strcmp(sampleData.Status, 'Tumour');
    
    [coeff,score,latent,tsquared,explained,mu] = pca(zscore(geneData_combat_p'));
    subplot(1, 2, 1)
    hold on
    scatter(score(~idx,1), score(~idx,2), 5+3*sampleData.Age(~idx), cb2(1,:), 'o', 'filled')
    scatter(score(idx,1), score(idx,2), 5+3*sampleData.Age(idx), cb2(2,:), 'o', 'filled')
    box on
    hold off

    legend({'Control', 'Tumour'}, 'Location', 'northwest', 'FontSize', 20)
    xlabel(['PC1 (', num2str(round(explained(1), 1)), '%)'], 'FontSize', 20, 'FontWeight', 'bold')
    ylabel(['PC2 (', num2str(round(explained(2), 1)), '%)'], 'FontSize', 20, 'FontWeight', 'bold')
    title('Gene expression', 'FontSize', 20)
    
    [coeff,score,latent,tsquared,explained,mu] = pca(zscore(fluxes));
    subplot(1, 2, 2)
    hold on
    scatter(score(~idx,1), score(~idx,2), 5+3*sampleData.Age(~idx), cb2(1,:), 'o', 'filled')
    scatter(score(idx,1), score(idx,2), 5+3*sampleData.Age(idx), cb2(2,:), 'o', 'filled')
    box on
    hold off

    legend({'Control', 'Tumour'}, 'Location', 'southwest', 'FontSize', 20)
    xlabel(['PC1 (', num2str(round(explained(1), 1)), '%)'], 'FontSize', 20, 'FontWeight', 'bold')
    ylabel(['PC2 (', num2str(round(explained(2), 1)), '%)'], 'FontSize', 20, 'FontWeight', 'bold')
    title('Metabolic fluxes', 'FontSize', 20)
    
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperSize = [5 10.5];
    fig.PaperPosition = [-.85 .01 12 5];
    fig.PaperOrientation = 'landscape';
    print('tumour groups PCA', '-dpdf', '-r600');

end
