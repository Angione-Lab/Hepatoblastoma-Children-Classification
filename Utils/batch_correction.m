
addpath('path\BatchRemove.build'); % from https://github.com/weikanggong/BatchRemove.build

geneData1 = readtable([path1 'gene_exp_GSE75271.csv'], 'ReadRowNames', true);
geneData2 = readtable([path2 'gene_exp_GSE131329.csv'], 'ReadRowNames', true);
geneData3 = readtable([path3 'gene_exp_E_MEXP_1851.csv'], 'ReadRowNames', true);
sampleData1 = readtable([path1 'metadata_GSE75271.csv']);
sampleData2 = readtable([path2 'metadata_GSE131329.csv']);
sampleData3 = readtable([path3 'metadata_E_MEXP_1851.csv']);

sampleData1.Sex(strcmp(sampleData1.Sex, 'M')) = {'Male'};
sampleData1.Sex(strcmp(sampleData1.Sex, 'F')) = {'Female'};
for i = 1:length(sampleData1.Age)
    sampleData1.Age{i} = str2num(sampleData1.Age{i}(8:end));
end
sampleData1.Age{cellfun(@isempty, sampleData1.Age)} = NaN;
sampleData1.Age = cell2mat(sampleData1.Age);
sampleData2.Status(strcmp(sampleData2.Status, 'noncancerous liver tissue')) = {'Normal Liver'};
sampleData2.Status(strcmp(sampleData2.Status, 'tumor tissue')) = {'Tumor'};
for i = 1:length(sampleData2.Age)
    sampleData2.Age{i} = str2num(sampleData2.Age{i}(10:end)) / 12;
end
sampleData2.Age = cell2mat(sampleData2.Age);
sampleData3.Properties.VariableNames = {'ID', 'Age', 'Unit', 'Sex', 'Status', 'ClinicalCourse'};
sampleData3.Sex(strcmp(sampleData3.Sex, 'male')) = {'Male'};
sampleData3.Sex(strcmp(sampleData3.Sex, 'female')) = {'Female'};
sampleData3.Status(strcmp(sampleData3.Status, 'normal')) = {'Normal Liver'};
sampleData3.Status(strcmp(sampleData3.Status, 'hepatoblastoma')) = {'Tumor'};
sampleData3.Age = sampleData3.Age / 12;
sampleData3.Unit = [];

geneData = [geneData1{:, :} geneData2{:, :} geneData3{:, :}];
sampleData = outerjoin(sampleData1, sampleData2, 'MergeKeys', true);
sampleData = outerjoin(sampleData, sampleData3, 'MergeKeys', true);

sampleData_binary = dummyvar({sampleData.Status, sampleData.Sex});
sampleData_binary = sampleData_binary(:, [1, 3]);
batch = [ones(size(sampleData1, 1), 1); 2*ones(size(sampleData2, 1), 1); 3*ones(size(sampleData3, 1), 1)];
batch_labels = [repmat({'GSE75271'}, size(sampleData1, 1), 1); repmat({'GSE131329'}, size(sampleData2, 1), 1); repmat({'Third'}, size(sampleData3, 1), 1)];

[coeff,score,latent,tsquared,explained,mu] = pca(zscore(geneData'));
subplot(1, 3, 1)
gscatter(score(:,1), score(:,2), batch_labels)

% parametric version of ComBat
geneData_combat_p = ComBat(geneData, batch, sampleData_binary, 1);

[coeff,score,latent,tsquared,explained,mu] = pca(zscore(geneData_combat_p'));
subplot(1, 3, 2)
gscatter(score(:,1), score(:,2), batch_labels)

save('batch_correction')
save('combined_dataset', 'geneData_combat_p', 'sampleData')
