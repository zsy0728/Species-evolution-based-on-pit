clc
clear
warning off

addpath(genpath('/home/songyao/matlab_path'))
temp_fs = vtkSurfRead('/home/songyao/data/cross-species/pit_count_others_fs6.vtk');
%% schaefer100
label_fs6 = temp_fs.Pdata{1,5}.val;
label_schaefer100 = zeros(1,50);
schaefer_100_fs6 = fetch_parcellation('fsaverage6', 'schaefer', 100);
for i = 1:50
    area = find(schaefer_100_fs6==i);
    label_schaefer100(i) = mode(label_fs6(area));
end

schaefer_100_fs5 = fetch_parcellation('fsaverage5', 'schaefer', 100);
label_fs5 = zeros(1,20484);
for i = 1:50
    area = find(schaefer_100_fs5==i);
    label_fs5(area) = label_schaefer100(i);
end

schaefer_1000 = fetch_parcellation('fsaverage5' ,'schaefer', 1000);
label_fs5_sch1000 = full2parcel(label_fs5,schaefer_1000);
label_lh = label_fs5_sch1000(1:500); 
label_shared = zeros(size(label_lh));  
label_unique = zeros(size(label_lh)); 
label_shared(label_lh > 0) = label_lh(label_lh > 0);
label_unique(label_lh < 0) = label_lh(label_lh < 0);

%% Genetics
%%%%%%%%%%%%%%%%% load gene data %%%%%%%%%%%%%%%%% 
[expression, gene_names] = surface_genetic_expression('schaefer', 1000);
expression_lh = expression(1:500,:);

[M, N] = size(expression_lh);
nan_rows = find(all(isnan(expression_lh), 2));
valid_rows = find(~all(isnan(expression_lh), 2));
for i = 1:length(nan_rows)
    row_idx = nan_rows(i);
    [~, nearest_idx] = min(abs(valid_rows - row_idx));
    nearest_row = valid_rows(nearest_idx);
    expression_lh(row_idx, :) = expression_lh(nearest_row, :);
end
%%%%%%%%%%%%%%%%% unique %%%%%%%%%%%%%%%%% 
X = zscore(expression_lh);  
Y = zscore(label_unique);   
% CV
[XL, YL, XS, YS, beta, PCTVAR, MSE, stats] = plsregress(X, Y, 20, 'CV', 10); 
figure('Color','w', 'Position', [100, 100, 600, 400]);
plot(0:20, MSE(2,:), '-o', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 6, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', [0.2 0.6 1]);  
xlabel('Number of Components', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('MSE (Specific Region)', 'FontSize', 14, 'FontWeight', 'bold');
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.2;
ax.Box = 'off';
ax.XLim = [0 20];
ax.XTick = 0:2:20;
grid on
saveas(gcf, '/home/songyao/data/cross-species/plsr_cv_unique.png')
% PLSR
nComp = 10;
[XL, YL, XS, YS, beta, PCTVAR, MSE, stats] = plsregress(X, Y, nComp);
comp_weights = PCTVAR(2, :) / sum(PCTVAR(2, :)); 
weighted_W = stats.W .* comp_weights;
gene_contributions = sum(abs(weighted_W), 2);
% gene_contributions = sum(abs(stats.W), 2);
[sorted_gene_contributions, sorted_idx] = sort(gene_contributions, 'descend');
top_2000_genes = sorted_idx(1:2000);
top_2000_gene_names = gene_names(top_2000_genes)';
writecell(top_2000_gene_names, '/home/songyao/data/cross-species/unique_top_2000_genes.csv');
T = array2table(top_2000_genes, 'VariableNames', {'entrez_id'});
writetable(T, '/home/songyao/data/cross-species/unique_top_2000_geneID.csv');

%%%%%%%%%%%%%%%%% shared %%%%%%%%%%%%%%%%% 
X = zscore(expression_lh);  
Y = zscore(label_shared);
% CV
[XL, YL, XS, YS, beta, PCTVAR, MSE, stats] = plsregress(X, Y, 20, 'CV', 10); 
figure('Color','w', 'Position', [100, 100, 600, 400]);
plot(0:20, MSE(2,:), '-o', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 6, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', [0.2 0.6 1]);  % 蓝色系填充
xlabel('Number of Components', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('MSE (Shared Region)', 'FontSize', 14, 'FontWeight', 'bold');
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.2;
ax.Box = 'off';
ax.XLim = [0 20];
ax.XTick = 0:2:20;
grid on
saveas(gcf, '/home/songyao/data/cross-species/plsr_cv_shared.png')
% PLSR
nComp = 10;
[XL, YL, XS, YS, beta, PCTVAR, MSE, stats] = plsregress(X, Y, nComp);
comp_weights = PCTVAR(2, :) / sum(PCTVAR(2, :)); 
weighted_W = stats.W .* comp_weights;
gene_contributions = sum(abs(weighted_W), 2);
% gene_contributions = sum(abs(stats.W), 2);
[sorted_gene_contributions, sorted_idx] = sort(gene_contributions, 'descend');
top_2000_genes = sorted_idx(1:2000);
top_2000_gene_names = gene_names(top_2000_genes)';
writecell(top_2000_gene_names, '/home/songyao/data/cross-species/shared_top_2000_genes.csv');
T = array2table(top_2000_genes, 'VariableNames', {'entrez_id'});
writetable(T, '/home/songyao/data/cross-species/shared_top_2000_geneID.csv');



