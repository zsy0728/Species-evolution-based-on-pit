% Species kinship and cortical similarity v4
% 20250423
% There is no need for matrix automatic clustering, and both rows and columns are sorted in the order of the phylogenetic tree
clc
clear
warning off

infodir = '/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/';
load('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/evolutionary_path_of_all_species.mat')

atlasdir = '/mnt/sda/songyao/results/Evolution_cortical_shape/atlas/';
surfl = vtkSurfRead([atlasdir,'sub-020_species-Homo+sapiens_hemi-L_topo-Homo.sapiens.surf_altas.vtk']);
surfr = vtkSurfRead([atlasdir,'sub-020_species-Homo+sapiens_hemi-R_topo-Homo.sapiens.surf_altas.vtk']);
yeo7_lh = surfl.Pdata{1,1}.val;
yeo7_rh = surfr.Pdata{1,1}.val;
Schaefer100_lh = surfl.Pdata{1,6}.val;
Schaefer100_rh = surfr.Pdata{1,6}.val;
Schaefer300_lh = surfl.Pdata{1,12}.val;
Schaefer300_rh = surfr.Pdata{1,12}.val;
load('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic_landmarks/pit_id_all_LH.mat')
load('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/evolutionary_path_of_all_species.mat')
load('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/kinship/Species_kinship_by_phyTree.mat')

label_l = Schaefer100_lh;
label_r = Schaefer100_rh;

ratio = zeros(length(unique(label_l))-1,90);
for species = 1:90
    evolution_pit_lh = pit_id_all{all_paths{1,species}(1)};
    region_lh = label_l(evolution_pit_lh);
    region_lh = region_lh(region_lh~=0);
    to = tabulate([region_lh]);
    if size(to, 1) < length(unique(label_l))-1
        all_classes = (1:length(unique(label_l))-1)';
        missing_classes = setdiff(all_classes, to(:, 1));
        to = [to; [missing_classes, zeros(length(missing_classes), 2)]];
        to = sortrows(to, 1);
    end
    ratio(:,species) = to(:,3)/100;
end
pit_corr = corr(ratio);
ratio(:,44) = 0.0001;

%% method 1
% 1. 强制对称化
pit_corr = (pit_corr + pit_corr') / 2;
% 2. 处理非负性（方法1：线性缩放）
pit_corr_scaled = (pit_corr + 1) / 2; % [-1,1] -> [0,1]
distance = 1 - pit_corr_scaled;       % 距离矩阵
% 3. 层次聚类
Z = linkage(squareform(distance), 'ward'); % 使用 Ward 方法
dendrogram(Z); % 绘制树状图

%% method 2
distance_matrix = 1-pit_corr;
distance_matrix = (distance_matrix + distance_matrix') / 2;  % 确保对称
Z = linkage(distance_matrix, 'average');
optimal_order = optimalleaforder(Z, distance_matrix);


save('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/kinship/Hierarchical_clustering_result.mat','optimal_order')
