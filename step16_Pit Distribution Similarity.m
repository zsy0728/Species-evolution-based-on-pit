% pit similarity by N-ring
%% Species_kinship v.s. pit density in YEO7/Schaefer100
clc
clear
% atlasdir = '/mnt/sda/songyao/results/Evolution_cortical_shape/atlas/';
% surfl = vtkSurfRead([atlasdir,'sub-020_species-Homo+sapiens_hemi-L_topo-Homo.sapiens.surf_altas.vtk']);
% surfr = vtkSurfRead([atlasdir,'sub-020_species-Homo+sapiens_hemi-R_topo-Homo.sapiens.surf_altas.vtk']);
surface_dir = '/mnt/sda/songyao/results/Evolution_cortical_shape_v1/surface_vtk/';
homol = ReadSurf_2([surface_dir,'sub-020_species-Homo+sapiens_hemi-L_topo-Homo.sapiens_c.vtk'],{},1);
load('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic_landmarks/pit_id_all_LH.mat')
load('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/evolutionary_path_of_all_species.mat')
% load('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/kinship/Species_kinship_by_phyTree.mat')
% load('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/kinship/DIST_normalized.mat')
% sp = readtable('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/Species_info.xlsx');
% % na = readtable('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/node_ages.xlsx');
% % ages = cellfun(@str2double, na{:,3});
max_search_ring = 10; % 可以根据实际情况调整
distance_ring = zeros(90, 90);
R = max_search_ring; % 距离阈值

for species_1 = 1:90
    species_1
    pit_lh_1 = pit_id_all{species_1};
    % 从 species_1 + 1 开始，只计算上三角部分
    for species_2 = species_1 + 1:90 
        pit_lh_2 = pit_id_all{species_2};
        % 计算两个物种所有Pit之间的距离
        distance_matrix = zeros(length(pit_lh_1), length(pit_lh_2));
        for i = 1:length(pit_lh_1)
            v1 = pit_lh_1(i);
            for j = 1:length(pit_lh_2)
                v2 = pit_lh_2(j);
                
                ring_size = 1; % 定义 ring_size
                found = false;
                ring_distance = 0;
                for current_ring = 1:max_search_ring
                    [Neighbor_ID,Ring_ID] = Search_Neighbor_ID(homol,v1-1,current_ring,Inf);
                    % 检查 vertex2 是否在当前圈的邻居中
                    if any(Neighbor_ID == v2-1)
                        ring_distance = current_ring;
                        found = true;
                        break;
                    end
                end
                % 如果未找到连接，将距离设为 Inf
                if ~found
                    ring_distance = Inf;
                end
                % 将计算得到的距离存入矩阵
                distance_matrix(i, j) = ring_distance;
            end
        end
        %% 
        row_min = min(distance_matrix, [], 2);
        col_min = min(distance_matrix, [], 1);
        
        similar = (sum(row_min<Inf)/length(row_min)+sum(col_min<Inf)/length(col_min))/2;
         % 存储 dice 相似度
        distance_ring(species_1, species_2) = similar;
        distance_ring(species_2, species_1) = similar;
        
    end
    % 对角线元素设为 1，表示自身与自身的相似度为 1
    distance_ring(species_1, species_1) = 1; 
end
load('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/phytree_indices.mat')
pit_corr_order = distance_ring(phytree_indices, phytree_indices);
heatmap(pit_corr_order)
save('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/kinship/pit_corr_Nring.mat','distance_ring')
save('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/kinship/pit_corr_Nring_order.mat','pit_corr_order')

load('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/kinship/pit_corr_Nring_order.mat')
load('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/kinship/DIST_normalized.mat')
% Level 1
A = DIST_normalized;
B = pit_corr_order;
% Level 2
A = DIST_normalized(1:58,1:58);
B = pit_corr_order(1:58,1:58);

A = DIST_normalized(61:86,61:86);
B = pit_corr_order(61:86,61:86);

% 余弦相似度（Cosine Similarity）
vec1 = A(:);
vec2 = B(:);
similarity = (vec1' * vec2) / (norm(vec1) * norm(vec2));
disp(['Cosine Similarity: ', num2str(similarity)]);
% 相关系数（Pearson Correlation Coefficient）
cc = corr2(A,B);
disp(['PCC: ', num2str(cc)]);
