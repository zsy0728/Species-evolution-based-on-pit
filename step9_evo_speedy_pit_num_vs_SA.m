% Compare surface area expansion speedy with sulcal pits number speedy
clc
clear
warning off

folder = '/mnt/sda/songyao/matlab_path/';
addpath(genpath(folder))

atlasdir = '/mnt/sda/songyao/results/Evolution_cortical_shape/atlas/';
temp_surfl = vtkSurfRead([atlasdir,'sub-020_species-Homo+sapiens_hemi-L_topo-Homo.sapiens.surf_altas.vtk']);
temp_surfr = vtkSurfRead([atlasdir,'sub-020_species-Homo+sapiens_hemi-R_topo-Homo.sapiens.surf_altas.vtk']);
yeo7_lh = temp_surfl.Pdata{1,1}.val;
yeo7_rh = temp_surfr.Pdata{1,1}.val;
sp = readtable('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/Species_info.xlsx');

load('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/id_primate.mat')
load('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/id_non-primate.mat')

surfdir = '/mnt/sda/songyao/results/Evolution_cortical_shape/surface_vtk/';
files = dir(fullfile(surfdir,'*-L_topo-Homo.sapiens.surf.vtk'));
filenames = {files.name}';

num_surfaces = 90;               % 表面数量
num_classes = 7;                  % 分类数量
area_matrix = zeros(num_classes, num_surfaces);  % 存储面积的矩阵 (7x90)

for subid = 1:90
    surfnamel = [surfdir,filenames{subid}];
    surfl = vtkSurfRead(surfnamel);

    Vl = surfl.Vtx';
    Fl = surfl.Face';
    v1 = Vl(Fl(:, 1), :);  % Coordinates of the first vertex in each triangle
    v2 = Vl(Fl(:, 2), :);  % Coordinates of the second vertex in each triangle
    v3 = Vl(Fl(:, 3), :);  % Coordinates of the third vertex in each triangle
    edge1 = v2 - v1;
    edge2 = v3 - v1;
    crossProd = cross(edge1, edge2, 2);  % Cross product for each triangle
    triangleArea = 0.5 * sqrt(sum(crossProd.^2, 2));  % Magnitude of the cross product
    area_sum = zeros(1, num_classes);   % 存储每类的表面积和

    for i = 1:size(Fl, 1)
        % 获取面片的顶点标签
        face_labels = yeo7_lh(Fl(i, :));
        most_common_label = mode(face_labels);
        if most_common_label ~= 0
            area_sum(most_common_label) = area_sum(most_common_label) + triangleArea(i);
        end
    end
    area_matrix(:, subid) = area_sum';
    fprintf('表面 %d 完成！\n', subid);  % 进度提示
end

total_area_per_surface = sum(area_matrix, 1);   % 1x90，每个表面的总面积
% 按列归一化
normalized_area_yeo7 = area_matrix ./ total_area_per_surface;
save('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/cortical_features/surface_area_yeo7.mat','normalized_area_yeo7')

% ancestor surface area
subid = 179;
ancestor_SA = zeros(1,7);
surfnamel = [surfdir,filenames{subid}];
surfl = vtkSurfRead(surfnamel);
Vl = surfl.Vtx';
Fl = surfl.Face';
v1 = Vl(Fl(:, 1), :);  % Coordinates of the first vertex in each triangle
v2 = Vl(Fl(:, 2), :);  % Coordinates of the second vertex in each triangle
v3 = Vl(Fl(:, 3), :);  % Coordinates of the third vertex in each triangle
edge1 = v2 - v1;
edge2 = v3 - v1;
crossProd = cross(edge1, edge2, 2);  % Cross product for each triangle
triangleArea = 0.5 * sqrt(sum(crossProd.^2, 2));  % Magnitude of the cross product

for i = 1:size(Fl, 1)
    % 获取面片的顶点标签
    face_labels = yeo7_lh(Fl(i, :));
    most_common_label = mode(face_labels);
    if most_common_label ~= 0
        ancestor_SA(most_common_label) = ancestor_SA(most_common_label) + triangleArea(i);
    end
end

area_expansion = area_matrix-ancestor_SA';
save('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/cortical_features/surface_area_expansion_yeo7.mat','area_expansion')

%% sa vs pit_num speed
clc
clear

load('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/cortical_features/surface_area_expansion_yeo7.mat')
load('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/evo_speedy/all_pit_evo_speed_7x90.mat')

A = area_expansion;
B = pit_growth_ratio;

%% all species
corr_coeff = corr(A(:), B(:));
fprintf('皮尔逊相关系数：%.4f\n', corr_coeff);
ssim_value = ssim(A, B);
fprintf('结构相似性 (SSIM)：%.4f\n', ssim_value);
A_vec = A(:);
B_vec = B(:);
cosine_sim = dot(A_vec, B_vec) / (norm(A_vec) * norm(B_vec));
fprintf('余弦相似度：%.4f\n', cosine_sim);

% in each region
for roi = 1:7
    corr_roi(roi) = corr(A(roi, :)', B(roi, :)');
end
disp('按 ROI 计算的皮尔逊相关系数:');
disp(corr_roi);

%% pri and non-pri
load('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/id_primate.mat')
load('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/id_non-primate.mat')
A_pri = A(:,pri_id);
A_non = A(:,nonpri_id);
B_pri = B(:,pri_id);
B_non = B(:,nonpri_id);

pcc_pri = zeros(7, 1);    
pcc_non = zeros(7, 1);    
for i = 1:7
    pcc_pri(i) = corr(A_pri(i, :)', B_pri(i, :)', 'type', 'Pearson');
    pcc_non(i) = corr(A_non(i, :)', B_non(i, :)', 'type', 'Pearson');
end

% 显示结果
disp('Pri 类别 7 个脑区的 PCC：');
disp(pcc_pri);
disp('Nonpri 类别 7 个脑区的 PCC：');
disp(pcc_non);

%% visualize results on fsaverage6
surfl_fsaverage6 = vtkSurfRead('/mnt/sda/songyao/results/Evolution_cortical_shape/show_pits/homo_sulcal_depth_fsaverage6_lh.vtk');
atlasdir = '/mnt/sda/songyao/results/Evolution_cortical_shape/atlas/';
surfl = vtkSurfRead([atlasdir,'sub-020_species-Homo+sapiens_hemi-L_topo-Homo.sapiens.surf_altas.vtk']);
yeo7_lh = surfl.Pdata{1,1}.val;

pitnum_coupled_sa_lh = zeros(1,40962);
pitnum_coupled_sa_pri_lh = zeros(1,40962);
pitnum_coupled_sa_non_lh = zeros(1,40962);
for i = 1:7
    v = find(yeo7_lh == i);
    pitnum_coupled_sa_lh(v) = corr_roi(i);
    pitnum_coupled_sa_pri_lh(v) = pcc_pri(i);
    pitnum_coupled_sa_non_lh(v) = pcc_non(i);
end
surfl_fsaverage6.Pdata = [];
surfl_fsaverage6.Pdata{1,1}.val = pitnum_coupled_sa_lh;
surfl_fsaverage6.Pdata{1,1}.name = 'pitnum_coupled_sa_lh';
surfl_fsaverage6.Pdata{1,2}.val = pitnum_coupled_sa_pri_lh;
surfl_fsaverage6.Pdata{1,2}.name = 'pitnum_coupled_sa_pri_lh';
surfl_fsaverage6.Pdata{1,3}.val = pitnum_coupled_sa_non_lh;
surfl_fsaverage6.Pdata{1,3}.name = 'pitnum_coupled_sa_non_lh';
surfl_fsaverage6.Face = surfl_fsaverage6.Face-1;
vtkSurfWrite(['/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/evo_speedy/pitnum_growth_coupled_SA_expansion_lh.vtk'],surfl_fsaverage6)

save('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/evo_speedy/pitnum_growth_coupled_SA_expansion_primate_lh.mat','pcc_pri')
save('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/evo_speedy/pitnum_growth_coupled_SA_expansion_nonprimate_lh.mat','pcc_non')







