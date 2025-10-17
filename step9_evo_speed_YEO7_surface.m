% Comparison of the development of 7 brain networks between different species
clc
clear
warning off

folder = '/mnt/sda/songyao/matlab_path/';
addpath(genpath(folder))

atlasdir = '/mnt/sda/songyao/results/Evolution_cortical_shape/atlas/';
surfl = vtkSurfRead([atlasdir,'sub-020_species-Homo+sapiens_hemi-L_topo-Homo.sapiens.surf_altas.vtk']);
surfr = vtkSurfRead([atlasdir,'sub-020_species-Homo+sapiens_hemi-R_topo-Homo.sapiens.surf_altas.vtk']);
yeo7_lh = surfl.Pdata{1,1}.val;
yeo7_rh = surfr.Pdata{1,1}.val;
load('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/evolutionary_path_of_all_species.mat')
load('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic_landmarks/pit_id_all_LH.mat')
sp = readtable('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/Species_info.xlsx');

pit_growth_ratio = zeros(7,90);
for sub = 1:90
    pits_id = pit_id_all{sub};
    region_lh = yeo7_lh(pits_id);
    region_lh = region_lh(region_lh~=0);
    
    to = tabulate(region_lh);
    to_complete = [(1:7)', zeros(7, 2)];
    to_complete(ismember(to_complete(:, 1), to(:, 1)), 2:3) = to(:, 2:3);
    number = to_complete(:,2);
    pit_growth_ratio(:,sub) = number - [1;1;1;0;0;0;0];
end

% area for yeo7
s = zeros(1,7);
for j=1:7
    s(j) = length(find(yeo7_lh == j));
end
pit_growth_ratio_norm = pit_growth_ratio./s';

save('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/evo_speedy/all_pit_evo_speed_7x90.mat','pit_growth_ratio')
save('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/evo_speedy/all_pit_evo_speed_7x90_normSA.mat','pit_growth_ratio_norm')


%% visualize atlas on fsaverage6
surfl_fsaverage6 = vtkSurfRead(['/mnt/sda/songyao/results/Evolution_cortical_shape/show_pits/homo_sulcal_depth_fsaverage6_lh.vtk']);
surfr_fsaverage6 = vtkSurfRead(['/mnt/sda/songyao/results/Evolution_cortical_shape/show_pits/homo_sulcal_depth_fsaverage6_rh.vtk']);
surfl_fsaverage6.Pdata{1,32} = surfl_fsaverage6.Pdata{1,4};
surfr_fsaverage6.Pdata{1,32} = surfr_fsaverage6.Pdata{1,4};
for n=1:31
    surfl_fsaverage6.Pdata{1,n} = surfl.Pdata{1,n};
    surfr_fsaverage6.Pdata{1,n} = surfr.Pdata{1,n};
end
surfl_fsaverage6.Face = surfl_fsaverage6.Face-1;
surfr_fsaverage6.Face = surfr_fsaverage6.Face-1;

vtkSurfWrite([atlasdir,'homo_atlas_fsaverage6_lh.vtk'],surfl_fsaverage6)
vtkSurfWrite([atlasdir,'homo_atlas_fsaverage6_rh.vtk'],surfr_fsaverage6)

%%
pri_id = find(strcmp(sp.Order, 'Primata'));
nonpri_id = find(~strcmp(sp.Order, 'Primata'));
save('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/id_primate.mat','pri_id')
save('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/id_non-primate.mat','nonpri_id')

%% 统计每个脑网络的进化程度
pri = mean(pit_growth_ratio(:,pri_id),2);
non_pri = mean(pit_growth_ratio(:,nonpri_id),2);
all = mean(pit_growth_ratio,2);
save('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/evo_speedy/pri_pit_evo_speed.mat','pri')
save('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/evo_speedy/nonpri_pit_evo_speed.mat','non_pri')
save('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/evo_speedy/all_pit_evo_speed.mat','all')

pri_ratio = mean(pit_growth_ratio_norm(:,pri_id),2);
non_pri_ratio = mean(pit_growth_ratio_norm(:,nonpri_id),2);
all_ratio = mean(pit_growth_ratio_norm,2);
save('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/evo_speedy/pri_pit_evo_speed_normSA.mat','pri_ratio')
save('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/evo_speedy/nonpri_pit_evo_speed_normSA.mat','non_pri_ratio')
save('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/evo_speedy/all_pit_evo_speed_normSA.mat','all_ratio')

%% 统计每个脑网络作为最快进化的次数
expansion_speed_lh_pri = zeros(1,40962);
expansion_speed_lh_nonpri = zeros(1,40962);
expansion_speed_lh_all = zeros(1,40962);
expansion_speed_lh_pri_ratio = zeros(1,40962);
expansion_speed_lh_nonpri_ratio = zeros(1,40962);
expansion_speed_lh_all_ratio = zeros(1,40962);
for i = 1:7
    v = find(yeo7_lh == i);
    expansion_speed_lh_pri(v) = pri(i);
    expansion_speed_lh_nonpri(v) = non_pri(i);
    expansion_speed_lh_all(v) = all(i);
    expansion_speed_lh_pri_ratio(v) = pri_ratio(i);
    expansion_speed_lh_nonpri_ratio(v) = non_pri_ratio(i);
    expansion_speed_lh_all_ratio(v) = all_ratio(i);
end
surfl_fsaverage6.Pdata = [];
surfl_fsaverage6.Pdata{1,1}.val = expansion_speed_lh_pri;
surfl_fsaverage6.Pdata{1,1}.name = 'expansion_speed_lh_pri';
surfl_fsaverage6.Pdata{1,2}.val = expansion_speed_lh_nonpri;
surfl_fsaverage6.Pdata{1,2}.name = 'expansion_speed_lh_nonpri';
surfl_fsaverage6.Pdata{1,3}.val = expansion_speed_lh_all;
surfl_fsaverage6.Pdata{1,3}.name = 'expansion_speed_lh_all';
surfl_fsaverage6.Pdata{1,4}.val = expansion_speed_lh_pri_ratio;
surfl_fsaverage6.Pdata{1,4}.name = 'expansion_speed_lh_pri_ratio';
surfl_fsaverage6.Pdata{1,5}.val = expansion_speed_lh_nonpri_ratio;
surfl_fsaverage6.Pdata{1,5}.name = 'expansion_speed_lh_nonpri_ratio';
surfl_fsaverage6.Pdata{1,6}.val = expansion_speed_lh_all_ratio*100000;
surfl_fsaverage6.Pdata{1,6}.name = 'expansion_speed_lh_all_ratio';
surfl_fsaverage6.Face = surfl_fsaverage6.Face-1;
vtkSurfWrite(['/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/evo_speedy/homo_expansion_speed_lh.vtk'],surfl_fsaverage6)

% 遍历每一列，计算最大值及其索引
% ratio
maxIndices_all = nan(1, 90);
for col = 1:90
    colData = pit_growth_ratio_norm(:, col);
    maxVal = max(colData); % 找到当前列的最大值
    maxIndices = find(colData == maxVal); % 找到最大值的所有索引
    if numel(maxIndices) == 1
        maxIndices_all(col) = maxIndices; % 如果最大值唯一，则记录最大值
    end
end
tabulate(maxIndices_all)

% count
maxIndices_all = nan(1, 90);
for col = 1:90
    colData = pit_growth_ratio(:, col);
    maxVal = max(colData); % 找到当前列的最大值
    maxIndices = find(colData == maxVal); % 找到最大值的所有索引
    if numel(maxIndices) == 1
        maxIndices_all(col) = maxIndices; % 如果最大值唯一，则记录最大值
    end
end
tabulate(maxIndices_all)

% count: DMN>Limbic>VAN>Visual>SMN>DAN=FPN 
% ratio: VAN>DMN>Limbic>Visual>FPN>DAN>SMN


tabulate(maxIndices_all(pri_id))
tabulate(maxIndices_all(nonpri_id))

pri_yeo7 = maxIndices_all(pri_id);
nonpri_yeo7 = maxIndices_all(nonpri_id);
save('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/evo_speedy/pri_yeo7.mat','pri_yeo7')
save('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/evo_speedy/nonpri_yeo7.mat','nonpri_yeo7')

% a = find(isnan(maxIndices_all)~=1);
% indices = find(strcmp(sp.Order, 'Primata'));
% intersect(a,indices)



