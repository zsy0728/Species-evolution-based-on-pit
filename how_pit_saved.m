% step14_how_pit_saved.m
% no use 0324
clc
clear
warning off

atlasdir = '/mnt/sda/songyao/results/Evolution_cortical_shape/atlas/fsaverage6/';
surfl = vtkSurfRead([atlasdir,'sub-020_species-Homo+sapiens_hemi-L_topo-Homo.sapiens.surf_altas.vtk']);
surfr = vtkSurfRead([atlasdir,'sub-020_species-Homo+sapiens_hemi-R_topo-Homo.sapiens.surf_altas.vtk']);
yeo7_lh = surfl.Pdata{1,1}.val;
yeo7_rh = surfr.Pdata{1,1}.val;
sp = readtable('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/Species_info.xlsx');
na = readtable('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/node_ages.xlsx');

%% Part1: From the beginning to the end of evolution
% 演化终点 v.s. 演化起点
% no use 0324
load('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic_landmarks/pit_id_all.mat')
load('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/evolutionary_path_of_all_species.mat')

ancestor_pit_lh = yeo7_lh(pit_id_all{357});
ancestor_pit_rh = yeo7_rh(pit_id_all{358});
inter_lh = zeros(1,90);
inter_rh = zeros(1,90);

for path_num = 1:90 %surface number
    pit_lh = pit_id_all{path_num*2-1};
    pit_rh = pit_id_all{path_num*2};
    
    region_lh = yeo7_lh(pit_lh);
    region_lh = region_lh(region_lh~=0);
    region_rh = yeo7_rh(pit_rh);
    region_rh = region_rh(region_rh~=0);
    
    inter_lh(path_num) = length(intersect(region_lh,ancestor_pit_lh)); % save的pit数量
    inter_rh(path_num) = length(intersect(region_rh,ancestor_pit_rh));
end

% is saved pit number correlated with species type

Primata_rows = find(sp{:, 4} == "Primata");
Rodentia_rows = find(sp{:, 4} == "Rodentia");
other_rows = setdiff([1:90],[Primata_rows',Rodentia_rows']);
% LH
Primata_saved_lh = inter_lh(Primata_rows);
Rodentia_saved_lh = inter_lh(Rodentia_rows);
other_saved_lh = inter_lh(other_rows);

sum(Primata_saved_lh==2)/length(Primata_saved_lh)
sum(Rodentia_saved_lh==2)/length(Rodentia_saved_lh)
sum(other_saved_lh==2)/length(other_saved_lh)

% RH
Primata_saved_rh = inter_rh(Primata_rows);
Rodentia_saved_rh = inter_rh(Rodentia_rows);
other_saved_rh = inter_rh(other_rows);

sum(Primata_saved_rh==2)/length(Primata_saved_rh)
sum(Rodentia_saved_rh==2)/length(Rodentia_saved_rh)
sum(other_saved_rh==2)/length(other_saved_rh)

% WB 1
(sum(Primata_saved_rh==2)+sum(Primata_saved_lh==2))/(length(Primata_saved_rh)+length(Primata_saved_lh))
(sum(Rodentia_saved_rh==2)+sum(Rodentia_saved_lh==2))/(length(Rodentia_saved_rh)+length(Rodentia_saved_lh))
(sum(other_saved_rh==2)+sum(other_saved_lh==2))/(length(other_saved_rh)+length(other_saved_lh))

% WB2
length(find(inter_rh(Primata_rows)==2 & inter_lh(Primata_rows)==2))/58
length(find(inter_rh(Rodentia_rows)==2 & inter_lh(Rodentia_rows)==2))/28
length(find(inter_rh(other_rows)==2 & inter_lh(other_rows)==2))/4



%% Part2: from now to 78 mya ago
% 演化路径中所有的Pit和演化终点计算重合度
% no use 0324
clc
clear
warning off

atlasdir = '/mnt/sda/songyao/results/Evolution_cortical_shape/atlas/';
surfl = vtkSurfRead([atlasdir,'sub-020_species-Homo+sapiens_hemi-L_topo-Homo.sapiens.surf_altas.vtk']);
surfr = vtkSurfRead([atlasdir,'sub-020_species-Homo+sapiens_hemi-R_topo-Homo.sapiens.surf_altas.vtk']);
Parcels100_7Networks_lh = surfl.Pdata{1,6}.val;
Parcels100_7Networks_rh = surfr.Pdata{1,6}.val;
sp = readtable('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/Species_info.xlsx');
na = readtable('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/node_ages.xlsx');
load('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic_landmarks/pit_id_all_LH.mat')
load('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/evolutionary_path_of_all_species.mat')
ages = cellfun(@str2double, na{:,3});

figure('Position', [100, 100, 3000, 1500]);
subplot_idx = 1;
for path_num = 1:90
    path = flip(all_paths{path_num});
    age = ages(path);
    dice = [];
    end_pit_lh = Parcels100_7Networks_lh(pit_id_all{path_num*2-1});
    end_pit_rh = Parcels100_7Networks_rh(pit_id_all{path_num*2});
    end_pit_lh = end_pit_lh(end_pit_lh~=0);
    end_pit_rh = end_pit_rh(end_pit_rh~=0);
    
    for j = 1:length(all_paths{1,path_num})
        evolution_pit_lh = pit_id_all{all_paths{1,path_num}(j)*2-1};
        evolution_pit_rh = pit_id_all{all_paths{1,path_num}(j)*2};
        region_lh = Parcels100_7Networks_lh(evolution_pit_lh);
        region_lh = region_lh(region_lh~=0);
        region_rh = Parcels100_7Networks_rh(evolution_pit_rh);
        region_rh = region_rh(region_rh~=0);
        Cl = intersect(region_lh,end_pit_lh);
        indices_Al = find(ismember(region_lh, Cl));
        indices_Bl = find(ismember(end_pit_lh, Cl));
        Cr = intersect(region_rh,end_pit_rh);
        indices_Ar = find(ismember(region_rh, Cr));
        indices_Br = find(ismember(end_pit_rh, Cr));
        dice(length(all_paths{1,path_num})+1-j) = (length(indices_Al)+length(indices_Bl)+length(indices_Ar)+length(indices_Br))/(length(region_lh)+length(end_pit_lh)+length(region_rh)+length(end_pit_rh));
    end
    if sp.Order{path(end)} == "Primata"
        subplot(10, 9, subplot_idx);
        plot(age, dice, 'r*-');
        title(['Evo ' num2str(path_num)], 'Color', 'r');  % 红色标题
    else
        % 如果条件不满足，线条为蓝色，标题为蓝色
        subplot(10, 9, subplot_idx);
        plot(age, dice, 'b*-');  % 蓝色线条
        title(['Evo ' num2str(path_num)], 'Color', 'b');  % 蓝色标题
    end
    % 设置横坐标范围并反转
    xlim([0 80]);
    set(gca, 'XDir', 'reverse');  % 将X轴反转，确保从大到小
    xticks(sort(age, 'ascend'));  % 将xticks设置为升序的年龄值
    xticklabels(string(sort(round(age), 'ascend')));  % 标注每个年龄值
    subplot_idx = subplot_idx + 1;
end
% 全局设置X轴和Y轴标签
han = axes(gcf, 'visible', 'off');
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
xlabel(han, 'Time');
ylabel(han, 'Pit Dice with Terminal');
saveas(gcf, '/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/saved/proportion_coincides_with_endpoint.png');


%% Part3: The evolutionary order of different brain networks
% 计算不同脑网络中出现pit的比例占全脑的比值，反映脑功能网络的进化顺序
% no use 0324
clc
clear
warning off

atlasdir = '/mnt/sda/songyao/results/Evolution_cortical_shape/atlas/fsaverage6/';
surfl = vtkSurfRead([atlasdir,'sub-020_species-Homo+sapiens_hemi-L_topo-Homo.sapiens.surf_altas.vtk']);
surfr = vtkSurfRead([atlasdir,'sub-020_species-Homo+sapiens_hemi-R_topo-Homo.sapiens.surf_altas.vtk']);
yeo7_lh = surfl.Pdata{1,1}.val;
yeo7_rh = surfr.Pdata{1,1}.val;
sp = readtable('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/Species_info.xlsx');
na = readtable('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/node_ages.xlsx');
load('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic_landmarks/pit_id_all.mat')
load('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/evolutionary_path_of_all_species.mat')
ages = cellfun(@str2double, na{:,3});

for path_num = 1%:90
    path = flip(all_paths{path_num});
    age = round(ages(path));
    name = sp.Species{path_num};
    ratio = zeros(7,length(age));
    for j = 1:length(all_paths{1,path_num})
        evolution_pit_lh = pit_id_all{all_paths{1,path_num}(j)*2-1};
        evolution_pit_rh = pit_id_all{all_paths{1,path_num}(j)*2};
        region_lh = yeo7_lh(evolution_pit_lh);
        region_lh = region_lh(region_lh~=0);
        region_rh = yeo7_lh(evolution_pit_rh);
        region_rh = region_rh(region_rh~=0);
        to = tabulate([region_lh;region_rh]);
        ratio(:, length(all_paths{1,path_num})+1-j) = to(:,3)/100;
    end
    imagesc(ratio);
    colormap(flipud(cbrewer2('seq','YlGnBu')))
    c = colorbar;
    if sp.Order{path(end)} == "Primata"
        title_color = 'r';
    else
        title_color = 'b';
    end
    title(name, 'Color', title_color, ...
        'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold');
    
    xlabel('Evolution Time');
    %     ylabel('Yeo 7');
    set(gca, 'XTick', 1:length(all_paths{1, path_num}), 'XTickLabel', age(1:1:end));
    set(gca, 'YTick', 1:7, 'YTickLabel', {'Visual', 'SMN', 'DAN','VAN', 'Limbic', 'FPN', 'DMN'});
    axis equal;
    axis tight;
    saveas(gcf, ['/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/saved/proportion_changes_in_yeo7_Species',num2str(path_num),'_',name,'.png']);
end

%% Part4: The evolutionary order of different brain networks
% 计算所有演化路径中所有surface pit和相邻surface计算重合度
% no use 0324
clc
clear
warning off

atlasdir = '/mnt/sda/songyao/results/Evolution_cortical_shape/atlas/fsaverage6/';
surfl = vtkSurfRead([atlasdir,'sub-020_species-Homo+sapiens_hemi-L_topo-Homo.sapiens.surf_altas.vtk']);
surfr = vtkSurfRead([atlasdir,'sub-020_species-Homo+sapiens_hemi-R_topo-Homo.sapiens.surf_altas.vtk']);
% yeo7_lh = surfl.Pdata{1,1}.val;
% yeo7_rh = surfr.Pdata{1,1}.val;
Parcels100_7Networks_lh = surfl.Pdata{1,6}.val;
Parcels100_7Networks_rh = surfr.Pdata{1,6}.val;
sp = readtable('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/Species_info.xlsx');
na = readtable('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/node_ages.xlsx');
load('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic_landmarks/pit_id_all.mat')
load('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/evolutionary_path_of_all_species.mat')
ages = cellfun(@str2double, na{:,3});
figure('Position', [100, 100, 3000, 1500]);
subplot_idx = 1;
% p_all = zeros(1,90);
r_primata=[];
r_other = [];
for path_num = 1:90
    path = flip(all_paths{path_num});
%     age = round(ages(path));
    age = ages(path);
    name = sp.Species{path_num};
    dice = zeros(1,length(path)-1);
    pit_region_lh={};pit_region_rh={};
    for j = 1:length(path)
        evolution_pit_lh = pit_id_all{path(j)*2-1};
        evolution_pit_rh = pit_id_all{path(j)*2};
        region_lh = Parcels100_7Networks_lh(evolution_pit_lh);
        region_lh = region_lh(region_lh~=0);
        region_rh = Parcels100_7Networks_rh(evolution_pit_rh);
        region_rh = region_rh(region_rh~=0);
        pit_region_lh{j} = region_lh;
        pit_region_rh{j} = region_rh;
    end
    for j = 1:length(path)-1
        t1_lh = pit_region_lh{j};
        t2_lh = pit_region_lh{j+1};
        t1_rh = pit_region_rh{j};
        t2_rh = pit_region_rh{j+1};
        Cl = intersect(t1_lh,t2_lh);
        indices_Al = find(ismember(t1_lh, Cl));
        indices_Bl = find(ismember(t2_lh, Cl));
        Cr = intersect(t1_rh,t2_rh);
        indices_Ar = find(ismember(t1_rh, Cr));
        indices_Br = find(ismember(t2_rh, Cr));
        dice(j) = (length(indices_Al)+length(indices_Bl)+length(indices_Ar)+length(indices_Br))/(length(t1_lh)+length(t2_lh)+length(t1_rh)+length(t2_rh));
    end    
    
    age=age';
    mid_time_points = (age(1:end-1) + age(2:end)) / 2;
    subplot(10, 9, path_num);  
    scatter(mid_time_points, dice, 'filled', 'r');
    hold on;
    p = polyfit(mid_time_points, dice, 1);
    fitted_line = polyval(p, mid_time_points);
    plot(mid_time_points, fitted_line, 'b', 'LineWidth', 1.5);
    % 计算 R^2 值
    y_residuals = dice - fitted_line;  % 计算残差
    SS_residual = sum(y_residuals.^2);  % 残差平方和
    SS_total = sum((dice - mean(dice)).^2);  % 总平方和
    R2 = 1 - (SS_residual / SS_total);  % 拟合优度 R^2
%     % 计算F检验的p值
% n = length(mid_time_points);  % 数据点个数
% SS_regression = SS_total - SS_residual;  % 回归平方和
% F = (SS_regression / 1) / (SS_residual / (n - 2));  % F 统计量 (1 自由度)
% p_value = 1 - fcdf(F, 1, n-2);  % 根据 F 分布计算 p 值
% p_all(path_num) = p_value;
    grid on;
    % 设置横坐标范围并反转
    xlim([0 80]);
    set(gca, 'XDir', 'reverse');  % 将X轴反转，确保从大到小
    xticks(sort(age, 'ascend'));

%     if sp.Order{path(end)} == "Primata"
%         plot(mid_time_points, dice, 'r*');  % 红色线条
%         title(['Evo ' num2str(path_num),' : ',name,' | R^2 = ' num2str(R2, '%.2f')], 'Color', 'r');  % 红色标题 'Evo ' num2str(path_num)
%     else
%         plot(mid_time_points, dice, 'b*');  % 蓝色线条
%         title(['Evo ' num2str(path_num),' : ',name,' | R^2 = ' num2str(R2, '%.2f')], 'Color', 'b');  % 蓝色标题
%     end
    if sp.Order{path(end)} == "Primata"
        plot(mid_time_points, dice, 'r*');  % 红色线条
%         title([num2str(path_num),' : ',name,'|p=' num2str(p_value, '%.2f')], 'Color', 'r');
        title([num2str(path_num),' : ',name], 'Color', 'r');
%         p_primata = [p_primata,p_value];
        r_primata = [r_primata,R2];
    else
        plot(mid_time_points, dice, 'b*');  % 蓝色线条
%         title([num2str(path_num),' : ',name,'|p=' num2str(p_value, '%.2f')], 'Color', 'b');  % 蓝色标题
        title([num2str(path_num),' : ',name], 'Color', 'b');  % 蓝色标题
        
%         p_other = [p_other,p_value];
        r_other = [r_other,R2];
    end
    subplot_idx = subplot_idx + 1;    
end
han = axes(gcf, 'visible', 'off');
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
xlabel(han, 'Time');
ylabel(han, 'Pit Dice with adjacent surfaces');
saveas(gcf, '/mnt/sda/songyao/results/Evolution_cortical_shape/statistic/saved/proportion_coincides_with_adjacent_surf.png');

mean(r_primata)
mean(r_other)

