%% 
clc
clear
warning off

load('/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/evolutionary_path_of_all_species.mat')
load('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic_landmarks/number_of_pits_species_LH.mat')

% 创建一个新图形窗口并设置大小
figure('Units', 'inches', 'Position', [1, 1, 36, 20]); % 设置宽12英寸，高10英寸
for j = 1:90
    path = all_paths{j}; 
    pits_num = pit_num_LH(path); % 假设该函数返回与 path 相关的 pit 数量

    % 创建一个新的子图
    subplot(10, 9, j); % 10 行 9 列的子图布局

    % 绘制散点图
    scatter(1:length(pits_num), pits_num, 'filled', 'MarkerEdgeColor', 'k');
    hold on; % 保持当前图形以绘制折线图

    % 绘制折线图
    plot(1:length(pits_num), pits_num, 'LineWidth', 1.5);

    % 设置标题和轴标签
    title(['Path ' num2str(j)]);
    xlabel('Index');
    ylabel('Pit Count');

    % 设置 Y 轴的范围
    ylim([0 max(pits_num) * 1.1]); % 适当调整 Y 轴范围
    grid on; % 添加网格线
end
% 调整整个图形的布局
sgtitle('Pit Count for 90 Paths'); % 添加总标题
saveas(gcf, '/mnt/sda/songyao/results/Evolution_cortical_shape/figure/pit_num_evo_lh.png'); 

