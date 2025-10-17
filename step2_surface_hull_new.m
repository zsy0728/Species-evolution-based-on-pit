clc
clear
warning off

folder = '/media/songyao/songyao/matlab_path/';
addpath(genpath(folder))

folder_path = '/mnt/sda/songyao/results/Evolution_cortical_shape/cortical_vtk/';
surf_files = dir(fullfile(folder_path, '*_cerebral.vtk'));
outdir = '/mnt/sda/songyao/results/Evolution_cortical_shape/hull_vtk_new/';
for subid = 71%1:numel(surf_files)
    %%
    surfname = fullfile(folder_path, surf_files(subid).name);
    surf = vtkSurfRead(surfname);
    [~, base_name, ~] = fileparts(surf_files(subid).name);
    % generate hull surface method1
    vertex = surf.Vtx;
    [K,v] = convhulln(vertex');

    used_vertices = unique(K(:));
    % 删除未使用的顶点
    new_vertex = vertex(:,used_vertices);
    % 更新面片列表
    [~, new_K] = ismember(K, used_vertices);
    [new_K, new_v] = convhulln(new_vertex');

    hull_surf.Vtx = new_vertex;
    hull_surf.Face = new_K'-1;
    vtkSurfWrite([outdir,base_name(1:end-8), 'hull.vtk'],hull_surf);
end

%% xi fen
clc
clear

folder_path = '/mnt/sda/songyao/results/Evolution_cortical_shape/cortical_vtk/';
surf_files = dir(fullfile(folder_path, '*_cerebral.vtk'));
outdir = '/mnt/sda/songyao/results/Evolution_cortical_shape/hull_vtk_new/';

area_thr = 0.2;
for subid = 71%1:numel(surf_files)    
    %%
    surfname = fullfile(folder_path, surf_files(subid).name);
    [~, base_name, ~] = fileparts(surf_files(subid).name);
%     if exist([outdir,base_name(1:end-4), 'hull.upsample.vtk'],'file')
%         continue;
%     end
    base_name(1:7)
    surf = vtkSurfRead([outdir,base_name(1:end-8), 'hull.vtk']);
    vertex = surf.Vtx';
    triangles = surf.Face';
    
    % 初始化标志位为 false
    has_large_area = true;
    round = 0;
    new_vertex = vertex;
    
    while has_large_area
        disp(['Round ', num2str(round)]);
        has_large_area = false; % 重置标志位
        % 创建一个空的顶点列表和三角形列表
        new_triangles = [];
        round = round+1;
        % 对每个三角形进行细分
        for i = 1:size(triangles, 1)
            % 获取当前三角形的顶点索引
            tri = triangles(i, :);
            % 计算当前三角形的三个顶点的坐标
            v1 = vertex(tri(1), :);
            v2 = vertex(tri(2), :);
            v3 = vertex(tri(3), :);
            % 计算三角形的面积
            area = 0.5 * norm(cross(v2 - v1, v3 - v1));
            % 如果三角形的面积大于1，进行细分
            if area > area_thr
                len1 = norm(v2 - v1);
                len2 = norm(v3 - v2);
                len3 = norm(v1 - v3);
                alpha = acos((len2^2 + len3^2 - len1^2) / (2 * len2 * len3));
                beta = acos((len1^2 + len3^2 - len2^2) / (2 * len1 * len3));
                gamma = acos((len1^2 + len2^2 - len3^2) / (2 * len1 * len2));
                angle1_deg = rad2deg(alpha);
                angle2_deg = rad2deg(beta);
                angle3_deg = rad2deg(gamma);
                if angle1_deg>90||angle2_deg>90||angle3_deg>90
                    % 找出最长边的长度和对应的中点
                    [max_len, idx] = max([len1, len2, len3]);
                    % 将最长边的中点作为新的顶点
                    midpoints = (vertex(tri(idx), :) + vertex(tri(mod(idx, 3)+1), :)) / 2;
                else
                     % 计算三个顶点的中点
                    midpoints = (v1 + v2 + v3) / 3;
                end
                % 将中点添加到顶点列表中
                new_vertex = [new_vertex; midpoints];
                % 获取中点的索引
                midpoint_idx = size(new_vertex, 1);
                % 创建新的三角形，将当前三角形分割为三个子三角形
                new_triangles = [new_triangles; tri(1)-1, tri(2)-1, midpoint_idx-1];
                new_triangles = [new_triangles; tri(2)-1, tri(3)-1, midpoint_idx-1];
                new_triangles = [new_triangles; tri(3)-1, tri(1)-1, midpoint_idx-1];
            else
                % 如果三角形的面积不大于2，保持原样
                new_triangles = [new_triangles; tri-1];
            end
        end
        %% 判断三角形面积和角度
        area_all = zeros(1,size(new_triangles, 1));
        for i = 1:size(new_triangles, 1)
            % 获取当前三角形的顶点索引
            tri = new_triangles(i, :)+1;
            % 计算当前三角形的三个顶点的坐标
            v1 = new_vertex(tri(1), :);
            v2 = new_vertex(tri(2), :);
            v3 = new_vertex(tri(3), :);
            % 计算三角形的面积
            area = 0.5 * norm(cross(v2 - v1, v3 - v1));
            area_all(i) = area;
        end
        if length(find(area_all>area_thr))>1
            has_large_area = true;
            triangles = new_triangles+1;
            vertex = new_vertex;
        else
            has_large_area = false;
        end
    end


    hull_surf.Vtx = new_vertex';
    hull_surf.Face = new_triangles';
    vtkSurfWrite([outdir,base_name(1:end-4), 'hull.upsample.vtk'],hull_surf);
end



