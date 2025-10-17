clc
clear
warning off

folder = '/media/songyao/songyao/matlab_path/';
addpath(genpath(folder))
folder_path = '/media/songyao/songyao/send_2_tuo/cross90/orig_files/_surfaces';
gii_files = dir(fullfile(folder_path, '*.gii'));
% 循环处理每个 .gii 文件
for i = 1:numel(gii_files)
    % 读取当前 .gii 文件
    gii_file_path = fullfile(folder_path, gii_files(i).name);
    surf = gifti(gii_file_path);
    % 获取文件名（不含扩展名）
    [~, base_name, ~] = fileparts(gii_files(i).name);
    output_surf.Vtx = surf.vertices';
    output_surf.Face = surf.faces'-1;
    vtkSurfWrite(['/media/songyao/songyao/send_2_tuo/cross90/orig_files/surface_vtk/',base_name, '.vtk'],output_surf); 
end


