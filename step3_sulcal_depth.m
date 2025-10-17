% 去掉嗅球部分的surface上面的depth
clc
clear
warning off

folder_path = '/mnt/sda/songyao/matlab_path/';
addpath(genpath(folder_path));

surfdir = '/mnt/sda/songyao/results/Evolution_cortical_shape/cortical_vtk/';
hulldir = '/mnt/sda/songyao/results/Evolution_cortical_shape/hull_vtk_new/';
outdir = '/mnt/sda/songyao/results/Evolution_cortical_shape/pits_v3/';
files = dir(fullfile(surfdir,'*cerebral.vtk'));
filenames = {files.name}';
for subid = 71%:length(filenames)
    surfname = [surfdir,filenames{subid}];
    surf = vtkSurfRead(surfname);
    hull = vtkSurfRead([hulldir,filenames{subid}(1:end-8),'hull.upsample.vtk']);
    depth = zeros(1,size(surf.Vtx,2));
    for i = 1:size(surf.Vtx,2)
        surf_coord = surf.Vtx(:,i);
        distance = dist(surf_coord,hull.Vtx);
        [mind,minx] = min(distance);
        depth(i) = mind;
    end
    surf.Pdata{1,1}.val = depth;  
    surf.Pdata{1,1}.name = 'depth';
    surf.Face = surf.Face-1;
    vtkSurfWrite([outdir,filenames{subid}(1:end-12),'depth.vtk'],surf)
end

%% smooth
clc
clear
warning off

folder_path = '/mnt/sda/songyao/matlab_path/';
addpath(genpath(folder_path));

surfdir = '/mnt/sda/songyao/results/Evolution_cortical_shape/pits_v3/';
outdir = surfdir;

files = dir(fullfile(surfdir,'*_topo_depth.vtk'));
filenames = {files.name}';

com = table2array(Complexity);
for subid = 71%1:length(filenames)
    surfname = [surfdir,filenames{subid}];
%     surfname
    surf1 = vtkSurfRead(surfname);    
    depth = surf1.Pdata{1,1}.val;
    curvname = ['/mnt/sda/songyao/results/Evolution_cortical_shape/cortical_vtk/',filenames{subid}(1:end-9),'c.vtk'];
    
    surf_c = ReadSurf_2(curvname,{},1);
    if com(round(subid/2))==1
        depth_new = smooth_surf(surf_c,depth,5);
    else
        depth_new = smooth_surf(surf_c,depth,1);
    end
    surf1.Pdata{1,2}.val = depth_new;
    surf1.Pdata{1,2}.name = 'depth_smooth';
    surf1.Face = surf1.Face -1;
    vtkSurfWrite([surfdir,filenames{subid}(1:end-4),'.vtk'],surf1);
end









