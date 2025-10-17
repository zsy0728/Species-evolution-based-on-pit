clc
clear
warning off

folder_path = '/mnt/sda/songyao/matlab_path/';
addpath(genpath(folder_path));

surfdir = '/mnt/sda/songyao/results/Evolution_cortical_shape/pits_v3/';
outdir = surfdir;

files = dir(fullfile(surfdir,'*_topo_depth.vtk'));
filenames = {files.name}';
ring_size = 4;
thr = 0.2;

for subid = 39% 1:length(filenames)
    surfname = [surfdir,filenames{subid}];
    surf1 = vtkSurfRead(surfname);
    surf = ReadSurf(surfname,[],0);
    depth = surf1.Pdata{1,2}.val;
    
    curvname = ['/mnt/sda/songyao/results/Evolution_cortical_shape/cortical_vtk/',filenames{subid}(1:end-9),'c.vtk'];
%     vtk2ply_ascii(surfname,'lh.white.vtk.ply');
%     system(['/mnt/sda/songyao/matlab_path/ComputeCurvatures ','lh.white.vtk.ply ',curvname]);
%     delete('lh.white.vtk.ply');
    surf_c = ReadSurf_2(curvname,{},1);
    surf.Neighbors = surf_c.Neighbors;
    
    %% pit
    pit_thr = thr*max(depth);
    % single pit
    local_maximum_vtxID_sulc=[];
    indx=1;
    
    for i = 1:size(surf.vertice,2)
        if depth(i)<pit_thr
            continue;
        end
        [Neighbor_ID,Ring_ID] = Search_Neighbor_ID(surf,i-1,ring_size,Inf);
        tmp_label = depth(Neighbor_ID+1);
        if max(tmp_label) == tmp_label(1)
            local_maximum_vtxID_sulc(indx) = i;
            indx = indx + 1;
        else
        end
    end
    
    % shaixuan 
    delete_id=[];
    for n = 1:size(local_maximum_vtxID_sulc,2)
         [Neighbor_ID,Ring_ID] = Search_Neighbor_ID(surf,local_maximum_vtxID_sulc(n)-1,4,Inf);
         if sum(depth(Neighbor_ID+1)<pit_thr)>=sum(depth(Neighbor_ID+1)>pit_thr)
             %shanchu
             delete_id = [delete_id,local_maximum_vtxID_sulc(n)];
         else
         end
    end
    local_maximum_vtxID_sulc_new = setdiff(local_maximum_vtxID_sulc,delete_id);
    vtx={};
    outfname = [outdir,filenames{subid}(1:end-9),'single_pit_new.vtk'];
    vtx.vertice = surf.vertice(:,local_maximum_vtxID_sulc_new);
    VertexWrite(outfname,vtx,{},'float');
%     pit_num{subid} = local_maximum_vtxID_sulc_new;
    
    
    
%     % pits
%     indx=1;
%     local_maximum_vtxID_ring=[];
%     ring_size2 = 2;
%     for i = 1:size(surf.vertice,2)
%         if depth(i)<pit_thr
%             continue;
%         end
%         [Neighbor_ID,Ring_ID] = Search_Neighbor_ID(surf_c,i-1,ring_size,Inf);
%         tmp_label = depth(Neighbor_ID+1);
%         if max(tmp_label) == tmp_label(1)
%             [Max_Neighbor_ID,Max_Ring_ID] = Search_Neighbor_ID(surf_c,i-1,ring_size2,Inf);
%             local_maximum_vtxID_ring = [local_maximum_vtxID_ring,(Max_Neighbor_ID+1)];
%             indx = indx + 1;
%         else
%         end
%     end
%     vtx={};
%     outfname = [outdir,filenames{subid}(1:end-14),'pits.vtk'];
%     vtx.vertice = surf.vertice(:,local_maximum_vtxID_ring);
%     VertexWrite(outfname,vtx,{},'float');
    
%     %% peak
%     peak_thr = 0.9*max(depth);
%     % single peak
%     local_maximum_vtxID_sulc=[];
%     indx=1;
%     for i = 1:size(surf.vertice,2)
%         if depth(i)>peak_thr
%             continue;
%         end
%         [Neighbor_ID,Ring_ID] = Search_Neighbor_ID(surf,i-1,ring_size,Inf);
%         tmp_label = depth(Neighbor_ID+1);
%         if min(tmp_label) == tmp_label(1)
%             local_maximum_vtxID_sulc(indx) = i;
%             indx = indx + 1;
%         else
%         end
%     end
%     vtx={};
%     outfname = [outdir,filenames{subid}(1:end-9),'single_peak.vtk'];
%     vtx.vertice = surf.vertice(:,local_maximum_vtxID_sulc);
%     VertexWrite(outfname,vtx,{},'float');
%     
%     % peaks
%     indx=1;
%     local_maximum_vtxID_ring=[];
%     ring_size2 = 2;
%     for i = 1:size(surf.vertice,2)
%         if depth(i)>pit_thr
%             continue;
%         end
%         [Neighbor_ID,Ring_ID] = Search_Neighbor_ID(surf_c,i-1,ring_size,Inf);
%         tmp_label = depth(Neighbor_ID+1);
%         if max(tmp_label) == tmp_label(1)
%             [Max_Neighbor_ID,Max_Ring_ID] = Search_Neighbor_ID(surf_c,i-1,ring_size2,Inf);
%             local_maximum_vtxID_ring = [local_maximum_vtxID_ring,(Max_Neighbor_ID+1)];
%             indx = indx + 1;
%         else
%         end
%     end
%     vtx={};
%     outfname = [outdir,filenames{subid}(1:end-9),'peaks.vtk'];
%     vtx.vertice = surf.vertice(:,local_maximum_vtxID_ring);
%     VertexWrite(outfname,vtx,{},'float');
end

% save('/mnt/sda/songyao/results/Evolution_cortical_shape/statistic_landmarks/pit_num_all_species_new.mat','pit_num')
