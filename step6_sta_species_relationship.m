clc
clear
warning off

surfdir = '/mnt/sda/songyao/results/Evolution_cortical_shape/surface_vtk/';
outdir = '/mnt/sda/songyao/results/Evolution_cortical_shape/data_info/';

files = dir(fullfile(surfdir,'*L_topo-Homo.sapiens.surf.vtk'));
filenames = {files.name}';

for subid = 91:length(filenames)
    surfname = filenames{subid};
    Branch1_name = surfname(13:15);
    if Branch1_name(1) == '0'
        Branch1_name = Branch1_name(2:end);
    end
    Branch2_name = surfname(17:19);
    if Branch2_name(1) == '0'
        Branch2_name = Branch2_name(2:end);
    end 
    Root_name = surfname(5:7);
    if Root_name(1) == '0'
        Root_name = Root_name(2:end);
    end
    relation_table{subid-90,1} = Root_name;
    relation_table{subid-90,2} = Branch1_name;
    relation_table{subid-90,3} = Branch2_name;
end
save([outdir,'phyTree_table.mat'],'relation_table')

