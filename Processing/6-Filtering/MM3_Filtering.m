clear
close all

load_sub_neuron_data = false;

use_log_bounding = true;

%sort_order = [2,8,6,1,3,4,7,5];

sort_order = [2,8,7,1,3,4,6,5];

% color_codes = ["15A3C7","FCDE9C","045275","089099","F0746E","7CCBA2","7C1D6F","DC3977"]';

% color_codes = ["13C4DC","5F5F5F","0855B1","BBE4EB","88C4D8","18D2C4","ED7014","F7D42F"]';

color_codes = ["2EAA90","E3242B","1F5A66","18D2C4","0855B1","88C4D8","F7D42F","ED7014"]';

label_strings = ["data   FLYWIRE (F.F.)","data   H01 (Human)","data   Hemibrain (F.F.)","data   Antennal Lobe (F.F.)","data   LARVA (F.F.)","data   MANC (F.F.)","data   ML (Mouse)","data   MM (Mouse)"];
label_strings = label_strings';

length_scale_vector = [1,1000,125,125,1000,125,1,1];
length_scale_vector = length_scale_vector';

directories = ["","","","","","","","D:\MM3_New_Full"];
directories = directories';

brain_volume_human = [103022.3,103022.3,103022.3];
brain_volume_mouse = [5535.4,5535.4,5535.4];
brain_volume_adult_fly = [590,340,120];

colors = [];
for i=1:length(color_codes)
    str = char(color_codes(i));
    colors = [colors;sscanf(str(1:end),'%2x%2x%2x',[1 3])/255];
end

inclusion_lists = cell(size(directories,1),1);

min_x = Inf;
min_y = Inf;
min_z = Inf;
max_x = 0;
max_y = 0;
max_z = 0;

cd ..
cd ..
cd Datasets\MM3\Skeletons

for i=8:8

    disp(i)

    cd(directories(i));

    files = dir('*.txt');
    
    disp(size(files,1))

    for j=1:size(files,1)
        
        tic
        
        filename = files(j).name;
        filename = string(filename(1:end-4));
        
        swc_string = strcat(filename,".txt");
        
        disp(j)

        if i==1
            swc_data = dlmread(swc_string);
            swc_points = swc_data(:,3:6);
            swc_points = swc_points ./ length_scale_vector(i);
        elseif i==2
            swc_data = dlmread(swc_string, ' ', 14, 0);
            swc_points = swc_data(:,3:6);
            swc_points = swc_points ./ length_scale_vector(i);
        elseif i==3 || i==4
            swc_data = dlmread(swc_string, ' ', 7, 0);
            swc_points = swc_data(:,3:6);
            swc_points = swc_points ./ length_scale_vector(i);
        elseif i==5
            swc_data = dlmread(swc_string);
            swc_points = swc_data(:,3:6);
            swc_points = swc_points ./ length_scale_vector(i);
        elseif i==6
            swc_data = dlmread(swc_string, ' ', 7, 0);
            swc_points = swc_data(:,3:6);
            swc_points = swc_points ./ length_scale_vector(i);
        elseif i==7
            swc_data = dlmread(swc_string, ' ', 1, 0);
            swc_points = swc_data(:,3:6);
            swc_points = swc_points ./ length_scale_vector(i);
        elseif i==8
            swc_data = dlmread(swc_string);
            swc_points = swc_data(:,3:6);
            swc_points = swc_points ./ length_scale_vector(i);
        end
        
        if min(swc_data(:,3)) < min_x
            min_x = min(swc_data(:,3));
        end
        
        if min(swc_data(:,4)) < min_y
            min_y = min(swc_data(:,4));
        end
        
        if min(swc_data(:,5)) < min_z
            min_z = min(swc_data(:,5));
        end
        
        if max(swc_data(:,3)) > max_x
            max_x = max(swc_data(:,3));
        end
        
        if max(swc_data(:,4)) > max_y
            max_y = max(swc_data(:,4));
        end
        
        if max(swc_data(:,5)) > max_z
            max_z = max(swc_data(:,5));
        end
        
        toc
        
    end
    
end

bounding_volumes = [min_x,max_x,min_y,max_y,min_z,max_z];

for i=8:8

    disp(i)

    cd(directories(i));

    files = dir('*.txt');
    
    inclusion_list = [];
    
    disp(size(files,1))

    for j=1:size(files,1)
        
        tic
        
        filename = files(j).name;
        filename = string(filename(1:end-4));
        
        swc_string = strcat(filename,".txt");
        
        disp(j)

        if i==1
            swc_data = dlmread(swc_string);
            swc_points = swc_data(:,3:6);
            swc_points = swc_points ./ length_scale_vector(i);
        elseif i==2
            swc_data = dlmread(swc_string, ' ', 14, 0);
            swc_points = swc_data(:,3:6);
            swc_points = swc_points ./ length_scale_vector(i);
        elseif i==3 || i==4
            swc_data = dlmread(swc_string, ' ', 7, 0);
            swc_points = swc_data(:,3:6);
            swc_points = swc_points ./ length_scale_vector(i);
        elseif i==5
            swc_data = dlmread(swc_string);
            swc_points = swc_data(:,3:6);
            swc_points = swc_points ./ length_scale_vector(i);
        elseif i==6
            swc_data = dlmread(swc_string, ' ', 7, 0);
            swc_points = swc_data(:,3:6);
            swc_points = swc_points ./ length_scale_vector(i);
        elseif i==7
            swc_data = dlmread(swc_string, ' ', 1, 0);
            swc_points = swc_data(:,3:6);
            swc_points = swc_points ./ length_scale_vector(i);
        elseif i==8
            swc_data = dlmread(swc_string);
            swc_points = swc_data(:,3:6);
            swc_points = swc_points ./ length_scale_vector(i);
        end
        
        if min(swc_points(:,1)) > bounding_volumes(1) + 0.05*(bounding_volumes(2) - bounding_volumes(1))
            if max(swc_points(:,1)) < bounding_volumes(2) - 0.05*(bounding_volumes(2) - bounding_volumes(1))
                if min(swc_points(:,2)) > bounding_volumes(3) + 0.05*(bounding_volumes(4) - bounding_volumes(3))
                    if max(swc_points(:,2)) < bounding_volumes(4) - 0.05*(bounding_volumes(4) - bounding_volumes(3))
                        if min(swc_points(:,3)) > bounding_volumes(5) + 0.05*(bounding_volumes(6) - bounding_volumes(5))
                            if max(swc_points(:,3)) < bounding_volumes(6) - 0.05*(bounding_volumes(6) - bounding_volumes(5))
                                inclusion_list = [inclusion_list;filename];
                            end
                        end
                    end
                end
            end
        end
        
        toc
        
    end
    
end

cd ..
cd ..
cd ..
cd ProcessedData\MM3

save("MM3_Inclusion_List.mat",'inclusion_list');

