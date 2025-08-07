clear
close all

min_x = Inf;
min_y = Inf;
min_z = Inf;
max_x = 0;
max_y = 0;
max_z = 0;

cd ..
cd ..
cd Datasets\MM3\Skeletons

files = dir('*.txt');

disp(size(files,1))

for j=1:size(files,1)

    tic

    filename = files(j).name;
    filename = string(filename(1:end-4));

    swc_string = strcat(filename,".txt");

    disp(j)

    swc_data = dlmread(swc_string);
    swc_points = swc_data(:,3:6);

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

bounding_volumes = [min_x,max_x,min_y,max_y,min_z,max_z];

files = dir('*.txt');

inclusion_list = [];

disp(size(files,1))

for j=1:size(files,1)

    tic

    filename = files(j).name;
    filename = string(filename(1:end-4));

    swc_string = strcat(filename,".txt");

    disp(j)

    swc_data = dlmread(swc_string);
    swc_points = swc_data(:,3:6);

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

cd ..
cd ..
cd ..
cd ProcessedData\MM3

save("MM3_Inclusion_List.mat",'inclusion_list');

