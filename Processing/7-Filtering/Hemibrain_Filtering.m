clear

cd ..
cd ..
cd Datasets\Hemibrain

%% Filtering out low quality neurons first

hemibrain_labels = readtable('Hemibrain_ROIs.csv');

Hemibrain_Exclusion = [];

Body_ID_List = string(hemibrain_labels{:,2});

[row,col] = find(string(hemibrain_labels{:,13})=="Leaves");

Hemibrain_Exclusion = [Hemibrain_Exclusion;Body_ID_List(row)];

[row,col] = find(string(hemibrain_labels{:,13})=="Orphan");

Hemibrain_Exclusion = [Hemibrain_Exclusion;Body_ID_List(row)];

[row,col] = find(string(hemibrain_labels{:,13})=="Orphan hotknife");

Hemibrain_Exclusion = [Hemibrain_Exclusion;Body_ID_List(row)];

[row,col] = find(string(hemibrain_labels{:,13})=="Orphan-artifact");

Hemibrain_Exclusion = [Hemibrain_Exclusion;Body_ID_List(row)];

cd ..
cd ..
cd Datasets\Hemibrain\Skeletons

%% Then filtering out neurons near the boundary

min_x = Inf;
min_y = Inf;
min_z = Inf;
max_x = 0;
max_y = 0;
max_z = 0;

files = dir('*.swc');

disp(size(files,1))

for j=1:size(files,1)

    tic

    filename = files(j).name;
    filename = string(filename(1:end-4));

    swc_string = strcat(filename,".swc");

    disp(j)

    swc_data = dlmread(swc_string, ' ', 7, 0);
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

files = dir('*.swc');

exclusion_list = [];

disp(size(files,1))

for j=1:size(files,1)

    tic

    filename = files(j).name;
    filename = string(filename(1:end-4));

    swc_string = strcat(filename,".swc");

    disp(j)

    swc_data = dlmread(swc_string, ' ', 7, 0);
    swc_points = swc_data(:,3:6);

    if min(swc_points(:,1)) > bounding_volumes(1) + 0.05*(bounding_volumes(2) - bounding_volumes(1))
        if max(swc_points(:,1)) < bounding_volumes(2) - 0.05*(bounding_volumes(2) - bounding_volumes(1))
            if min(swc_points(:,2)) > bounding_volumes(3) + 0.05*(bounding_volumes(4) - bounding_volumes(3))
                if max(swc_points(:,2)) < bounding_volumes(4) - 0.05*(bounding_volumes(4) - bounding_volumes(3))
                    if min(swc_points(:,3)) > bounding_volumes(5) + 0.05*(bounding_volumes(6) - bounding_volumes(5))
                        if max(swc_points(:,3)) < bounding_volumes(6) - 0.05*(bounding_volumes(6) - bounding_volumes(5))
                        else
                            exclusion_list = [exclusion_list;filename];
                        end
                    else
                        exclusion_list = [exclusion_list;filename];
                    end
                else
                    exclusion_list = [exclusion_list;filename];
                end
            else
                exclusion_list = [exclusion_list;filename];
            end
        else
            exclusion_list = [exclusion_list;filename];
        end
    else
        exclusion_list = [exclusion_list;filename];
    end

    toc

end

cd ..
cd ..
cd ..
cd ProcessedData\Hemibrain

Hemibrain_Exclusion = unique(vertcat(Hemibrain_Exclusion,exclusion_list));

writematrix(Hemibrain_Exclusion,"Hemibrain_Exclusion_List.csv");