clear

cd ..
cd ..
cd Datasets\MANC

%% First filtering out low quality neurons

MANC_labels = readtable('MANC_ROIs.csv');

Body_ID_List = string(MANC_labels{:,2});

all_row = [];

[row,col] = find(string(MANC_labels{:,12})=="Prelim Roughly traced");

all_row = [all_row;row];

[row,col] = find(string(MANC_labels{:,12})=="0.5assign");

all_row = [all_row;row];

[row,col] = find(string(MANC_labels{:,12})=="Orphan");

all_row = [all_row;row];

[row,col] = find(string(MANC_labels{:,12})=="Sensory Anchor");

all_row = [all_row;row];

[row,col] = find(string(MANC_labels{:,12})=="RT Orphan");

all_row = [all_row;row];

[row,col] = find(string(MANC_labels{:,12})=="PRT Orphan");

all_row = [all_row;row];

Body_ID_List(all_row,:) = [];

cd ..
cd ..
cd Datasets\MANC\Skeletons

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

inclusion_list = [];

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
cd ProcessedData\MANC

Body_ID_List = intersect(Body_ID_List,inclusion_list);

writematrix(Body_ID_List,"MANC_Inclusion_List.csv");