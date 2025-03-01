clear

collect_lengths = true;

collect_figure_data = false;

% synapse_locations(:,1:2) = fliplr(synapse_locations(:,1:2));

if collect_lengths

    cd ..
    cd ..
    cd ..
    cd ProcessedData/Zebrafish/Cleaned_Data

    files = dir('*.mat');
    
    bytes_array = zeros(size(files,1),1);
    for i=1:length(files)
        bytes_array(i) = files(i).bytes;
    end
    
    [B,I_Bytes] = sort(bytes_array,'ascend');

    neuron_lengths = [];
    
    for i=1:3
        
        disp(i)

%         try

            if collect_lengths
                segment_lengths = [];
                segment_lengths_neurons = [];
                neuron_lengths = [];
                neuron_lengths_neurons = [];
                segment_curvatures = [];
                segment_curvatures_neurons = [];
                neuron_curvatures = [];
                neuron_curvatures_neurons = [];
                segment_torsions = [];
                segment_torsions_neurons = [];
                neuron_torsions = [];
                neuron_torsions_neurons = [];
            end
            filename = files(I_Bytes(i)).name;
            dataset = load(string(filename));
            dataset = dataset.skeleton_cell;
            for j=1:size(dataset,1)
                if iscell(dataset{j,2}) && dataset{j,4} ~= "NaN"
                    if collect_lengths
                        segment_subset_lengths = [];
                        segment_subset_curvatures = [];
                        segment_subset_torsions = [];
                        path_set = dataset{j,2};
                        point_set = dataset{j,1};
                        point_set = point_set;
                        for k=1:length(path_set)
                            point_subset = point_set(path_set{k},:);
                            segment_length = sum(sqrt(sum((point_subset(2:end,:) - point_subset(1:end-1,:)).^2,2)));
                            segment_subset_lengths = [segment_subset_lengths;segment_length];
                            if size(point_subset,1) >= 7
                                path1_points = point_subset;
                                torsion_segments = NaN(size(path1_points,1),1);
                                curvature_segments = NaN(size(path1_points,1),1);
                                for iter=4:1:size(path1_points)-3
                                    %% Method 1: Torsion
                                    path_length = sqrt(sum((path1_points(iter,:) - path1_points(iter-1,:)).^2));
                                    x_prime = (-1*path1_points(iter+2,1) + 8*path1_points(iter+1,1) - 8*path1_points(iter-1,1) + 1*path1_points(iter-2,1))/(12*path_length);
                                    x_primeprime = (-1*path1_points(iter+2,1) + 16*path1_points(iter+1,1) - 30*path1_points(iter,1) + 16*path1_points(iter-1,1) - 1*path1_points(iter-2,1))/(12*path_length^2);
                                    x_primeprimeprime = (-1*path1_points(iter+3,1) + 8*path1_points(iter+2,1) - 13*path1_points(iter+1,1) + 13*path1_points(iter-1,1) - 8*path1_points(iter-2,1) + 1*path1_points(iter-3,1))/(8*path_length^3);
                                    y_prime = (-1*path1_points(iter+2,2) + 8*path1_points(iter+1,2) - 8*path1_points(iter-1,2) + 1*path1_points(iter-2,2))/(12*path_length);
                                    y_primeprime = (-1*path1_points(iter+2,2) + 16*path1_points(iter+1,2) - 30*path1_points(iter,2) + 16*path1_points(iter-1,2) - 1*path1_points(iter-2,2))/(12*path_length^2);
                                    y_primeprimeprime = (-1*path1_points(iter+3,2) + 8*path1_points(iter+2,2) - 13*path1_points(iter+1,2) + 13*path1_points(iter-1,2) - 8*path1_points(iter-2,2) + 1*path1_points(iter-3,2))/(8*path_length^3);
                                    z_prime = (-1*path1_points(iter+2,3) + 8*path1_points(iter+1,3) - 8*path1_points(iter-1,3) + 1*path1_points(iter-2,3))/(12*path_length);
                                    z_primeprime = (-1*path1_points(iter+2,3) + 16*path1_points(iter+1,3) - 30*path1_points(iter,3) + 16*path1_points(iter-1,3) - 1*path1_points(iter-2,3))/(12*path_length^2);
                                    z_primeprimeprime = (-1*path1_points(iter+3,3) + 8*path1_points(iter+2,3) - 13*path1_points(iter+1,3) + 13*path1_points(iter-1,3) - 8*path1_points(iter-2,3) + 1*path1_points(iter-3,3))/(8*path_length^3);
                                    %% Method 2: Torsion Energy
                                    torsion_energy_density = x_primeprimeprime^2 + y_primeprimeprime^2 + z_primeprimeprime^2;
                                    torsion_segments(iter,1) = torsion_energy_density;
                                    %% Method 3: Curvature
                                    curvature = sqrt((z_primeprime*y_prime - y_primeprime*z_prime)^2 + (x_primeprime*z_prime - z_primeprime*x_prime)^2 + (y_primeprime*x_prime - x_primeprime*y_prime)^2) / ((x_prime^2 + y_prime^2 + z_prime^2)^(3/2));
                                    curvature_segments(iter,1) = curvature;
                                end
                                total_path_length = sum(sqrt(sum((path1_points(1:end-1,:) - path1_points(2:end,:)).^2)));
                                segment_curvature = nanmean(curvature_segments .* total_path_length);
                                segment_torsion = nanmean(torsion_segments .* total_path_length.^4);
                            else
                                segment_curvature = NaN;
                                segment_torsion = NaN;
                            end
                            segment_subset_curvatures = [segment_subset_curvatures;segment_curvature];
                            segment_subset_torsions = [segment_subset_torsions;segment_torsion];
                        end
                        segment_lengths = [segment_lengths;segment_subset_lengths];
                        segment_lengths_neurons = [segment_lengths_neurons;repelem(dataset{j,4},length(segment_subset_lengths))'];
                        neuron_lengths = [neuron_lengths;sum(segment_subset_lengths)];
                        neuron_lengths_neurons = [neuron_lengths_neurons;dataset{j,4}];
                        segment_curvatures = [segment_curvatures;segment_subset_curvatures];
                        segment_curvatures_neurons = [segment_curvatures_neurons;repelem(dataset{j,4},length(segment_subset_lengths))'];
                        neuron_curvatures = [neuron_curvatures;nanmean(segment_subset_curvatures).*length(segment_subset_lengths)];
                        neuron_curvatures_neurons = [neuron_curvatures_neurons;dataset{j,4}];
                        segment_torsions = [segment_torsions;segment_subset_torsions];
                        segment_torsions_neurons = [segment_torsions_neurons;repelem(dataset{j,4},length(segment_subset_lengths))'];
                        neuron_torsions = [neuron_torsions;nanmean(segment_subset_torsions).*length(segment_subset_lengths)];
                        neuron_torsions_neurons = [neuron_torsions_neurons;dataset{j,4}];
                    end
                end
            end
                
            if collect_lengths

                cd ..

                cd Cleaned_Lengths

                length_data = horzcat(segment_lengths,segment_lengths_neurons);

                save(strcat(string(filename(1:end-4)),"_Lengths.mat"),'length_data')

                cd ..

                cd Cleaned_Lengths_Neurons

                length_data = horzcat(neuron_lengths,neuron_lengths_neurons);

                save(strcat(string(filename(1:end-4)),"_Lengths.mat"),'length_data')

                cd ..

                cd Cleaned_Curvatures

                curvature_data = horzcat(segment_curvatures,segment_curvatures_neurons);

                save(strcat(string(filename(1:end-4)),"_Curvatures.mat"),'curvature_data')

                cd ..

                cd Cleaned_Curvatures_Neurons

                curvature_data = horzcat(neuron_curvatures,neuron_curvatures_neurons);

                save(strcat(string(filename(1:end-4)),"_Curvatures.mat"),'curvature_data')

                cd ..

                cd Cleaned_Torsions

                torsion_data = horzcat(segment_torsions,segment_torsions_neurons);

                save(strcat(string(filename(1:end-4)),"_Torsions.mat"),'torsion_data')

                cd ..

                cd Cleaned_Torsions_Neurons

                torsion_data = horzcat(neuron_torsions,neuron_torsions_neurons);

                save(strcat(string(filename(1:end-4)),"_Torsions.mat"),'torsion_data')

                cd ..

                cd Cleaned_Data

            end

%         end  

    end

    cd ..

end

cd Cleaned_Synapses_Neurons
syn_files = dir('*.mat');
neuron_synapses = [];
neuron_synapse_ids = [];
for i=1:length(syn_files)
    filename = syn_files(i).name;
    dataset = load(string(filename));
    dataset = dataset.synapse_data;
    neuron_synapses = [neuron_synapses;(double(dataset(:,1)))];
    neuron_synapse_ids = [neuron_synapse_ids;dataset(:,2)];
end

cd ..

cd Cleaned_Synapses
syn_files = dir('*.mat');
segment_synapses = [];
segment_synapse_ids = [];
for i=1:length(syn_files)
    filename = syn_files(i).name;
    dataset = load(string(filename));
    dataset = dataset.synapse_data;
    segment_synapses = [segment_synapses;(double(dataset(:,1)))];
    segment_synapse_ids = [segment_synapse_ids;dataset(:,2)];
end

cd ..

cd Cleaned_Lengths_Neurons
syn_files = dir('*.mat');
neuron_lengths = [];
neuron_length_ids = [];
for i=1:length(syn_files)
    filename = syn_files(i).name;
    dataset = load(string(filename));
    dataset = dataset.length_data;
    neuron_lengths = [neuron_lengths;(double(dataset(:,1)))];
    neuron_length_ids = [neuron_length_ids;dataset(:,2)];
end

cd ..

cd Cleaned_Lengths
syn_files = dir('*.mat');
segment_lengths = [];
segment_length_ids = [];
for i=1:length(syn_files)
    filename = syn_files(i).name;
    dataset = load(string(filename));
    dataset = dataset.length_data;
    segment_lengths = [segment_lengths;(double(dataset(:,1)))];
    segment_length_ids = [segment_length_ids;dataset(:,2)];
end

cd ..

cd Cleaned_Curvatures_Neurons
syn_files = dir('*.mat');
neuron_curvatures = [];
neuron_curvature_ids = [];
for i=1:length(syn_files)
    filename = syn_files(i).name;
    dataset = load(string(filename));
    dataset = dataset.curvature_data;
    neuron_curvatures = [neuron_curvatures;(double(dataset(:,1)))];
    neuron_curvature_ids = [neuron_curvature_ids;dataset(:,2)];
end

cd ..

cd Cleaned_Curvatures
syn_files = dir('*.mat');
segment_curvatures = [];
segment_curvature_ids = [];
for i=1:length(syn_files)
    filename = syn_files(i).name;
    dataset = load(string(filename));
    dataset = dataset.curvature_data;
    segment_curvatures = [segment_curvatures;(double(dataset(:,1)))];
    segment_curvature_ids = [segment_curvature_ids;dataset(:,2)];
end

cd ..

cd Cleaned_Torsions_Neurons
syn_files = dir('*.mat');
neuron_torsions = [];
neuron_torsion_ids = [];
for i=1:length(syn_files)
    filename = syn_files(i).name;
    dataset = load(string(filename));
    dataset = dataset.torsion_data;
    neuron_torsions = [neuron_torsions;(double(dataset(:,1)))];
    neuron_torsion_ids = [neuron_torsion_ids;dataset(:,2)];
end

cd ..

cd Cleaned_Torsions
syn_files = dir('*.mat');
segment_torsions = [];
segment_torsion_ids = [];
for i=1:length(syn_files)
    filename = syn_files(i).name;
    dataset = load(string(filename));
    dataset = dataset.torsion_data;
    segment_torsions = [segment_torsions;(double(dataset(:,1)))];
    segment_torsion_ids = [segment_torsion_ids;dataset(:,2)];
end

cd ..

neuron_array = horzcat(neuron_length_ids,neuron_lengths);

%save("Fig3Data1.mat",'neuron_array','-v7.3');

clear neuron_length_ids neuron_lengths

segment_synapse_array = horzcat(segment_synapse_ids,segment_synapses);

%save("Fig3Data2.mat",'segment_synapse_array','-v7.3');

clear segment_synapses segment_synapse_ids

%save("Fig3Data3.mat",'neuron_synapse_ids','-v7.3');

%clear neuron_synapse_ids

segment_array = cell(1,6);

segment_array{1} = segment_length_ids;

clear segment_length_ids

segment_array{2} = segment_lengths;

clear segment_lengths

segment_array{3} = segment_curvature_ids;

clear segment_curvature_ids

segment_array{4} = segment_curvatures;

clear segment_curvatures

segment_array{5} = segment_torsion_ids;

clear segment_torsion_ids

segment_array{6} = segment_torsions;

clear segment_torsions

segment_array = horzcat(segment_array{:});

%save("Fig3Data4.mat",'segment_array','-v7.3');

%clear segment_array
    
clearvars -except neuron_array segment_array segment_synapse_array neuron_synapse_ids neuron_synapses

neuron_dataset = [];
segment_dataset = [];

for i=1:length(neuron_synapse_ids)
    disp(i)
    tic
    [row,col] = find(neuron_array(:,1) == string(neuron_synapse_ids(i)));
    [row1,col1] = find(segment_array(:,1) == string(neuron_synapse_ids(i)));
    [row2,col2] = find(segment_synapse_array(:,1) == string(neuron_synapse_ids(i)));
    %[row3,col3] = find(segment_curvature_ids == neuron_synapse_ids(i));
    %[row4,col4] = find(segment_torsion_ids == neuron_synapse_ids(i));
    if length(row1) == length(row2)
        segment_dataset = [segment_dataset;horzcat(segment_array(row1,2),segment_synapse_array(row2,2),repelem(neuron_synapse_ids(i),length(row1))',segment_array(row1,4),segment_array(row1,6))];
        if ~isempty(row)
            if length(row) == 1
                neuron_dataset = [neuron_dataset;neuron_array(row,2),neuron_synapses(i),mean(double(segment_synapse_array(row2,2))),neuron_synapse_ids(i),nanmean(double(segment_array(row1,4))),nanmean(double(segment_array(row1,6)))];
            end
        end
        %neuron_array(row,:) = [];
        %segment_array(row1,:) = [];
        %segment_synapse_array(row2,:) = [];
    end
    toc
end

synapse_density_data = cell(1,2);

synapse_density_data{1} = neuron_dataset;

synapse_density_data{2} = segment_dataset;

save("Geometry_Data.mat",'synapse_density_data');

