clear

collect_synapses = true;

cd ..
cd ..
cd ProcessedData/FlyWire/Cleaned_Data

files = dir('*.mat');

bytes_array = zeros(size(files,1),1);
for i=1:length(files)
    bytes_array(i) = files(i).bytes;
end

[B,I_Bytes] = sort(bytes_array,'ascend');

if collect_synapses
    cd ..
    cd ..
    cd ..
    cd Datasets\FlyWire\Synapses
    opts = detectImportOptions('synapse_coordinates.csv');
    opts = setvartype(opts, 'pre_root_id', 'string');
    opts = setvartype(opts, 'post_root_id', 'string');%or 'char' if you prefer
    synapse_locations = readtable('synapse_coordinates.csv',opts);
    synapse_pre_id = string(synapse_locations{:,1});
    synapse_post_id = string(synapse_locations{:,2});
    synapse_locations = table2array(synapse_locations(:,3:5)) ./ 1000;
    cd ..
    cd ..
    cd ..
    cd ProcessedData/FlyWire/Cleaned_Data
    for i=2:size(synapse_pre_id)
        if ismissing(synapse_pre_id(i))
            synapse_pre_id(i) = synapse_pre_id(i-1);
        end
        if ismissing(synapse_post_id(i))
            synapse_post_id(i) = synapse_post_id(i-1);
        end
    end
end

neuron_lengths = [];
% 85 to 129 when collect next
for i=1:129
    disp(i)

    if collect_synapses
        segment_synapses = [];
        segment_synapse_neurons = [];
        neuron_synapses = [];
        neuron_synapse_neurons = [];
    end
    filename = files(I_Bytes(i)).name;
    dataset = load(string(filename));
    dataset = dataset.skeleton_cell;
    
    for j=1:size(dataset,1)
        
        disp(j)

        try

            if iscell(dataset{j,2})
                
                point_set = dataset{j,1};
                point_set = point_set ./ 1000;

                if collect_synapses
                    
                    [row,col] = find(synapse_pre_id==dataset{j,4});
                    [row2,col2] = find(synapse_post_id==dataset{j,4});
                    
                    synapse_subset = synapse_locations(union(row,row2),:);
                    
                    [Idx,D] = knnsearch(point_set,synapse_subset,'K',1);
                    
                    neuron_synapses = [neuron_synapses;sum(D < 10)];
                    neuron_synapse_neurons = [neuron_synapse_neurons;dataset{j,4}];
                    Idx_List = Idx(D < 10);
                    counts_list = zeros(length(point_set),1);
                    for k=1:length(Idx_List)
                        counts_list(Idx_List(k)) = counts_list(Idx_List(k)) + 1;
                    end
                    segments = dataset{j,2};
                    data_points = dataset{j,1};
                    for k=1:length(segments)
                        segment_synapses = [segment_synapses;sum(counts_list(segments{1,k}))];
                        segment_synapse_neurons = [segment_synapse_neurons;dataset{j,4}];
                    end

                end
                
                neuron_lengths = [neuron_lengths;length(dataset{j,2})];
                
            end

        end

    end
    
    if collect_synapses

        cd ..

        cd Cleaned_Synapses

        synapse_data = horzcat(segment_synapses,segment_synapse_neurons);

        save(strcat(string(filename(1:end-4)),"_Synapses.mat"),'synapse_data')

        cd ..

        cd Cleaned_Synapses_Neurons

        synapse_data = horzcat(neuron_synapses,neuron_synapse_neurons);

        save(strcat(string(filename(1:end-4)),"_Synapses.mat"),'synapse_data')

        cd ..

        cd Cleaned_Data

    end
    
end