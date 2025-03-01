clear

collect_synapses = true;

cd ..
cd ..
cd Datasets\H01_104\Synapses


if collect_synapses

    synapse_locations = readmatrix("full_synapse_locations.txt");

end

cd ..
cd ..
cd ..
cd ProcessedData/H01_104/Cleaned_Data

files = dir('*.mat');

neuron_lengths = [];
for i=1:1
    disp(i)

    if collect_synapses
        segment_synapses = [];
        segment_synapse_neurons = [];
        neuron_synapses = [];
        neuron_synapse_neurons = [];
    end
    filename = files(i).name;
    dataset = load(string(filename));
    dataset = dataset.skeleton_cell;
    for j=1:size(dataset,1)
        disp(j)

        try

            if iscell(dataset{j,2})
                
                point_set = dataset{j,1};
                dataset{j,1} = point_set .* [32 32 33];

                if collect_synapses

                    [Idx,D] = knnsearch(dataset{j,1},synapse_locations,'K',1);
                    neuron_synapses = [neuron_synapses;sum(D < 1000)];
                    neuron_synapse_neurons = [neuron_synapse_neurons;dataset{j,4}];
                    Idx_List = Idx(D < 1000);
                    counts_list = zeros(length(dataset{j,1}),1);
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
