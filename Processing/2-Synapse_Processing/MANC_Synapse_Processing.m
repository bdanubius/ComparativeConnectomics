clear

collect_synapses = true;

cd ..
cd ..
cd ProcessedData/MANC/Cleaned_Data

files = dir('*.mat');

neuron_lengths = [];
for i=1:22
    disp(i)
    try
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
            if iscell(dataset{j,2})
                neuron_lengths = [neuron_lengths;length(dataset{j,2})];
                if collect_synapses 
                    
                    cd ..
                    cd ..
                    cd ..
                    cd Datasets\MANC\Synapses

                    synapse_data = readtable(strcat(string(dataset{j,4}),".csv"));
                    synapse_locations = table2array(synapse_data(:,5:7));
                    
                    cd ..
                    cd ..
                    cd ..
                    cd ProcessedData/MANC/Cleaned_Data

                    [Idx,D] = knnsearch(dataset{j,1},synapse_locations,'K',1);
                    neuron_synapses = [neuron_synapses;sum(D < 300)];
                    neuron_synapse_neurons = [neuron_synapse_neurons;dataset{j,4}];
                    Idx_List = Idx(D < 300);
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
    
end
