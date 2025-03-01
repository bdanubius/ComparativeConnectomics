clear

collect_synapses = true;

cd ..
cd ..
cd ProcessedData/MM3/Cleaned_Data

files = dir('*.mat');

neuron_lengths = [];
for i=1:53
    disp(i)
%     try
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
                    cd Datasets\MM3\Synapses

                    synapse_data = load(strcat(string(dataset{j,4}),".csv"));
                    synapse_locations = synapse_data(:,1:3) ./ [250 250 25];
                    
                    [row,col] = find(synapse_data(:,4)==0);
                    
                    synapse_locations = synapse_locations(row,:);
                    
                    cd ..
                    cd ..
                    cd ..
                    cd ProcessedData/MM3/Cleaned_Data

                    [Idx,D] = knnsearch(dataset{j,1},synapse_locations,'K',1);
                    
                    disp( 100 * length(D(D < 5)) / length(D) )
                    
                    neuron_synapses = [neuron_synapses;sum(D < 5)];
                    neuron_synapse_neurons = [neuron_synapse_neurons;dataset{j,4}];
                    Idx_List = Idx(D < 5);
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
        
%     end
    
end