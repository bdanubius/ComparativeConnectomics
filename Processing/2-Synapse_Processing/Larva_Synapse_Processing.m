clear

% Flag to determine whether to collect synapse data
collect_synapses = true;

% Navigate to the directory containing cleaned skeleton data
cd ..
cd ..
cd ProcessedData/Larva/Cleaned_Data

% List all processed neuron skeleton data files (.mat)
files = dir('*.mat');

% Navigate to the directory containing synapse data
cd ..
cd ..
cd ..
cd Datasets\Larva\Synapses

% List all synapse data files (.csv) and store their names without extensions
filelist = dir('*.csv');
filestrings = string(zeros(size(filelist,1),1));
for i=1:size(filestrings,1)
    filestring = char(filelist(i).name);
    filestring = filestring(1:end-4); % Remove .csv extension
    filestrings(i) = string(filestring);
end

% Return to the cleaned skeleton data directory
cd ..
cd ..
cd ..
cd ProcessedData/Larva/Cleaned_Data

% Initialize an array to store neuron lengths
neuron_lengths = [];

% Loop through the first four cleaned skeleton data files
for i=1:4
    disp(i) % Display the current iteration
    
    % Initialize synapse storage variables if collecting synapses
    if collect_synapses
        segment_synapses = [];
        segment_synapse_neurons = [];
        neuron_synapses = [];
        neuron_synapse_neurons = [];
    end
    
    % Load the processed skeleton data
    filename = files(i).name;
    dataset = load(string(filename));
    dataset = dataset.skeleton_cell;
    
    % Loop through the dataset, processing each neuron's structure
    for j=1:size(dataset,1)
        disp(j) % Display the current sub-iteration

        try
            % Check if the current neuron entry contains valid segment data
            if iscell(dataset{j,2})
                
                if collect_synapses
                    % Navigate to the synapse data directory
                    cd ..
                    cd ..
                    cd ..
                    cd Datasets\Larva\Synapses
                    
                    % Load synapse location data for the current neuron
                    synapse_locations = readmatrix(strcat(string(dataset{j,4}),".csv"));
                    
                    % Return to the cleaned skeleton data directory
                    cd ..
                    cd ..
                    cd ..
                    cd ProcessedData/Larva/Cleaned_Data
                    
                    % Perform nearest neighbor search to match synapses to neuron segments
                    [Idx,D] = knnsearch(dataset{j,1},synapse_locations,'K',1);
                    neuron_synapses = [neuron_synapses;sum(D < 3000)]; % Count synapses within threshold distance (excluding synapses unreasonably far from neuron)
                    neuron_synapse_neurons = [neuron_synapse_neurons;dataset{j,4}];
                    
                    % Filter indices for synapses within the distance threshold
                    Idx_List = Idx(D < 3000);
                    counts_list = zeros(length(dataset{j,1}),1);
                    for k=1:length(Idx_List)
                        counts_list(Idx_List(k)) = counts_list(Idx_List(k)) + 1;
                    end
                    
                    % Extract segment structures
                    segments = dataset{j,2};
                    data_points = dataset{j,1};
                    
                    % Count synapses per segment
                    for k=1:length(segments)
                        segment_synapses = [segment_synapses;sum(counts_list(segments{1,k}))];
                        segment_synapse_neurons = [segment_synapse_neurons;dataset{j,4}];
                    end
                end
                
                % Check if the neuron has corresponding synapse data
                if sum(double(contains(filestrings,string(dataset{j,4})))) > 0
                    neuron_lengths = [neuron_lengths;length(dataset{j,2})];
                end
            end
        end
    end
    
    % Save synapse data if collection is enabled
    if collect_synapses
        cd ..
        
        % Save segment-level synapse data
        cd Cleaned_Synapses
        synapse_data = horzcat(segment_synapses, segment_synapse_neurons);
        save(strcat(string(filename(1:end-4)),"_Synapses.mat"),'synapse_data')
        cd ..
        
        % Save neuron-level synapse data
        cd Cleaned_Synapses_Neurons
        synapse_data = horzcat(neuron_synapses, neuron_synapse_neurons);
        save(strcat(string(filename(1:end-4)),"_Synapses.mat"),'synapse_data')
        cd ..
        
        % Return to cleaned data directory
        cd Cleaned_Data
    end
end
