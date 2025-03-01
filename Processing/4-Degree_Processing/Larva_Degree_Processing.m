clear  % Clear workspace variables

cd ..  % Move up one directory
cd ..  % Move up another directory
cd Datasets\Larva\Synapses  % Navigate to the folder containing synapse data

% Load synapse connection data from a text file
synapse_locations = readmatrix("synapse_locations.txt");

% Extract presynaptic and postsynaptic neuron IDs as strings
presyn_id = string(synapse_locations(:,1));
postsyn_id = string(synapse_locations(:,2));

% Get a list of unique neurons involved in synapses
unique_neurons = unique(vertcat(presyn_id, postsyn_id));

% Initialize degree lists for each neuron
full_degree_list = zeros(length(unique_neurons),1); % Total connections
degree_list = zeros(length(unique_neurons),1); % Unique connections

% Map presynaptic and postsynaptic IDs to indices in the unique neuron list
[Lia, Lib] = ismember(presyn_id, unique_neurons);
[Lia2, Lib2] = ismember(postsyn_id, unique_neurons);

% Initialize a sparse adjacency matrix for synaptic connections
A_sparse = sparse(length(unique_neurons), length(unique_neurons));

% Iterate through all synapses and update degree counts and adjacency matrix
for i = 1:length(Lib)
    disp(i) % Display iteration number (for debugging)
    full_degree_list(Lib(i)) = full_degree_list(Lib(i)) + 1; % Increase count for presynaptic neuron
    full_degree_list(Lib2(i)) = full_degree_list(Lib2(i)) + 1; % Increase count for postsynaptic neuron
    A_sparse(Lib(i), Lib2(i)) = A_sparse(Lib(i), Lib2(i)) + 1; % Store connection in adjacency matrix
end

% Store neuron index pairs corresponding to synaptic connections
Lib_Pairs = horzcat(Lib, Lib2);

% Clear temporary variables to free up memory
clear Lia Lia2 Lib Lib2

% Ensure unique connections are recorded bidirectionally
Lib_Pairs = unique(Lib_Pairs, 'rows'); % Remove duplicate rows
Lib_Pairs = vertcat(Lib_Pairs, fliplr(Lib_Pairs)); % Add reversed pairs to ensure undirected representation
Lib_Pairs = unique(Lib_Pairs, 'rows'); % Remove duplicates again

% Initialize a second sparse adjacency matrix for bidirectional connections
A_sparse2 = sparse(length(unique_neurons), length(unique_neurons));

% Iterate through the unique neuron connections to compute final degree counts
for i = 1:length(Lib_Pairs)
    disp(i) % Display iteration number (for debugging)
    if Lib_Pairs(i,1) < Lib_Pairs(i,2) % Prevent double counting
        degree_list(Lib_Pairs(i,1)) = degree_list(Lib_Pairs(i,1)) + 1;
        degree_list(Lib_Pairs(i,2)) = degree_list(Lib_Pairs(i,2)) + 1;
    end
    A_sparse2(Lib_Pairs(i,1), Lib_Pairs(i,2)) = A_sparse2(Lib_Pairs(i,1), Lib_Pairs(i,2)) + 1; 
end

% Combine neuron IDs with degree lists for final output
Larva_Degrees = horzcat(unique_neurons, degree_list, full_degree_list);

% Navigate to the processed data directory
cd ..
cd ..
cd ..
cd ProcessedData\Larva

% Save the degree data to a .mat file
save("Larva_Degree_Data.mat", 'Larva_Degrees');
