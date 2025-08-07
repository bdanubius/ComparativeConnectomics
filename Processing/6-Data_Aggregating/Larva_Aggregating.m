clear  % Clear workspace

% Navigate to the directory containing processed larval data
cd ..
cd ..
cd ProcessedData\Larva

% Load geometry data containing synapse density information
Larva_Data = load("Geometry_Data.mat");
Larva_Data = Larva_Data.synapse_density_data;
Larva_Data = Larva_Data{1};  % Extract data from the cell array

% Load degree data containing connectivity information
Larva_Degree_Data = load("Larva_Degree_Data.mat");
Larva_Degree_Data = Larva_Degree_Data.Larva_Degrees;

% Match neuron IDs in the geometry data with those in the degree data
[Lia,Lib] = ismember(double(Larva_Data(:,4)), double(Larva_Degree_Data(:,1)));

% Extract corresponding degree data based on matched indices
Larva_Degree_Data = double(Larva_Degree_Data(Lib,:));
Larva_Data = double(Larva_Data);  % Ensure numerical format

% Aggregate data into a single matrix:
% Column 1: Neuron ID
% Column 2: Normalized length (converted to mm)
% Column 3: Number of synapses
% Column 4: Synaptic density (synapses per mm)
% Column 5: Degree (number of unique connections)
% Column 6: Weighted degree (total number of connections including multiple synapses)
Larva_Total_Data = horzcat(Larva_Data(:,4), Larva_Data(:,1)./1000, Larva_Data(:,2), 1000.*Larva_Data(:,2)./Larva_Data(:,1), Larva_Degree_Data(:,2:3));

% Convert numeric data to cell array of strings for table formatting
A = cellstr(horzcat(string(int64(Larva_Total_Data(:,1))), string(Larva_Total_Data(:,2:end))));

% Create a table with appropriate column names
tC = cell2table(A, 'VariableNames', {'CellID', 'Length', 'Synapses', 'Density', 'Degree', 'WeightedDegree'});

% Save the aggregated data as a CSV file
writetable(tC, "Larva_AggregatedData.csv");

