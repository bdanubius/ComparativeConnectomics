clear

cd ..
cd ..
cd Datasets\Celegans\Synapses

Synapse_Info = readtable("Subgraph_Synaptic_Connections_Summary.csv");

cd ..
cd Skeletons

neuron_lengths = zeros(size(Synapse_Info,1),1);
for i=1:size(Synapse_Info,1)
    disp(i)
    segment_data = readtable(strcat(string(Synapse_Info{i,1}),"_segments.csv"));
    ID_List1 = double(segment_data{:,1});
    ID_List2 = double(segment_data{2:end,2});
    [Lia,connecting_IDs] = ismember(ID_List2,ID_List1);
    %connecting_IDs = double(segment_data{2:end,2}+1);
    all_points = double(segment_data{:,3:5});
    neuron_length = sum(sqrt(sum((all_points(2:end,:) - all_points(connecting_IDs,:)).^2,2)));
    neuron_lengths(i) = neuron_length;
end
   
cd ..
cd ..
cd ..

cd ProcessedData\Celegans

CELEGANS_Total_Data = horzcat(string(Synapse_Info{:,1}),string(neuron_lengths),string(Synapse_Info{:,7}),string(double(Synapse_Info{:,7}) ./ neuron_lengths),string(Synapse_Info{:,6}),string(Synapse_Info{:,7}));

CELEGANS_Total_Data = num2cell(CELEGANS_Total_Data);

tC=cell2table(CELEGANS_Total_Data,'VariableNames',{'CellID','Length','Synapses','Density','Degree','Weighted_Degree'});

writetable(tC,"Celegans_AggregatedData.csv");