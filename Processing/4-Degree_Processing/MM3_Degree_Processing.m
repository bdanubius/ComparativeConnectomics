clear

close all

read_synapses = true;

%% Data Collection

if ~read_synapses
    
    cd ..
    cd ..
    cd Datasets\MM3\Synapses
    
    files = dir('*.csv');
    
    synapse_connectivity = cell(size(files,1),1);
    for i=1:size(files,1)
        disp(i)
        filename = files(i).name;
        s = string(filename(1:end-4));
        if isfile(string(filename))
            opts = detectImportOptions(string(filename));
            opts = setvartype(opts, "Var5", 'string');
            synapse_locations = readtable(string(filename),opts);
            [row,col] = find(synapse_locations{:,4} == 0);
            synapse_locations = synapse_locations(row,:);
            synapse_connectivity{i,1} = horzcat(int64(double(repelem(s,size(synapse_locations,1))')),int64(synapse_locations{:,1:3}),int64(double(string(synapse_locations{:,5}))),int64(synapse_locations{:,6:8}));
        end
    end
    
    cd ..
    cd ..
    cd ..
    cd ProcessedData/MM3
    
    save("Synapse_Connectivity.mat",'synapse_connectivity','-v7.3');
    
else
    
    cd ..
    cd ..
    cd Datasets\MM3\Synapses
    
    cd ..
    cd ..
    cd ..
    cd ProcessedData/MM3
    
    synapse_connectivity = load("Synapse_Connectivity.mat");
    synapse_connectivity = synapse_connectivity.synapse_connectivity;
    
end

synapse_connectivity(:,2:4) = [];
synapse_connectivity(:,3:5) = [];

presyn_id = synapse_connectivity(:,1);
postsyn_id = synapse_connectivity(:,2);

clear synapse_connectivity

cd ..
cd ..
cd Datasets/MM3/Skeletons

files = dir("*.txt");
filenames = zeros(size(files,1),1);

for i=1:length(files)
    filename = files(i).name;
    filenames(i) = int64(double(string(filename(1:end-4))));
end

unique_neurons = filenames;

full_degree_list = zeros(length(unique_neurons),1);
degree_list = zeros(length(unique_neurons),1);

[Lia,Lib] = ismember(presyn_id,unique_neurons);

[Lia2,Lib2] = ismember(postsyn_id,unique_neurons);

A_sparse = sparse(length(unique_neurons),length(unique_neurons));

Lib_Pairs = horzcat(Lib,Lib2);

clear Lia Lia2 Lib Lib2

Lib_Pairs = unique(Lib_Pairs,'rows');
Lib_Pairs = vertcat(Lib_Pairs,fliplr(Lib_Pairs));
Lib_Pairs = unique(Lib_Pairs,'rows');

A_sparse2 = sparse(length(unique_neurons),length(unique_neurons));

for i=1:length(Lib_Pairs)
    disp(i)
    if Lib_Pairs(i,1) < Lib_Pairs(i,2) && Lib_Pairs(i,1)~=0 && Lib_Pairs(i,2)~=0
        degree_list(Lib_Pairs(i,1)) = degree_list(Lib_Pairs(i,1)) + 1;
        degree_list(Lib_Pairs(i,2)) = degree_list(Lib_Pairs(i,2)) + 1;
        A_sparse2(Lib_Pairs(i,1),Lib_Pairs(i,2)) = A_sparse2(Lib_Pairs(i,1),Lib_Pairs(i,2)) + 1;
    end
end


MM3_Degrees = horzcat(unique_neurons,degree_list);

MM3_Degrees(:,2) = string(degree_list);

save("MM3_Degree_Data.mat",'MM3_Degrees');

