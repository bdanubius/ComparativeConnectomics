clear

% Note: Neurons have fragments, can't use unique!

cd ..
cd ..
cd Datasets\Zebrafish\Synapses

synapse_locations = readtable("zebrafish_synapses.csv");

presyn_id = double(synapse_locations{:,20});
postsyn_id = double(synapse_locations{:,19});
synapse_size = double(synapse_locations{:,18});

presyn_id(isnan(postsyn_id)) = [];
synapse_size(isnan(postsyn_id)) = [];
postsyn_id(isnan(postsyn_id)) = [];
postsyn_id(isnan(presyn_id)) = [];
synapse_size(isnan(presyn_id)) = [];
presyn_id(isnan(presyn_id)) = [];

cd ..
cd Skeletons

files = dir("*.swc");
filenames = zeros(size(files,1),1);

for i=1:length(files)
    filename = files(i).name;
    filenames(i) = double(string(filename(1:end-4)));
end

unique_neurons = filenames;

full_degree_list = zeros(length(unique_neurons),1);
degree_list = zeros(length(unique_neurons),1);

[Lia,Lib] = ismember(presyn_id,unique_neurons);

[Lia2,Lib2] = ismember(postsyn_id,unique_neurons);

A_sparse = sparse(length(unique_neurons),length(unique_neurons));

for i=1:length(Lib)
    disp(i)
    if Lib(i) ~= 0 
        full_degree_list(Lib(i)) = full_degree_list(Lib(i)) + 1;
    end
    if Lib2(i) ~= 0
        full_degree_list(Lib2(i)) = full_degree_list(Lib2(i)) + 1;
    end
    if Lib(i) ~= 0 && Lib2(i) ~= 0
        A_sparse(Lib(i),Lib2(i)) = A_sparse(Lib(i),Lib2(i)) + 1;
    end
end

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

Zebrafish_Degrees = horzcat(unique_neurons,degree_list,full_degree_list);

cd ..
cd ..
cd ..

cd ProcessedData/Zebrafish

save("Zebrafish_Degree_Data.mat",'Zebrafish_Degrees');