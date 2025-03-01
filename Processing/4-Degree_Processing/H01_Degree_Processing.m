clear

cd ..
cd ..

cd Datasets\H01_104\Skeletons

files = dir("*.swc");
filenames = zeros(size(files,1),1);

for i=1:length(files)
    filename = files(i).name;
    filenames(i) = double(string(filename(1:end-4)));
end

unique_neurons = filenames;

full_degree_list = zeros(length(unique_neurons),1);
degree_list = zeros(length(unique_neurons),1);

cd ..
cd Synapses

files = dir("*.csv");

for iter=1:size(files,1)
    
    synapse_connectivity = load(files(iter).name);
    presyn_id = synapse_connectivity(:,1);
    postsyn_id = synapse_connectivity(:,5);

    [Lia,Lib] = ismember(presyn_id,unique_neurons);

    [Lia2,Lib2] = ismember(postsyn_id,unique_neurons);

    for i=1:length(Lib)
        disp(i)
        if Lib(i) ~= 0 
            full_degree_list(Lib(i)) = full_degree_list(Lib(i)) + 1;
        end
        if Lib2(i) ~= 0
            full_degree_list(Lib2(i)) = full_degree_list(Lib2(i)) + 1;
        end
    end

    Lib_Pairs = horzcat(Lib,Lib2);

    clear Lia Lia2 Lib Lib2

    Lib_Pairs = unique(Lib_Pairs,'rows');

    Lib_Pairs = vertcat(Lib_Pairs,fliplr(Lib_Pairs));

    Lib_Pairs = unique(Lib_Pairs,'rows');

    for i=1:length(Lib_Pairs)
        disp(i)
        if Lib_Pairs(i,1) ~= 0
            degree_list(Lib_Pairs(i,1)) = degree_list(Lib_Pairs(i,1)) + 1;
        end
        if Lib_Pairs(i,2) ~= 0
            degree_list(Lib_Pairs(i,2)) = degree_list(Lib_Pairs(i,2)) + 1;
        end
    end
    
end

degree_list = degree_list ./ 2;

H01_Degrees = horzcat(unique_neurons,degree_list,full_degree_list);

cd ..
cd ..
cd ..
cd ProcessedData\H01_104

save("H01_Degree_Data.mat",'H01_Degrees');