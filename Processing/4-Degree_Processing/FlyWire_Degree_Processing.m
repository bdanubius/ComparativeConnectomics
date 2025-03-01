clear

cd ..
cd ..
cd Datasets\FlyWire\Synapses

opts = detectImportOptions('synapse_coordinates.csv');
opts = setvartype(opts, 'pre_root_id', 'string');
opts = setvartype(opts, 'post_root_id', 'string');%or 'char' if you prefer
presyn_pts = readtable('synapse_coordinates.csv',opts);
presyn_id = string(presyn_pts{:,1});
postsyn_id = string(presyn_pts{:,2});
presyn_pts = table2array(presyn_pts(:,3:5)) ./ 1000;
postsyn_pts = presyn_pts;

cd ..
cd ..
cd ..
cd ProcessedData/FlyWire/Cleaned_Data

presyn_id = fillmissing(presyn_id,'previous');
postsyn_id = fillmissing(postsyn_id,'previous');

unique_neurons = unique(vertcat(presyn_id,postsyn_id));

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

% for i=1:length(unique_neurons)
%     disp(i)
%     tic
%     neuron_id = unique_neurons(i);
%     [row,col] = find(presyn_id==neuron_id);
%     [row2,col] = find(postsyn_id==neuron_id);
%     full_degree_list(i) = length(row) + length(row2);
%     degree_list(i) = length(unique(vertcat(postsyn_id(row),presyn_id(row2))));
%     toc
% end

FlyWire_Degrees = horzcat(unique_neurons,degree_list,full_degree_list);

cd ..

save("FlyWire_Degree_List.mat",'FlyWire_Degrees');