clear

%% Mode

use_soma = true;

%% Continue

save_visualizations = false;

% test_index = 1;

%1:1:46378

skeleton_cell = {};

counter = 0;

bytes_sort = 1;

cd ..
cd Datasets\FlyWire\Skeletons

if bytes_sort

    files = dir('*.swc');

    bytes_array = zeros(size(files,1),1);
    for i=1:length(files)
        bytes_array(i) = files(i).bytes;
    end

    [B,I_Bytes] = sort(bytes_array,'ascend');

end

% Width_Ratios = [];

iter_range = 1:1:99000;

% 1:1:129000

first_iter = iter_range(1);

All_Chain_Data = cell(1000,8);

Chain_Names = string(zeros(length(files),1));
for i=1:length(Chain_Names)
    Chain_Names(i) = string(files(I_Bytes(i)).name);
end

cd ..
cd Synapses

opts = detectImportOptions('synapse_coordinates.csv');
opts = setvartype(opts, 'pre_root_id', 'string');
opts = setvartype(opts, 'post_root_id', 'string');%or 'char' if you prefer
presyn_pts = readtable('synapse_coordinates.csv',opts);
presyn_id = string(presyn_pts{:,1});
postsyn_id = string(presyn_pts{:,2});
presyn_pts = table2array(presyn_pts(:,3:5)) ./ 1000;

presyn_id = fillmissing(presyn_id,'previous');
postsyn_id = fillmissing(postsyn_id,'previous');

cd ..
cd Skeletons

% size(files,1)
for iter=iter_range

    counter = counter + 1;

    disp(iter)
    
    try

        swc_string = string(files(I_Bytes(iter)).name);

        swc_ID = char(swc_string);
        swc_ID = string(swc_ID(1:end-4));

    %         [Lia,Lib] = ismember(Inclusion_List,swc_ID);
    %         if sum(double(Lia)) == 0
    %             error("Not in inclusion List.")
    %         end

        swc_data = dlmread(files(I_Bytes(iter)).name, ' ', 7, 0);

    %     [Lia,Lib] = ismember(Inclusion_List,swc_ID);
    %     if sum(double(Lia)) == 0 || swc_data(1,2) ~= 1
    %         
    %     else

        %[Lia,Lib] = ismember(Inclusion_List,swc_ID);

        s = swc_data(:,1);
        t = swc_data(:,7);
        x = swc_data(:,3);
        y = swc_data(:,4);
        z = swc_data(:,5);
        widths = swc_data(:,6);

        [s_sorted,I] = sort(s);

        x = x(I) ./ 1000;
        y = y(I) ./ 1000;
        z = z(I) ./ 1000;
        widths = widths ./ 1000;

        [row,col] = find(t==-1);

        s(row,:) = [];
        t(row,:) = [];

    %     [row,col] = find(s > length(x));
    %     
    %     s(row,:) = [];
    %     t(row,:) = [];

        w = sqrt((x(s) - x(t)).^2 + (y(s) - y(t)).^2 + (z(s) - z(t)).^2);

        Node_List = sort(swc_data(1:end,1));

        G = graph(s,t,w);

        adjacency_matrix = adjacency(G);

        % Find all degree two nodes
        degree_list = degree(G);
        degree_1 = find(degree_list == 1);
        degree_3 = find(degree_list >= 3);

        % Initialize a cell array to store connected chains
        connected_chains = {};

    %     for i = 1:length(degree_1)
    %         for j = 1:length(degree_1)
    %             Node1 = Node_List(degree_1(i));
    %             Node2 = Node_List(degree_1(j));
    %             path = shortestpath(G,Node1,Node2);
    %             if sum(ismember(Node_List(degree_1),path)) == 2
    %                 if length(path) > 2
    %                     connected_chains{end+1} = path(2:end-1);
    %                 end
    %             end
    %         end
    %     end

        dist_array = distances(G,Node_List(degree_1),Node_List(degree_3));

        [M,I] = min(dist_array,[],2);

        for i = 1:length(M)
            if ~isinf(M(i))
                Node1 = Node_List(degree_1(i));
                Node2 = Node_List(degree_3(I(i)));
                path = shortestpath(G,Node1,Node2);
                connected_chains{end+1} = path;
            end
        end      

    %     for i = 1:length(degree_1)
    %         for j = 1:length(degree_3)
    %             Node1 = Node_List(degree_1(i));
    %             Node2 = Node_List(degree_3(j));
    %             path = shortestpath(G,Node1,Node2);
    %             if sum(horzcat(ismember(Node_List(degree_1),path)',ismember(Node_List(degree_3),path)')) == 2
    %                 if length(path) >= 2
    %                     connected_chains{end+1} = path;
    %                 end
    %             end
    %         end
    %     end

        dist_array = distances(G,Node_List(degree_3),Node_List(degree_3));

        dist_array(isinf(dist_array)) = 0;

        G_MST = graph(dist_array,'upper');

        T = minspantree(G_MST,'Type','forest');

        EL = T.Edges;

        EL = table2array(EL(:,1));

    %     for i = 1:length(degree_3)
    %         for j = i+1:length(degree_3)
    %             Node1 = Node_List(degree_3(i));
    %             Node2 = Node_List(degree_3(j));
    %             path = shortestpath(G,Node1,Node2);
    %             if sum(ismember(Node_List(degree_3),path)) == 2
    %                 if length(path) >= 2
    %                     connected_chains{end+1} = path;
    %                 end
    %             end
    %         end
    %     end

        for i=1:size(EL,1)
            Node1 = Node_List(degree_3(EL(i,1)));
            Node2 = Node_List(degree_3(EL(i,2)));
            path = shortestpath(G,Node1,Node2);
            connected_chains{end+1} = path;
        end

        if length(connected_chains) < 1

        else

            sub_chains = connected_chains;

    %         cd Meshes
    % 
    %         mesh_string = char(swc_string);
    %         mesh_string = mesh_string(1:end-4);
    %         mesh_string = strcat(string(mesh_string),".csv");
    % 
    %         mesh_points = readmatrix(mesh_string);
    % 
    %         disp(length(sub_chains))
    %         
    %         best_D = Inf(length(mesh_points),1);
    %         best_Idx = zeros(length(mesh_points),1);
    %         
    %         for j=1:length(sub_chains)
    %             disp(j)
    %             [Idx,D] = knnsearch(horzcat(x(sub_chains{j},1),y(sub_chains{j},1),z(sub_chains{j},1)),mesh_points);
    %             [row,col] = find(D < best_D);
    %             if ~isempty(row)
    %                 best_D(row) = D(row);
    %                 best_Idx(row) = j;
    %             end
    %         end

    %         min_values = best_D;
    %         min_indices = best_Idx;
    % 
    %         sub_mesh_points = cell(1,length(sub_chains));
    %         for i=1:length(sub_chains)
    %             [row,col] = find(min_indices == i);
    %             sub_mesh_points{i} = mesh_points(row,:);
    %         end

            linking_chain = [];
            sub_swc_points = cell(1,length(sub_chains));
            for i=1:length(sub_chains)
                sub_chain = sub_chains{i};
                if any(ismember(sub_chain,1))
                    linking_chain = [linking_chain;i];
                end
                sub_swc_points{i} = horzcat(x(sub_chain,1),y(sub_chain,1),z(sub_chain,1));
            end

            if length(linking_chain) >= 1

        %         cd Classifications
        % 
                class_string = char(swc_string);
                class_string = class_string(1:end-4);

                final_string = class_string;
        % 
        %         class_string = strcat("Classification_",string(class_string),".csv");
        % 
        %         classifications = readmatrix(class_string);
        % 
        %         sub_classifications = cell(1,length(sub_chains));
        %         for i=1:length(sub_chains)
        %             sub_chain = sub_chains{i};
        %             sub_classifications{i} = classifications(1,sub_chain)';
        %         end
        % 
        %         cd ..

                all_data = cell(1,4);

        %         all_data{1} = sub_classifications;
        %         all_data{2} = sub_mesh_points;
                all_data{1} = horzcat(x,y,z);
                all_data{2} = sub_chains;

                A_chains = zeros(length(sub_chains),length(sub_chains));

                chain_widths = zeros(length(sub_chains),1);
                chain_lengths = zeros(length(sub_chains),1);
                for chain_iter_i=1:1:length(sub_chains)
                    chain_i_width = mean(widths(sub_chains{chain_iter_i}));
                    sub_chain = sub_chains{chain_iter_i};
                    chain_points = horzcat(x(sub_chain),y(sub_chain),z(sub_chain));
                    chain_i_length = sum(sqrt(sum((chain_points(2:end,:) - chain_points(1:end-1,:)).^2,2)));
                    chain_widths(chain_iter_i) = chain_i_width;
                    chain_lengths(chain_iter_i) = chain_i_length;
                    for chain_iter_j=chain_iter_i+1:1:length(sub_chains)
                        chain_j_width = mean(widths(sub_chains{chain_iter_j}));
                        TF = any(ismember(sub_chains{chain_iter_i},sub_chains{chain_iter_j})); 
                        if TF
                            A_chains(chain_iter_i,chain_iter_j) = 1;
                            A_chains(chain_iter_j,chain_iter_i) = 1;
                        end
                    end
                end

                chain_betweenness = centrality(graph(A_chains),'betweenness');

                for chain_iter_i=1:1:length(sub_chains)
                    starting_index = 0;
                    if ismember(1,sub_chains{chain_iter_i})
                        starting_index = chain_iter_i;
                        break
                    end
                    for chain_iter_j=chain_iter_i+1:1:length(sub_chains)
                        TF = any(ismember(sub_chains{chain_iter_i},sub_chains{chain_iter_j})); 
                        if TF
                            weight1 = chain_betweenness(chain_iter_i) * (chain_widths(chain_iter_i)).^(2/1);
                            %weight1 = chain_betweenness(chain_iter_i);
                            weight2 = chain_betweenness(chain_iter_j) * (chain_widths(chain_iter_j)).^(2/1);
                            %weight2 = chain_betweenness(chain_iter_j);
                            A_chains(chain_iter_i,chain_iter_j) = weight1*weight2+1;
                            A_chains(chain_iter_j,chain_iter_i) = weight1*weight2+1;
                        end
                    end
                    if sum(length(nonzeros(A_chains(chain_iter_i,:)))) <= 1
                        A_chains(chain_iter_i,:) = zeros(1,length(chain_betweenness)); 
                        A_chains(:,chain_iter_i) = zeros(length(chain_betweenness),1); 
                    end
                end

                if starting_index~=0

                    if ~use_soma
                        [M,I] = max(chain_betweenness.*chain_widths.^2);
                    else
                        %[M,max_Idx] = max(chain_widths(linking_chain));
                        %I = linking_chain(max_Idx);
                        I = starting_index;
                    end

                    A_chains2 = A_chains;
                    A_chains2(A_chains2~=0) = 1;

    %                     priority_vector = chain_betweenness.*chain_widths.^2;

    %                     [B,I_sort] = sort(priority_vector,'descend');
    % 
    %                     priority_vector = priority_vector(I_sort);
    %                     A_chains2(:,:) = A_chains2(I_sort,I_sort);
    %                     chain_betweenness = chain_betweenness(I_sort);
    %                     chain_widths = chain_widths(I_sort);
    %                     chain_lengths = chain_lengths(I_sort);
    %                     sub_chains = sub_chains(I_sort);
    %                     sub_swc_points = sub_swc_points(I_sort);

    %                     [M,I2] = max(A_chains2(:,1).*(priority_vector));

    %                     A_row = zeros(1,size(A_chains2,1));
    %                     A_row(I2) = 1;
    %                     A_row(starting_index) = 1;
    %                     A_chains2(1,:) = A_row;
    %                     A_chains2(:,1) = A_row';

                    A_chains2 = circshift(A_chains2,[-starting_index+1,-starting_index+1]);
                    chain_betweenness = circshift(chain_betweenness,-starting_index+1);
                    chain_widths = circshift(chain_widths,-starting_index+1);
                    chain_lengths = circshift(chain_lengths,-starting_index+1);
                    sub_chains = circshift(sub_chains,-starting_index+1);
                    sub_swc_points = circshift(sub_swc_points,-starting_index+1);

                    segment_synapse_vector = zeros(length(sub_swc_points),1);
                    
                    [row,col] = find(presyn_id==swc_ID);
                    [row2,col2] = find(postsyn_id==swc_ID);
                    
                    if ~isempty(row) || ~isempty(row2)
                        
                        synapse_subset = presyn_pts(union(row,row2),:);

                        for i=1:length(sub_swc_points)
                            [Idx,D] = knnsearch(sub_swc_points{i},synapse_subset,'K',1);
                            if ~isempty(D)
                                segment_synapse_vector(i) = sum(D<=1);
                            end
                        end
                        
                    end
                    
                    [T,E] = bfsearch(graph(A_chains2),1,'edgetonew');

                    if ~isempty(T)
                        All_Chain_Data{mod(iter-first_iter+1,1000)+1,1} = chain_betweenness;
                        All_Chain_Data{mod(iter-first_iter+1,1000)+1,2} = chain_widths;
                        All_Chain_Data{mod(iter-first_iter+1,1000)+1,3} = chain_lengths;
                        All_Chain_Data{mod(iter-first_iter+1,1000)+1,4} = sub_chains;
                        All_Chain_Data{mod(iter-first_iter+1,1000)+1,5} = sub_swc_points;
                        All_Chain_Data{mod(iter-first_iter+1,1000)+1,6} = T;
                        All_Chain_Data{mod(iter-first_iter+1,1000)+1,7} = A_chains2;
                        All_Chain_Data{mod(iter-first_iter+1,1000)+1,8} = segment_synapse_vector;
                    end

                end
            end

        end
    
    end
    
    if mod(iter,1000) == 0
        
        cd ..\..\..
        cd Figure45_Data\Soma_Tree_Data
        
        if use_soma
            save(strcat(string(iter),"Soma_Tree_Data.mat"),'All_Chain_Data','-v7.3');
        else
            save(strcat(string(iter),"Dendrite_Tree_Data.mat"),'All_Chain_Data','-v7.3');
        end
        
        cd ..\..
        cd Datasets\FlyWire\Skeletons
        
        All_Chain_Data = cell(1000,8);
        
    end

end
