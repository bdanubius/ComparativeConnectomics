clear

skeleton_cell = {};

counter = 0;

bytes_sort = 1;

cd ..
cd ..
cd Datasets/FlyWire/Skeletons

if bytes_sort

    files = dir('*.swc');

    bytes_array = zeros(size(files,1),1);
    for i=1:length(files)
        bytes_array(i) = files(i).bytes;
    end

    [B,I_Bytes] = sort(bytes_array,'ascend');

end

for iter=1:1:129278

    counter = counter + 1;

    disp(iter)
    
    try

        swc_string = string(files(I_Bytes(iter)).name);
        
        swc_ID = char(swc_string);
        swc_ID = string(swc_ID(1:end-4));

        swc_data = dlmread(files(I_Bytes(iter)).name, ' ', 7, 0);

        s = swc_data(:,1);
        t = swc_data(:,7);
        x = swc_data(:,3);
        y = swc_data(:,4);
        z = swc_data(:,5);

        [s_sorted,I] = sort(s);

        x = x(I);
        y = y(I);
        z = z(I);

        [row,col] = find(t==-1);

        s(row,:) = [];
        t(row,:) = [];

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

        dist_array = distances(G,Node_List(degree_3),Node_List(degree_3));

        dist_array(isinf(dist_array)) = 0;

        G_MST = graph(dist_array,'upper');

        T = minspantree(G_MST,'Type','forest');

        EL = T.Edges;

        EL = table2array(EL(:,1));

        for i=1:size(EL,1)
            Node1 = Node_List(degree_3(EL(i,1)));
            Node2 = Node_List(degree_3(EL(i,2)));
            path = shortestpath(G,Node1,Node2);
            connected_chains{end+1} = path;
        end

        if length(connected_chains) < 1
            
        else

            sub_swc_points = cell(1,length(sub_chains));
            for i=1:length(sub_chains)
                sub_chain = sub_chains{i};
                sub_swc_points{i} = horzcat(x(sub_chain,1),y(sub_chain,1),z(sub_chain,1));
            end

            class_string = char(swc_string);
            class_string = class_string(1:end-4);

            final_string = class_string;

            all_data = cell(1,4);
            
            all_data{1} = horzcat(x,y,z);
            all_data{2} = sub_chains;

            A_chains = zeros(length(sub_chains),length(sub_chains));

            for chain_iter_i=1:1:length(sub_chains)
                for chain_iter_j=chain_iter_i+1:1:length(sub_chains)
                    TF = any(ismember(sub_chains{chain_iter_i},sub_chains{chain_iter_j})); 
                    if TF
                        A_chains(chain_iter_i,chain_iter_j) = 1;
                        A_chains(chain_iter_j,chain_iter_i) = 1;
                    end
                end
            end

            all_data{3} = sparse(A_chains);

            all_data{4} = string(final_string);

            skeleton_cell = [skeleton_cell;all_data];

        end
        
    catch
        
        skeleton_cell = [skeleton_cell;num2cell(NaN(1,4))];
    
    end
     
    if mod(iter,1000) == 0 || iter == size(files,1)
        
        if counter > 1

            cd ..
            cd ..
            cd ..

            cd ProcessedData/FlyWire/Cleaned_Data

            save(strcat("Skeletons_",string(iter),".mat"),'skeleton_cell');
            
            skeleton_cell = {};
            
            cd ..
            cd .. 
            cd ..
            cd Datasets/FlyWire/Skeletons
            
        end
        
    end

end
        
close all