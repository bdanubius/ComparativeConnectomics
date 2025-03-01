skeleton_cell = {};

counter = 0;

cd ..
cd ..
cd Datasets/H01_104/Skeletons

files = dir('*.swc');

filenames = [];
IDs_104 = [];
for i=1:length(files)
    filenames = [filenames;string(files(i).name)];
    IDs_104 = [IDs_104;extractBefore(string(files(i).name),".")];
end

IDs_104 = unique(IDs_104);
filename_sets = cell(length(IDs_104),1);

for i=1:length(IDs_104)
    filename_set = [];
    for j=1:length(filenames)
        idx = strfind(filenames(j),IDs_104(i));
        if ~isempty(idx)
            filename_set = [filename_set;filenames(j)];
        end
    end
    filename_sets{i} = filename_set;
end

cd ..

for iter=1:length(IDs_104)

    counter = counter + 1;

    disp(iter)
    
    filename_set = filename_sets{iter};
    
    all_data = cell(1,4);

    all_sub_chains = {};

    all_points = [];
    
    try
        
        disp(length(filename_set))
        
        for iter2=1:length(filename_set)
            
            disp(iter2)

            swc_string = filename_set(iter2);

            swc_data = importdata(swc_string);

            s = swc_data(:,1);
            t = swc_data(:,7);
            x = swc_data(:,3) .* 32;
            y = swc_data(:,4) .* 32;
            z = swc_data(:,5) .* 33;

            [s_sorted,I] = sort(s);

            x = x(I);
            y = y(I);
            z = z(I);

            [row,col] = find(t==-1);

            s(row,:) = [];
            t(row,:) = [];
            
            s = s+1;
            t = t+1;

            w = sqrt((x(s) - x(t)).^2 + (y(s) - y(t)).^2 + (z(s) - z(t)).^2);

            Node_List = sort(swc_data(1:end,1)) + 1;

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

            if isempty(connected_chains)

            else

                sub_chains = connected_chains;

                sub_swc_points = cell(1,length(sub_chains));
                for i=1:length(sub_chains)
                    sub_chain = sub_chains{i};
                    sub_swc_points{i} = horzcat(x(sub_chain,1),y(sub_chain,1),z(sub_chain,1));
                end
                
                class_string = char(swc_string);
                class_string = class_string(1:end-4);

                final_string = class_string;
                
                if ~isempty(all_points)
                    
                    for iter3 = 1:length(sub_chains)
                        sub_chains{iter3} = sub_chains{iter3} + size(all_points,1);
                    end
                    
                    all_chains = horzcat(all_chains,sub_chains);
                    
                else
                    
                    all_chains = sub_chains;
                    
                end
                
                all_points = [all_points;horzcat(x,y,z)];

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

            end
            
        end
        
        all_data{1} = all_points;
        
        all_data{2} = all_chains;
        
        all_data{3} = [];
        
        all_data{4} = IDs_104(iter);
        
        skeleton_cell = [skeleton_cell;all_data];
        
    catch
        
        skeleton_cell = [skeleton_cell;num2cell(NaN(1,4))];
    
    end
     
    if iter == length(IDs_104)
        
        cd ..
        cd ..
        cd ..

        cd ProcessedData/H01_104/Cleaned_Data

        save(strcat("Skeletons_1000.mat"),'skeleton_cell');
        
        cd ..
        cd .. 
        cd ..
        cd Datasets/H01_104/Skeletons
        
    end

end
        
close all