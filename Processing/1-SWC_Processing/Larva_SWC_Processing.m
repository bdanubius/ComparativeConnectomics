clear

% Define range of neuron skeletons to process
skeletons_idx = 1:1:3066; % Total 3066 entries

skeleton_cell = {}; % Cell array to store processed skeleton data
counter = 0; % Counter for processed files

bytes_sort = 1; % Flag to sort files by size

% Change directory to dataset folder
cd ..
cd ..
cd Datasets/Larva/Skeletons

% Get list of all .swc files and sort them by file size
files = dir('*.swc');
if bytes_sort
    bytes_array = [files.bytes];
    [B, I_Bytes] = sort(bytes_array);
end

for iter=skeletons_idx
    counter = counter + 1;
    disp(iter) % Display current index
    
    try

        % Get filename of current .swc file
        if bytes_sort
            swc_string = string(files(I_Bytes(iter)).name);
        else
            swc_string = string(files(I_Bytes(iter)).name);
        end
            
        % Load SWC data (Neuron morphology)
        swc_data = load(swc_string);
        
        % Extract node IDs, parent connections, and 3D coordinates
        s = swc_data(:,1); % Node indices
        t = swc_data(:,7); % Parent indices
        x = swc_data(:,3); % X-coordinates
        y = swc_data(:,4); % Y-coordinates
        z = swc_data(:,5); % Z-coordinates
        
        % Map node indices to their respective parent indices
        [Lia, Lib] = ismember(s, t);
        
        % Adjust indices to a sequential format
        s = (1:length(s))';
        t = Lib;
        
        s_original = s; % Store original node indices
        
        % Remove root node connections (parent = 0)
        [row, ~] = find(t == 0);
        if ~isempty(row)
            s(row) = [];
            t(row) = [];
        end
        
        % Compute edge weights (Euclidean distance between nodes)
        w = sqrt((x(s) - x(t)).^2 + (y(s) - y(t)).^2 + (z(s) - z(t)).^2);
        
        Node_List = s_original; % Store node list
        
        % Construct a weighted graph
        G = graph(s, t, w);
        
        % Compute adjacency matrix
        adjacency_matrix = adjacency(G);
        
        % Identify terminal (degree-1) and branch (degree>=3) nodes
        degree_list = degree(G);
        degree_1 = find(degree_list == 1);
        degree_3 = find(degree_list >= 3);
        
        % Initialize a cell array for connected chains
        connected_chains = {};
        
        % Compute shortest paths from terminals to branch points
        dist_array = distances(G, Node_List(degree_1), Node_List(degree_3));
        [M, I] = min(dist_array, [], 2);
        
        % Extract paths from terminal to nearest branch node
        for i = 1:length(M)
            if ~isinf(M(i))
                path = shortestpath(G, Node_List(degree_1(i)), Node_List(degree_3(I(i))));
                connected_chains{end+1} = path;
            end
        end
        
        % Find paths between all terminal and branch nodes
        for i = 1:length(degree_1)
            for j = 1:length(degree_3)
                path = shortestpath(G, Node_List(degree_1(i)), Node_List(degree_3(j)));
                if sum(ismember(Node_List(degree_1), path)) + sum(ismember(Node_List(degree_3), path)) == 2
                    if length(path) >= 2
                        connected_chains{end+1} = path;
                    end
                end
            end
        end
        
        % Find paths between all branch nodes
        for i = 1:length(degree_3)
            for j = i+1:length(degree_3)
                path = shortestpath(G, Node_List(degree_3(i)), Node_List(degree_3(j)));
                if sum(ismember(Node_List(degree_3), path)) == 2 && length(path) >= 2
                    connected_chains{end+1} = path;
                end
            end
        end
        
        % Compute distances between all terminal nodes
        D_matrix = distances(G);
        for i = 1:length(degree_1)
            for j = i+1:length(degree_1)
                if D_matrix(Node_List(degree_1(i)), Node_List(degree_1(j)))
                    path = shortestpath(G, Node_List(degree_1(i)), Node_List(degree_1(j)));
                    if sum(ismember(Node_List(degree_1), path)) == 2 && length(path) >= 2
                        connected_chains{end+1} = path;
                    end
                end
            end
        end
        
        % Compute minimum spanning tree of branch nodes
        dist_array = distances(G, Node_List(degree_3), Node_List(degree_3));
        dist_array(isinf(dist_array)) = 0;
        G_MST = graph(dist_array, 'upper');
        T = minspantree(G_MST, 'Type', 'forest');
        
        % Extract edges of the spanning tree and store paths
        EL = table2array(T.Edges(:,1));
        for i = 1:size(EL,1)
            path = shortestpath(G, Node_List(degree_3(EL(i,1))), Node_List(degree_3(EL(i,2))));
            connected_chains{end+1} = path;
        end
        
        % Store processed data
        if ~isempty(connected_chains)
            sub_swc_points = cellfun(@(p) [x(p), y(p), z(p)], connected_chains, 'UniformOutput', false);
            all_data = {horzcat(x, y, z), connected_chains, sparse(zeros(length(connected_chains))), string(swc_string(1:end-4))};
            skeleton_cell = [skeleton_cell; all_data];
        else
            skeleton_cell = [skeleton_cell; num2cell(NaN(1,4))];
        end
        
    catch
        skeleton_cell = [skeleton_cell; num2cell(NaN(1,4))]; % Store NaN for failed files
    end
    
    % Save results every 1000 iterations
    if mod(iter, 1000) == 0 || iter == size(files, 1)
        if counter > 1
            cd ../../..
            cd ProcessedData/Larva/Cleaned_Data
            save(strcat("Skeletons_", string(iter), ".mat"), 'skeleton_cell');
            skeleton_cell = {}; % Reset data storage
            cd ../../..
            cd Datasets/Larva/Skeletons
        end
    end
end

close all