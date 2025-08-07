clear

Micron_Tests = 1000;

Sample_Length = 10;

N_neurons = 2500;

cd ../..
cd Datasets/Larva/Skeletons

files = dir('*.swc');
filenames = string(zeros(size(files,1),1));
for i=1:size(files,1)
    disp(i)
    filename = files(i).name;
    filename = filename(1:end-4);
    filenames(i) = string(filename);
end

outdata = cell(5,1);

Synapse_Start = true;

all_walk_synapses = cell(N_neurons,1);
all_walk_lengths = cell(N_neurons,1);
for i=1:N_neurons
    disp(i)
    swc_data = load(strcat(filenames(i),".swc"));
    swc_data(:,3:5) = swc_data(:,3:5) ./ 1000;
    swc_keys = containers.Map(swc_data(:,1),linspace(1,size(swc_data,1),size(swc_data,1)));
    cd ..
    cd Synapses
    synapse_data = load(strcat(filenames(i),".csv"));
    synapse_data = synapse_data ./ 1000;
    if ~isempty(synapse_data)
        [Idx,D] = knnsearch(swc_data(:,3:5),synapse_data);
        synapse_idx = zeros(size(swc_data,1),1);
        synapse_idx = accumarray(Idx, 1, [length(synapse_idx), 1]);
        [synapse_row,col] = find(synapse_idx>=1);
        synapse_test_counts = NaN(Micron_Tests,1);
        length_test_counts = NaN(Micron_Tests,1);
        for j=1:Micron_Tests
            if ~Synapse_Start
                current_point = randi(size(swc_data,1));
            else
                if ~isempty(synapse_row)
                    current_point = synapse_row(randi(length(synapse_row)));
                else
                    break
                end
            end
            current_length = 0;
            current_synapses = 0;
            synapse_found = false;
            while current_length < Sample_Length
                current_synapses = current_synapses + synapse_idx(current_point);
                next_point = swc_data(current_point,7);
                if next_point == -1
                    break
                elseif synapse_idx(swc_keys(next_point)) == 1
                    current_length = current_length + sqrt(sum((swc_data(current_point,3:5) - swc_data(swc_keys(next_point),3:5)).^2,2));
                    current_point = swc_keys(next_point);
                else
                    current_length = current_length + sqrt(sum((swc_data(current_point,3:5) - swc_data(swc_keys(next_point),3:5)).^2,2));
                    current_point = swc_keys(next_point);
                end    
            end
            synapse_test_counts(j) = current_synapses;
            length_test_counts(j) = current_length;
        end
        all_walk_synapses{i} = synapse_test_counts;
        all_walk_lengths{i} = length_test_counts;
    end
    cd ..
    cd Skeletons
end

outdata{1,1} = all_walk_synapses;

outdata{2,1} = all_walk_lengths;

Synapse_Start = true;

all_walk_synapses = cell(N_neurons,1);
all_walk_lengths = cell(N_neurons,1);
for i=1:N_neurons
    disp(i)
    swc_data = load(strcat(filenames(i),".swc"));
    swc_data(:,3:5) = swc_data(:,3:5) ./ 1000;
    swc_keys = containers.Map(swc_data(:,1),linspace(1,size(swc_data,1),size(swc_data,1)));
    cd ..
    cd Synapses
    synapse_data = load(strcat(filenames(i),".csv"));
    synapse_data = synapse_data ./ 1000;
    if ~isempty(synapse_data)
        [Idx,D] = knnsearch(swc_data(:,3:5),synapse_data);
        synapse_idx = zeros(size(swc_data,1),1);
        synapse_idx = accumarray(Idx, 1, [length(synapse_idx), 1]);
        [synapse_row,col] = find(synapse_idx>=1);
        synapse_test_counts = NaN(Micron_Tests,1);
        length_test_counts = NaN(Micron_Tests,1);
        for j=1:Micron_Tests
            if ~Synapse_Start
                current_point = randi(size(swc_data,1));
            else
                if ~isempty(synapse_row)
                    current_point = synapse_row(randi(length(synapse_row)));
                else
                    break
                end
            end
            current_length = 0;
            current_synapses = 0;
            synapse_found = false;
            while ~synapse_found && current_length < Sample_Length
                current_synapses = current_synapses + synapse_idx(current_point);
                next_point = swc_data(current_point,7);
                if next_point == -1
                    break
                elseif synapse_idx(swc_keys(next_point)) == 1
                    synapse_found = true;
                    break
                else
                    current_length = current_length + sqrt(sum((swc_data(current_point,3:5) - swc_data(swc_keys(next_point),3:5)).^2,2));
                    current_point = swc_keys(next_point);
                end    
            end
            if synapse_found
                synapse_test_counts(j) = current_synapses;
                length_test_counts(j) = current_length;
            end
        end
        all_walk_synapses{i} = synapse_test_counts;
        all_walk_lengths{i} = length_test_counts;
    end
    cd ..
    cd Skeletons
end

outdata{3,1} = all_walk_synapses;

outdata{4,1} = all_walk_lengths;

outdata{5,1} =  filenames(1:N_neurons);

cd ..\..\..
cd ProcessedData\Larva

save("LARVA_Local_Density.mat",'outdata');

