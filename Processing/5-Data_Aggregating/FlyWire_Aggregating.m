clear

cd ..
cd ..
cd ProcessedData\FlyWire

Chain_Names = load("Filenames.mat");
Chain_Names = Chain_Names.Chain_Names;

cd ..
cd ..
cd Datasets\FlyWire\Skeletons

neuron_lengths = zeros(length(Chain_Names),1);

for i=1:length(Chain_Names)
    
    disp(i)

    swc_data = dlmread(Chain_Names(i), ' ', 7, 0);

    swc_data(:,3:5) = swc_data(:,3:5) ./ 1000;

    neuron_length = sum(sqrt(sum((swc_data(2:end,3:5) - swc_data(swc_data(2:end,7),3:5)).^2,2)));

    neuron_lengths(i) = neuron_length;
    
end

cd ..
cd ..
cd ..
cd ProcessedData\FlyWire

Neuron_Lengths = neuron_lengths;
Chain_Names = extractBefore(Chain_Names,".swc");

HEMI_Data = load("Geometry_Data.mat");
HEMI_Data = HEMI_Data.synapse_density_data;
HEMI_Data = HEMI_Data{1};

HEMI_Degree_Data = load("FlyWire_Degree_Data.mat");
HEMI_Degree_Data = HEMI_Degree_Data.FlyWire_Degrees;

[Lia,Lib] = ismember(string(HEMI_Data(:,4)),Chain_Names);

HEMI_Data2 = string(NaN(length(Chain_Names),6));

for i=1:length(Lib)
    if Lib(i)==0
    else
        HEMI_Data2(Lib(i),:) = HEMI_Data(i,:);
    end
end

[Lia,Lib] = ismember(string(HEMI_Degree_Data(:,1)),Chain_Names);

HEMI_Degree_Data2 = string(NaN(length(Chain_Names),3));

for i=1:length(Lib)
    if Lib(i)==0
    else
        HEMI_Degree_Data2(Lib(i),:) = HEMI_Degree_Data(i,:);
    end
end

HEMI_Total_Data = horzcat(Chain_Names,string(Neuron_Lengths),string(HEMI_Degree_Data2(:,3)),string(double(HEMI_Degree_Data2(:,3))./double(Neuron_Lengths)),string(HEMI_Degree_Data2(:,2:3)));

A = cellstr(string(HEMI_Total_Data));

tC=cell2table(A,'VariableNames',{'CellID','Length','Synapses','Density','Degree','WeightedDegree'});

writetable(tC,"FlyWire_AggregatedData.csv");
