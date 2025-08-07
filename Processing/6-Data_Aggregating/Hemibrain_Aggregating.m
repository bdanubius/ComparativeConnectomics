clear

cd ..
cd ..
cd ProcessedData\Hemibrain

HEMI_Data = load("Geometry_Data.mat");
HEMI_Data = HEMI_Data.synapse_density_data;
HEMI_Data = HEMI_Data{1};

HEMI_Data = double(HEMI_Data);

HEMI_Total_Data = horzcat(HEMI_Data(:,4),HEMI_Data(:,1)./125,HEMI_Data(:,2),125.*HEMI_Data(:,2)./HEMI_Data(:,1));

HEMI_Total_Data = horzcat(HEMI_Total_Data,zeros(size(HEMI_Total_Data,1),2));

cd ..
cd ..
cd Datasets\Hemibrain\Synapses

for i=1:size(HEMI_Total_Data)
    disp(i)
    try
        subtable = readtable(strcat(string(HEMI_Total_Data(i,1)),".txt"));
        HEMI_Total_Data(i,6) = size(subtable,1);
        HEMI_Total_Data(i,5) = length(unique(double(subtable{:,3})));
    end
end

cd ..
cd ..
cd ..
cd ProcessedData\Hemibrain

HEMI_Total_Data = num2cell(HEMI_Total_Data);

tC=cell2table(HEMI_Total_Data,'VariableNames',{'CellID','Length','Synapses','Density','Degree','Weighted_Degree'});

writetable(tC,"Hemibrain_AggregatedData.csv");