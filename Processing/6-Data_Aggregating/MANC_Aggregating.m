clear

cd ..
cd ..
cd ProcessedData\MANC

MANC_Data = load("Geometry_Data.mat");
MANC_Data = MANC_Data.synapse_density_data;
MANC_Data = MANC_Data{1};

MANC_Data = double(MANC_Data);

MANC_Total_Data = horzcat(MANC_Data(:,4),MANC_Data(:,1)./125,MANC_Data(:,2),125.*MANC_Data(:,2)./MANC_Data(:,1));

MANC_Total_Data = horzcat(MANC_Total_Data,zeros(size(MANC_Total_Data,1),2));

cd ..
cd ..
cd Datasets\MANC\Synapses

for i=1:size(MANC_Total_Data)
    disp(i)
    try
        subtable = readtable(strcat(string(MANC_Total_Data(i,1)),".txt"));
        MANC_Total_Data(i,6) = size(subtable,1);
        MANC_Total_Data(i,5) = length(unique(double(subtable{:,3})));
    end
end

cd ..
cd ..
cd ..
cd ProcessedData\MANC

MANC_Total_Data = num2cell(MANC_Total_Data);

tC=cell2table(MANC_Total_Data,'VariableNames',{'CellID','Length','Synapses','Density','Degree','Weighted_Degree'});

writetable(tC,"MANC_AggregatedData.csv");
