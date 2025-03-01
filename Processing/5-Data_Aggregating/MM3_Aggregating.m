clear

cd ..
cd ..
cd ProcessedData\MM3

MM3_Data = load("Geometry_Data.mat");
MM3_Data = MM3_Data.synapse_density_data;
MM3_Data = MM3_Data{1};

MM3_Degree_Data = load("MM3_Degree_Data.mat");
MM3_Degree_Data = MM3_Degree_Data.MM3_Degrees;

[Lia,Lib] = ismember(double(MM3_Data(:,4)),MM3_Degree_Data(:,1));

MM3_Degree_Data = double(MM3_Degree_Data(Lib,:));
MM3_Data = double(MM3_Data);

MM3_Total_Data = horzcat(MM3_Data(:,4),MM3_Data(:,1:2),MM3_Data(:,2)./MM3_Data(:,1),MM3_Degree_Data(:,2:3));

A = cellstr(horzcat(string(int64(MM3_Total_Data(:,1))),string(MM3_Total_Data(:,2:end))));

tC=cell2table(A,'VariableNames',{'CellID','Length','Synapses','Density','Degree','WeightedDegree'});

writetable(tC,"MM3_AggregatedData.csv");
