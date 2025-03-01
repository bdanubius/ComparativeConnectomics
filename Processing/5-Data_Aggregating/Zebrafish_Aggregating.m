clear

cd ..
cd ..
cd ProcessedData\Zebrafish

Zebrafish_Data = load("Geometry_Data.mat");
Zebrafish_Data = Zebrafish_Data.synapse_density_data;
Zebrafish_Data = Zebrafish_Data{1};

Zebrafish_Degree_Data = load("Zebrafish_Degree_Data.mat");
Zebrafish_Degree_Data = Zebrafish_Degree_Data.Zebrafish_Degrees;

[Lia,Lib] = ismember(double(Zebrafish_Data(:,4)),Zebrafish_Degree_Data(:,1));

Zebrafish_Degree_Data = double(Zebrafish_Degree_Data(Lib,:));
Zebrafish_Data = double(Zebrafish_Data);

Zebrafish_Total_Data = horzcat(Zebrafish_Data(:,4),Zebrafish_Data(:,1:2),Zebrafish_Data(:,2)./Zebrafish_Data(:,1),Zebrafish_Degree_Data(:,2:3));

A = cellstr(horzcat(string(int64(Zebrafish_Total_Data(:,1))),string(Zebrafish_Total_Data(:,2:end))));

tC=cell2table(A,'VariableNames',{'CellID','Length','Synapses','Density','Degree','WeightedDegree'});

writetable(tC,"Zebrafish_AggregatedData.csv");
