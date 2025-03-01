clear

cd ..
cd ..
cd ProcessedData\H01_104

H01_Data = load("Geometry_Data.mat");
H01_Data = H01_Data.synapse_density_data;
H01_Data = H01_Data{1};

H01_Degree_Data = load("H01_Degree_Data.mat");
H01_Degree_Data = H01_Degree_Data.H01_Degrees;

[Lia,Lib] = ismember(double(H01_Data(:,4)),H01_Degree_Data(:,1));

H01_Total_Data = horzcat(H01_Data(:,4),double(H01_Data(:,1))./1000,H01_Data(:,2),1000.*double(H01_Data(:,2))./double(H01_Data(:,1)),string(NaN(size(H01_Data,1),2)));

for i=1:length(Lib)
    
    if Lib(i) ~= 0
        
        H01_Total_Data(i,5:6) = H01_Degree_Data(Lib(i),2:3);
        
    end
    
end

A = cellstr(H01_Total_Data);

tC=cell2table(A,'VariableNames',{'CellID','Length','Synapses','Density','Degree','WeightedDegree'});

writetable(tC,"H01_AggregatedData.csv");
