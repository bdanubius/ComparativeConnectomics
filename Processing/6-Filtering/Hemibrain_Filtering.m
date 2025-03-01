clear

cd ..
cd ..
cd Datasets\Hemibrain

hemibrain_labels = readtable('Hemibrain_ROIs.csv');

Hemibrain_Exclusion = [];

Body_ID_List = string(hemibrain_labels{:,2});

[row,col] = find(string(hemibrain_labels{:,13})=="Leaves");

Hemibrain_Exclusion = [Hemibrain_Exclusion;Body_ID_List(row)];

[row,col] = find(string(hemibrain_labels{:,13})=="Orphan");

Hemibrain_Exclusion = [Hemibrain_Exclusion;Body_ID_List(row)];

[row,col] = find(string(hemibrain_labels{:,13})=="Orphan hotknife");

Hemibrain_Exclusion = [Hemibrain_Exclusion;Body_ID_List(row)];

[row,col] = find(string(hemibrain_labels{:,13})=="Orphan-artifact");

Hemibrain_Exclusion = [Hemibrain_Exclusion;Body_ID_List(row)];

cd ..
cd ..
cd ProcessedData\Hemibrain

writematrix(Hemibrain_Exclusion,"Hemibrain_Exclusion_List.csv");