clear

cd ..
cd ..
cd Datasets\MANC

MANC_labels = readtable('MANC_ROIs.csv');

Body_ID_List = string(MANC_labels{:,2});

all_row = [];

[row,col] = find(string(MANC_labels{:,12})=="Prelim Roughly traced");

all_row = [all_row;row];

[row,col] = find(string(MANC_labels{:,12})=="0.5assign");

all_row = [all_row;row];

[row,col] = find(string(MANC_labels{:,12})=="Orphan");

all_row = [all_row;row];

[row,col] = find(string(MANC_labels{:,12})=="Sensory Anchor");

all_row = [all_row;row];

[row,col] = find(string(MANC_labels{:,12})=="RT Orphan");

all_row = [all_row;row];

[row,col] = find(string(MANC_labels{:,12})=="PRT Orphan");

all_row = [all_row;row];

Body_ID_List(all_row,:) = [];

cd ..
cd ..
cd ProcessedData\MANC

writematrix(Body_ID_List,"MANC_Inclusion_List.csv");