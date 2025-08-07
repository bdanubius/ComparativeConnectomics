clear
close all

%% Table Data

table_data = cell(7,1);

cd ..\..\..
cd ProcessedData\H01_104

H01_table = readtable("H01_AggregatedData.csv");
table_data{1} = H01_table;

cd ..
cd MM3

MM3_table = readtable("MM3_AggregatedData.csv");

% Setting synapses as weighted degree

MM3_table{:,3} = MM3_table{:,6};

%% MM3 Filtering

MM3_Inclusion_List = load("MM3_Inclusion_List.mat");
MM3_Inclusion_List = MM3_Inclusion_List.inclusion_list;

[Lia,Lib] = ismember(double(MM3_table{:,1}),double(MM3_Inclusion_List));

MM3_table = MM3_table(Lia,:);

%% Continue

table_data{2} = MM3_table;

cd ..
cd FlyWire

FlyWire_table = readtable("FlyWire_AggregatedData.csv");
table_data{3} = FlyWire_table;

cd ..
cd MANC

MANC_table = readtable("MANC_AggregatedData.csv");

MANC_Inclusion = load("MANC_Inclusion_List.csv");

[Lia,Lib] = ismember(double(MANC_table{:,1}),MANC_Inclusion);

MANC_table = MANC_table(Lia,:);

table_data{4} = MANC_table;

cd ..
cd Hemibrain

HEMI_table = readtable("Hemibrain_AggregatedData.csv");

HEMI_Exclusion = load("Hemibrain_Exclusion_List.csv");

[Lia,Lib] = ismember(double(HEMI_table{:,1}),HEMI_Exclusion);

HEMI_table = HEMI_table(~Lia,:);

table_data{5} = HEMI_table;

cd ..
cd Larva

Larva_table = readtable("Larva_AggregatedData.csv");
table_data{6} = Larva_table;

cd ..
cd Zebrafish

Zebrafish_table = readtable("Zebrafish_AggregatedData.csv");

Zebrafish_synapses = Zebrafish_table{:,3};

% Setting Synapses as weighted degree

Zebrafish_table{:,3} = Zebrafish_table{:,6};

table_data{7} = Zebrafish_table;

%% Continue

cd ..
cd H01_104

All_Density = cell(7,1);

Data = load("H01_Local_Density2.mat");
Data = Data.outdata;

All_Density{1} = Data;

cd ..
cd MM3

Data = load("MM3_Local_Density.mat");
Data = Data.outdata;

All_Density{2} = Data;

cd ..
cd FlyWire

Data = load("FLYWIRE_Local_Density.mat");
Data = Data.outdata;

All_Density{3} = Data;

cd ..
cd MANC

Data = load("MANC_Local_Density2.mat");
Data = Data.outdata;

All_Density{4} = Data;

cd ..
cd Hemibrain

Data = load("HEMI_Local_Density2.mat");
Data = Data.outdata;

All_Density{5} = Data;

cd ..
cd Larva

Data = load("LARVA_Local_Density.mat");
Data = Data.outdata;

All_Density{6} = Data;

cd ..
cd Zebrafish

Data = load("Zebrafish_Local_Density.mat");
Data = Data.outdata;

All_Density{7} = Data;

Density_Vectors = cell(7,2);

LD_corr = zeros(7,1);

Artificial_Dens_Data = cell(7,3);

for i=1:length(All_Density)
    
    Subdata = All_Density{i};
    SynapseDist_data = Subdata{3};
    LengthDist_data = Subdata{4};
    Synapse10_data = Subdata{1};
    Length10_data = Subdata{2};
    SubFilenames = Subdata{5};
    
    if i==1
        
        SubFilenames = extractBefore(SubFilenames,".");
        
        UniqueFilenames = unique(SubFilenames);
        
        Density10_vector = NaN(length(UniqueFilenames),1);
        
        for j=1:length(UniqueFilenames)
            
            disp(j)
            
            [Lia,Lib] = ismember(SubFilenames,UniqueFilenames(j));
            
            [row,col] = find(double(Lia));
            
            SubSubS = 0;
            SubSubL = 0;
            
            for k=1:length(row)
                
                if ~isnan(nanmean(Synapse10_data{row(k)})) && ~isnan(nanmean(Length10_data{row(k)}))
                
                    SubSubS = SubSubS + nanmean(Synapse10_data{row(k)});
                    SubSubL = SubSubL + nanmean(Length10_data{row(k)});
                    
                end
                
            end
            
            SubSubD = SubSubS ./ SubSubL;
            Density10_vector(j) = nanmean(SubSubD);
        end
        
        SubFilenames = UniqueFilenames;
        
    else
        
        Density10_vector = NaN(length(Synapse10_data),1);

        for j=1:length(Length10_data)
            disp(j)
            SubSubS = Synapse10_data{j};
            SubSubL = Length10_data{j};
            [row,col] = find(SubSubL > 10);
            SubSubD = (SubSubS(row)) ./ SubSubL(row);
            Density10_vector(j) = nanmean(SubSubD);
        end
    
    end
    
    SubFilenames(isnan(Density10_vector)) = [];
    Density10_vector(isnan(Density10_vector)) = [];
    
    Density_Vectors{i,1} = Density10_vector;
    
    DensityDist_vector = NaN(length(SynapseDist_data),1);
    
    for j=1:length(LengthDist_data)
        disp(j)
        SubSubL = LengthDist_data{j};
        [row,col] = find(SubSubL>0);
        DensityDist_vector(j) = nanmean(1 ./ SubSubL(row));
    end
    DensityDist_vector(isnan(DensityDist_vector)) = [];
    
    if i==1
        subtable = table_data{i};
        AllFileNames = double(subtable{:,1});
        AllLength = double(subtable{:,2});
        %SubFilenames = extractBefore(SubFilenames,".");
        [Lia,Lib] = ismember(double(SubFilenames),AllFileNames);
        SubDensCorr = Density10_vector(Lia);
        SubLength = AllLength(nonzeros(Lib));
        OutFileNames = AllFileNames(Lia);
        D_subdata = log10(SubDensCorr);
        L_subdata = log10(SubLength);
        [R,P] = corrcoef(D_subdata,L_subdata);
        LD_corr(i) = R(2,1);
    else
        subtable = table_data{i};
        AllFileNames = double(subtable{:,1});
        AllLength = double(subtable{:,2});
        [Lia,Lib] = ismember(double(SubFilenames),AllFileNames);
        SubDensCorr = Density10_vector(Lia);
        SubLength = AllLength(nonzeros(Lib));
        OutFileNames = AllFileNames(Lia);
        D_subdata = log10(SubDensCorr);
        L_subdata = log10(SubLength);
        [R,P] = corrcoef(D_subdata(~isnan(L_subdata)),L_subdata(~isnan(L_subdata)));
        LD_corr(i) = R(2,1);
        
    end
    
    Artificial_Dens_Data{i,1} = 10.^D_subdata(~isnan(L_subdata));
    Artificial_Dens_Data{i,2} = 10.^L_subdata(~isnan(L_subdata)) .* 10.^D_subdata(~isnan(L_subdata));
    Artificial_Dens_Data{i,3} = OutFileNames;
    
    Density_Vectors{i,2} = DensityDist_vector;
    
end

cd ..\..
cd FigureGeneration

save("Local_Dens_Data.mat",'Artificial_Dens_Data');
