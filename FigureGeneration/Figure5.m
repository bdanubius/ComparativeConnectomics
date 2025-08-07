clear
close all

% Set default font for axes
set(groot, 'DefaultAxesFontName', 'Helvetica');
% Set default font for text in figures
set(groot, 'DefaultTextFontName', 'Helvetica');
% Set default font for legends
set(groot, 'DefaultLegendFontName', 'Helvetica');

% Set default font for axes
set(groot, 'DefaultAxesFontSize', 20);
set(groot, 'DefaultAxesTickLabelInterpreter', 'tex');
% Set default font for text in figures
set(groot, 'DefaultTextFontSize', 20);
set(groot, 'DefaultTextInterpreter', 'tex');
% Set default font for legends
set(groot, 'DefaultLegendFontSize', 20);
set(groot, 'DefaultLegendInterpreter', 'tex');

% Set default font size for axes titles (relative to axes font size)
set(groot, 'DefaultAxesLabelFontSizeMultiplier', 1.5)
set(groot, 'DefaultAxesTitleFontSizeMultiplier', 1.5);

set(groot, 'DefaultAxesLineWidth', 0.5)

%%  Figure 5D

cd ..
cd Figure45_Data

Chain_Names = load("Chain_Names.mat");
Chain_Names = Chain_Names.Chain_Names;
Chain_Names = extractBefore(Chain_Names,".swc");

outdata = load("Xn_Tree_Data2.mat");
outdata = outdata.outdata;

tree_sizes = outdata{1};
survival_probability = outdata{2};
tree_depth = outdata{3};
max_tree_size = outdata{4};
xn = outdata{5};
mu = outdata{6};
survival_probability2 = outdata{7};
cell_labels = outdata{8};

opts = detectImportOptions("FlyWire_CellTypeData.csv");
opts = setvartype(opts,'CellID', 'string');
FlyWire_Table = readtable("FlyWire_CellTypeData.csv",opts);

FlyWire_DepthLength = readtable("FlyWireDepthLengthTable4.csv");

FlyWire_Labels = string(FlyWire_Table{:,1});
FlyWire_CellLabels = string(FlyWire_Table{:,10});

%Kenyon_Cell
%mAL
%DAN
%LHCENT

% This part is reading through the depth / length / p data stored for each
% cell type, and storing it in a cell structure

Cell_Type_Array = unique(FlyWire_CellLabels);
Cell_Type_Array = Cell_Type_Array(2:end);

Cell_Type_Array = ["ALL","visual","central","Kenyon_Cell",'optic_lobes',"gustatory"];

Cell_Type_Data = cell(1,length(Cell_Type_Array));

for iter = 1:length(Cell_Type_Array)
    
    disp(iter)
    
    Cell_Type = Cell_Type_Array(iter);

    if Cell_Type == "ALL"

        subtable = FlyWire_Table;

        FlyWire_Sublabels = string(subtable{:,1});

    else

        [row,col] = find(FlyWire_CellLabels == Cell_Type);

        subtable = FlyWire_Table(row,:);

        FlyWire_Sublabels = string(subtable{:,1});

    end

    [Lia,Lib] = ismember(FlyWire_Sublabels,cell_labels);

    xn_subdata = xn(:,Lib(Lib>0));

    p_subdata = cell(1,length(xn_subdata));

    for i=1:length(xn_subdata)
        xn_segment = xn_subdata{i};
        p_subdata{i} = xn_segment(2:end) ./ (xn_segment(1:end-1).*2);
    end

    super_p_vector = horzcat(p_subdata{:})';
    p_idx_vector = zeros(length(super_p_vector),1);

    counter = 0;
    for i=1:length(p_subdata)
        subidx = linspace(1,length(p_subdata{i}),length(p_subdata{i}))';
        p_idx_vector(1+counter:length(p_subdata{i})+counter,:) = subidx;
        counter = counter + length(p_subdata{i}); 
    end

    layer_mean = [];
    layer_std = [];
    layer_counts = [];
    depth_histogram = [];
    maxidx_histogram = [];
    maxwidth_histogram = [];

    for i=1:length(p_subdata)
        depth_histogram = [depth_histogram;length(p_subdata{i})];
        [M,I] = max(xn_subdata{i});
        maxidx_histogram = [maxidx_histogram;I(1)];
        maxwidth_histogram = [maxwidth_histogram;M(1)];
    end

    for i=1:max(p_idx_vector)
        [row,col] = find(p_idx_vector==i);
        if ~isempty(row)
            layer_mean = [layer_mean;mean(super_p_vector(row))];
            layer_std = [layer_std;std(super_p_vector(row))];
            layer_counts = [layer_counts;length(row)];
        end
    end
    
    [row,col] = find(FlyWire_CellLabels==Cell_Type);
    
    length_histogram = FlyWire_Table{row,2};
    
    if ~isempty(length_histogram)
        M = median(length_histogram);
        [subidx_row,col] = find(round(length_histogram)==round(M));
        %disp(string(FlyWire_Table{row(subidx_row(1)),1}))
    end
    
    length_histogram = double(string(length_histogram(~ismissing(length_histogram))));
    
    length_histogram = log10(length_histogram);
    
    length_histogram(isinf(length_histogram),:) = [];
    
    outdata = cell(6,1);
    
    outdata{1} = layer_mean;
    outdata{2} = layer_std;
    outdata{3} = layer_counts;
    outdata{4} = depth_histogram;
    outdata{5} = maxidx_histogram;
    outdata{6} = length_histogram;
    outdata{7} = maxwidth_histogram;
    
    Cell_Type_Data{iter} = outdata;
    
end

Colors =          [0 0 0;
    0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840];

dataset_shuffle_vector = [8,4,7,3,5,1,6,2];

color_codes = ["2EAA90","E3242B","1F5A66","734F96","88C4D8","F7D42F","0855B1","D341D1"]';

colors = [];
for i=1:length(color_codes)
    str = char(color_codes(i));
    colors = [colors;sscanf(str(1:end),'%2x%2x%2x',[1 3])/255];
end

colors = colors(dataset_shuffle_vector,:);

colors = colors([1,2,3,4,5,6,7,8],:);

colors(1:5,:) = ([115, 80, 151; 238, 135, 40; 124, 203, 221; 134, 199, 108; 92, 107, 178] ./ 255);

p_STD = [0.695535429977907,0.252369877977962,0.232751151497085,0.230493016454119,0.210304719106992,0.197644394583968,0.190165411468545,0.178989864028709,0.171309997463044,0.169098270273740,0.164467973191630,0.160601196822332,0.158000888224647,0.160994697593523,0.159155566182279,0.159322243080546,0.160539479982327,0.163667752847260,0.165796622040120,0.162293456305564,0.168008039550561,0.164836863441451,0.166929460505272,0.167149846848826,0.167481573267498,0.167358148839206,0.166212151750443,0.168959990612990,0.168362514065006,0.171431308245077,0.172142518608431,0.175157737770696,0.171850902013596,0.176836367455948,0.177543075804106,0.179406264887942,0.178651772175485,0.183146817275179,0.182202404021374,0.184298654051057,0.188396290248425,0.186005005736119,0.184604779940529,0.188148607357319,0.187412940777900,0.191713330994320,0.188728633336184,0.192798670190682,0.191915523791876,0.189036332698515,0.202217072883747,0.192261028966874,0.197006098530451,0.194843870498071,0.197845803160546,0.200774165102795,0.199661300216890,0.199197939951150,0.196654008552513,0.190605443921816,0.197917394920559,0.192370846745879,0.191011047339573,0.192295003883379,0.194257295374166,0.175149948644691,0.193267487895270,0.191590639223255,0.182995617703511,0.194194654807931,0.161635164232321,0.188704100591559,0.198081304911207,0.194632595649882,0.194534990519121,0.209029817694213,0.185893298844860,0.173916993687399,0.195048352969420,0.197653154648074,0.222186396334869,0.163759561542535,0.165083865612491,0.131937557733853,0.173095722679070,0.157664456069401,0.0994825969381712,0.171866620432047,0.182611023948511,0.191243426235958,0.197223187657615,0.253599489659238,0.132779733445909,0.288385710802836,0.192239190026210,0.272540453648098,0.177895082305971,0.247494628577456,0.210184335479819,0.282270465443184];

figure('Renderer', 'painters', 'Position', [10 10 900 600])
hold on
DL_labels = string(FlyWire_DepthLength{:,1});
N_array = zeros(length(Cell_Type_Data),1);
for i=1:length(Cell_Type_Data)
    subdata = Cell_Type_Data{i};
    counts = subdata{3};
    N_array(i) = counts(1);
    [row_idx,col] = find(DL_labels==Cell_Type_Array(i));
    Depth_Subdata = double(FlyWire_DepthLength{row_idx,2:end});
    if i==1
        xline(Depth_Subdata(2)+0.2,'-.','Characteristic Layer Number','Color',[0 0 0],'HandleVisibility','off','LineWidth',1,'LabelVerticalAlignment','middle','FontSize',8);
        %xline(Depth_Subdata(1),'-','Mean Max Depth','Color',colors(1,:),'HandleVisibility','off','LineWidth',1.5);
    else
        xline(Depth_Subdata(2),'-.','Characteristic Layer Number','Color',colors(i,:),'HandleVisibility','off','LineWidth',1,'FontSize',8);
        %xline(Depth_Subdata(1),'-','Mean Max Depth','Color',colors(i-1,:),'HandleVisibility','off','LineWidth',1.5);
    end
    layer_counts = subdata{1};
    layer_counts = layer_counts(1:round(Depth_Subdata(1)));
    if ~isempty(layer_counts)
        if i==0
            %errorbar(linspace(1,45,45),layer_counts(1:45),p_STD(1:45),'-*k','LineWidth',2)
            plot(linspace(1,length(layer_counts),length(layer_counts)),layer_counts,'-*','Color',[0 0 0],'LineWidth',2.5)
        elseif i==1
            plot(linspace(1,length(layer_counts),length(layer_counts)),layer_counts,'Color',[0 0 0],'LineWidth',2.5)
        else
            plot(linspace(1,length(layer_counts),length(layer_counts)),layer_counts,'Color',colors(i,:),'LineWidth',1.5)
        end
    end
end

Cell_Type_Array(4) = "Kenyon cell";
Cell_Type_Array(5) = "optic lobes";

%plot([0 45],[0.5 0.5],'-.k','LineWidth',2)
yline(0.5,'-.k',"   p_s^c",'FontSize',18,'LineWidth',2,'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom')

hold off
ylabel("$p_s^j$",'Interpreter','Latex')
xlabel("Layer # (j)")
%title("Branching Probability vs. Layer by Cell Type (Vertical Lines for Mean Max Depth and Widest Tree Layer)")
legend(vertcat(strcat(Cell_Type_Array'," (N=",string(N_array),")")),"Location",'northeast')
xlim([0 45])
ylim([0.45 0.7])

cd ..
cd Figures\Figure5

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('5D.pdf', '-dpdf', '-painters', '-bestfit');

cd ..\..
cd Figure45_Data

%% 5H

dataset_shuffle_vector = [8,4,7,3,5,1,6,2];

color_codes = ["2EAA90","E3242B","1F5A66","734F96","88C4D8","F7D42F","0855B1","D341D1"]';

colors = [];
for i=1:length(color_codes)
    str = char(color_codes(i));
    colors = [colors;sscanf(str(1:end),'%2x%2x%2x',[1 3])/255];
end

colors = colors(dataset_shuffle_vector,:);

colors(1:5,:) = ([115, 80, 151; 238, 135, 40; 124, 203, 221; 134, 199, 108; 92, 107, 178] ./ 255);

CellTable = readtable("FlyWireDepthLengthTable4.csv");

%Lengths = double(CellTable{1:end-1,4});
Widths = double(CellTable{1:end-1,5});
PeakDepth = double(CellTable{1:end-1,3});
MaxDepth = double(CellTable{1:end-1,2});
Labels = string(CellTable{1:end-1,1});
P_Vector = double(CellTable{1:end-1,7});

Widths = double(CellTable{1:end-1,3});
Weights = double(CellTable{1:end-1,end});
Weights([24,27],:) = [];

% These are the points to label and color
UseLabels = [7,10,13,15,18,24,20,23,27,30];

%Lengths = log(10.^(Lengths));

Lengths = double(CellTable{1:end-1,10});

figure('Renderer', 'painters', 'Position', [10 10 900 600])
hold on
for i=1:length(P_Vector)
    if i==30
        scatter(log(P_Vector(i))+log(2),Lengths(i),80,'filled','MarkerFaceColor',colors(1,:),"MarkerEdgeColor",colors(1,:),'LineWidth',1.5,'HandleVisibility','on')
    elseif i==24
        scatter(log(P_Vector(i))+log(2),Lengths(i),80,'filled','MarkerFaceColor',colors(4,:),"MarkerEdgeColor",colors(4,:),'LineWidth',1.5,'HandleVisibility','off')
    elseif i==13
        scatter(log(P_Vector(i))+log(2),Lengths(i),80,'filled','MarkerFaceColor',colors(2,:),"MarkerEdgeColor",colors(2,:),'LineWidth',1.5,'HandleVisibility','off')
    elseif i==7
        scatter(log(P_Vector(i))+log(2),Lengths(i),80,'filled','MarkerFaceColor',colors(3,:),"MarkerEdgeColor",colors(3,:),'LineWidth',1.5,'HandleVisibility','off')
    elseif i==15
        scatter(log(P_Vector(i))+log(2),Lengths(i),80,'filled','MarkerFaceColor',colors(5,:),"MarkerEdgeColor",colors(5,:),'LineWidth',1.5,'HandleVisibility','off')
    else
        scatter(log(P_Vector(i))+log(2),Lengths(i),20,'filled','MarkerFaceColor',[0.4 0.4 0.4],'HandleVisibility','off')
    end
    if ismember(i,UseLabels)
        if i==7
            text(log(P_Vector(i))+log(2)+0.0003,Lengths(i)-0.03,"Kenyon cell","FontSize",12)
        elseif i==23
            text(log(P_Vector(i))+log(2)+0.0003,Lengths(i)+0.03,Labels(i),"FontSize",12)
        elseif i==24
            text(log(P_Vector(i))+log(2)+0.0003,Lengths(i)+0.03,"optic lobes","FontSize",12)
        else
            text(log(P_Vector(i))+log(2)+0.0003,Lengths(i)-0.03,Labels(i),"FontSize",12)
        end
    end
end

%xlim([0.09 0.19])

x = log(P_Vector)+log(2);
x([24,27],:) = [];

y = Lengths;
y([24,27],:) = [];

mdl = fitlm(x,y,'Weights',Weights);

coeffs = mdl.Coefficients;
slope = double(coeffs{2,1});

slope_mu = slope;

% h = plot(mdl,'HandleVisibility','off');
% delete(h(1));

legend('off');

y_pred = predict(mdl,x);

plot(x,y_pred,'r','LineWidth',0.75);

%plot(log([nanmin(P_Vector),nanmax(P_Vector)])+log(2),(log([nanmin(P_Vector),nanmax(P_Vector)])+log(2)).*a,'-.k')

text(0.11,4.85,strcat("n = ",string(round(slope,2))," ± ",string(round(double(coeffs{2,2}),2))),'Color','red','FontSize',12,'Rotation',30)

hold off
%xlim([10 70])
%ylabel("$\mu_{L}$ measured, $\langle log (L) \rangle$",'Interpreter','Latex')
ylabel("$\langle \log{\frac{L}{\langle l \rangle}} \rangle$",'Interpreter','Latex')
xlabel("$\langle \log{p_{s}^j} \rangle + log(2)$",'Interpreter','Latex')
title("$\mu_L$ dependence on $p_{surv}$",'Interpreter','Latex')
title("")
%legend("Cell Type",strcat("Fit y \propto ",string(round(slope,2)),"x"),'Location','northwest')

xlim([0.09 0.185])
ylim([4.4 6.6])

cd ..
cd Figures\Figure5

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('5H.pdf', '-dpdf', '-painters', '-bestfit');

cd ..\..
cd Figure45_Data

P_sigma = double(CellTable{1:end-1,9}).^2;
%L_sigma = double(CellTable{1:end-1,8}).^2;
L_sigma = double(CellTable{1:end-1,11}).^2;

figure('Renderer', 'painters', 'Position', [10 10 900 600])
hold on
for i=1:length(P_Vector)
    if i==30
        scatter(P_sigma(i),L_sigma(i),80,'filled','MarkerFaceColor',colors(1,:),"MarkerEdgeColor",colors(1,:),'LineWidth',1.5,'HandleVisibility','on')
    elseif i==24
        scatter(P_sigma(i),L_sigma(i),80,'filled','MarkerFaceColor',colors(4,:),"MarkerEdgeColor",colors(4,:),'LineWidth',1.5,'HandleVisibility','off')
    elseif i==13
        scatter(P_sigma(i),L_sigma(i),80,'filled','MarkerFaceColor',colors(2,:),"MarkerEdgeColor",colors(2,:),'LineWidth',1.5,'HandleVisibility','off')
    elseif i==7
        scatter(P_sigma(i),L_sigma(i),80,'filled','MarkerFaceColor',colors(3,:),"MarkerEdgeColor",colors(3,:),'LineWidth',1.5,'HandleVisibility','off')
    elseif i==15
        scatter(P_sigma(i),L_sigma(i),80,'filled','MarkerFaceColor',colors(5,:),"MarkerEdgeColor",colors(5,:),'LineWidth',1.5,'HandleVisibility','off')
    else
        scatter(P_sigma(i),L_sigma(i),20,'filled','MarkerFaceColor',[0.4 0.4 0.4],'HandleVisibility','off')
    end
    if ismember(i,UseLabels)
        if i==7
            text(P_sigma(i)+0.0003,L_sigma(i)-0.003,"Kenyon cell","FontSize",12)
        elseif i==23
            text(P_sigma(i)+log(2)+0.0003,L_sigma(i)+0.003,Labels(i),"FontSize",12)
        elseif i==24
            text(P_sigma(i)+0.0003,L_sigma(i)+0.003,"optic lobes","FontSize",12)
        else
            text(P_sigma(i)+0.0003,L_sigma(i)-0.003,Labels(i),"FontSize",12)
        end
    end
end

x = P_sigma;
x([24,27],:) = [];

y = L_sigma;
y([24,27],:) = [];

mdl = fitlm(x,y,'Weights',Weights);

coeffs = mdl.Coefficients;
slope = double(coeffs{2,1});

intercept = mean(y - slope_mu * x);

% h = plot(mdl,'HandleVisibility','off');
% delete(h(1));

%xlim([0.06 Inf])
%ylim([0 Inf])

legend('off');

%plot(log([nanmin(P_Vector),nanmax(P_Vector)])+log(2),(log([nanmin(P_Vector),nanmax(P_Vector)])+log(2)).*a,'-.k')

%plot(x,slope_mu * x + intercept);

y_pred = predict(mdl,x);

plot(x,y_pred,'r','LineWidth',0.75);

text(0.015,-0.255,strcat("n = ",string(round(slope,2))," ± ",string(round(double(coeffs{2,2}),2))),'Color','red','FontSize',12,'Rotation',25)

hold off
%xlim([10 70])
ylabel("$\sigma_{\frac{L}{\langle l \rangle}}^2$",'Interpreter','Latex')
xlabel("$\sigma_{p}^2$",'Interpreter','Latex')
title("$\mu_L$ dependence on $p_{surv}$",'Interpreter','Latex')
title("")
%legend("Cell Type",strcat("Fit y \propto ",string(round(slope,2)),"x"),'Location','northwest')

xlim([0.01 0.054])
ylim([-0.7 0.72])

cd ..
cd Figures\Figure5

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('sigma_Lp.pdf', '-dpdf', '-painters', '-bestfit');

%% Data Loading

Use_Fits = true;
normalize = true;
balanced_sampling = false;
sample_counts = 10000;

cd ..\..
cd FigureGeneration

%color_codes = ["2EAA90","E3242B","1F5A66","0855B1","88C4D8","F7D42F","734F96","D341D1"]';

color_codes = ["2EAA90","E3242B","1F5A66","734F96","88C4D8","F7D42F","0855B1","D341D1"]';

label_strings = ["data   FlyWire (FF)","data   H01 (Human)","data   Hemibrain (FF)","data   Larva (FF)","data   MANC (FF)","data   MM3 (Mouse)","data   Zebrafish","data   C.Elegans"];

% Unused now, but the length scale vector is something I've used in
% previous versions to make sure everything is in microns (right now all
% the length data in the spreadsheets already has this correction)

length_scale_vector = [1,1000,125,1000,125,1,1,1];

colors = [];
for i=1:length(color_codes)
    str = char(color_codes(i));
    colors = [colors;sscanf(str(1:end),'%2x%2x%2x',[1 3])/255];
end

% Loading in neuron counts / brain volumes / bounding volumes as lists

neuron_counts = [2*10^5,86*10^9,2*10^5,10^4,2*10^5,70*10^6,10^5,302];
brain_volumes = [80*10^6,1.0934*10^15,80*10^6,2.6*10^6,80*10^6,1.696*10^11,1*10^8,1];
bounding_volumes = [7.8*10^7,1.07*10^9,2.04*10^7,2.61*10^6,8.17*10^7,5.09*10^8,2.4*10^6,1];

dataset_shuffle_vector = [8,4,7,3,5,1,6,2];

colors = colors(dataset_shuffle_vector,:);

label_strings = extractAfter(label_strings(dataset_shuffle_vector),"data   ");

length_scale_vector = length_scale_vector(dataset_shuffle_vector);
neuron_counts = neuron_counts(dataset_shuffle_vector);
bounding_volumes = bounding_volumes(dataset_shuffle_vector);
brain_volumes = brain_volumes(dataset_shuffle_vector);

sub_shuffle_vector = [6,7,5,4,3,2,1];

Local_Dens_Data = load("Local_Dens_Data.mat");
Local_Dens_Data = Local_Dens_Data.Artificial_Dens_Data;

Local_Dens_Data = Local_Dens_Data(sub_shuffle_vector,:);

table_data = cell(8,1);

cd ..
cd ProcessedData\H01_104

H01_table = readtable("H01_AggregatedData.csv");
table_data{8} = H01_table;

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

table_data{7} = MM3_table;

cd ..
cd FlyWire

FlyWire_table = readtable("FlyWire_AggregatedData.csv");
table_data{6} = FlyWire_table;

cd ..
cd MANC

MANC_table = readtable("MANC_AggregatedData.csv");

MANC_Inclusion = load("MANC_Inclusion_List.csv");

[Lia,Lib] = ismember(double(MANC_table{:,1}),MANC_Inclusion);

MANC_table = MANC_table(Lia,:);

table_data{5} = MANC_table;

cd ..
cd Hemibrain

HEMI_table = readtable("Hemibrain_AggregatedData.csv");

HEMI_Exclusion = load("Hemibrain_Exclusion_List.csv");

[Lia,Lib] = ismember(double(HEMI_table{:,1}),HEMI_Exclusion);

HEMI_table = HEMI_table(~Lia,:);

table_data{4} = HEMI_table;

cd ..
cd Larva

Larva_table = readtable("Larva_AggregatedData.csv");
table_data{2} = Larva_table;

cd ..
cd Zebrafish

Zebrafish_table = readtable("Zebrafish_AggregatedData.csv");

Zebrafish_synapses = Zebrafish_table{:,3};

% Setting Synapses as weighted degree

Zebrafish_table{:,3} = Zebrafish_table{:,6};

table_data{3} = Zebrafish_table;

cd ..
cd Celegans

CELEGANs_table = readtable("CELEGANS_AggregatedData.csv");
table_data{1} = CELEGANs_table;

cd ..
cd ..
cd FigureGeneration

Moments_table = readtable("Datasets_All_Properties.csv");

cd ..
cd Figure45_Data

sub_shuffle_vector = [6,7,5,4,3,2,1];

% This section loads the local density data, which is currently stored in a
% .mat file (not the .csvs)

FlyWire_CellType = load("FlyWire_CellTypeData_Giulia.mat");
FlyWire_CellType = FlyWire_CellType.Cell_Type_Data;

FlyWire_CellType = FlyWire_CellType(1,[35,18,10,29,20]);

label_strings = ["visual","central","Kenyon Cell","optic lobes","gustatory"];

opts = detectImportOptions("FlyWire_CellTypeData.csv");
opts = setvartype(opts,'CellID', 'string');
FlyWire_Table = readtable("FlyWire_CellTypeData.csv",opts);

table_data2 = cell(5,1);

for i=1:length(table_data2)
    subdata_cell = FlyWire_CellType{i};
    [Lia,Lib] = ismember(FlyWire_Table{:,1},subdata_cell{3});
    table_data2{i} = FlyWire_Table(Lia,:);
end

%% FlyWire 5E

fig = figure('Renderer', 'painters', 'Position', [10 10 1200 300]);

mu = [];
sigma = [];
deviation = [];
mean_val = [];
mse = [];

subplot1(1,5,'Gap',[0.03 0.16],'Min',[0.1 0.2],'Max',[0.9,0.9])

legend_strings = [];
for i=1:length(table_data2)
    
    subplot1(i)
    
    hold on
    
    subtable = table_data2{i};
    
    xdata = double(subtable{:,2});
    
    plot_data = xdata;
    
    plot_data(plot_data==0) = [];
    plot_data(isnan(plot_data)) = [];
    plot_data(isinf(plot_data)) = [];
    
    pd = fitdist(plot_data,'Exponential');
    
    disp(pd.mu);
    
    %disp(mean(plot_data.^2));
    
    [xBinMeans,yBinMeans,y_vals,mu_val,sigma_val] = fitlogncdf(plot_data,30);
    mu = [mu;mu_val];
    sigma = [sigma;sigma_val];
    deviation = [deviation;std(log(plot_data))];
    mean_val = [mean_val;mean(log(plot_data))];
    pd = makedist('Lognormal','mu',mu_val,'sigma',sigma_val);
    
    [N,edges] = histcounts(log(plot_data),20);
    
    % You'll see this section repeated frequently, this is an adjustment to
    % represent density accurately

    edge_difs = exp(edges(2:end)) - exp(edges(1:end-1));
    edges = exp((edges(2:end) + edges(1:end-1)) ./ 2);
    x = edges;
    y = pdf(pd,x);
    y_fit = y.*(edge_difs.*sum(N));
    y_empirical = N;    
    y_max = max(y_fit);
    y_fit = y_fit ./ y_max;
    y_empirical = y_empirical ./ y_max;
    %plot(x,y_fit,'-','LineWidth',1,'MarkerFaceColor',colors(i,:),'Color',colors(i,:));
    %plot(x,y_fit,'-','LineWidth',0.75,'MarkerFaceColor',colors(i,:),'Color',colors(i,:),'HandleVisibility','off');
    %scatter(edges,y_empirical,20,colors(i,:),'filled')
    scatter(edges(y_empirical~=0),y_empirical(y_empirical~=0),10,colors(i,:),'o','LineWidth',1.5);
    mse = [mse;goodnessOfFit(y_empirical,y_fit,'MSE')];
    legend_strings = [legend_strings,label_strings(i)];
    if Use_Fits
        %legend_strings = [legend_strings,strcat("Fit: \mu=",string(round(exp(mu_val),2)),", \sigma=",string(round(sigma_val,2)))];
    end
    set(gca, 'XScale', 'log')
    xlim([10^(1) 10^4])
    ylim([0 1.4])
    title(strcat(label_strings(i)," \mu=",string(round(mu_val,2))," \sigma=",string(round(sigma_val,2))),'FontSize',8)
    
    hold off
    
    xticks(logspace(1,4,4));
    %xticklabels(["","10^{1}  ","","10^{3}  ","","10^{5}  "]);
    
    if i==2 || i==3 || i==4 || i==6 || i==7 || i==8
        set(gca, 'YTick', []);
    else
        set(gca, 'YTick', [0 0.5 1]);
    end
    
    ax = gca;
    
    ax.TickLength = [0.06, 0.025];
    ax.LineWidth = 0.6;
    
    set(ax, 'Box', 'off');
    
    ax.XTick = ax.XTick;
    ax.YTick = ax.YTick;
    
end
han=axes(fig,'visible','off');
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,"P(L)",'FontSize',24);
xlabel(han,"$L (\mu m)$",'Interpreter','Latex','FontSize',20);
%sgtitle("Neuron Length PDF - $P(L)=\frac{1}{L\sigma\sqrt{2\pi}}e^{\frac{-(ln(L)-\mu_{*})^2}{2\sigma^2}} $, $\mu_{*}=e^{\mu}$",'Interpreter','latex', 'FontSize', 16);
sgtitle("Neuron Length",'Interpreter','latex', 'FontSize', 24);

ax = gca;
ax.LooseInset = [0.08 0 0 0];

cd ..
cd Figures\Figure5

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('5E.pdf', '-dpdf', '-painters', '-bestfit');

cd ..\..
cd Figure45_Data

%% 5A

fig = figure('Renderer', 'painters', 'Position', [10 10 1200 300]);

mu = [];
sigma = [];
deviation = [];
mean_val = [];
mse = [];

subplot1(1,5,'Gap',[0.03 0.16],'Min',[0.1 0.2],'Max',[0.9,0.9])

legend_strings = [];
for i=1:length(table_data2)
    
    subplot1(i)
    
    hold on
    
    subtable = table_data2{i};
    
    xdata = double(subtable{:,3});
    
    plot_data = xdata;
    
    plot_data(plot_data==0) = [];
    plot_data(isnan(plot_data)) = [];
    plot_data(isinf(plot_data)) = [];
    
    pd = fitdist(plot_data,'Exponential');
    
    disp(pd.mu);
    
    %disp(mean(plot_data.^2));
    
    [xBinMeans,yBinMeans,y_vals,mu_val,sigma_val] = fitlogncdf(plot_data,30);
    mu = [mu;mu_val];
    sigma = [sigma;sigma_val];
    deviation = [deviation;std(log(plot_data))];
    mean_val = [mean_val;mean(log(plot_data))];
    pd = makedist('Lognormal','mu',mu_val,'sigma',sigma_val);
    
    [N,edges] = histcounts(log(plot_data),20);
    
    % You'll see this section repeated frequently, this is an adjustment to
    % represent density accurately

    edge_difs = exp(edges(2:end)) - exp(edges(1:end-1));
    edges = exp((edges(2:end) + edges(1:end-1)) ./ 2);
    x = edges;
    y = pdf(pd,x);
    y_fit = y.*(edge_difs.*sum(N));
    y_empirical = N;    
    y_max = max(y_fit);
    y_fit = y_fit ./ y_max;
    y_empirical = y_empirical ./ y_max;
    %plot(x,y_fit,'-','LineWidth',1,'MarkerFaceColor',colors(i,:),'Color',colors(i,:));
    %plot(x,y_fit,'-','LineWidth',0.75,'MarkerFaceColor',colors(i,:),'Color',colors(i,:),'HandleVisibility','off');
    %scatter(edges,y_empirical,20,colors(i,:),'filled')
    scatter(edges(y_empirical~=0),y_empirical(y_empirical~=0),10,colors(i,:),'o','LineWidth',1.5);
    mse = [mse;goodnessOfFit(y_empirical,y_fit,'MSE')];
    legend_strings = [legend_strings,label_strings(i)];
    if Use_Fits
        %legend_strings = [legend_strings,strcat("Fit: \mu=",string(round(exp(mu_val),2)),", \sigma=",string(round(sigma_val,2)))];
    end
    set(gca, 'XScale', 'log')
    xlim([10^(0) 10^3])
    ylim([0 1.4])
    title(strcat(label_strings(i)," \mu=",string(round(mu_val,2))," \sigma=",string(round(sigma_val,2))),'FontSize',8)
    
    hold off
    
    xticks(logspace(0,3,4));
    %xticklabels(["","10^{1}  ","","10^{3}  ","","10^{5}  "]);
    
    if i==2 || i==3 || i==4 || i==6 || i==7 || i==8
        set(gca, 'YTick', []);
    else
        set(gca, 'YTick', [0 0.5 1]);
    end
    
    ax = gca;
    
    ax.TickLength = [0.06, 0.025];
    ax.LineWidth = 0.6;
    
    set(ax, 'Box', 'off');
    
    ax.XTick = ax.XTick;
    ax.YTick = ax.YTick;
    
end
han=axes(fig,'visible','off');
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,"P(S)",'FontSize',24);
xlabel(han,"$S$",'Interpreter','Latex','FontSize',20);
%sgtitle("Neuron Length PDF - $P(L)=\frac{1}{L\sigma\sqrt{2\pi}}e^{\frac{-(ln(L)-\mu_{*})^2}{2\sigma^2}} $, $\mu_{*}=e^{\mu}$",'Interpreter','latex', 'FontSize', 16);
sgtitle("Strength (Synapses)",'Interpreter','latex', 'FontSize', 24);

ax = gca;
ax.LooseInset = [0.08 0 0 0];

cd ..
cd Figures\Figure5

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('5A.pdf', '-dpdf', '-painters', '-bestfit');

cd ..\..
cd Figure45_Data

%% 5F

figure('Renderer', 'painters', 'Position', [10 10 900 600])

mu = [];
sigma = [];
mse = [];

all_plot_data = [];

legend_strings = [];
hold on
for i=1:length(table_data2)
    
    subtable = table_data2{i};
    
    xdata = double(subtable{:,2});
    
    plot_data = xdata;
    plot_data(plot_data==0) = [];
    plot_data(isnan(plot_data)) = [];
    plot_data(isinf(plot_data)) = [];
    
    % This section is where the rescaling actually happens
    
    plot_data = log(plot_data);
    plot_data = plot_data - mean(plot_data);
    plot_data = plot_data ./ std(plot_data);
    plot_data = exp(plot_data);
    
    if balanced_sampling
        subsample = datasample(plot_data,sample_counts);
        all_plot_data = [all_plot_data;subsample];
    else
        all_plot_data = [all_plot_data;plot_data];
    end

    [xBinMeans,yBinMeans,y_vals,mu_val,sigma_val] = fitlogncdf(plot_data,30);
    mu = [mu;mu_val];
    sigma = [sigma;sigma_val];
    pd = makedist('Lognormal','mu',mu_val,'sigma',sigma_val);

    [N,edges] = histcounts(log(plot_data),20);
    
    edge_difs = exp(edges(2:end)) - exp(edges(1:end-1));
    edges = exp((edges(2:end) + edges(1:end-1)) ./ 2);
    x = edges;
    y = pdf(pd,x);
    y_fit = y.*(edge_difs.*sum(N));
    y_empirical = N;    
    y_max = max(y_fit);
    y_fit = y_fit ./ y_max;
    y_empirical = y_empirical ./ y_max;
    
    scatter(edges,y_empirical,40,colors(i,:),'filled')
    
    if i==4
        [N,edges] = histcounts(log(plot_data),100);
    else
        [N,edges] = histcounts(log(plot_data),20);
    end
    
    edge_difs = exp(edges(2:end)) - exp(edges(1:end-1));
    edges = exp((edges(2:end) + edges(1:end-1)) ./ 2);
    x = edges;
    y = pdf(pd,x);
    y_fit = y.*(edge_difs.*sum(N));
    y_empirical = N;    
    y_max = max(y_fit);
    y_fit = y_fit ./ y_max;
    y_empirical = y_empirical ./ y_max;
    
    % Draw one line, using fit parameters from Hemibrain
    
    if i==4
        plot(x,y_fit,'-','LineWidth',1,'MarkerFaceColor',[0 0 0],'Color',[0 0 0],'HandleVisibility','off');
    end
    
    mse = [mse;goodnessOfFit(y_empirical,y_fit,'MSE')];
    %legend_strings = [legend_strings,label_strings(i)];
    %legend_strings = [legend_strings,strcat("Fit: \mu=",string(round(exp(mu_val),2)),", \sigma=",string(round(sigma_val,2)))];
end
hold off
set(gca, 'XScale', 'log')
%legend(legend_strings,"Location",'best')
xlabel("$\overline{L} = exp(\frac{\ln(L) - \mu_{L}}{\sigma_{L}})$", 'Interpreter','Latex','FontSize',24);
ylabel("$P(\overline{L})$", 'Interpreter','Latex','FontSize',24);
%title("Neuron Lengths Rescaled - $P(\overline{L})=\frac{1}{\overline{L}\sigma\sqrt{2\pi}}e^{\frac{-(ln(\overline{L})-\mu)^2}{2\sigma^2}}$",'Interpreter','latex', 'FontSize', 16);
xlim([10^(-2) 10^2])
ylim([0 1.4])
legend(label_strings,'FontSize',16)

ax = gca;
ax.FontSize = 20;
ax.TickLength = [0.03, 0.025];
ax.LineWidth = 0.6;

cd ..
cd Figures\Figure5

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('5F.pdf', '-dpdf', '-painters', '-bestfit');

cd ..\..
cd Figure45_Data

%% 5B

figure('Renderer', 'painters', 'Position', [10 10 900 600])

mu = [];
sigma = [];
mse = [];

all_plot_data = [];

legend_strings = [];
hold on
for i=1:length(table_data2)
    
    subtable = table_data2{i};
    
    xdata = double(subtable{:,3});
    
    plot_data = xdata;
    plot_data(plot_data==0) = [];
    plot_data(isnan(plot_data)) = [];
    plot_data(isinf(plot_data)) = [];
    
    % This section is where the rescaling actually happens
    
    plot_data = log(plot_data);
    plot_data = plot_data - mean(plot_data);
    plot_data = plot_data ./ std(plot_data);
    plot_data = exp(plot_data);
    
    if balanced_sampling
        subsample = datasample(plot_data,sample_counts);
        all_plot_data = [all_plot_data;subsample];
    else
        all_plot_data = [all_plot_data;plot_data];
    end

    [xBinMeans,yBinMeans,y_vals,mu_val,sigma_val] = fitlogncdf(plot_data,30);
    mu = [mu;mu_val];
    sigma = [sigma;sigma_val];
    pd = makedist('Lognormal','mu',mu_val,'sigma',sigma_val);

    [N,edges] = histcounts(log(plot_data),20);
    
    edge_difs = exp(edges(2:end)) - exp(edges(1:end-1));
    edges = exp((edges(2:end) + edges(1:end-1)) ./ 2);
    x = edges;
    y = pdf(pd,x);
    y_fit = y.*(edge_difs.*sum(N));
    y_empirical = N;    
    y_max = max(y_fit);
    y_fit = y_fit ./ y_max;
    y_empirical = y_empirical ./ y_max;
    
    scatter(edges,y_empirical,40,colors(i,:),'filled')
    
    if i==4
        [N,edges] = histcounts(log(plot_data),100);
    else
        [N,edges] = histcounts(log(plot_data),20);
    end
    
    edge_difs = exp(edges(2:end)) - exp(edges(1:end-1));
    edges = exp((edges(2:end) + edges(1:end-1)) ./ 2);
    x = edges;
    y = pdf(pd,x);
    y_fit = y.*(edge_difs.*sum(N));
    y_empirical = N;    
    y_max = max(y_fit);
    y_fit = y_fit ./ y_max;
    y_empirical = y_empirical ./ y_max;
    
    % Draw one line, using fit parameters from Hemibrain
    
    if i==4
        plot(x,y_fit,'-','LineWidth',1,'MarkerFaceColor',[0 0 0],'Color',[0 0 0],'HandleVisibility','off');
    end
    
    mse = [mse;goodnessOfFit(y_empirical,y_fit,'MSE')];
    %legend_strings = [legend_strings,label_strings(i)];
    %legend_strings = [legend_strings,strcat("Fit: \mu=",string(round(exp(mu_val),2)),", \sigma=",string(round(sigma_val,2)))];
end
hold off
set(gca, 'XScale', 'log')
%legend(legend_strings,"Location",'best')
xlabel("$\overline{S} = exp(\frac{\ln(S) - \mu_{S}}{\sigma_{S}})$", 'Interpreter','Latex','FontSize',24);
ylabel("$P(\overline{S})$", 'Interpreter','Latex','FontSize',24);
%title("Neuron Lengths Rescaled - $P(\overline{L})=\frac{1}{\overline{L}\sigma\sqrt{2\pi}}e^{\frac{-(ln(\overline{L})-\mu)^2}{2\sigma^2}}$",'Interpreter','latex', 'FontSize', 16);
xlim([10^(-2) 10^2])
ylim([0 1.4])
legend(label_strings,'FontSize',16)

ax = gca;
ax.FontSize = 20;
ax.TickLength = [0.03, 0.025];
ax.LineWidth = 0.6;

cd ..
cd Figures\Figure5

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('5B.pdf', '-dpdf', '-painters', '-bestfit');

cd ..\..
cd Figure45_Data


%% KSSTAT Calculations

% This part calculates KS statistics for each of the quantities, and stores
% that data for the quantity plots later

KSSTATs_All = cell(length(table_data2),1);

for i=1:length(table_data2)
    
    subtable = table_data2{i};
    
    KSSTATs_Array = [];
    
    xdata = double(subtable{:,2});
    
    plot_data = xdata;
    
    plot_data(plot_data==0) = [];
    plot_data(isnan(plot_data)) = [];
    plot_data(isinf(plot_data)) = [];
    
    KSSTATs = KSSTAT_Lognormal3(plot_data);
    
    KSSTATs_Array = [KSSTATs_Array,KSSTATs];
    
    xdata = double(subtable{:,3});
    
    plot_data = xdata;
    
    plot_data(plot_data==0) = [];
    plot_data(isnan(plot_data)) = [];
    plot_data(isinf(plot_data)) = [];
    
    KSSTATs = KSSTAT_Lognormal3(plot_data);
    
    KSSTATs_Array = [KSSTATs_Array,KSSTATs];
    
    xdata = double(subtable{:,4});
    
    plot_data = xdata;
    
    plot_data(plot_data==0) = [];
    plot_data(isnan(plot_data)) = [];
    plot_data(isinf(plot_data)) = [];
    
    KSSTATs = KSSTAT_Lognormal3(plot_data);
    
    KSSTATs_Array = [KSSTATs_Array,KSSTATs];
    
    xdata = double(subtable{:,5});
    
    plot_data = xdata;
    
    plot_data(plot_data==0) = [];
    plot_data(isnan(plot_data)) = [];
    plot_data(isinf(plot_data)) = [];
    
    KSSTATs = KSSTAT_Lognormal3(plot_data);
    
    KSSTATs_Array = [KSSTATs_Array,KSSTATs];
    
    xdata = double(subtable{:,6});
    
    plot_data = xdata;
    
    plot_data(plot_data==0) = [];
    plot_data(isnan(plot_data)) = [];
    plot_data(isinf(plot_data)) = [];
    
    KSSTATs = KSSTAT_Lognormal3(plot_data);
    
    KSSTATs_Array = [KSSTATs_Array,KSSTATs];
    
    KSSTATs_All{i} = KSSTATs_Array;
    
end

Quantity_KSSTATs = cell(size(KSSTATs_Array,2),1);

for j=1:size(KSSTATs_Array,2)
    
    Quantity_Data = [];
    
    for i=1:length(KSSTATs_All)

        substats = KSSTATs_All{i};
        
        Quantity_Data = [Quantity_Data;substats(:,j)'];

    end
    
    Quantity_KSSTATs{j,1} = Quantity_Data;
    
end

Quantity_Labels = ["L","S","\rho","K","K_{w}"];

%% 5G

close all

fig = figure('Renderer', 'painters', 'Position', [10 10 900 600]);

for i=1:1
    hold on
    substats = Quantity_KSSTATs{i};
    
    plot(substats(:,1:end-1),':.','LineWidth',2,'MarkerSize',25)
    
    plot(substats(:,end),'.-','LineWidth',2.5,'MarkerSize',25)
    
    %plot(KSSTATs_All{i}',':','LineWidth',2)
    hold off
    %legend(["Poisson","Exponential","Power Law","LN Fit"],"Location",'northeast')
    %ylabel("KS Statistic (Lower Value Indicates Better Fit)",'FontSize',12)
    %xlabel("Measurement",'FontSize',12)
    xticks([1,2,3,4,5])
    
    label_strings2 = strcat(label_strings,"(");
    label_strings2 = extractBefore(label_strings2,"(");
    
    xticklabels(label_strings2)
    
    ax = gca;
    ax.FontSize = 16;
    %title(Quantity_Labels(i),'FontSize',16)
end

ylim([0 0.8])

legend(["Poisson","Exponential","Power Law","LN Fit"],"Location",'northeast')

han=axes(fig,'visible','off');
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,{'KS Statistic (L)';''},'FontSize',20);
%xlabel(han,'Dataset');

cd ..
cd Figures\Figure5

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('5G.pdf', '-dpdf', '-painters', '-bestfit');

cd ..\..
cd Figure45_Data

%% 5C

close all

fig = figure('Renderer', 'painters', 'Position', [10 10 900 600]);

for i=2:2
    hold on
    substats = Quantity_KSSTATs{i};
    
    plot(substats(:,1:end-1),':.','LineWidth',2,'MarkerSize',25)
    
    plot(substats(:,end),'.-','LineWidth',2.5,'MarkerSize',25)
    
    %plot(KSSTATs_All{i}',':','LineWidth',2)
    hold off
    %legend(["Poisson","Exponential","Power Law","LN Fit"],"Location",'northeast')
    %ylabel("KS Statistic (Lower Value Indicates Better Fit)",'FontSize',12)
    %xlabel("Measurement",'FontSize',12)
    xticks([1,2,3,4,5])
    
    label_strings2 = strcat(label_strings,"(");
    label_strings2 = extractBefore(label_strings2,"(");
    
    xticklabels(label_strings2)
    
    ax = gca;
    ax.FontSize = 16;
    %title(Quantity_Labels(i),'FontSize',16)
end

ylim([0 0.8])

legend(["Poisson","Exponential","Power Law","LN Fit"],"Location",'northeast')

han=axes(fig,'visible','off');
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,{'KS Statistic (S)';''},'FontSize',20);
%xlabel(han,'Dataset');

cd ..
cd Figures\Figure5

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('5C.pdf', '-dpdf', '-painters', '-bestfit');