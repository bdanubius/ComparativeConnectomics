clear

Use_Fits = true;
normalize = true;
balanced_sampling = false;
sample_counts = 10000;

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

%% Data Loading

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
cd FigureGeneration\HelperFunctions

%% Figure 3A

figure('Renderer', 'painters', 'Position', [10 10 1200 600])

a = [];
b = [];
r2 = [];

legend_strings = [];
hold on
for i=2:7
    
    subtable = table_data{i};
    
    xdata = double(subtable{:,2});
    
    subtable = table_data{i};
    
    ydata = double(subtable{:,3});
    
    xdata(ydata==0) = [];
    xdata(isnan(ydata)) = [];
    xdata(isinf(ydata)) = [];
    ydata(ydata==0) = [];
    ydata(isnan(ydata)) = [];
    ydata(isinf(ydata)) = [];
    ydata(xdata==0) = [];
    ydata(isnan(xdata)) = [];
    ydata(isinf(xdata)) = [];
    xdata(xdata==0) = [];
    xdata(isnan(xdata)) = [];
    xdata(isinf(xdata)) = [];
    
    [xBinMeans,yBinMeans,bincounts] = LinBinsCounts(xdata,ydata,20);
    bincounts(isnan(xBinMeans)) = [];
    yBinMeans(isnan(xBinMeans)) = [];
    xBinMeans(isnan(xBinMeans)) = [];
    if i==1
        scatter(linspace(0,7000,10),zeros(10,1)+6.24,100*sqrt(10)+1,colors(i,:),'filled')
    else
        scatter(xBinMeans,yBinMeans,30,colors(i,:),'filled')
    end
    P = polyfit(xdata,ydata,1);
    yfit = P(1)*xBinMeans+P(2);
    a = [a;P(1)];
    b = [b;P(2)];
    mdl = fitlm(xdata,ydata);
    r2 = [r2;mdl.Rsquared.Ordinary];
    if i==1
        %plot(linspace(0,7000,10),zeros(10,1)+6.24,'LineWidth',1,'MarkerFaceColor',colors(i,:),'Color',colors(i,:),'HandleVisibility','off');
    else
        %plot(xBinMeans,yfit,'LineWidth',1,'MarkerFaceColor',colors(i,:),'Color',colors(i,:),'HandleVisibility','off');
    end
    legend_strings = [legend_strings,label_strings(i)];
    %legend_strings = [legend_strings,strcat("Fit: a=",string(round(P(1),2)),", b=",string(round(P(2),2)))];
    %legend_strings = [legend_strings,strcat("Fit: a=",string(round(P(1),2)))];
end
plot(logspace(1,5,10),logspace(1,5,10),'-.k','LineWidth',2)
hold off
legend(horzcat(legend_strings,"y=x"),"Location",'northeast')
xlabel("L (\mum)");
ylabel("S");
%title("Neuron Synapses vs. Lengths - $S = aL + b$",'Interpreter','latex');
xlim([10^1 10^5])
ylim([10^1 10^5])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

ax = gca; % Get the current axis handle
xTicks = ax.XTick; % Get current x-axis tick positions
yTicks = ax.YTick; % Get current y-axis tick positions

% Set ticks explicitly to make them fixed
ax.XTick = xTicks;
ax.YTick = yTicks;

ax.TickLength = [0.04, 0.025];

f = gcf;
f.Position = [10, 10, 600, 700];

cd ..\..
cd Figures/Figure3

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('4A.pdf', '-dpdf', '-painters', '-bestfit');

cd ..\..
cd FigureGeneration/HelperFunctions

colors2 = colors;
label_strings2 = label_strings;

%% Figure 3B

close all

table_data(1:1:7,1) = table_data(8:-1:2,1);
table_data(8,:) = [];

cd ..\..
cd ProcessedData\H01_104
All_Density = cell(7,1);

Data = load("H01_Local_Density.mat");
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

Data = load("MANC_Local_Density.mat");
Data = Data.outdata;

All_Density{4} = Data;

cd ..
cd Hemibrain

Data = load("HEMI_Local_Density.mat");
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

Artificial_Dens_Data = cell(7,2);

% As I did many random walks per neuron, I need to aggregate that data,
% that's what this part is doing.

for i=1:length(All_Density)
    
    Subdata = All_Density{i};
    SynapseDist_data = Subdata{3};
    LengthDist_data = Subdata{4};
    Synapse10_data = Subdata{1};
    Length10_data = Subdata{2};
    SubFilenames = Subdata{5};
    
    Density10_vector = NaN(length(Synapse10_data),1);
    
    for j=1:length(Length10_data)
        disp(j)
        SubSubS = Synapse10_data{j};
        SubSubL = Length10_data{j};
        [row,col] = find(SubSubL > 10);
        SubSubD = (SubSubS(row)) ./ SubSubL(row);
        Density10_vector(j) = exp(nanmean(log(nonzeros(SubSubD))));
    end
    
    SubFilenames(isnan(Density10_vector)) = [];
    Density10_vector(isnan(Density10_vector)) = [];
    
    Density_Vectors{i,1} = Density10_vector;
    
    DensityDist_vector = NaN(length(SynapseDist_data),1);
    
    for j=1:length(LengthDist_data)
        disp(j)
        SubSubL = LengthDist_data{j};
        [row,col] = find(SubSubL>0);
        DensityDist_vector(j) = exp(nanmean(log(1 ./ SubSubL(row))));
    end
    DensityDist_vector(isnan(DensityDist_vector)) = [];
    
    if i==1
        subtable = table_data{i};
        AllFileNames = double(subtable{:,1});
        AllLength = double(subtable{:,2});
        SubFilenames = extractBefore(SubFilenames,".");
        [Lia,Lib] = ismember(double(SubFilenames),AllFileNames);
        SubDensCorr = Density10_vector(Lia);
        SubLength = AllLength(nonzeros(Lib));
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
        D_subdata = log10(SubDensCorr);
        L_subdata = log10(SubLength);
        [R,P] = corrcoef(D_subdata(~isnan(L_subdata)),L_subdata(~isnan(L_subdata)));
        LD_corr(i) = R(2,1);
    end
    
    Artificial_Dens_Data{i,1} = 10.^D_subdata(~isnan(L_subdata));
    Artificial_Dens_Data{i,2} = 10.^L_subdata(~isnan(L_subdata)) .* 10.^D_subdata(~isnan(L_subdata));
    
    Density_Vectors{i,2} = DensityDist_vector;
end

cd ..\..
cd FigureGeneration

D_col = [];
D_dev = [];
D_avg = [];

label_strings = ["data   H01 (Human)","data   MM3 (Mouse)","data   FLYWIRE (F.F.)","data   MANC (F.F.)","data   Hemibrain (F.F.)","data   LARVA (F.F.)","data   Zebrafish"];

label_strings = extractAfter(label_strings,"data   ");

figure('Renderer', 'painters', 'Position', [10 10 1200 600])
for i=1:size(Density_Vectors,1)
    D_col = [D_col;mean(log10(nonzeros(Density_Vectors{i,1})))];
    D_avg = [D_avg;mean(nonzeros(Density_Vectors{i,1}))];
    D_dev = [D_dev;std(log10(nonzeros(Density_Vectors{i,1})))];
    subplot(2,4,i)
    nhist(log10(nonzeros(Density_Vectors{i,1})))
    xlabel("")
    title(label_strings(i))
    xlabel("Log <\rho_{loc}>")
    ylabel("PDF")
end

colors = [0.890196078431373,0.141176470588235,0.168627450980392;0.968627450980392,0.831372549019608,0.184313725490196;0.180392156862745,0.666666666666667,0.564705882352941;0.533333333333333,0.768627450980392,0.847058823529412;0.121568627450980,0.352941176470588,0.400000000000000;0.450980392156863,0.309803921568627,0.588235294117647;0.0313725490196078,0.333333333333333,0.694117647058824];

all_data = readtable("Datasets_All_Properties.csv");

% Loading moments

L_col = double(all_data{:,2});
%D_col = double(all_data{:,4});
S_col = double(all_data{:,6});
K_col = double(all_data{:,8});
Kw_col = double(all_data{:,10});
PSK_col = log10(double(all_data{:,15})./double(all_data{:,14}));

S_mean_pred = L_col(1:7) + D_col;
K_mean_pred = S_col + PSK_col;

L_dev = double(all_data{:,3});
%D_dev = double(all_data{:,5});
S_dev = double(all_data{:,7});
K_dev = double(all_data{:,9});
Kw_dev = double(all_data{:,11});
PSK_dev = double(all_data{:,16});
%LD_corr = double(all_data{:,12});
PSKS_corr = double(all_data{:,17});

% Use weighted degree as synapses hemibrain and MANC

S_dev(5) = Kw_dev(5);
S_col(5) = Kw_col(5);

S_dev(4) = Kw_dev(4);
S_col(4) = Kw_col(4);

S_dev_pred = sqrt(L_dev(1:7).^2 + D_dev.^2 + 2.*L_dev(1:7).*D_dev.*LD_corr(1:7));
%S_dev_pred = sqrt(L_dev(1:7).^2 + D_dev.^2);
K_dev_pred = sqrt(S_dev.^2 + PSK_dev.^2 + 2.*S_dev.*PSK_dev.*PSKS_corr);

global_min = min(vertcat(L_col(1:7)./log10(exp(1)),S_col(1:7)./log10(exp(1)),L_col(1:7)./log10(exp(1))+D_col./log10(exp(1))));
global_max = max(vertcat(L_col(1:7)./log10(exp(1)),S_col(1:7)./log10(exp(1)),L_col(1:7)./log10(exp(1))+D_col./log10(exp(1))));

figure('Renderer', 'painters', 'Position', [10 10 1200 600])

hold on
plot(linspace(global_min-0.5,global_max+0.5,10),linspace(global_min-0.5,global_max+0.5,10),'-.k','LineWidth',1)

for i=1:7
    if i==1
        %errorbar(S_col(i),L_col(i)+log10(D_avg(i)),L_dev(i),L_dev(i),S_dev(i),S_dev(i),'+','LineWidth',0.5,'Color',colors(i,:).*0.8,'HandleVisibility','off')
        %scatter(S_col(i),L_col(i)+log10(D_avg(i)),60,'o','MarkerEdgeColor',colors(i,:).*0.8,'LineWidth',1)
        
        %errorbar(K_mean,K_mean_partial_pred,S_std,S_std,K_std,K_std,'g*','LineWidth',0.5,'HandleVisibility','off')
        %scatter(K_mean,K_mean_partial_pred,100,'go','LineWidth',2)
        errorbar(S_col(i)./log10(exp(1)),S_mean_pred(i)./log10(exp(1)),S_dev_pred(i)./log10(exp(1)),S_dev_pred(i)./log10(exp(1)),S_dev(i)./log10(exp(1)),S_dev(i)./log10(exp(1)),'*','LineWidth',0.5,'Color',colors(i,:),'HandleVisibility','off')
        scatter(S_col(i)./log10(exp(1)),S_mean_pred(i)./log10(exp(1)),100,'filled','MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:))
    else
        %errorbar(S_col(i),L_col(i)+log10(D_avg(i)),L_dev(i),L_dev(i),S_dev(i),S_dev(i),'+','LineWidth',0.5,'Color',colors(i,:).*0.8,'HandleVisibility','off')
        %scatter(S_col(i),L_col(i)+log10(D_avg(i)),60,'o','MarkerEdgeColor',colors(i,:).*0.8,'LineWidth',1,'HandleVisibility','off')
        
        %errorbar(K_mean,K_mean_partial_pred,S_std,S_std,K_std,K_std,'g*','LineWidth',0.5,'HandleVisibility','off')
        %scatter(K_mean,K_mean_partial_pred,100,'go','LineWidth',2)
        errorbar(S_col(i)./log10(exp(1)),S_mean_pred(i)./log10(exp(1)),S_dev_pred(i)./log10(exp(1)),S_dev_pred(i)./log10(exp(1)),S_dev(i)./log10(exp(1)),S_dev(i)./log10(exp(1)),'*','LineWidth',0.5,'Color',colors(i,:),'HandleVisibility','off')
        scatter(S_col(i)./log10(exp(1)),S_mean_pred(i)./log10(exp(1)),100,'filled','MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:),'HandleVisibility','off')
    end
end

hold off
xlim([global_min-0.5 global_max+0.5])
ylim([global_min-0.5 global_max+0.5])
ylabel("\mu_{S*} (Predicted)")
xlabel("\mu_{S} (Measured)")
%legend([" y=x line","\mu_{S*} = \mu_{L} + ln(<\rho_{loc}>) (Length)","\mu_{S*} = \mu_{L} + \mu_{\rho_{loc}} (Length x Density)"],"Location",'northeast','FontSize',8)
%legend([" y=x line","\mu_{S*} = \mu_{L} + \mu_{\rho_{loc}} (Length x Density)"],"Location",'northeast')
%title("Predicted Synapses vs. Actual Synapses (per each dataset)")
xlim([1.5./log10(exp(1)) 3.8./log10(exp(1))])
ylim([1.5./log10(exp(1)) 4.5./log10(exp(1))])

ax = gca; % Get the current axis handle
xTicks = ax.XTick; % Get current x-axis tick positions
yTicks = ax.YTick; % Get current y-axis tick positions

% Set ticks explicitly to make them fixed
ax.XTick = xTicks;
ax.YTick = yTicks;

ax.TickLength = [0.04, 0.025];

f = gcf;
f.Position = [10, 10, 600, 700];

cd ..
cd Figures\Figure3

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('4B.pdf', '-dpdf', '-painters', '-bestfit');

global_min = min(vertcat(L_dev(1:7)./log10(exp(1)),S_dev(1:7)./log10(exp(1)),S_dev_pred(1:7)./log10(exp(1))));
global_max = max(vertcat(L_dev(1:7)./log10(exp(1)),S_dev(1:7)./log10(exp(1)),S_dev_pred(1:7)./log10(exp(1))));

%% Figure 3C

figure('Renderer', 'painters', 'Position', [10 10 1200 600])

hold on

for i=1:7
    if i==1
        %scatter(S_dev(i),L_dev(i),60,'o','MarkerEdgeColor',colors(i,:).*0.8,'LineWidth',1)
        scatter(S_dev(i)./log10(exp(1)),S_dev_pred(i)./log10(exp(1)),100,'filled','MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:))
    else
        %scatter(S_dev(i),L_dev(i),60,'o','MarkerEdgeColor',colors(i,:).*0.8,'LineWidth',1,'HandleVisibility','off')
        scatter(S_dev(i)./log10(exp(1)),S_dev_pred(i)./log10(exp(1)),100,'filled','MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:),'HandleVisibility','off')
    end
end 

plot(linspace(global_min-0.5,global_max+0.5,10),linspace(global_min-0.5,global_max+0.5,10),'-.k','LineWidth',1)
hold off
% xlim([global_min-0.5 global_max+0.5])
% ylim([global_min-0.5 global_max+0.5])
ylabel("\sigma_{S*} (Predicted)")
xlabel("\sigma_{S} (Measured)")
%legend(["\sigma_{S*} = \sigma_{L}","\sigma_{S*} = (\sigma^2_{L}+\sigma^2_{\rho_{loc}}+2\sigma_{L}\sigma_{\rho_{loc}}corr(L,\rho_{loc}))^{1/2}"," y=x line"],"Location",'northeast','FontSize',7)
%legend(["\sigma_{S*} = (\sigma^2_{L}+\sigma^2_{\rho_{loc}}+2\sigma_{L}\sigma_{\rho_{loc}}corr(L,\rho_{loc}))^{1/2}"," y=x line"],"Location",'northeast')
%title("Predicted Synapse Deviation vs. Actual Synapse Deviation")
xlim([0.2./log10(exp(1)) 0.55./log10(exp(1))])
ylim([0.2./log10(exp(1)) 0.55./log10(exp(1))])
xticks(linspace(0.5,1.1,4))

ax = gca; % Get the current axis handle
xTicks = ax.XTick; % Get current x-axis tick positions
yTicks = ax.YTick; % Get current y-axis tick positions

% Set ticks explicitly to make them fixed
ax.XTick = xTicks;
ax.YTick = yTicks;

ax.TickLength = [0.04, 0.025];

f = gcf;
f.Position = [10, 10, 600, 700];

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('4C.pdf', '-dpdf', '-painters', '-bestfit');

% figure('Renderer', 'painters', 'Position', [10 10 1200 600]);
% 
% hold on
% 
% for i=1:7
%     if i==1
%         %scatter(exp(S_dev(i)),exp(L_dev(i)),60,'o','MarkerEdgeColor',colors(i,:).*0.8,'LineWidth',1)
%         scatter(exp(S_dev(i)),exp(S_dev_pred(i)),100,'filled','MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:))
%     else
%         %scatter(exp(S_dev(i)),exp(L_dev(i)),60,'o','MarkerEdgeColor',colors(i,:).*0.8,'LineWidth',1,'HandleVisibility','off')
%         scatter(exp(S_dev(i)),exp(S_dev_pred(i)),100,'filled','MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:),'HandleVisibility','off')
%     end
% end 
% 
% plot(linspace(0,2,10),linspace(0,2,10),'-.k','LineWidth',1)
% hold off
% % xlim([global_min-0.5 global_max+0.5])
% % ylim([global_min-0.5 global_max+0.5])
% ylabel("exp(\sigma_{S*}) (Predicted)")
% xlabel("exp(\sigma_{S}) (Measured)")
% %legend(["\sigma_{S*} = \sigma_{L}","\sigma_{S*} = (\sigma^2_{L}+\sigma^2_{\rho_{loc}}+2\sigma_{L}\sigma_{\rho_{loc}}corr(L,\rho_{loc}))^{1/2}"," y=x line"],'Location','northeast','FontSize',12)
% legend(["\sigma_{S*} = (\sigma^2_{L}+\sigma^2_{\rho_{loc}}+2\sigma_{L}\sigma_{\rho_{loc}}corr(L,\rho_{loc}))^{1/2}"," y=x line"],'Location','northeast')
% title("Predicted Synapse Deviation vs. Actual Synapse Deviation")
% % xlim([0.2 0.55])
% % ylim([0.2 0.55])
% 
% cd C:\Users\bspia\OneDrive\Desktop\Matlab_Scripts\Daniel_Table2Figure\Figures
% 
% ax = gca; % Get the current axis handle
% xTicks = ax.XTick; % Get current x-axis tick positions
% yTicks = ax.YTick; % Get current y-axis tick positions
% 
% % Set ticks explicitly to make them fixed
% ax.XTick = xTicks;
% ax.YTick = yTicks;
% 
% ax.TickLength = [0.04, 0.025];
% 
% f = gcf;
% f.Position = [10, 10, 600, 700];
% 
% cd Figure4
% 
% pos = get(gcf, 'Position');
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
% set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
% set(gcf, 'PaperPositionMode', 'auto');
% print('PredSynDev10_Linear.pdf', '-dpdf', '-painters', '-bestfit');

%% Figure 4D-F

% Note: Have to reinitialize tables

Use_Fits = true;
normalize = true;
balanced_sampling = false;
sample_counts = 10000;

cd ..\..
cd FigureGeneration

%% Data Loading

Moments_table = readtable("Datasets_All_Properties.csv");

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
cd FigureGeneration\HelperFunctions

% L and S correctly refer to the corresponding quantities

L_col = [3.89000000000000;3.41000000000000;2.84000000000000;3.06000000000000;3.01000000000000;2.41000000000000;2.61000000000000];

S_pred = [3.29678899068438;3.35702370804304;2.50627700381322;3.09434941284233;2.75731011737809;2.19930344681399;2.62907108048627];

S_col = [3.34000000000000;3.34000000000000;2.36000000000000;2.98000000000000;2.70000000000000;2.15000000000000;2.59000000000000];

L_dev = [0.360000000000000;0.220000000000000;0.370000000000000;0.460000000000000;0.440000000000000;0.250000000000000;0.360000000000000];

S_dev_pred = [0.426675997125073;0.274899314619755;0.449870009499464;0.579742379878768;0.473918788322101;0.402234924687398;0.396710837334093];

S_dev = [0.380000000000000;0.280000000000000;0.540000000000000;0.500000000000000;0.400000000000000;0.520000000000000;0.400000000000000;0.320000000000000];

LD_corr = [-0.186869090722982;0.387270905953448;0.0555151177424913;0.414559647673023;0.0365190774967295;0.240245685957107;0.0595324981266745];

shuffle_vector = [6,7,5,4,3,2,1];

L_col = L_col(shuffle_vector);
S_col = S_col(shuffle_vector);
L_dev = L_dev(shuffle_vector);
S_dev = S_dev(shuffle_vector);
S_dev_pred = S_dev_pred(shuffle_vector);
S_pred = S_pred(shuffle_vector);
LD_corr = LD_corr(shuffle_vector);

sub_shuffle_vector = [6,7,5,4,3,2,1];

Dens_mu_vector = [];
Dens_sigma_vector = [];

for i=1:size(Local_Dens_Data,1)
    
     Density_Vec = Local_Dens_Data{i,1};
     
     plot_data = Density_Vec;
     
     plot_data(plot_data==0) = [];
     plot_data(isnan(plot_data)) = [];
     plot_data(isinf(plot_data)) = [];
     
     [xBinMeans,yBinMeans,y_vals,mu_val,sigma_val] = fitlogncdf(plot_data,30);
     
     Dens_mu_vector = [Dens_mu_vector;mu_val];
     Dens_sigma_vector = [Dens_sigma_vector;sigma_val];
     
end

Dens_mu_vector = Dens_mu_vector(sub_shuffle_vector);
Dens_sigma_vector = Dens_sigma_vector(sub_shuffle_vector);

xlim_range = [0,5];

xlim_array = [0,4;1,4;2,4;2,4;1,4;2,5;2,5];

LN_S_mu = [];
LN_S_sigma = [];

LN_L_mu = [];
LN_L_sigma = [];

LN_S_pred_mu = [];
LN_S_pred_sigma = [];

for i=2:length(table_data)
    
    figure('Renderer', 'painters', 'Position', [10 10 1200 600])

    mu = [];
    sigma = [];
    deviation = [];
    mean_val = [];
    mse = [];

    legend_strings = [];
    hold on
        
    subtable = table_data{i};

    if i==3
        xdata = Zebrafish_synapses;
    elseif i==4 || i==5
        xdata = double(subtable{:,6});
    else
        xdata = double(subtable{:,3});
    end
    
    plot_data = xdata;
    
    plot_data(plot_data==0) = [];
    plot_data(isnan(plot_data)) = [];
    plot_data(isinf(plot_data)) = [];
    
    [xBinMeans,yBinMeans,y_vals,mu_val,sigma_val] = fitlogncdf(plot_data,30);
    mu = [mu;mu_val];
    sigma = [sigma;sigma_val];
    deviation = [deviation;std(log(plot_data))];
    mean_val = [mean_val;mean(log(plot_data))];
    
    %pd = makedist('Lognormal','mu',mu_val,'sigma',sigma_val);
    pd = makedist('Lognormal','mu',mu_val,'sigma',sigma_val);
    
    LN_S_mu = [LN_S_mu;mu_val];
    LN_S_sigma = [LN_S_sigma;sigma_val];
    
    if i==2 || i==3 || i==8
        [N,edges] = histcounts(log(plot_data),20);
    elseif i==6
        [N,edges] = histcounts(log(plot_data),20);
    else
        [N,edges] = histcounts(log(plot_data),30);
    end
    
    %[N,edges] = histcounts(log(plot_data),20);
    
    disp(mean(log10(plot_data)))
    disp(std(log10(plot_data)))

    edge_difs = exp(edges(2:end)) - exp(edges(1:end-1));
    edges = exp((edges(2:end) + edges(1:end-1)) ./ 2);
    x = edges;
    y = pdf(pd,x);
    y_fit = y.*(edge_difs.*sum(N));
    y_empirical = N;    
    y_max = max(y_fit);
    y_fit = y_fit ./ y_max;
    y_empirical = y_empirical ./ y_max;
    if Use_Fits
        %plot(x,y_fit,'-','LineWidth',1,'MarkerFaceColor',colors(i,:),'Color',colors(i,:));
        %plot(x,y_fit,'-','LineWidth',1.5,'MarkerFaceColor',colors(i,:),'Color',colors(i,:),'HandleVisibility','off');
        %scatter(edges,y_empirical,20,colors(i,:),'filled')
        scatter(edges,y_empirical,60,colors(i,:),'filled');
    else
        scatter(edges,y_empirical,60,colors(i,:),'filled');
    end
    
    x_spline = logspace(log10(min(edges)/5),log10(max(edges)*5),1000);
    
    s = spline(x,y_fit,x_spline);
    
    %plot(x_spline,s,'-','LineWidth',3,'MarkerFaceColor',colors(i,:),'Color',colors(i,:));

    mse = [mse;goodnessOfFit(y_empirical,y_fit,'MSE')];
    legend_strings = [legend_strings,label_strings(i)];
    if Use_Fits
        %legend_strings = [legend_strings,strcat("Fit: \mu=",string(round(exp(mu_val),2)),", \sigma=",string(round(sigma_val,2)))];
    end
    
    subtable = table_data{i};
    
    xdata = double(subtable{:,2});
    
    plot_data = xdata;
    
    plot_data(plot_data==0) = [];
    plot_data(isnan(plot_data)) = [];
    plot_data(isinf(plot_data)) = [];
    
    [xBinMeans,yBinMeans,y_vals,mu_val,sigma_val] = fitlogncdf(plot_data,30);
    mu = [mu;mu_val];
    sigma = [sigma;sigma_val];
    deviation = [deviation;std(log(plot_data))];
    mean_val = [mean_val;mean(log(plot_data))];
    
    pd = makedist('Lognormal','mu',mu_val,'sigma',sigma_val);
    
    LN_L_mu = [LN_L_mu;mu_val];
    LN_L_sigma = [LN_L_sigma;sigma_val];
    
    [N,edges] = histcounts(log(plot_data),1000);
    
%     y_fit = pdf(pd,logspace(xlim_range(1),xlim_range(2),1000));
%     plot(logspace(xlim_range(1),xlim_range(2),1000),y_fit./max(y_fit),':','LineWidth',1.5,'MarkerFaceColor',colors(i,:),'Color',colors(i,:),'HandleVisibility','off');
        
    edge_difs = exp(edges(2:end)) - exp(edges(1:end-1));
    edges = exp((edges(2:end) + edges(1:end-1)) ./ 2);
    x = edges;
    y = pdf(pd,x);
    y_fit = y.*(edge_difs.*sum(N));
    y_empirical = N;    
    y_max = max(y_fit);
    y_fit = y_fit ./ y_max;
    y_empirical = y_empirical ./ y_max;
    
    %plot(x,y_fit,':','LineWidth',1.5,'MarkerFaceColor',colors(i,:),'Color',colors(i,:),'HandleVisibility','off');
    
    x_spline = logspace(log10(min(edges)/5),log10(max(edges)*5),1000);
    
    s = spline(x,y_fit,x_spline);
    
    %plot(x_spline,s,':','LineWidth',1.5,'MarkerFaceColor',colors(i,:),'Color',colors(i,:));

    %pd = makedist('Lognormal','mu',mu_val+Dens_mu_vector(i-1),'sigma',sqrt(sigma_val.^2 + Dens_sigma_vector(i-1).^2 + 2*sigma_val*Dens_sigma_vector(i-1)*LD_corr(i-1)));
    pd = makedist('Lognormal','mu',log(10.^S_pred(i-1)),'sigma',log(10.^S_dev_pred(i-1)));
    
    LN_S_pred_mu = [LN_S_pred_mu;log(10.^S_pred(i-1))];
    LN_S_pred_sigma = [LN_S_pred_sigma;log(10.^S_dev_pred(i-1))];
    
    if i==3
        xdata = Zebrafish_synapses;
    else
        xdata = double(subtable{:,3});
    end
    
    plot_data = xdata;
    
    plot_data(plot_data==0) = [];
    plot_data(isnan(plot_data)) = [];
    plot_data(isinf(plot_data)) = [];
    
    [N,edges] = histcounts(log(plot_data),20);
    
%     y_fit = pdf(pd,logspace(xlim_range(1),xlim_range(2),1000));
%     plot(logspace(xlim_range(1),xlim_range(2),1000),y_fit./max(y_fit),'-','LineWidth',1.5,'MarkerFaceColor',colors(i,:),'Color',colors(i,:),'HandleVisibility','off');
    
    edge_difs = exp(edges(2:end)) - exp(edges(1:end-1));
    edges = exp((edges(2:end) + edges(1:end-1)) ./ 2);
    x = edges;
    y = pdf(pd,x);
    y_fit = y.*(edge_difs.*sum(N));
    y_empirical = N;    
    y_max = max(y_fit);
    y_fit = y_fit ./ y_max;
    y_empirical = y_empirical ./ y_max;
    
    %plot(x,y_fit,'-','LineWidth',1.5,'MarkerFaceColor',colors(i,:),'Color',colors(i,:),'HandleVisibility','off');
    
    if i==4
        x_spline = logspace(log10(min(edges)/2),log10(max(edges)*2),1000);
    else
        x_spline = logspace(log10(min(edges)/5),log10(max(edges)*5),1000);
    end
    s = spline(x,y_fit,x_spline);
    
    plot(x_spline,s,'-','LineWidth',1.5,'MarkerFaceColor',colors(i,:),'Color',colors(i,:));

    hold off
    
    set(gca, 'XScale', 'log')
    %legend([strcat("Synapse Distribution Data (\mu=",string(round(log10(exp(LN_S_mu(i-1))),2)),",\sigma=",string(round(log10(exp(LN_S_sigma(i-1))),2)),")"),strcat("Predicted fit for synapses (L) (\mu=",string(round(log10(exp(LN_L_mu(i-1))),2)),",\sigma=",string(round(log10(exp(LN_L_sigma(i-1))),2)),")"),strcat("Predicted fit for synapses (L,\rho) (\mu=",string(round(log10(exp(LN_S_pred_mu(i-1))),2)),",\sigma=",string(round(log10(exp(LN_S_pred_sigma(i-1))),2)),")")],"Location",'northeast','FontSize',12)
    legend([strcat("Data (\mu=",string(round(log10(exp(LN_S_mu(i-1)))./log10(exp(1)),2)),",\sigma=",string(round(log10(exp(LN_S_sigma(i-1)))./log10(exp(1)),2)),")"),strcat("Predicted (\mu=",string(round(log10(exp(LN_S_pred_mu(i-1)))./log10(exp(1)),2)),",\sigma=",string(round(log10(exp(LN_S_pred_sigma(i-1)))./log10(exp(1)),2)),")")],"Location",'northeast')
    xlabel("S");
    ylabel("P(S)");
    title(strcat(label_strings(i)," Neuron Predicted Synapses vs. Length Fits"),'Interpreter','latex');
    xlim([10^xlim_array(i-1,1) 10^xlim_array(i-1,2)])
    ylim([0 1.4])
    
    ax = gca; % Get the current axis handle
    xTicks = ax.XTick; % Get current x-axis tick positions
    yTicks = ax.YTick; % Get current y-axis tick positions

    % Set ticks explicitly to make them fixed
    ax.XTick = xTicks;
    ax.YTick = yTicks;

    ax.TickLength = [0.04, 0.025];
    
    f = gcf;
    f.Position = [10, 10, 600, 700];
    
    cd ..
    cd ..
    cd Figures\Figure3
    
    pos = get(gcf, 'Position');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
    set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
    set(gcf, 'PaperPositionMode', 'auto');
    print(strcat(label_strings(i),'PredSyn.pdf'), '-dpdf', '-painters', '-bestfit');
    
    cd ..
    cd ..
    cd FigureGeneration\HelperFunctions
    
end

%% Figure 3G

cd ..\..
cd GraphData

% All of this is loading the shuffled versions of the weighted adjacencies,
% as we needed those for a figure Giulia asked for.

load_shuffled_data = true;

Graph_Data = load("Larva_Graphs");
Graph_Data = Graph_Data.outgraphs;

Weighted_A = Graph_Data{1};

Degree_Data = load("Larva_Degree_List.mat");
Degree_Data = Degree_Data.Larva_Degrees;

if ~load_shuffled_data

    [row,col,v] = find(Weighted_A);

    v = v(randperm(numel(v)));
    
    Shuffled_A = sparse(row, col, v, size(Weighted_A, 1), size(Weighted_A, 2));
    
    save("LarvaShuffled.mat",'Shuffled_A');
    
else
    
    Shuffled_A = load("LarvaShuffled.mat");
    Shuffled_A = Shuffled_A.Shuffled_A;
    
end

Larva_Shuffled_A = Shuffled_A;
Larva_Normal_A = Weighted_A;
Larva_Degree_Data = Degree_Data;

Graph_Data = load("Zebrafish_Graphs");
Graph_Data = Graph_Data.outgraphs;

Weighted_A = Graph_Data{1};

Degree_Data = load("Zebrafish_Degree_List.mat");
Degree_Data = Degree_Data.Zebrafish_Degrees;

if ~load_shuffled_data

    [row,col,v] = find(Weighted_A);

    v = v(randperm(numel(v)));
    
    Shuffled_A = sparse(row, col, v, size(Weighted_A, 1), size(Weighted_A, 2));
    
    save("ZebrafishShuffled.mat",'Shuffled_A');
    
else
    
    Shuffled_A = load("ZebrafishShuffled.mat");
    Shuffled_A = Shuffled_A.Shuffled_A;
    
end

Zebrafish_Shuffled_A = Shuffled_A;
Zebrafish_Normal_A = Weighted_A;
Zebrafish_Degree_Data = Degree_Data;

Graph_Data = load("FlyWire_Graphs");
Graph_Data = Graph_Data.outgraphs;

Weighted_A = Graph_Data{1};

Degree_Data = load("FlyWire_Degree_List.mat");
Degree_Data = Degree_Data.FlyWire_Degrees;

if ~load_shuffled_data

    [row,col,v] = find(Weighted_A);

    v = v(randperm(numel(v)));
    
    Shuffled_A = sparse(row, col, v, size(Weighted_A, 1), size(Weighted_A, 2));
    
    save("FlyWireShuffled.mat",'Shuffled_A');
    
else
    
    Shuffled_A = load("FlyWireShuffled.mat");
    Shuffled_A = Shuffled_A.Shuffled_A;
    
end

FlyWire_Shuffled_A = Shuffled_A;
FlyWire_Normal_A = Weighted_A;
FlyWire_Degree_Data = Degree_Data;

Graph_Data = load("MM3_Graphs");
Graph_Data = Graph_Data.outgraphs;

Weighted_A = Graph_Data{1};

Degree_Data = load("MM3_Degree_List2.mat");
Degree_Data = Degree_Data.MM3_Degrees;

if ~load_shuffled_data

    [row,col,v] = find(Weighted_A);

    v = v(randperm(numel(v)));
    
    Shuffled_A = sparse(row, col, v, size(Weighted_A, 1), size(Weighted_A, 2));
    
    save("MM3Shuffled.mat",'Shuffled_A');
    
else
    
    Shuffled_A = load("MM3Shuffled.mat");
    Shuffled_A = Shuffled_A.Shuffled_A;
    
end

MM3_Shuffled_A = Shuffled_A;
MM3_Normal_A = Weighted_A;
MM3_Degree_Data = Degree_Data;

Degree_Data = load("H01_Degree_List.mat");
Degree_Data = Degree_Data.H01_Degrees;
H01_Degree_Data = Degree_Data;

HEMI_Shuffled_A = load("HEMIShuffled.mat");
HEMI_Shuffled_A = HEMI_Shuffled_A.Shuffled_A;
HEMI_Normal_A = load("HEMI_Weighted_A2.mat");
HEMI_Normal_A = HEMI_Normal_A.Weighted_A2;
HEMI_Degree_Data = load("HEMIDegreeList.mat");
HEMI_Degree_Data = HEMI_Degree_Data.Degree_List;

MANC_Shuffled_A = load("MANCShuffled.mat");
MANC_Shuffled_A = MANC_Shuffled_A.Shuffled_A;
MANC_Normal_A = load("MANC_Weighted_A2.mat");
MANC_Normal_A = MANC_Normal_A.Weighted_A2;
MANC_Degree_Data = load("MANCDegreeList.mat");
MANC_Degree_Data = MANC_Degree_Data.Degree_List;

Degree_Data_All = cell(8,1);
Shuffled_Data = cell(8,1);
Normal_Data = cell(8,1);

Shuffled_Data{7} = Zebrafish_Shuffled_A;
Normal_Data{7} = Zebrafish_Normal_A;
Degree_Data_All{7} = Zebrafish_Degree_Data;
Shuffled_Data{6} = Larva_Shuffled_A;
Normal_Data{6} = Larva_Normal_A;
Degree_Data_All{6} = Larva_Degree_Data;
Shuffled_Data{5} = HEMI_Shuffled_A;
Normal_Data{5} = HEMI_Normal_A;
Degree_Data_All{5} = HEMI_Degree_Data;
Shuffled_Data{4} = MANC_Shuffled_A;
Normal_Data{4} = MANC_Normal_A;
Degree_Data_All{4} = MANC_Degree_Data;
Shuffled_Data{3} = FlyWire_Shuffled_A;
Normal_Data{3} = FlyWire_Normal_A;
Degree_Data_All{3} = FlyWire_Degree_Data;
Shuffled_Data{2} = MM3_Shuffled_A;
Normal_Data{2} = MM3_Normal_A;
Degree_Data_All{2} = MM3_Degree_Data;
Degree_Data_All{1} = H01_Degree_Data;

%label_strings = ["MM3 (Mouse)","FLYWIRE (F.F.)","MANC (F.F.)","Hemibrain (F.F.)","LARVA (F.F.)","Zebrafish"];
% 
% colors = [0.968627450980392,0.831372549019608,0.184313725490196;
%     0.180392156862745,0.666666666666667,0.564705882352941;
%     0.533333333333333,0.768627450980392,0.847058823529412;
%     0.121568627450980,0.352941176470588,0.400000000000000;
%     0.450980392156863,0.309803921568627,0.588235294117647;
%     0.0313725490196078,0.333333333333333,0.694117647058824];

color_codes = ["2EAA90","E3242B","1F5A66","734F96","88C4D8","F7D42F","0855B1","D341D1"]';

label_strings = ["data   FLYWIRE (F.F.)","data   H01 (Human)","data   Hemibrain (F.F.)","data   LARVA (F.F.)","data   MANC (F.F.)","data   MM3 (Mouse)","data   Zebrafish","data   C.Elegans"];

colors = [];
for i=1:length(color_codes)
    str = char(color_codes(i));
    colors = [colors;sscanf(str(1:end),'%2x%2x%2x',[1 3])/255];
end

dataset_shuffle_vector = [8,4,7,3,5,1,6,2];

colors = colors(dataset_shuffle_vector,:);

label_strings = extractAfter(label_strings(dataset_shuffle_vector),"data   ");

% Actual plotting starts here

cd ..
cd FigureGeneration\HelperFunctions

close all

figure('Renderer', 'painters', 'Position', [10 10 1200 600])

a = [];
b = [];
legend_strings = [];
hold on

for i=1:7
    
    disp(i)
    
    if i==1
       subtable = table_data{1};
       xdata = double(subtable{:,3});
       ydata = double(subtable{:,5});
    elseif i==8
       subtable = table_data{8};
       xdata = double(subtable{:,3});
       ydata = double(subtable{:,5});
    else
        Degree_Data = Degree_Data_All{9-i};
        xdata = double(Degree_Data(:,3));
        ydata = double(Degree_Data(:,2));
    end
    
    xdata(ydata==0) = [];
    xdata(isnan(ydata)) = [];
    xdata(isinf(ydata)) = [];
    ydata(ydata==0) = [];
    ydata(isnan(ydata)) = [];
    ydata(isinf(ydata)) = [];
    ydata(xdata==0) = [];
    ydata(isnan(xdata)) = [];
    ydata(isinf(xdata)) = [];
    xdata(xdata==0) = [];
    xdata(isnan(xdata)) = [];
    xdata(isinf(xdata)) = [];
    
    [xBinMeans,yBinMeans,bincounts] = LinBinsCounts(log(xdata),log(ydata),100);
    
    bincounts(isnan(xBinMeans)) = [];
    yBinMeans(isnan(xBinMeans)) = [];
    xBinMeans(isnan(xBinMeans)) = [];
    
    scatter(exp(xBinMeans),exp(yBinMeans),30*sqrt(bincounts./mean(bincounts))+1,colors(i,:),'filled','HandleVisibility','off')
    
    P = polyfit(log(xdata),log(ydata),1);
    
    %[robust_coefs,stats] = robustfit(log(xdata),log(ydata),'bisquare',0.6);
    
    %stats = stats.se;
    
    %P = robust_coefs([2,1]);
    
    yfit = P(2) + P(1).*xBinMeans;
    a = [a;P(1)];
    b = [b;P(2)];
    mdl = fitlm(log(xdata),log(ydata));
    
    ci = coefCI(mdl);
    
    alpha_ci = (ci(2,2) - ci(2,1)) ./ 2;
    beta_ci = (ci(1,2) - ci(1,1)) ./ 2;
    
    plot(exp(xBinMeans),exp(yfit),'LineWidth',1,'MarkerFaceColor',colors(i,:),'Color',colors(i,:));

    %legend_strings = [legend_strings,strcat("Fit: \alpha=",string(round(P(1),2))," ± ",string(round(alpha_ci,3))," \beta=",string(round(exp(P(2)),2))," ± ",string(round(beta_ci,3)))];
    legend_strings = [legend_strings,strcat(label_strings(i)," fit: \alpha=",string(round(P(1),2))," ± ",string(round(alpha_ci,3)))];
    
end

%plot(logspace(0,5,10),logspace(0,5,10),'-.k','LineWidth',2)
hold off
%legend(horzcat(legend_strings,"y \propto x^{1}"),"Location",'northeast','FontSize',12)

%legend(legend_strings,"Location",'northeast')

% Get all plot handles (in the order they were plotted)
handles = findall(gca, 'Type', 'scatter');

%legend_strings([3,4]) = legend_strings([4,3]);
%legend_strings([5,6]) = legend_strings([6,5]);

%handles([1,2]) = handles([2,1]);flip

% Reverse handles and legend entries
legend(flipud(handles), legend_strings,"Location",'northeast')

xlabel("S");
ylabel("k");
%title("Neuron Degrees vs. Synapses - $k = \beta S^{\alpha}$",'Interpreter','latex');
xlim([10^0 10^5])
ylim([10^0 10^4])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

ax = gca; % Get the current axis handle
xTicks = ax.XTick; % Get current x-axis tick positions
yTicks = ax.YTick; % Get current y-axis tick positions

% Set ticks explicitly to make them fixed
ax.XTick = xTicks;
ax.YTick = yTicks;

ax.TickLength = [0.04, 0.025];

f = gcf;
f.Position = [10, 10, 600, 700];

cd ..\..
cd Figures\Figure3

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('3G.pdf', '-dpdf', '-painters', '-bestfit');

cd ..\..
cd GraphData

%% Figure 3H 

S_mean = [];
S_std = [];
K_mean = [];
K_std = [];
K_mean_pred = [];
K_mean_partial_pred = [];
K_std_pred = [];

for i=2:7
    
    Degree_Data = Degree_Data_All{i};
    
    S_submean = mean(log(nonzeros(double(Degree_Data(:,3)))));
    K_submean = mean(log(nonzeros(double(Degree_Data(:,2)))));
    S_substd = std(log(nonzeros(double(Degree_Data(:,3)))));
    K_substd = std(log(nonzeros(double(Degree_Data(:,2)))));
    
    S_mean = [S_mean;S_submean];
    K_mean = [K_mean;K_submean];
    %K_mean_pred = [K_mean_pred;b(2*i-1)+a(2*i-1)*S_submean];
    K_mean_pred = [K_mean_pred;b(9-i)+a(9-i)*S_submean];
    %K_mean_partial_pred = [K_mean_partial_pred;b(2*i-1)+S_submean];
    K_mean_partial_pred = [K_mean_partial_pred;b(9-i)+S_submean];
    
    S_std = [S_std;S_substd];
    K_std = [K_std;K_substd];
    %K_std_pred = [K_std_pred;(a(2*i-1))*(S_substd)];
    K_std_pred = [K_std_pred;(a(9-i))*(S_substd)];
    
end

% Plotting for H

close all

figure('Renderer', 'painters', 'Position', [10 10 1200 600])
hold on
plot(linspace(1,9,10),linspace(1,9,10),'-.k','LineWidth',1)
for i=1:length(K_mean)
    if i==1
        %errorbar(K_mean(i),S_mean(i),S_std(i),S_std(i),K_std(i),K_std(i),'*','LineWidth',0.5,'Color',colors(i,:).*0.8,'HandleVisibility','off')
        %scatter(K_mean(i),S_mean(i),100,'o','MarkerEdgeColor',colors(i,:).*0.8,'LineWidth',1)
        %errorbar(K_mean,K_mean_partial_pred,S_std,S_std,K_std,K_std,'g*','LineWidth',0.5,'HandleVisibility','off')
        %scatter(K_mean,K_mean_partial_pred,100,'go','LineWidth',2)
        errorbar(K_mean(i),K_mean_pred(i),K_std_pred(i),K_std_pred(i),K_std(i),K_std(i),'*','LineWidth',0.5,'Color',colors(8-i,:),'HandleVisibility','off')
        scatter(K_mean(i),K_mean_pred(i),100,'filled','MarkerFaceColor',colors(8-i,:),'MarkerEdgeColor',colors(8-i,:))
    else
        %errorbar(K_mean(i),S_mean(i),S_std(i),S_std(i),K_std(i),K_std(i),'*','LineWidth',0.5,'Color',colors(i,:).*0.8,'HandleVisibility','off')
        %scatter(K_mean(i),S_mean(i),100,'o','MarkerEdgeColor',colors(i,:).*0.8,'LineWidth',1,'HandleVisibility','off')
        %errorbar(K_mean,K_mean_partial_pred,S_std,S_std,K_std,K_std,'g*','LineWidth',0.5,'HandleVisibility','off')
        %scatter(K_mean,K_mean_partial_pred,100,'go','LineWidth',2)
        errorbar(K_mean(i),K_mean_pred(i),K_std_pred(i),K_std_pred(i),K_std(i),K_std(i),'*','LineWidth',0.5,'Color',colors(8-i,:),'HandleVisibility','off')
        scatter(K_mean(i),K_mean_pred(i),100,'filled','MarkerFaceColor',colors(8-i,:),'MarkerEdgeColor',colors(8-i,:),'HandleVisibility','off')
    end
end
hold off
ylabel("\mu_{k*} (Predicted)")
xlabel("\mu_{k} (Measured)")
%legend([" y=x line","\mu_{k*} = \mu_{S} (Synapses)","\mu_{k*} = ln(\beta) + \mu_{S} (Rescaled Synapses)","\mu_{k*} = ln(\beta) + \alpha x \mu_{S} (PL Transformed S)"],'Location','southeast','FontSize',8)
%legend([" y=x line","\mu_{k*} = \mu_{S} (Synapses)","\mu_{k*} = ln(\beta) + \alpha \mu_{S} (PL Transformed S)"],"Location",'northeast','FontSize',8)
%legend([" y=x line","\mu_{k*} = ln(\beta) + \alpha \mu_{S} (PL Transformed S)"],"Location",'northeast')
%title("Predicted Degree vs. Actual Degree (per each dataset)")
xlim([1 7])

ax = gca; % Get the current axis handle
xTicks = ax.XTick; % Get current x-axis tick positions
yTicks = ax.YTick; % Get current y-axis tick positions

% Set ticks explicitly to make them fixed
ax.XTick = xTicks;
ax.YTick = yTicks;

ax.TickLength = [0.04, 0.025];

f = gcf;
f.Position = [10, 10, 600, 700];

cd ..
cd Figures\Figure3

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('3H.pdf', '-dpdf', '-painters', '-bestfit');

cd ..\..
cd GraphData

% Plotting for I

figure('Renderer', 'painters', 'Position', [10 10 1200 600])
hold on
for i=1:length(K_std)
    if i==1
        %scatter(K_std(i),S_std(i),100,'o','MarkerEdgeColor',colors(i,:).*0.8,'LineWidth',1)
        scatter(K_std(i),K_std_pred(i),100,'filled','MarkerFaceColor',colors(8-i,:),'MarkerEdgeColor',colors(8-i,:))
    else
        %scatter(K_std(i),S_std(i),100,'o','MarkerEdgeColor',colors(i,:).*0.8,'LineWidth',1,'HandleVisibility','off')
        scatter(K_std(i),K_std_pred(i),100,'filled','MarkerFaceColor',colors(8-i,:),'MarkerEdgeColor',colors(8-i,:),'HandleVisibility','off')
    end
end 
plot(linspace(0,2,10),linspace(0,2,10),'-.k','LineWidth',1)
hold off
ylabel("\sigma_{k*} (Predicted)")
xlabel("\sigma_{k} (Measured)")
%legend(["\sigma_{k*} = \sigma_{S}","\sigma_{k*} = \alpha \sigma_{S}"],"Location",'northeast','FontSize',12)
%legend("\sigma_{k*} = \alpha \sigma_{S}","Location",'northeast')
%title("Predicted Degree Deviation vs. Actual Degree Deviation")
xlim([0.5 1.5])
ylim([0.5 1.5])
set(gca, 'XScale', 'linear')

ax = gca; % Get the current axis handle
xTicks = ax.XTick; % Get current x-axis tick positions
yTicks = ax.YTick; % Get current y-axis tick positions

% Set ticks explicitly to make them fixed
ax.XTick = xTicks;
ax.YTick = yTicks;

ax.TickLength = [0.04, 0.025];

f = gcf;
f.Position = [10, 10, 600, 700];

cd ..
cd Figures\Figure3

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('3I.pdf', '-dpdf', '-painters', '-bestfit');

cd ..\..
cd GraphData

cd ..
cd FigureGeneration

%% Figure 4J-L

clear
close all

% Have to reload table

Use_Fits = true;
normalize = true;
balanced_sampling = false;
sample_counts = 10000;

%% Data Loading

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
cd FigureGeneration\HelperFunctions

cd ..

Moments_table = readtable("Datasets_All_Properties.csv");

cd HelperFunctions

% Have to reinitialize tables again

% Here L and S map to S and K respectively, reusing code from earlier

L_col = [7.73418504191466;5.43399010702799;7.27684433194986;6.41466622438756;4.95146824022541;4.66796093161001];

S_pred = [4.70698050707048;3.20680852777293;5.67072060263290;5.24803574040497;3.96194992279263;2.94670511125058];

S_col = [4.66939281716240;3.07105280591794;5.57701296150031;5.11402048350142;3.83378870449483;2.81821848732995];

L_dev = [0.614863717531972;1.24322579612817;1.06706597959258;1.09560321636564;1.18818956356150;1.24463359873474];

S_dev_pred = [0.582266240110145;1.03538243516567;0.755937258197191;0.933109580017290;0.943292840072858;1.17292969560994];

S_dev = [0.822412579665220;1.07895847082380;0.818390914390835;0.862367924581365;0.955572229496762;1.26598377747912];

shuffle_vector = [5,6,4,3,2,1];

L_col = L_col(shuffle_vector);
S_col = S_col(shuffle_vector);
L_dev = L_dev(shuffle_vector);
S_dev = S_dev(shuffle_vector);
S_dev_pred = S_dev_pred(shuffle_vector);
S_pred = S_pred(shuffle_vector);

exp_beta = S_pred-L_col;

xlim_array = [0,3;0,3;1,3;1,4;0,4;0,4];

LN_S_mu = [];
LN_S_sigma = [];

LN_L_mu = [];
LN_L_sigma = [];

LN_S_pred_mu = [];
LN_S_pred_sigma = [];

for i=2:7
    
    figure('Renderer', 'painters', 'Position', [10 10 1200 600])

    mu = [];
    sigma = [];
    deviation = [];
    mean_val = [];
    mse = [];

    legend_strings = [];
    hold on
    
    subtable = table_data{i};
    
    xdata = double(subtable{:,5});
    
    plot_data = xdata;
    
    plot_data(plot_data==0) = [];
    plot_data(isnan(plot_data)) = [];
    plot_data(isinf(plot_data)) = [];
    
    [xBinMeans,yBinMeans,y_vals,mu_val,sigma_val] = fitlogncdf(plot_data,30);
    mu = [mu;mu_val];
    sigma = [sigma;sigma_val];
    deviation = [deviation;std(log(plot_data))];
    mean_val = [mean_val;mean(log(plot_data))];
    
    LN_S_mu = [LN_S_mu;mu_val];
    LN_S_sigma = [LN_S_sigma;sigma_val];
    
    mu_val = S_pred(i-1);
    sigma_val = S_dev_pred(i-1);
    pd = makedist('Lognormal','mu',mu_val,'sigma',sigma_val);
    
    LN_S_pred_mu = [LN_S_pred_mu;mu_val];
    LN_S_pred_sigma = [LN_S_pred_sigma;sigma_val];
    
    [N,edges] = histcounts(log(plot_data),100);

    edge_difs = exp(edges(2:end)) - exp(edges(1:end-1));
    edges = exp((edges(2:end) + edges(1:end-1)) ./ 2);
    x = edges;
    y = pdf(pd,x);
    y_fit = y.*(edge_difs.*sum(N));
    y_empirical = N;    
    y_max = max(y_fit);
    y_fit = y_fit ./ y_max;
    y_empirical = y_empirical ./ y_max;
    
    x_spline = logspace(log10(min(edges)/2),log10(max(edges)*2),1000);
    
    s = spline(log(x),y_fit,log(x_spline));
    
    %plot(x,y_fit,'-.','LineWidth',1,'MarkerFaceColor',0.7.*colors(i,:),'Color',0.7.*colors(i,:));
    h2 = plot(x_spline,s,'-','LineWidth',1,'MarkerFaceColor',1.*colors(i,:),'Color',1.*colors(i,:));

    
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
    
    if Use_Fits
        %plot(x,y_fit,'-','LineWidth',1.5,'MarkerFaceColor',colors(i,:),'Color',colors(i,:),'HandleVisibility','off');
        %scatter(edges,y_empirical,20,colors(i,:),'filled')
        h1 = scatter(edges,y_empirical,80,colors(i,:),'filled');
    else
        h1 = scatter(edges,y_empirical,80,colors(i,:),'filled');
    end
    mu_val = L_col(i-1);
    sigma_val = L_dev(i-1);
    
    LN_L_mu = [LN_L_mu;mu_val];
    LN_L_sigma = [LN_L_sigma;sigma_val];
    
    pd = makedist('Lognormal','mu',mu_val,'sigma',sigma_val);
    
    [N,edges] = histcounts(log(plot_data),100);

    edge_difs = exp(edges(2:end)) - exp(edges(1:end-1));
    edges = exp((edges(2:end) + edges(1:end-1)) ./ 2);
    x = edges;
    y = pdf(pd,x);
    y_fit = y.*(edge_difs.*sum(N));
    y_empirical = N;    
    y_max = max(y_fit);
    y_fit = y_fit ./ y_max;
    y_empirical = y_empirical ./ y_max;
    
    x_spline = logspace(log10(min(edges)/5),log10(max(edges)*5),1000);
    
    s = spline(log(x),y_fit,log(x_spline));
    
    %plot(x,y_fit,'-.','LineWidth',1,'MarkerFaceColor',0.7.*colors(i,:),'Color',0.7.*colors(i,:));
    %plot(x_spline,s,'-.','LineWidth',1,'MarkerFaceColor',0.7.*colors(i,:),'Color',0.7.*colors(i,:));

    mse = [mse;goodnessOfFit(y_empirical,y_fit,'MSE')];
    legend_strings = [legend_strings,label_strings(i)];
    if Use_Fits
        %legend_strings = [legend_strings,strcat("Fit: \mu=",string(round(exp(mu_val),2)),", \sigma=",string(round(sigma_val,2)))];
    end
    
    hold off
    set(gca, 'XScale', 'log')
    %legend([strcat("Predicted fit for degrees (\beta S^{\alpha}) (\mu=",string(round(LN_S_pred_mu(i-1),2)),",\sigma=",string(round(LN_S_pred_sigma(i-1),2)),")"),strcat("Degree Distribution Data (\mu=",string(round(LN_S_mu(i-1),2)),",\sigma=",string(round(LN_S_sigma(i-1),2)),")"),strcat("Predicted fit for degrees (S) (\mu=",string(round(LN_L_mu(i-1),2)),",\sigma=",string(round(LN_L_sigma(i-1),2)),")")],"Location",'northeast','FontSize',12)
    legend([h1,h2],[strcat("Data (\mu=",string(round(LN_S_mu(i-1),2)),",\sigma=",string(round(LN_S_sigma(i-1),2)),")"),strcat("Predicted (\mu=",string(round(LN_S_pred_mu(i-1),2)),",\sigma=",string(round(LN_S_pred_sigma(i-1),2)),")")],"Location",'northeast')
    xlabel("k");
    ylabel("P(k)");
    %title(strcat(label_strings(i)," Neuron Predicted Degrees vs. Synapses Fits"),'Interpreter','latex');
    xlim([10^xlim_array(i-1,1) 10^xlim_array(i-1,2)])
    ylim([0 1.4])
    
    ax = gca; % Get the current axis handle
    xTicks = ax.XTick; % Get current x-axis tick positions
    yTicks = ax.YTick; % Get current y-axis tick positions

    % Set ticks explicitly to make them fixed
    ax.XTick = xTicks;
    ax.YTick = yTicks;

    ax.TickLength = [0.04, 0.025];
    
    f = gcf;
    f.Position = [10, 10, 600, 700];

    cd ..\..
    cd Figures\Figure3
    
    pos = get(gcf, 'Position');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
    set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
    set(gcf, 'PaperPositionMode', 'auto');
    print(strcat(label_strings(i),'PredDeg.pdf'), '-dpdf', '-painters', '-bestfit');
    
    cd ..\..
    cd FigureGeneration\HelperFunctions
    
end