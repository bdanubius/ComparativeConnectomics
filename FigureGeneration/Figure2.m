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

%% Figure 2A

fig = figure('Renderer', 'painters', 'Position', [10 10 900 600]);

mu = [];
sigma = [];
deviation = [];
mean_val = [];
mse = [];

%subplot1(2,4,'Gap',[0.01 0.08])
subplot1(2,4,'Gap',[0.03 0.16],'Min',[0.1,0.15],'Max',[0.9 1])

legend_strings = [];
for i=1:length(table_data)
    
    subplot1(i)
    
    hold on
    
    subtable = table_data{i};
    
    xdata = double(subtable{:,5});
    
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
    
    
    %  || i==7 || i==6 || i==3
    [N,edges] = histcounts(log(plot_data),20);
    
    linspace_edges = exp(edges);

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
    xlim([10^(-0.9) 10^4])
    ylim([0 1.4])
    
    hold off
    
    xticks(logspace(0,4,5));
    xticklabels(["10^{0}  ","","10^{2}  ","","10^{4}  "]);
    
    if i==2 || i==3 || i==4 || i==6 || i==7 || i==8
        set(gca, 'YTick', []);
    else
        set(gca, 'YTick', [0 0.5 1]);
    end
    
    ax = gca;
    
    ax.FontSize = 20;
    
    ax.TickLength = [0.06, 0.025];
    ax.LineWidth = 0.6;
    
    set(ax, 'Box', 'off');
    
    ax.XTick = ax.XTick;
    ax.YTick = ax.YTick;
    
    title(label_strings(i),'FontSize',16)
    
end
han=axes(fig,'visible','off');
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,"P(k)");
xlabel(han,"k");
sgtitle("Neuron Degree",'Interpreter','latex', 'FontSize', 24);
%sgtitle("Neuron Degree PDF - $P(k)=\frac{1}{k\sigma\sqrt{2\pi}}e^{\frac{-(ln(k)-\mu_{*})^2}{2\sigma^2}} $, $\mu_{*}=e^{\mu}$",'Interpreter','latex', 'FontSize', 16);

ax = gca;
ax.LooseInset = [0.08 0 0 0];

cd ..
cd ..
cd Figures/Figure2

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('2A.pdf', '-dpdf', '-painters', '-bestfit');
%print('DegDistAll_NoLine.eps', '-depsc', '-painters');

cd ..
cd ..
cd FigureGeneration/HelperFunctions

%% Figure 2D

fig = figure('Renderer', 'painters', 'Position', [10 10 900 600]);

mu = [];
sigma = [];
deviation = [];
mean_val = [];
mse = [];

subplot1(2,4,'Gap',[0.03 0.16],'Min',[0.1,0.15],'Max',[0.9 1])

legend_strings = [];
for i=1:length(table_data)
    
    subplot1(i)
    
    hold on
    
    subtable = table_data{i};
    
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
    xlim([10^(0) 10^5])
    ylim([0 1.4])
    
    hold off
    
    xticks(logspace(0,5,6));
    xticklabels(["","10^{1}  ","","10^{3}  ","","10^{5}  "]);
    
    if i==2 || i==3 || i==4 || i==6 || i==7 || i==8
        set(gca, 'YTick', []);
    else
        set(gca, 'YTick', [0 0.5 1]);
    end
    
    ax = gca;
    
    ax.FontSize = 20;
    
    ax.TickLength = [0.06, 0.025];
    ax.LineWidth = 0.6;
    
    set(ax, 'Box', 'off');
    
    ax.XTick = ax.XTick;
    ax.YTick = ax.YTick;
    
    title(label_strings(i),'FontSize',16)
    
end
han=axes(fig,'visible','off');
han.XLabel.Visible='on';
han.YLabel.Visible='on';
han.Title.Visible='on';
xlabel(han,"S");
ylabel(han,"P(S)");
%title(han,"Neuron Synapses PDF - $P(S)=\frac{1}{S\sigma\sqrt{2\pi}}e^{\frac{-(ln(S)-\mu_{*})^2}{2\sigma^2}} $, $\mu_{*}=e^{\mu}$",'Interpreter','latex', 'FontSize', 16);
sgtitle("Strength (Synapses)",'Interpreter','latex', 'FontSize', 24);
%sgtitle("Neuron Synapses PDF - $P(S)=\frac{1}{S\sigma\sqrt{2\pi}}e^{\frac{-(ln(S)-\mu_{*})^2}{2\sigma^2}} $, $\mu_{*}=e^{\mu}$",'Interpreter','latex', 'FontSize', 16);

ax = gca;
ax.LooseInset = [0.08 0 0 0];

cd ..
cd ..
cd Figures/Figure2

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('2D.pdf', '-dpdf', '-painters', '-bestfit');

cd ..
cd ..
cd FigureGeneration/HelperFunctions

%% Figure 2G

fig = figure('Renderer', 'painters', 'Position', [10 10 900 600]);

mu = [];
sigma = [];
deviation = [];
mean_val = [];
mse = [];

subplot1(2,4,'Gap',[0.03 0.16],'Min',[0.1,0.15],'Max',[0.9 1])

legend_strings = [];
for i=1:length(table_data)
    
    subplot1(i)
    
    hold on
    
    subtable = table_data{i};
    
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
    xlim([10^(0) 10^5])
    ylim([0 1.4])
    
    hold off
    
    xticks(logspace(0,5,6));
    xticklabels(["","10^{1}  ","","10^{3}  ","","10^{5}  "]);
    
    if i==2 || i==3 || i==4 || i==6 || i==7 || i==8
        set(gca, 'YTick', []);
    else
        set(gca, 'YTick', [0 0.5 1]);
    end
    
    ax = gca;
    
    ax.FontSize = 20;
    
    ax.TickLength = [0.06, 0.025];
    ax.LineWidth = 0.6;
    
    set(ax, 'Box', 'off');
    
    ax.XTick = ax.XTick;
    ax.YTick = ax.YTick;
    
    title(label_strings(i),'FontSize',16)
    
end
han=axes(fig,'visible','off');
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,"P(L)");
%xlabel(han,"L (\mum)");
xlabel(han,"$L (\mu m)$",'Interpreter','Latex');
%sgtitle("Neuron Length PDF - $P(L)=\frac{1}{L\sigma\sqrt{2\pi}}e^{\frac{-(ln(L)-\mu_{*})^2}{2\sigma^2}} $, $\mu_{*}=e^{\mu}$",'Interpreter','latex', 'FontSize', 16);
sgtitle("Neuron Length",'Interpreter','latex', 'FontSize', 24);

ax = gca;
ax.LooseInset = [0.08 0 0 0];

cd ..
cd ..
cd Figures/Figure2

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('2G.pdf', '-dpdf', '-painters', '-bestfit');

cd ..
cd ..
cd FigureGeneration/HelperFunctions

%% Figure 2B

figure('Renderer', 'painters', 'Position', [10 10 900 600])

mu = [];
sigma = [];
mse = [];

all_plot_data = [];

legend_strings = [];
hold on
for i=1:length(table_data)
    
    subtable = table_data{i};
    
    xdata = double(subtable{:,5});
    
    plot_data = xdata;
    plot_data(plot_data==0) = [];
    plot_data(isnan(plot_data)) = [];
    plot_data(isinf(plot_data)) = [];
    
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
%     if i~=1 && i~=8
%         [N,edges] = histcounts(log(plot_data),100);
%     else
%         [N,edges] = histcounts(log(plot_data));
%     end
%     if  i==1 || i == 8 || i==6 || i==2 || i==3
%         [N,edges] = histcounts(log(plot_data),20);
%     else
%         [N,edges] = histcounts(log(plot_data),100);
%     end

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
    
    if i==5
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
    
    if i==5
        plot(x,y_fit,'-','LineWidth',1,'MarkerFaceColor',[0 0 0],'Color',[0 0 0],'HandleVisibility','off');
    end
    
    mse = [mse;goodnessOfFit(y_empirical,y_fit,'MSE')];
    legend_strings = [legend_strings,label_strings(i)];
    %legend_strings = [legend_strings,strcat("Fit: \mu=",string(round(exp(mu_val),2)),", \sigma=",string(round(sigma_val,2)))];
end
hold off
set(gca, 'XScale', 'log')
legend(legend_strings,"Location",'northeast')
xlabel("$\overline{k} = exp(\frac{\log(k) - \mu_{k}}{\sigma_{k}})$", 'Interpreter','Latex');
ylabel("$P(\overline{k})$", 'Interpreter','Latex');
%title("Neuron Degree Rescaled - $P(\overline{k})=\frac{1}{\overline{k}\sigma\sqrt{2\pi}}e^{\frac{-(ln(\overline{k})-\mu)^2}{2\sigma^2}}$",'Interpreter','latex', 'FontSize', 16);

xlim([10^(-2) 10^2])
ylim([0 1.4])

cd ..
cd ..
cd Figures/Figure2

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('2B.pdf', '-dpdf', '-painters', '-bestfit');

cd ..
cd ..
cd FigureGeneration/HelperFunctions

%% Figure 2E

figure('Renderer', 'painters', 'Position', [10 10 900 600])

mu = [];
sigma = [];
deviation = [];
mean_val = [];
mse = [];

legend_strings = [];
hold on
for i=1:length(table_data)
    
    subtable = table_data{i};
    
    xdata = double(subtable{:,3});
    
    plot_data = xdata;
    
    plot_data(plot_data==0) = [];
    plot_data(isnan(plot_data)) = [];
    plot_data(isinf(plot_data)) = [];
    
    plot_data = log(plot_data);
    plot_data = plot_data - mean(plot_data);
    plot_data = plot_data ./ std(plot_data);
    plot_data = exp(plot_data);
    
    [xBinMeans,yBinMeans,y_vals,mu_val,sigma_val] = fitlogncdf(plot_data,30);
    mu = [mu;mu_val];
    sigma = [sigma;sigma_val];
    deviation = [deviation;std(log(plot_data))];
    mean_val = [mean_val;mean(log(plot_data))];
    pd = makedist('Lognormal','mu',mu_val,'sigma',sigma_val);
    
%     if i==8 || i==1
%         [N,edges] = histcounts(log(plot_data),20);
%     else
%         [N,edges] = histcounts(log(plot_data),50);
%     end

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
        %plot(x,y_fit,'-','LineWidth',1,'MarkerFaceColor',colors(i,:),'Color',colors(i,:));
        %plot(x,y_fit,'-','LineWidth',1.5,'MarkerFaceColor',colors(i,:),'Color',colors(i,:),'HandleVisibility','off');
        %scatter(edges,y_empirical,20,colors(i,:),'filled')
        scatter(edges,y_empirical,40,colors(i,:),'filled');
    else
        scatter(edges,y_empirical,40,colors(i,:),'filled');
    end
    
    if i==7
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
    
    if i==7
        plot(x,y_fit,'-','LineWidth',1,'MarkerFaceColor',[0 0 0],'Color',[0 0 0],'HandleVisibility','off');
    end
    mse = [mse;goodnessOfFit(y_empirical,y_fit,'MSE')];
    legend_strings = [legend_strings,label_strings(i)];
    if Use_Fits
        %legend_strings = [legend_strings,strcat("Fit: \mu=",string(round(exp(mu_val),2)),", \sigma=",string(round(sigma_val,2)))];
    end
end
hold off
set(gca, 'XScale', 'log')
%legend(legend_strings,"Location",'best')
xlabel("$\overline{S} = exp(\frac{\log(S) - \mu_{S}}{\sigma_{S}})$", 'Interpreter','Latex');
ylabel("$P(\overline{S})$", 'Interpreter','Latex');
%title("Neuron Synapses Rescaled - $P(\overline{S})=\frac{1}{\overline{S}\sigma\sqrt{2\pi}}e^{\frac{-(ln(\overline{S})-\mu)^2}{2\sigma^2}}$",'Interpreter','latex', 'FontSize', 16);
xlim([10^(-2) 10^2])
ylim([0 1.4])

cd ..
cd ..
cd Figures/Figure2

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('2E.pdf', '-dpdf', '-painters', '-bestfit');

cd ..
cd ..
cd FigureGeneration/HelperFunctions

%% Figure 2H

figure('Renderer', 'painters', 'Position', [10 10 900 600])

mu = [];
sigma = [];
mse = [];

all_plot_data = [];

legend_strings = [];
hold on
for i=1:length(table_data)
    
    subtable = table_data{i};
    
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
xlabel("$\overline{L} = exp(\frac{\log(L) - \mu_{L}}{\sigma_{L}})$", 'Interpreter','Latex');
ylabel("$P(\overline{L})$", 'Interpreter','Latex');
%title("Neuron Lengths Rescaled - $P(\overline{L})=\frac{1}{\overline{L}\sigma\sqrt{2\pi}}e^{\frac{-(ln(\overline{L})-\mu)^2}{2\sigma^2}}$",'Interpreter','latex', 'FontSize', 16);
xlim([10^(-2) 10^2])
ylim([0 1.4])

cd ..
cd ..
cd Figures/Figure2

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('2H.pdf', '-dpdf', '-painters', '-bestfit');

cd ..
cd ..
cd FigureGeneration/HelperFunctions

%% 2C

KSSTATs_All = cell(length(table_data),1);

for i=1:length(table_data)
    
    subtable = table_data{i};
    
    KSSTATs_Array = [];
    
    xdata = double(subtable{:,2});
    
    plot_data = xdata;
    
    plot_data(plot_data==0) = [];
    plot_data(isnan(plot_data)) = [];
    plot_data(isinf(plot_data)) = [];
    
    KSSTATs = KSSTAT_Lognormal2(plot_data);
    
    KSSTATs_Array = [KSSTATs_Array,KSSTATs];
    
    xdata = double(subtable{:,3});
    
    plot_data = xdata;
    
    plot_data(plot_data==0) = [];
    plot_data(isnan(plot_data)) = [];
    plot_data(isinf(plot_data)) = [];
    
    KSSTATs = KSSTAT_Lognormal2(plot_data);
    
    KSSTATs_Array = [KSSTATs_Array,KSSTATs];
    
    xdata = double(subtable{:,4});
    
    plot_data = xdata;
    
    plot_data(plot_data==0) = [];
    plot_data(isnan(plot_data)) = [];
    plot_data(isinf(plot_data)) = [];
    
    KSSTATs = KSSTAT_Lognormal2(plot_data);
    
    KSSTATs_Array = [KSSTATs_Array,KSSTATs];
    
    xdata = double(subtable{:,5});
    
    plot_data = xdata;
    
    plot_data(plot_data==0) = [];
    plot_data(isnan(plot_data)) = [];
    plot_data(isinf(plot_data)) = [];
    
    KSSTATs = KSSTAT_Lognormal2(plot_data);
    
    KSSTATs_Array = [KSSTATs_Array,KSSTATs];
    
    xdata = double(subtable{:,6});
    
    plot_data = xdata;
    
    plot_data(plot_data==0) = [];
    plot_data(isnan(plot_data)) = [];
    plot_data(isinf(plot_data)) = [];
    
    KSSTATs = KSSTAT_Lognormal2(plot_data);
    
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

%% Deg KS

default_colors = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560];

fig = figure('Renderer', 'painters', 'Position', [10 10 900 600]);

for i=4:4
    
    hold on
    substats = Quantity_KSSTATs{i};
    
    substats = substats(2:end-1,:);
    
    plot(1:6,substats(:,1),':.','Color',default_colors(1,:),'LineWidth',2,'MarkerSize',25)
    %plot(1:6,KSSTATs_Array(:,1),':o','Color',default_colors(1,:),'LineWidth',1.5,'MarkerSize',10)
    plot(1:6,substats(:,2),':.','Color',default_colors(2,:),'LineWidth',2,'MarkerSize',25)
    %plot(1:6,KSSTATs_Array(:,2),':o','Color',default_colors(2,:),'LineWidth',1.5,'MarkerSize',10)
    plot(1:6,substats(:,3),':.','Color',default_colors(3,:),'LineWidth',2,'MarkerSize',25)
    %plot(1:6,KSSTATs_Array(:,3),':o','Color',default_colors(3,:),'LineWidth',1.5,'MarkerSize',10)
    plot(1:6,substats(:,4),'.-','Color',default_colors(4,:),'LineWidth',2.5,'MarkerSize',25)
    %plot(1:6,KSSTATs_Array(:,4),'-o','Color',default_colors(4,:),'LineWidth',2,'MarkerSize',10)
    %plot(KSSTATs_All{i}',':','LineWidth',2)
    hold off
    %legend(["Poisson","Poisson (K*)","Exponential","Exponential (K*)","Power Law","Power Law (K*)","LN Fit","LN Fit (K*)"],"Location","northwest","FontSize",8)
    legend(["Poisson","Exponential","Power Law","LN Fit"],"Location","northeast")
    %ylabel("KS Statistic (Lower Value Indicates Better Fit)",'FontSize',12)
    %xlabel("Measurement",'FontSize',12)
    xticks([1,2,3,4,5,6])
    ax = gca;
    ax.XAxis.MinorTick = 'off';
    
    %xticklabels(horzcat(label_strings(2),label_strings(3),label_strings(4),label_strings(5),label_strings(6),label_strings(7)))
    
    label_strings2 = strcat(label_strings,"(");
    label_strings2 = extractBefore(label_strings2,"(");
    
    xticklabels(label_strings2(2:7))
    
    ax.FontSize = 16;
    
    %title(Quantity_Labels(i),'FontSize',16)
end

ylim([0 0.8])

han=axes(fig,'visible','off');
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'KS Statistic');
%xlabel(han,'Dataset');

cd ..
cd ..
cd Figures/Figure2

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('2C.pdf', '-dpdf', '-painters', '-bestfit');

cd ..
cd ..
cd FigureGeneration/HelperFunctions

%% 2F

default_colors = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560];

fig = figure('Renderer', 'painters', 'Position', [10 10 900 600]);

for i=2:2
    
    hold on
    substats = Quantity_KSSTATs{i};
    plot(substats(:,1),':.','Color',default_colors(1,:),'LineWidth',2,'MarkerSize',25)
    %plot(2:8,KSSTATs_Array(:,1),':o','Color',default_colors(1,:),'LineWidth',1.5,'MarkerSize',10)
    plot(substats(:,2),':.','Color',default_colors(2,:),'LineWidth',2,'MarkerSize',25)
    %plot(2:8,KSSTATs_Array(:,2),':o','Color',default_colors(2,:),'LineWidth',1.5,'MarkerSize',10)
    plot(substats(:,3),':.','Color',default_colors(3,:),'LineWidth',2,'MarkerSize',25)
    %plot(2:8,KSSTATs_Array(:,3),':o','Color',default_colors(3,:),'LineWidth',1.5,'MarkerSize',10)
    plot(substats(:,4),'.-','Color',default_colors(4,:),'LineWidth',2.5,'MarkerSize',25)
    %plot(2:8,KSSTATs_Array(:,4),'-o','Color',default_colors(4,:),'LineWidth',2,'MarkerSize',10)
    %plot(KSSTATs_All{i}',':','LineWidth',2)
    hold off
    %legend(["Poisson","Exponential","Power Law","LN Fit"],"Location","northeast")
    %ylabel("KS Statistic (Lower Value Indicates Better Fit)",'FontSize',12)
    %xlabel("Measurement",'FontSize',12)
    xticks([1,2,3,4,5,6,7,8])
    %xticklabels(label_strings)
    
    label_strings2 = strcat(label_strings,"(");
    label_strings2 = extractBefore(label_strings2,"(");
    
    xticklabels(label_strings2)
    
    ax = gca;
    ax.FontSize = 16;
    %title(Quantity_Labels(i),'FontSize',16)
end

ylim([0 0.8])

han=axes(fig,'visible','off');
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'KS Statistic');
%xlabel(han,'Dataset');

cd ..
cd ..
cd Figures/Figure2

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('2F.pdf', '-dpdf', '-painters', '-bestfit');

cd ..
cd ..
cd FigureGeneration/HelperFunctions

%% Figure 2I

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
    xticks([1,2,3,4,5,6,7,8])
    
    label_strings2 = strcat(label_strings,"(");
    label_strings2 = extractBefore(label_strings2,"(");
    
    xticklabels(label_strings2)
    
    ax = gca;
    ax.FontSize = 16;
    %title(Quantity_Labels(i),'FontSize',16)
end

ylim([0 0.8])

han=axes(fig,'visible','off');
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'KS Statistic');
%xlabel(han,'Dataset');

cd ..
cd ..
cd Figures/Figure2

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('2I.pdf', '-dpdf', '-painters', '-bestfit');

cd ..
cd ..
cd FigureGeneration/HelperFunctions