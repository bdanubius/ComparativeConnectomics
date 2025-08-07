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

clear
close all

%% Figure 4D

cd ..

% This is repurposed from a script where I processed the Xn for each neuron
% tree. It looks more complicated than it is, I'm just calculating this for
% one neuron, and superimposing that tree visualization with the mesh

load_giulia_data = false;

N_colors = 5;
color_set = [linspace(0,1,N_colors)',linspace(1,0,N_colors)',linspace(1,0,N_colors)'];

mesh_number = 1;

% Visualization #s
% 3: 720575940638731427
% 5: 720575940625163651
% 10: 720575940619146337

cd Figure45_Data

Neuron_IDs = ["720575940638731427","720575940625163651","720575940619146337"];

[F,V] = stlread(strcat(Neuron_IDs(mesh_number),".stl"));

V = V ./ 1000;

if ~load_giulia_data

    cd Soma_Tree_Data

    files = dir("*.mat");

    tree_sizes = [];
    survival_probability = [];
    tree_depth = [];
    max_tree_size = [];
    xn = {};
    mu = [];
    survival_probability2 = [];

    for file_iter = 1:size(files,1)

        disp(file_iter)

        filename = files(file_iter).name;

        all_tree_data = load(filename);
        all_tree_data = all_tree_data.All_Chain_Data;

        for i=3:size(all_tree_data,1)
            if ~isempty(all_tree_data{i,1})
                
                st = all_tree_data{i,6};
                
%                 [row,col] = find(all_tree_data{i,7});
%                 
%                 st = horzcat(row(row<col),col(row<col));
                
                tree_sizes = [tree_sizes;length(unique(vertcat(st(:,1),st(:,2))))-1];
                G = graph(st(:,1),st(:,2));
                
                [bin,binsize] = conncomp(G);
                
                idx = binsize(bin) == max(binsize);
                SG = subgraph(G, idx);
                
                %% Load Sample Neuron Data
                
                sample_neuron_data = load("Sample_Neuron_Plot_Data.mat");
                sample_neuron_data = sample_neuron_data.outdata;
                
                x_data_sample_neuron = sample_neuron_data{1};
                y_data_sample_neuron = sample_neuron_data{2};
                x_fill = sample_neuron_data{3};
                y_fill = sample_neuron_data{4};
                x_mean = sample_neuron_data{5};
                y_mean = sample_neuron_data{6};
                
                %% Continue
                
                % This specifically is 4D
                
                fig = figure('Renderer', 'painters', 'Position', [10 10 400 600]);

                p = plot(SG,'Layout','layered');
                
                hold on
                
                y_list = p.XData';
                x_list = p.YData';
                x_list = 23 - x_list;
                p.XData = x_list';
                p.YData = y_list';
                scatter(x_list,y_list,10,[x_list/max(x_list),1-x_list/max(x_list),1-x_list/max(x_list)],'filled')
                
                %plot(x_data_sample_neuron,y_data_sample_neuron.*60,'-k');
                
                hold off
                title(strcat("FlyWire Neuron - ",Neuron_IDs(mesh_number)," -"),'FontSize',8)
                for textiter = 1:22
                    %xline(textiter,':k',strcat(string(textiter)),'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','FontSize',6);
                    %text(-6,textiter,strcat("-",string(23-textiter),"-"));
                    [row,col] = find(x_list == textiter);
                    if mod(textiter,2) == 0
                        %text(textiter,61,string(length(row)),'FontSize',8,'HorizontalAlignment','right');
                    else
                        %text(textiter,62.5,string(length(row)),'FontSize',8,'HorizontalAlignment','right');
                    end
                end
                %yline(60,'-k',"X_{m}",'LineWidth',1,'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','FontSize',12)
                %yline(60,'-k',"p_{surv}=1",'LineWidth',1,'LabelHorizontalAlignment','left','FontSize',6)
                %yline(30,'-.k',"p_{surv}=0.5",'LineWidth',0.5,'LabelHorizontalAlignment','left','FontSize',6)
                xlabel("Layer #")
                ylabel("p_{surv} (By Layer)")
                xlim([0 22])
                
                ax = gca; % Get the current axes
                ax.XAxis.Visible = 'on';
                
                xticks([1:1:22]);
                
                a = get(gca,'XTickLabel');
                set(gca,'XTickLabel',a)
                
                cd ..\..
                
                cd Figures\Figure4
                
                pos = get(gcf, 'Position');
                set(gcf, 'PaperUnits', 'inches');
                set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
                set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
                set(gcf, 'PaperPositionMode', 'auto');
                print('4D.pdf', '-dpdf', '-painters', '-bestfit');
                
                % Note that if you were to go beyond the break here, the
                % next figure and the code after it is where I color the
                % neuron mesh for 4C
                
                break
                
                degree_list = degree(SG);
                
                [row,col] = find(degree_list==1);
                
                degree_one = length(row);
                
                degree_three = length(degree_list) - degree_one;
                
                D = distances(G,1);
                
                figure(2)
                hold on
                subpts = all_tree_data{i,5};
                max_D = max(D);
                all_segments = [];
                all_colors = [];
                for subidx=1:length(D)
                    if ~isinf(D(subidx))
                        sub_seg_pts = subpts{subidx};
                        all_segments = [all_segments;sub_seg_pts(:,1:3)];
                        all_colors = [all_colors;repmat([D(subidx)/max_D,1-D(subidx)/max_D,1-D(subidx)/max_D],size(sub_seg_pts,1),1)];
                        %plot3(sub_seg_pts(:,1),sub_seg_pts(:,2),sub_seg_pts(:,3),'-*','Color',[D(subidx)/max_D,1-D(subidx)/max_D,1-D(subidx)/max_D],'LineWidth',3)
                    end
                end
                [Idx,D] = knnsearch(all_segments,V);
                V_colors = zeros(size(V,1),3);
                %scatter3(V(:,1),V(:,2),V(:,3),'.k','MarkerFaceAlpha',.01,'MarkerEdgeAlpha',.01)
                for iter=1:length(all_segments)
                    disp(iter)
                    [row,col] = find(Idx==iter);
                    color = all_colors(iter,:);
                    scatter3(V(row,1),V(row,2),V(row,3),'.','MarkerFaceColor',color,'MarkerEdgeColor',color,'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1)
                    V_colors(row,:) = repmat(color,length(row),1);
                end
                hold off
                title(strcat("FlyWire Neuron - ",Neuron_IDs(mesh_number)," - (Cyan-Red is tree depth distance from soma)"))
                xlabel("x (\mu m)")
                ylabel("y (\mu m)")
                zlabel("z (\mu m)")
                
                axis equal
                
                plywrite('ColoredNeuron.ply', F, V, round(255.*V_colors));
                
                D(isinf(D)) = 0;
                tree_depth = [tree_depth;max(D)];
                max_size = 1;
                for j=1:max(D)
                    max_size = max_size + 2^j;
                end
                max_tree_size = [max_tree_size;max_size];
                
                [cnt_unique,unique_values] = hist(D,unique(D));

                cnt_unique(1) = 1;

                xn{end+1} = cnt_unique;

                if length(cnt_unique) > 2

                    p = polyfit(unique_values,log(cnt_unique),1);
                    
                    mu = [mu;exp(p(1))];
                    
                    if mu < 0
                        break
                    end
                    
                    survival_probability = [survival_probability;exp(p(1))/2];
                    
                    survival_probability2 = [survival_probability2;degree_three / (degree_one + degree_three)];
                    
                end
                
                cd ..
                
                save(strcat(string(Neuron_IDs(3)),"_tree.mat"),'xn');
                
%                 degree_list = degree(G);
%                 [degree1,col1] = find(degree_list==1);
%                 [degree3,col1] = find(degree_list==3);
%                 if ~isempty(degree1) && ~isempty(degree3)
%                     survival_probability = [survival_probability;(length(degree3) ./ (length(degree3)+length(degree1)))];
%                 else
%                     survival_probability = [survival_probability;NaN];
%                 end
            end
        end

    end
    
end

cd ..\..

%% Figure 4E

cd Figure45_Data

clear

N_bins = 20;

dynamic_data = true;

Starting_Layer_MeanSTD = 1;

Simulation_Data = load("Tree_Data_Static_Updated.mat");
Simulation_Data = Simulation_Data.outdata;

Size_Data = Simulation_Data{1};
Final_Data = Simulation_Data{2};

if dynamic_data
    p_data = Simulation_Data{5};
    a_data = Simulation_Data{4};
end

error_matrix = zeros(length(Size_Data{1,1}),size(Size_Data,1));
mean_matrix = zeros(length(Size_Data{1,1}),size(Size_Data,1));
std_matrix = zeros(length(Size_Data{1,1}),size(Size_Data,1));
real_matrix = zeros(length(Size_Data{1,1}),size(Size_Data,1));
if dynamic_data
    p_matrix = zeros(length(Size_Data{1,1}),size(Size_Data,1));
    a_matrix = zeros(length(Size_Data{1,1}),size(Size_Data,1));
end

for i=1:size(Size_Data,1)
    sub_size_real = Size_Data{i,1};
    sub_size_artificial = Size_Data{i,2};
    if dynamic_data
        p_subdata = p_data{i,1};
        a_subdata = a_data{i,1};
    end
    
    error_vector = sqrt(mean((sub_size_artificial - sub_size_real).^2,2));
    
    error_vector(sub_size_real==0) = NaN;
    error_matrix(:,i) = error_vector;
    
    sub_size_artificial(sub_size_artificial==0) = NaN;
    
    mean_matrix(:,i) = nanmean(log(sub_size_artificial),2);
    std_matrix(:,i) = nanstd(log(sub_size_artificial),[],2);
    
    sub_size_real(sub_size_real==0) = NaN;
    
    real_matrix(:,i) = log(sub_size_real);
    
    if dynamic_data
    
        p_subdata(p_subdata==0) = NaN;

        p_matrix(:,i) = p_subdata;

        a_subdata(a_subdata==0) = NaN;

        a_matrix(:,i) = a_subdata;
        
    end
    
end

vector_size = 0;
for i=1:size(error_matrix,1)
    subrow = error_matrix(i,:);
    last_nonNaN = find(~isnan(subrow), 1, 'last');
    vector_size = vector_size + last_nonNaN;
end

percentile_vector = zeros(vector_size,1);
error_percentile_vector = zeros(vector_size,1);
mean_percentile_vector = zeros(vector_size,1);
std_percentile_vector = zeros(vector_size,1);
real_percentile_vector = zeros(vector_size,1);
if dynamic_data
    p_percentile_vector = zeros(vector_size,1);
    a_percentile_vector = zeros(vector_size,1);
end

counter = 0;
for i=1:size(error_matrix,1)
    
    subrow = error_matrix(i,:);
    last_nonNaN = find(~isnan(subrow), 1, 'last');
    subrow = subrow(1:last_nonNaN);
    percentile_vector(counter+1:counter+length(subrow)) = linspace(0,1,length(subrow))';
    error_percentile_vector(counter+1:counter+length(subrow)) = subrow';
    
    subrow = mean_matrix(i,:);
    last_nonNaN = find(~isnan(subrow), 1, 'last');
    subrow = subrow(1:last_nonNaN);
    percentile_vector(counter+1:counter+length(subrow)) = linspace(0,1,length(subrow))';
    mean_percentile_vector(counter+1:counter+length(subrow)) = subrow';
    
    subrow = std_matrix(i,:);
    last_nonNaN = find(~isnan(subrow), 1, 'last');
    subrow = subrow(1:last_nonNaN);
    percentile_vector(counter+1:counter+length(subrow)) = linspace(0,1,length(subrow))';
    std_percentile_vector(counter+1:counter+length(subrow)) = subrow';
    
    subrow = real_matrix(i,:);
    last_nonNaN = find(~isnan(subrow), 1, 'last');
    subrow = subrow(1:last_nonNaN);
    percentile_vector(counter+1:counter+length(subrow)) = linspace(0,1,length(subrow))';
    real_percentile_vector(counter+1:counter+length(subrow)) = subrow';
    
    if dynamic_data
        subrow = p_matrix(i,:);
        last_nonNaN = find(~isnan(subrow), 1, 'last');
        subrow = subrow(1:last_nonNaN);
        percentile_vector(counter+1:counter+length(subrow)) = linspace(0,1,length(subrow))';
        p_percentile_vector(counter+1:counter+length(subrow)) = subrow';

        subrow = a_matrix(i,:);
        last_nonNaN = find(~isnan(subrow), 1, 'last');
        subrow = subrow(1:last_nonNaN);
        percentile_vector(counter+1:counter+length(subrow)) = linspace(0,1,length(subrow))';
        a_percentile_vector(counter+1:counter+length(subrow)) = subrow';
    end
    
    counter = counter + length(subrow);
    
end

% Ignore these figures, these were sanity checks for the SI

figure('Renderer', 'painters', 'Position', [10 10 900 600])

Size_Data = Simulation_Data{1};
Size_Data = Size_Data{1,1};
nhist(log10(Size_Data))
title("FlyWire Log Tree Sizes Histogram");
ylabel("Counts");
xlabel("Log Tree Size");

[xBinMeans,yBinMeans] = LinBins(percentile_vector,error_percentile_vector,N_bins);

figure('Renderer', 'painters', 'Position', [10 10 900 600])
plot(100.*xBinMeans(~isnan(yBinMeans)),yBinMeans(~isnan(yBinMeans)),"-o",'LineWidth',2)
title("Tree Size RMSE Error vs. Layer Percentile (Simulated vs. Real Data)")
ylabel("Tree Size RMSE Error (N=100 Simulated Trees)")
xlabel("Starting Layer Percentile (Percentile Real and Simulated Trees Measured From)")

figure('Renderer', 'painters', 'Position', [10 10 900 600])
[xBinMeans,yBinMeans] = LinBins(percentile_vector,mean_percentile_vector,N_bins);
plot(100.*xBinMeans(~isnan(yBinMeans)),yBinMeans(~isnan(yBinMeans)),"-o",'LineWidth',2)
hold on
[xBinMeans,yBinMeans] = LinBins(percentile_vector,real_percentile_vector,N_bins);
plot(100.*xBinMeans(~isnan(yBinMeans)),yBinMeans(~isnan(yBinMeans)),"-o",'LineWidth',2)
hold off
title("Simulated and Real Tree Size vs. Layer Percentile")
ylabel("Mean Log Tree Size")
xlabel("Starting Layer Percentile (Percentile Real and Simulated Trees Measured From)")
legend(["Real Trees","Simulated Trees"],"Location",'northeast')

[xBinMeans,yBinMeans] = LinBins(percentile_vector,p_percentile_vector,N_bins);

figure('Renderer', 'painters', 'Position', [10 10 900 600])
plot(100.*xBinMeans(~isnan(yBinMeans)),yBinMeans(~isnan(yBinMeans)),"-o",'LineWidth',2)
title("Branching Probability vs. Layer Percentile (Simulated vs. Real Data)")
ylabel("p value (branching probability)")
xlabel("Starting Layer Percentile (Percentile Real and Simulated Trees Measured From)")

[xBinMeans,yBinMeans,xBinData,yBinData] = LinBinsData(percentile_vector,a_percentile_vector,N_bins);

colorset = horzcat(zeros(length(yBinData),1),ones(length(yBinData),1),ones(length(yBinData),1));

figure('Renderer', 'painters', 'Position', [10 10 1100 550])
daviolinplot(yBinData,'xtlabels',string(round(xBinMeans,2)),'colors',colorset,'violinwidth',4,'scatter',1,'jitter',0,'boxwidth',2);
hold on
plot(linspace(1,length(xBinMeans),length(xBinMeans)),yBinMeans,'-r','LineWidth',3)
ylim([0 Inf])
xlabel("Starting Layer Percentile")
ylabel("# of Segments in Layer")
title("FlyWire Tree Segments in Layer by Percentile")

[xBinMeans,yBinMeans,xBinData,yBinData] = LinBinsData(percentile_vector,p_percentile_vector,N_bins);

ySTD = [];
for i=1:length(yBinData)
   ySTD = [ySTD;std(yBinData{i})]; 
end

deviation_plot = cell(2,1);
deviation_plot{1} = 1:1:20;
deviation_plot{2} = ySTD;

figure('Renderer', 'painters', 'Position', [10 10 900 600])
plot(1:1:20,ySTD,'-o','LineWidth',3)
title("Standard Deviation of p per layer")
ylabel("Standard Deviation (p)")
xlabel("Layer Number")

colorset = horzcat(zeros(length(yBinData),1),ones(length(yBinData),1),ones(length(yBinData),1));

figure('Renderer', 'painters', 'Position', [10 10 900 600])
daviolinplot(yBinData,'xtlabels',string(round(xBinMeans,2)),'colors',colorset,'violinwidth',4,'scatter',1,'jitter',0,'boxwidth',2);
hold on
plot(linspace(1,length(xBinMeans),length(xBinMeans)),yBinMeans,'-r','LineWidth',3)
ylim([0 1])
xlabel("Starting Layer Percentile")
ylabel("P(p) ($\frac{X_{n}}{2X_{n-1}}$) per Layer",'Interpreter','Latex')
title("FlyWire P(p) by Layer Percentile")

p_matrix2 = p_matrix;
p_matrix2(p_matrix2==0) = NaN;

violin_data = cell(size(p_matrix2,2),1);
violin_means = zeros(size(p_matrix2,2),1);
violin_std = zeros(size(p_matrix2,2),1);
for i=1:size(p_matrix2,2)
    violin_data{i} = p_matrix2(~isnan(p_matrix2(:,i)),i);
    violin_means(i) = mean(p_matrix2(~isnan(p_matrix2(:,i)),i));
    violin_std(i) = std(p_matrix2(~isnan(p_matrix2(:,i)),i));
end

colorset = horzcat(zeros(size(p_matrix,2),1),ones(size(p_matrix,2),1),ones(size(p_matrix,2),1));

close all

figure('Renderer', 'painters', 'Position', [10 10 900 600])

hold on

p_STD = nanstd(p_matrix,1);

STD_rootN = zeros(1,length(p_STD)); 

for i=1:length(STD_rootN)
    
    STD_rootN(i) = sqrt(sum(~isnan(p_matrix(i,:))));
    
end

%errorbar(1:1:size(p_matrix,2),violin_means,p_STD,'-*k','LineWidth',2)
plot(1:1:size(p_matrix,2),violin_means,'-*k','LineWidth',2)

curve1 = violin_means' + p_STD;
curve2 = violin_means' - p_STD;

x_points = 1:1:size(p_matrix,2);

x2 = [x_points,fliplr(x_points)];
inBetween = [(curve1),fliplr(curve2)];
fill(x2, inBetween,[0 0 0],'FaceAlpha',0.2,'HandleVisibility','off');

plot([20.5 20.51], [0.34 0.67],'k-.','MarkerSize',15,'LineWidth', 1);
scatter(20.5,0.67,50,"^k",'filled')
scatter(20.5,0.34,50,"vk",'filled')
yline(0.5,':',"    p_s^c",'FontSize',18,'LineWidth',1,'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom')

% quiver(20.5,0.5,0,-0.16, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 200,'AutoScaleFactor', 2);
% quiver(20.5,0.5,0,0.16, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 200,'AutoScaleFactor', 2);

text(20.7,0.6,"\sigma_{s}^{j}",'FontSize',22)


tree_data = load("720575940619146337_tree.mat");
tree_data = tree_data.xn;
tree_data = tree_data{1};
%plot(1:1:length(tree_data)-1,tree_data(2:end) ./ (2.*tree_data(1:end-1)),'.-','LineWidth',2,'MarkerSize',20)

tree_data = load("720575940625163651_tree.mat");
tree_data = tree_data.xn;
tree_data = tree_data{1};
%plot(1:1:length(tree_data)-1,tree_data(2:end) ./ (2.*tree_data(1:end-1)),'.-','LineWidth',2,'MarkerSize',20)

tree_data = load("720575940638731427_tree.mat");
tree_data = tree_data.xn;
tree_data = tree_data{1};
%plot(1:1:length(tree_data)-1,tree_data(2:end) ./ (2.*tree_data(1:end-1)),'.-','LineWidth',2,'MarkerSize',20)

%legend(["Dataset Mean ± 1STD","Neuron 720575940619146337","Neuron 720575940625163651","Neuron 720575940638731427"],'FontSize',12)

%legend(["Dataset Mean ± 1STD","Sample Neuron"],"Location",'northeast')

ylim([0.1 1.2])
xlabel("Layer j")
ylabel("$p^j_s = \frac{1}{2}\frac{X_j}{X_{j-1}}$",'Interpreter','Latex')
xlim([0 22])

%title("FlyWire P(p) per Layer (Single Neurons and Dataset Mean)")
%title("FlyWire P(p) per Layer")

x_data_sample_neuron = 1:1:length(tree_data)-1;
y_data_sample_neuron = tree_data(2:end) ./ (2.*tree_data(1:end-1));

x_fill = x2;
y_fill = inBetween;

x_mean = 1:1:size(p_matrix,2);
y_mean = violin_means;

outdata = cell(6,1);

outdata{1} = x_data_sample_neuron;
outdata{2} = y_data_sample_neuron;
outdata{3} = x_fill;
outdata{4} = y_fill;
outdata{5} = x_mean;
outdata{6} = y_mean;

xticks(1:1:21)

cd ..
cd Figures\Figure4

pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
set(gcf, 'PaperPosition', [0 0 pos(3)/100 pos(4)/100]);
set(gcf, 'PaperPositionMode', 'auto');
print('4E.pdf', '-dpdf', '-painters', '-bestfit');
