function [xBinMeans,yBinMeans,y_vals,mu,sigma] = fitlogncdf(plot_data,bincounts)
%FITLOGNCDF Summary of this function goes here
%   Detailed explanation goes here

f1 = figure(100);

h = cdfplot(plot_data);

x_CDF = h.XData;
y_CDF = h.YData;

y_CDF(isinf(x_CDF)) = [];
x_CDF(isinf(x_CDF)) = [];
y_CDF(isnan(x_CDF)) = [];
x_CDF(isnan(x_CDF)) = [];
y_CDF(x_CDF==0) = [];
x_CDF(x_CDF==0) = [];

x_CDF = log(x_CDF);

[xBinMeans,yBinMeans,Counts] = LinBinsCounts(x_CDF,y_CDF,bincounts);

yBinMeans(isnan(xBinMeans)) = [];
xBinMeans(isnan(xBinMeans)) = [];

xBinMeans = exp(xBinMeans);

pHat = lognfit(plot_data);

y_vals = (1/2)*abs(1+erf((log(xBinMeans)-pHat(1))./(pHat(2)*sqrt(2))));

close(f1);

mu = pHat(1);

sigma = pHat(2);

end

