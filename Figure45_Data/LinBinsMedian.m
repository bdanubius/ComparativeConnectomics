function [xBinMeans,yBinMeans,y25,y75] = LinBinsMedian(x,y,numBins)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Compute bin edges for x-axis
xEdges = linspace(min(x), max(x), numBins + 1);

% Preallocate arrays for binned data
xBinMeans = zeros(numBins, 1);
yBinMeans = zeros(numBins, 1);
y25 = zeros(numBins, 1);
y75 = zeros(numBins, 1);

% Compute mean of y-values in each x-bin
xBinData = cell(numBins,1);
yBinData = cell(numBins,1);
for i = 1:numBins
    binIndices = x >= xEdges(i) & x < xEdges(i+1);
    xBinMeans(i) = nanmean(x(binIndices));
    yBinMeans(i) = nanmedian(y(binIndices));
    y25(i) = prctile(y(binIndices),25);
    y75(i) = prctile(y(binIndices),75);
    xBinData{i} = x(binIndices);
    yBinData{i} = y(binIndices);
    
end

end

