function [xBinMeans,yBinMeans,Counts] = LinBinsCounts(x,y,numBins)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Compute bin edges for x-axis
xEdges = linspace(min(x), max(x), numBins + 1);

% Preallocate arrays for binned data
xBinMeans = zeros(numBins, 1);
yBinMeans = zeros(numBins, 1);
Counts = zeros(numBins, 1);

% Compute mean of y-values in each x-bin
for i = 1:numBins
    binIndices = x >= xEdges(i) & x < xEdges(i+1);
    xBinMeans(i) = nanmean(x(binIndices));
    yBinMeans(i) = nanmean(y(binIndices));
    Counts(i) = length(nonzeros(binIndices));
end

end

