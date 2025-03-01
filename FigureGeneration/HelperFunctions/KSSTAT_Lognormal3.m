function [KSSTATs] = KSSTAT_Lognormal3(Data)
%KSSTAT_VECTOR Summary of this function goes here
%   Detailed explanation goes here
KSSTATs = [];

try

    pd = fitdist(Data,'Poisson');

    [h,p,ksstat,cv] = kstest(Data,'CDF',pd);

    KSSTATs = [KSSTATs;ksstat];

    pd = fitdist(Data,'Exponential');

    [h,p,ksstat,cv] = kstest(Data,'CDF',pd);

    KSSTATs = [KSSTATs;ksstat];
    
    [N,edges] = histcounts(log(Data));

    edges = (edges(2:end) + edges(1:end-1)) ./ 2;

    N = log(N);

    N = N(~isnan(edges) & ~isinf(edges));
    edges = edges(~isnan(edges) & ~isinf(edges));
    edges = edges(~isnan(N) & ~isinf(N));
    N = N(~isnan(N) & ~isinf(N));

    % Perform linear fit in log-space
    fit_result = polyfit(edges, N, 1); % Linear fit in log-space

    Data_Fit = polyval(fit_result,edges);

    b_fit = -fit_result(1); % Slope corresponds to -b
    b4 = b_fit;
    a_fit = exp(-fit_result(2) / b_fit); % Recover a (if shifted)
    a4 = a_fit;

    edges = exp(edges);
    N = exp(N);
    Data_Fit = exp(Data_Fit);

    numSamples = 250000;

    probabilities_data = N / sum(N);

    data_samples = randsample(edges, numSamples, true, probabilities_data);

    probabilities_fit = Data_Fit / sum(Data_Fit);

    fit_samples = randsample(edges, numSamples, true, probabilities_fit);

    [h,p,ks2stat] = kstest2(data_samples,fit_samples);

    KSSTATs = [KSSTATs;ks2stat];
    
    pd = fitdist(Data,'Lognormal');

    [h,p,ksstat,cv] = kstest(Data,'CDF',pd);

    KSSTATs = [KSSTATs;ksstat];

end

if length(KSSTATs) ~= 4
    KSSTATs = NaN(1,4);
end
    
end

