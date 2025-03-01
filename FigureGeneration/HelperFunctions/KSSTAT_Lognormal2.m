function [KSSTATs] = KSSTAT_Lognormal2(Data)
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
    
    [f,x] = ecdf(Data);
    
    x = x ./ max(x);
    
    x(1) = [];
    f(1) = [];
    
    Data = Data(4:end-3);
    
    x = x(4:end-3);
    f = f(4:end-3);
    
    x = log(x);
    f = log(f);
    
    pd = polyfit(x,f,1);
    
    Data_Fit = polyval(pd,x);
    
    f2 = f;
    
    [h,p,ks2stat] = kstest2(f2,Data_Fit);

    KSSTATs = [KSSTATs;ks2stat];
    
    pd = fitdist(Data,'Lognormal');

    [h,p,ksstat,cv] = kstest(Data,'CDF',pd);

    KSSTATs = [KSSTATs;ksstat];

end

if length(KSSTATs) ~= 4
    KSSTATs = NaN(1,4);
end
    
end

