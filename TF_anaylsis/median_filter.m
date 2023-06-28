function [trend]=median_filter(X,K)
%X = raw signal
%==========================================================================
%prepare data
trend = zeros(size(X));
nX = [zeros(10000,1) ;X ;zeros(10000,1)];

%%
for i = 10001:(length(X)+10000)
    m = median(nX(i-K:i+K));
    trend(i-10000) = m;
end


 





