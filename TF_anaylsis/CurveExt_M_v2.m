%
% Curve extraction. Written by J. Lu.
%
% version 2:
% fund = the reference harmonics, e.g. the fundamental.
% k = the kth multiple of the fund
% gamma = penalty of the similarities between answer and k*[fund]
% shift, tic2: Used to calculate the actual frequency(Hz) of the answer.
% tic1: Used to calculate the actual frequency of [fund].

function [c, FVal] = CurveExt_M_v2(P, lambda, k, shift, tic1, fund, gamma, tic2)
% k: kth harmonics; c: 1st harmonics; gamma: penalty
% flag: extract 1st harmonics or not.
gamma = gamma*100;

eps = 1e-8;

E = P;
E = E/sum(sum(E));
E = -log(E+eps);

[m,n] = size(E);
FVal = inf(m,n);	% m is time, n is freq
FVal(1,:) = E(1,:);
c = zeros(m,1);

for ii = 2:m		%% time
    for jj = 1:n	%% freq
	    % calculate the penalty term
        for kk = 1:n
            % FVal(ii,jj) = min(FVal(ii,jj), FVal(ii-1,kk)+lambda*(kk-jj).^2);
            FVal(ii,jj) = min(FVal(ii,jj), ...
                FVal(ii-1,kk)+lambda*(kk-jj).^2+gamma*(tic2(kk+shift-1)-k*tic1(fund(ii-1))).^2);
        end
	    % E(ii,jj) is the SST value at time ii and freq jj
        FVal(ii,jj) = FVal(ii,jj) + E(ii,jj);
    end
end

[~, c(m)] = min(FVal(m,:)); %c(m) = Index of the minimum frequency
for ii = m-1:-1:1
    val = FVal(ii+1,c(ii+1)) - E(ii+1,c(ii+1));
    for kk = 1:n
        if (abs(val-FVal(ii,kk)-lambda*(kk-c(ii+1)).^2-gamma*(tic2(kk+shift-1)-k*tic1(fund(ii))).^2)<eps)
            c(ii) = kk;
            break
        end
    end

    if c(ii) == 0
	    c(ii) = round(n/2);
    end
end

return
end
