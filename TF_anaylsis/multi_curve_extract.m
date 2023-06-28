function [harmonics] = multi_curve_extract(tfrsq, tfrsqtic, sampling_rate, freqRng, numHarmonics,lambda)

disp('Curve extraction...');
%lambda = 30.0; 
gamma = 9000.0;
harmonics = zeros(numHarmonics, size(tfrsq,2));

% fundamental frequency
idx0 = find(tfrsqtic*sampling_rate>0.8 & tfrsqtic*sampling_rate<3.5);
[f] = CurveExt_M(abs(tfrsq(idx0,:))', lambda);
f = f + idx0(1) - 1;
harmonics(1,:) = f;

% the 2nd harmonic
frng_low = 0.5; frng_high = 0.5;
rng_low = ceil(frng_low*size(tfrsq,1)/freqRng);
rng_high = ceil(frng_high*size(tfrsq,1)/freqRng);
tfrsq1 = tfrsq;
for t = 1:length(f)
    tfrsq1(f(t)+floor(prctile(f,90))+rng_high+1:end, t) = 0;
    tfrsq1(1:f(t)+floor(prctile(f,10))-rng_low, t) = 0;
end
idx = find(tfrsqtic*sampling_rate >= min(2*tfrsqtic(harmonics(1,:))*sampling_rate)-frng_low);
[c] = CurveExt_M_v2(abs(tfrsq1(idx,:))', lambda, 2, idx(1), tfrsqtic, f, gamma, tfrsqtic);
c = c + idx(1) - 1;
harmonics(2,:) = c;

for k = 3:numHarmonics
    frng = 0.5;
    tfrsq1 = tfrsq;
    rng = ceil(frng*size(tfrsq1,1)/freqRng);
    for t = 1:length(f)
        tfrsq1(c(t)+floor(prctile(f,75))+rng+1:end, t) = 0;
        tfrsq1(1:c(t)+floor(prctile(f,25))-rng, t) = 0;
    end
    idx = find(tfrsqtic*sampling_rate >= min(tfrsqtic(f)+tfrsqtic(c))*sampling_rate-frng);
%     [c] = CurveExt_M_v2(abs(tfrsq1(idx,:))', lambda, k/(k-1), idx(1), tfrsqtic, c, gamma, tfrsqtic);
    [c] = CurveExt(abs(tfrsq1(idx,:))', 10.0);
    c = c + idx(1) - 1;
    harmonics(k,:) = c;
end
% idx0 = find(tfrsqtic2*sampling_rate>4.5 & tfrsqtic2*sampling_rate<6.0);
% [c] = CurveExt_M(abs(tfrsq2(idx0,:))', lambda);
% c = c + idx0(1) - 1;
% harmonics(5,:) = c;
% for k = 6:numHarmonics
%     frng = 0.5;
%     tfrsq1 = tfrsq2;
%     rng = ceil(frng*size(tfrsq2,1)/freqRng);
%     for t = 1:length(f)
%         tfrsq1(c(t)+floor(prctile(f,75))+rng+1:end, t) = 0;
%         tfrsq1(1:c(t)+floor(prctile(f,25))-rng, t) = 0;
%     end
%     idx = find(tfrsqtic2*sampling_rate >= min(tfrsqtic(f)+tfrsqtic2(c))*sampling_rate-frng);
% %     [c] = CurveExt_M_v2(abs(tfrsq1(idx,:))', lambda, k/(k-1), idx(1), tfrsqtic, c, gamma, tfrsqtic);
%     [c] = CurveExt_M(abs(tfrsq1(idx,:))', 50.0);
%     c = c + idx(1) - 1;
%     harmonics(k,:) = c;
% end

end