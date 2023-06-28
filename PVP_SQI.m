function [SQI] = PVP_SQI(PVP_dtr,sampling_rate)
%PVP_SQI Summary of this function goes here
%Input:
% 'PVP_dtr': the detrend PVP signal.
% 'sampling_rate': the sampling rate of the PVP signal.
%Output:
% 'SQI': the proposed SQI index for each ten second epoch
addpath('TF_anaylsis');

basicTF.win = sampling_rate*10+1;% window length of SST
basicTF.hop = sampling_rate/20;%window hopping frequency
basicTF.fs = sampling_rate ;
basicTF.fr = 0.01; % frequency resolution
% advanced parameters for STFT
advTF.HighFreq = 10/basicTF.fs; % highest frequency/sampling freq
advTF.LowFreq = 0.01/basicTF.fs;
%% Syncro Squeezing Transform of the PVP signal
[~, ~, PVP_tfrsq, ~, tfrsqtic] = ConceFT_sqSTFT_C(PVP_dtr, 0, advTF.HighFreq, basicTF.fr/basicTF.fs,basicTF.hop, basicTF.win, 1, 6, 1, 1, 0) ;
%% extract the first three harmonics of the PVP signal
tic;
[harmonics] = multi_curve_extract(PVP_tfrsq, tfrsqtic, sampling_rate, 10, 3,10);
toc;
%% reconstruct cardiac signal;
card_tfrsq = zeros(size(PVP_tfrsq)); % the tfr used for reconstruction

band_Hz = 0.2/2;% bandwith of cardiac component in Hz.

for l = 1:3
    for k = 1:size(PVP_tfrsq,2)
        card_tfrsq(harmonics(l,k)-band_Hz/basicTF.fr+1:harmonics(l,k)+band_Hz/basicTF.fr,k) = PVP_tfrsq(harmonics(l,k)-band_Hz/basicTF.fr+1:harmonics(l,k)+band_Hz/basicTF.fr,k);
    end
end

alpha = tfrsqtic(2)-tfrsqtic(1);
[h, ~] = hermf(basicTF.win, 1, 6) ;
coeff = h(0.5*(basicTF.win+1)) ;%
C = 2 * alpha / coeff ; % the normalizing coeffient for the reconstruction of SST

recon = zeros(1,size(card_tfrsq,2));% reconstruct cardiac signal with sampling rate same as the hopping frequency = basicTF.hop
for kk = 1: size(card_tfrsq,2)
    recon(kk) = C * sum(card_tfrsq(:,kk),1) ;
end
% resample the reconstructed cardiac signal into the samme sampling rate as
% the PVP signal.
recon_cardiac = interp1([1:basicTF.hop:length(PVP_dtr)]/sampling_rate,real(recon),[1:1:length(PVP_dtr)]/sampling_rate,'pchip')';
%% Get SQI
PVP_noise = PVP_dtr - recon_cardiac;

SQI_tmp = zeros(size(PVP_tfrsq,2),1);% the SQI index with sampling rate = hopping frequency.
for k = 1:size(PVP_tfrsq,2)-1
    cardiac_power = zeros(1,3);
    for l = 1:3
        cardiac_power(l)= sum(PVP_tfrsq(harmonics(l,k)-band_Hz/basicTF.fr+1:harmonics(l,k)+band_Hz/basicTF.fr,k));
    end
    SQI_tmp(k) = sum((abs(C*cardiac_power)).^2)/(sum((abs(C*cardiac_power)).^2)+sum(PVP_noise((k-1)*basicTF.hop+1:k*basicTF.hop).^2)/basicTF.hop);
end
epoch_len = 10;% time length of PVP epochs, unit=second.
epoch_num = floor(length(PVP_dtr)./(sampling_rate*epoch_len));
SQI = zeros(epoch_num,1);
for i = 1:epoch_num
    s = (i-1)*sampling_rate/basicTF.hop*epoch_len;
    t = i*sampling_rate/basicTF.hop*epoch_len;
    SQI(i)= median(SQI_tmp(s+1:t));
end

