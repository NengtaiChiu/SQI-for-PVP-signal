clear all;close all;
addpath('TF_anaylsis');
%% generate the simulated data
fs = 100;% sampling_rate
time = [1/fs:1/fs:120]' ;
N = length(time) ;
%% the instantaneous frequency of the simulated signal
if1 = smooth(cumsum(randn(N,1)) ./ fs, 400, 'loess') ;
if1 = 1.2 + 0.2 * if1 ./ max(abs(if1)) ;
if2 = 2.4 + 0.2 * if1 ./ max(abs(if1)) ;
if3 = 0.2 + 0.05 * if1 ./ max(abs(if1)) ;
phi1 = cumsum(if1) / fs ; 
phi2 = cumsum(if2) / fs ; 
phi3 = cumsum(if3) / fs ; 

%% the simulated signal.
trend = smooth(cumsum(randn(N,1))+10, 1000, 'loess') ;
s1 = 1* cos(2*pi*phi1) ; 
s2 = 0.4* cos(2*pi*phi1) ; 
s3 = 3 * cos(2*pi*phi1) ; 
clean = s1+s2+s3+trend;

%% add noise (Gaussian white noise)
sigma = 0.5;
noise = random('T',4,N,1) ;
noise = sigma * noise ; 
%% simulated observed time series
phantom_PVP = clean + noise ;
%% Preprocess signal
[b, a] = butter(6, [0.5, 15] / (100 / 2));% bandpass filter setting.
y = median_filter(phantom_PVP, 6*fs); % trend
y(1) = y(2);
y = smooth(y, 12*fs, 'loess'); % smooth trend
PVP_dtr = filtfilt(b, a, phantom_PVP-y);% bandpass filter
%% Get SQI index
[SQI] = PVP_SQI(PVP_dtr,fs);





