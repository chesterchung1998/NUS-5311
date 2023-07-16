%% Derivation of FFT plots and Gabor Transform Spectrogram .
clc ; clear ; close all
% dimensions of data
n_samples = 56825;
n_latitudes = 20;
n_longitudes = 128;
% load data
ncfile = ' DalZilioetal2020sim .nc ';
% visualize dataset content
ncinfo ( ncfile )
ncdisp ( ncfile )
% get SPR values
x = ncread ( ncfile , 'SPR ') ;
% get time snapshots
t = ncread ( ncfile , 'time ') ;
% get longitude values
lon = ncread ( ncfile , 'lon ') ;
% get latitude values
lat = ncread ( ncfile , 'lat ') ;
% Get the longitudes and latitudes of the maximum slip potency rates .
spr_mat = []; % Longitude , Latitude , Maximum SPR.
for i = 1: n_samples +1
data = x (1: end , 1: end , i ) ;
max_spr = max ( data (:) ) ;
[ ii , jj ] = find ( data == max_spr ) ; % Row and column number .
spr_mat (i ,1) = ii ;
spr_mat (i ,2) = jj ;
spr_mat (i ,3) = max_spr ;
end
spr_max = spr_mat (: ,3) ;
f_noisy = spr_max ;
dt = 1;

%% Compute the Fast Fourier Transform FFT
% Signal length .
n = length ( t ) ;
% Fast Fourier Transform .
fhat = fft ( f_noisy , n ) ;
% Power spectrum ( power per each frequency ).
PSD = fhat .* conj ( fhat ) / n ;
% Create x- axis of frequencies in Hz.
freq = 1 / ( dt * n ) * (0: n ) ;
% Plot only the first half of frequencies .
L = 1: floor ( n /2) ;

%% Use the PSD to filter out noise within signal .
% Find all frequencies with large power above a certain threshold .
indices = PSD >125;
% Zero out all other frequencies .
PSD_clean = PSD .* indices ;
% Zero out small Fourier coefficients .
fhat = indices .* fhat ;
% Apply Inverse FFT for filtered time signal .
ffilt = ifft ( fhat ) ;
% Adjustments ( Removing the sharp peak at 0Hz for FFT plot ).
freq_adj = freq (1 ,2: end) ;
PSD_adj = PSD (2: end ,1) ;
PSD_clean_adj = PSD_clean (2: end ,1) ;

%% Plots for raw & clean signals , and FFT for raw & clean signals .
figure ;
plot (t , f_noisy ,'r','LineWidth ' ,1) , hold on
plot (t , ffilt ,'b','LineWidth ' ,1.5)
ylim ([ -0.25 1])
xlabel ('Epoch (day)')
ylabel ('Slip Potency Rate ')
legend ('Raw ','Filtered ', 'Location ','southeast ')
title ('Plots for raw and filtered signals ')
figure ;
plot ( freq_adj ( L ) , PSD_adj ( L ) ,'r','LineWidth ' ,1.5) , hold on
plot ( freq_adj ( L ) , PSD_clean_adj ( L ) ,'-b','LineWidth ' ,1.2)
xlim ([0 0.005])
xlabel ('Frequency (1/ day)')
ylabel ('Power Spectral Density ')
legend ('Raw ','Filtered ')
title ('FFT plots for raw and filtered signals')

%% Spectrogram from Gabor Transform .
figure ;
% window = 365 (1 year ), noverlap = 292 (80% of window ), nfft = 365.
% fs = 1/(24 * 3600) ( sampling rate ).
spectrogram ( f_noisy ,365 ,292 ,365 ,1/(24*3600) ,'yaxis ') ;
title ('Spectrogram of raw signal')