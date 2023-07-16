%% Coherent spatial structure and Ritz eigenvalues from DMD .
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
% Define the window size for the moving average filter .
% window_size = 10;
% Define filtering neighbourhood m x n.
m = 5;
n = 5;
% Using median filtering to smoothen data for all time epochs .
for i = 1: n_samples
smooth_data (: ,: , i ) = medfilt2 ( x (: ,: , i ) , [ m n ]) ;
end
% Flattened matrix .
X0 = reshape ( smooth_data ,128*20 ,[]) ;
[r , c ] = size ( X0 ) ;
X1 = X0 (: , 1: c -1) ;
X2 = X0 (: , 2: c ) ;
% X1 = X0 (: , 1:c -100) ;
% X2 = X0 (: , 101: c);
% SVD
[U ,S , V ] = svd ( X1 ,'econ ') ;
% View plot to see which modes are dominant . (10 modes are dominant )
figure (1)
plot ( diag ( S ) / sum ( diag ( S ) ) ,'ko ','Linewidth ' ,1.2)
hold on
plot ( diag ( S (1:200 ,1:200) ) / sum ( diag ( S ) ) ,'ro ','Linewidth ' ,1.2)
hold off
xlabel ('Number of modes ') ;
ylabel ('Energy ') ;
title ('Energy vs Number of modes ') ;
% Compute DMD (Phi are eigenvectors )
r = 200; % Truncate at 200 modes .
U = U (: ,1: r ) ;
S = S (1: r ,1: r ) ;
V = V (: ,1: r ) ;
Atilde = U' * X2 * V * inv ( S ) ; % Best -fit linear model .
[W , eigs ] = eig ( Atilde ) ;
Phi = X2 * V * inv ( S ) * W ; % Reconstruct high -dim full state DMD e- vectors .
alpha1 = S * V (1 ,:) ';
b = ( W * eigs ) \ alpha1 ;

%% Plot the dominant DMD modes .
num_modes = 10; % Choose the number of dominant modes to visualize .
real_dominant_modes = [];
imag_dominant_modes = [];
% Discard negative slip potency rates .
for i = 1: num_modes
    for j = 1:2560
        real_dominant_modes (j , i ) = max ([0 real ( Phi (j , i ) ) ]) ;
        imag_dominant_modes (j , i ) = max ([0 imag ( Phi (j , i ) ) ]) ;
    end
end

for i =1: num_modes
    figure ;
    data_1 = reshape ( real_dominant_modes (: , i ) ,128 ,20) ;
    surf ( data_1 )
    xlabel ('Latitude ') ;
    ylabel ('Longitude ') ;
    zlabel ('Relative Slip Potency Rate ') ;
    title ('Spatial Structure ( Dominant DMD Mode )') ;
    colorbar
    figure ;
    data_2 = reshape ( imag_dominant_modes (: , i ) ,128 ,20) ;
    surf ( data_2 )
    colorbar
end

%% Ritz Values.
figure
theta = (0:1:100) *2* pi /100;
plot ( cos ( theta ) , sin ( theta ) ,'r- -') % plot unit circle
hold on , grid on
scatter ( real ( diag ( eigs ) ) , imag ( diag ( eigs ) ) ,'ok ')
axis ([ -1.1 1.1 -1.1 1.1]) ;
title ('DMD Spectrum ( Eigenvalues / Ritz Values )') ;
xlabel ('Re ') ;
ylabel ('Im ') ;

