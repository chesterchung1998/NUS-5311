%% Future state prediction using Dynamic Mode Decomposition (DMD).
clc ; clear ; close all

%% This code is used to predict the future states using DMD
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

%% Smoothening of original data
X_filt = x ;
X_original = reshape (x ,128*20 ,[]) ;
X_original = X_original (: ,1:56826) ;
X_filt_flat = reshape ( X_filt ,128*20 ,[]) ;
X_filt_flat = X_filt_flat (: ,1:56826) ;
[ row , c ] = size ( X_filt_flat ) ;
X1 = X_filt_flat (: , 1: c -1) ;
X2 = X_filt_flat (: , 2: c ) ;
count = -1;

%% Plot to compare original with smoothened data
figure (1)
for i = [100 ,28000 ,56000]
count = count +2;
subplot (3 ,2 , count )
plot ( X_original (: , i ) ,'g')
hold on
plot ( X1 (: , i ) ,'b')
hold off
xlim ([0 2560])
title ( strcat (" Flatened slip potenccy rate values at epoch =" , num2str ( i ) ) ) ;
legend ('True values ','filtered values ') ;
end

%% Compute DMD
% SVD
[U ,S , V ] = svd ( X1 , 'econ ') ;
r = 1000; % Truncate at 10 modes . (Try 500 modes to see difference )
U = U (: ,1: r ) ;
S = S (1: r ,1: r ) ;
V = V (: ,1: r ) ;
Atilde = U '* X2 * V * inv ( S ) ; % Best -fit linear model .
[W , eigs ] = eig ( Atilde ) ;
Phi = X2 * V * inv ( S ) * W ; % Reconstruct high -dim full state DMD e- vectors .
alpha1 = S * V (1 ,:) ';
b = ( W * eigs ) \ alpha1 ;

% Prediction of slip potency rates for epoch = 2.
alpha2 = Atilde * alpha1 ;
xt2 = U * alpha2 ;

%% Prediction of future state ( epoch = 56826)
alpha = alpha1 ;
for i = 2:56826
alpha_pre (: ,i -1) = S * V (i -1 ,:) ';
alpha (: , i ) = Atilde * alpha_pre (: ,i -1) ;
x_pre (: , i ) = U * alpha (: , i ) ;
end
x_pre = max ( x_pre ,0) ;
X_predected = reshape ( x_pre ,128 ,20 ,56826) ;

%% visualisation
figure (1)
count = 0;
for i = [100 ,28000 ,56000]
count = count +2;
subplot (3 ,2 , count )
plot ( X1 (: , i ) ,'b')
hold on
plot ( x_pre (: , i ) ,'r')
hold off
xlim ([0 2560])
title ( strcat (" Predicted value vs Filtered values at epoch =" , num2str ( i ) ) ) ;
legend ('Filtered values ','Predicted values ') ;
end

% figure (2)
levels = [0 ,0.01 ,0.03 ,0.05 ,0.08 ,20 ,100 ,660];
count = 0;
fig = 2;
for i = [1200 ,23000 ,48480 ,100 ,40000 ,50000]
figure ( fig )
fig = fig +1;
count = count +1;
subplot (1 ,2 , count )
contourf ( X_filt (: ,: , i ) , levels )
title ( strcat (" True Values at epoch =" , num2str ( i ) ) )
count = count +1;
subplot (1 ,2 , count )
contourf ( X_predected (: ,: , i ) , levels )
title ( strcat (" Predicted values at epoch =" , num2str ( i ) ) ) ;
% legend (' Filtered values ',' Predicted values ');
% xlim ([0 2560])
count =0;
hold off
pause (0.1)
end

