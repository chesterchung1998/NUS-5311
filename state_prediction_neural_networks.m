%% Designing Neural Network for future state prediction .
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
%% Smoothening the data .

% Define filtering neighbourhood m x n.
m = 5;
n = 5;

% Using median filtering to smoothen data for all time epochs .
for i = 1: n_samples
    smooth_data (: ,: , i ) = medfilt2 ( x (: ,: , i ) , [ m n ]) ;
end
future_time = 5000;
for k = 1:50: future_time
    X0 = reshape ( smooth_data ,128*20 ,[]) ;
    [r , c ] = size ( X0 ) ;
    X1 = X0 (: , 1: c - k ) ;
    X2 = X0 (: , 1+ k : c ) ;
    % SVD
    [U ,S , V ] = svd ( X1 ,'econ ') ;
    % Compute DMD (Phi are eigenvectors )
    r = 50; % Truncate at 10 modes . (Try 500 modes to see difference )
    U1 = U1 (: ,1: r ) ;
    S1 = S1 (1: r ,1: r ) ;
    V1 = V1 (: ,1: r ) ;
    U2 = U2 (: ,1: r ) ;
    S2 = S2 (1: r ,1: r ) ;
    V2 = V2 (: ,1: r ) ;
    alpha1 = S1 * V1 (: ,:) ';
    alpha2 = S2 * V2 (: ,:) ';
    x_norm = normalize ( alpha1 ) ;
    % Neural network configuration and training
    net = feedforwardnet ([15 15 15 15]) ;
    net . layers {1}. transferFcn = 'tansig ';
    net . layers {2}. transferFcn = 'radbas ';
    net . layers {3}. transferFcn = 'radbas ';
    net . layers {4}. transferFcn = 'tansig ';
    net . layers {5}. transferFcn = 'purelin ';
    net . divideFcn = 'dividerand ';
    net . divideMode = 'sample ';
    net . divideParam . trainRatio = 0.8;
    net . divideParam . valRatio = 0.1;
    net . divideParam . testRatio = 0.1;
    net . trainFcn = 'trainlm ';
    net . trainParam . max_epochs = 1000;
    net . trainParam . goal = 0.01;
    net = train ( net , x_norm , alpha2 ) ;

    %% Removing negative slip potency rates .
    for i = 1:2560
        xt2 (i ,1) = max ([0 xt2 (i ,1) ]) ;
    end

    %% Plotting contour plots .
    % Predicted states (For epoch = last + 1000) .
    contourf ( reshape ( xt2 (1:2560 ,1) ,128 ,20) )
    colorbar
    title ( strcat (" Predicted state ( Epoch " , num2str ( n_samples + k ) , ") ") ) ;
end
