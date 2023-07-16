%% Computation of Largest Lyapunov Exponents
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
for i = 1: n_samples
    data = x (1: end , 1: end , i ) ;
    max_spr = max ( data (:) ) ;
    [ ii , jj ] = find ( data == max_spr ) ; % Row and column number .
    spr_mat (i ,1) = ii ;
    spr_mat (i ,2) = jj ;
    spr_mat (i ,3) = max_spr ;
end

% Using MATLAB lyapunovExponent function .
x_data = spr_mat (: ,1) ; % Use value = 1 for x- data and value = 2 for y- data .
[~ , elag , eDim ] = phaseSpaceReconstruction ( x_data ) ; % Calculating lag and dimension .
fs = 1/(24*60*60) ; % Sampling frequency .
eRange = 1000;
lyapunovExponent ( x_data , fs , elag , eDim , ' ExpansionRange ', eRange )