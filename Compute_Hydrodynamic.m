function [siteOpts,buoy,siteName]=Compute_Hydrodynamic(Opt)

%% WAVE CLIMATE
% read netcdf file into matlab to get the probablity diagram of sea states

%%%%%%%%%%%%%
% Information on the netcdf file, put ncdisp(fname)
%%%%%%%%%%%%%

location = Opt.WaveModel;

switch location
    
    case 1  % Perth
        fname = 'HCtimeseries_aus_4m_114.80E_33.50S.nc';
        siteName = 'Perth';
        disp('Perth wave model')
    case 2  % Adelaide
        fname = 'HCtimeseries_aus_4m_136.62E_36.07S.nc';
        siteName = 'Adelaide';
        disp('Adelaide wave model')
    case 3  % Tasmania
        fname = 'Tasmania_145.0372E_42.8175S.nc';
        siteName = 'Tasmania';
        disp('Tasmania wave model')
    case 4  % Sydney
        fname = 'HCtimeseries_aus_4m_152.50E_34.00S.nc';
        siteName = 'Sydney';
        disp('Sydney wave model')
end


time    = ncread(fname,'time');
hs      = ncread(fname,'hs');
fp      = ncread(fname,'fp');
tm0m1   = ncread(fname,'tm0m1');
dirWave = ncread(fname,'dir');

n1 = length(find(time==-32767));
n2 = length(find(hs==-32767));
n3 = length(find(fp==-32767));
n4 = length(find(tm0m1==-32767));
n5 = length(find(dirWave==-32767));

nMax = max([n1, n2, n3, n4, n5]);

time(1:nMax)    = [];
hs(1:nMax)      = [];
fp(1:nMax)      = [];
tm0m1(1:nMax)   = [];
dirWave(1:nMax) = [];

dirWave(dirWave == 360) = 0;

%% 4D Histogram (wave height, wave period, wave direction)
xbins_hs = 0.5:1:[ceil(max(hs))+0.5];
xbins_tp = round(1/max(fp)):1:round(1/min(fp));
xbins_wa = 0:15:360;

[counts, edges, mid, loc] = histcn([hs 1./fp dirWave], xbins_hs, xbins_tp, xbins_wa);

%% Detect significant wave directions
N_hist_wa = squeeze(sum(sum(counts, 1), 2));

[N_hist_wa_sort, I_wa_sort] = sort(N_hist_wa, 'descend');

n_wa = 0;
ii = 0;
while n_wa < 0.90
    ii = ii + 1;
    n_wa = n_wa + N_hist_wa_sort(ii)/sum(N_hist_wa_sort);
end

I_wa = sort(I_wa_sort(1:ii));

%% Detect significant wave heights
N_hist_hs = squeeze(sum(sum(counts, 3), 2));

[N_hist_hs_sort, I_hs_sort] = sort(N_hist_hs, 'descend');

n_hs = 0;
ii = 0;
while n_hs < 0.90
    ii = ii + 1;
    n_hs = n_hs + N_hist_hs_sort(ii)/sum(N_hist_hs_sort);
end

I_hs = sort(I_hs_sort(1:ii));

%% Detect significant peak wave periods
N_hist_tp = squeeze(sum(sum(counts, 1), 3));

[N_hist_tp_sort, I_tp_sort] = sort(N_hist_tp, 'descend');

n_tp = 0;
ii = 0;
while n_tp < 0.90
    ii = ii + 1;
    n_tp = n_tp + N_hist_tp_sort(ii)/sum(N_hist_tp_sort);
end

I_tp = sort(I_tp_sort(1:ii));

%% Data reduction
N_hist          = counts(I_hs, I_tp, I_wa);
Hs_hist         = mid{1}(I_hs);
Tp_hist         = mid{2}(I_tp);
waveAngle_hist	= mid{3}(I_wa);

%% Define the parameters you want to run the model
siteOpts.waterDepth         = 30;                % Water depth
siteOpts.submergenceDepth   = 3;                 % Submergence depth of the buoy (from the water level to the top of the buoy)
siteOpts.waveFreqs          = linspace(0.2, 2, 50);
siteOpts.location.siteName  = siteName;
siteOpts.location.Hs        = Hs_hist;
siteOpts.location.Tp        = Tp_hist;
siteOpts.location.waveAngle = waveAngle_hist*pi/180;
siteOpts.location.N_hist    = N_hist;
%%%%%%%%
Num_Buoy                  = Opt.Buoy_Num;
array.radius = 5*ones(1,Num_Buoy );
% array.number = length(array.radius);
%
% array.sphereCoordinate(1,:) = x - min(x);
% array.sphereCoordinate(2,:) = y - min(y);
% array.sphereCoordinate(3,:) = -(siteOpts.submergenceDepth + array.radius);

% Parameters of the indivudual buoy
buoy.mass           = 3.7568e5;     % kg
buoy.volume         = 523.5988;     % m^3
buoy.tetherAngle	= 0.9553;       % rad
buoy.kPTO           = 407510;       % N/m
buoy.dPTO           = 97412;        % N/(m/s)
%% Calculate power for one buoy with given parameters
rho     = 1025;
g       = 9.80665;

wave.waterDensity	= rho;
wave.angle          = 0;
numApprox           = 4;

mass	= buoy.mass;
V       = buoy.volume;
Ft      = (rho*g*V - mass*g);

PowerIsolated = zeros(length(siteOpts.waveFreqs), 1);

for count_w = 1:length(siteOpts.waveFreqs)
    
    w = siteOpts.waveFreqs(count_w);
    
    K = w^2/g;
    
    % Structure of the isolated sphere
    sphere.radius           = array.radius(1);
    sphere.massMatrix       = eye(3)*mass;
    sphere.tetherPretention = Ft;
    sphere.number           = 1;
    sphere.sphereCoordinate = [0; 0; -(siteOpts.submergenceDepth + array.radius(1))];
    sphere.oceanDepth       = siteOpts.waterDepth;
    sphere.submergenceDepth = siteOpts.submergenceDepth + array.radius(1);
    sphere.waveAngle        = 0;
    sphere.tetherAngle      = buoy.tetherAngle;
    
    [A, Bd, X] = arraySubmergedSphereParfor(sphere, wave, w, K, numApprox, 1);
    
    sphere.addedMass = A;
    sphere.damping = Bd;
    sphere.excitationForce = X;
    sphere.waveFrequency = w;
    sphere.waveNumber = K;
    
    [PowerIsolated(count_w), ~, ~, ~] = powerFunSphere([buoy.kPTO buoy.dPTO], sphere);
    
end

buoy.PowerIsolated = PowerIsolated;
end