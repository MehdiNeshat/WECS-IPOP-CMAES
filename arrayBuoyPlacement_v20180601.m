% Array
% Submerged spheres, 3-tether arrangement, individual PTOs for each tether

% 26/10/2015    q-factor is changed to the RCW
% 07/01/2016    added calculation of power for each buoy
% 16/02/2016    added parameters - water depth, submergence depth, wave angle, PTO - use precalculated or calculate inside,
%               tether angle is changing according to the submergence depth
% 06/04/2016    added different PTO settings for each frequency
%               added q-factor
%               TetherForceWeighted instead of the Power/TetherForce
%               PTO parameters are included in the lookUpTable
% 19/04/2016    added possibility of one frequency
% 11/07/2016    added the range of wave angles
%               location is a subject of optimisation (NOT staggered array)
% 22/03/2017    corrected a bug associated with numFreqs ~= 50
% 08/10/2017    changed calculation routine for the total power (similar to spectral-domain in irregular waves)
% 01/06/2018    changed 'lookUpTable' to 'buoy' structure (PTO parameters are fixed for all frequencies)
%               changed calculation routine for the total power - now it
%               calculates the annual average power based on the
%               methodology proposed by Babarit

function [Parray_AAP, ParrayBuoy_AAP, q_AAP] = arrayBuoyPlacement_v20180601(array, siteOpts, buoy)
%% Input
% array:
%     radius            - radius of each sphere
%     number            - number of bodies within array
%     sphereCoordinates	- coordinates of all bodies (x,y,z), z is up from the water level
% siteOpts:
%     subDepthT         - submergence depth, distance from the top of the body to the water surface
%     oceanDepth        - water depth
%     siteOpts.location:
%           siteName
%           Hs          - significant wave heights
%           Tp          - peak wave periods
%           beta        - wave directions
% buoy:
%     mass
%     volume
%     tetherAngle
%     kPTO
%     dPTO

%% Output
% Parray_AAP              - annual average power absorbed by the array
% ParrayBuoy_AAP          - annual average power absorbed by the array

%%
subDepthT       = siteOpts.submergenceDepth;
oceanDepth      = siteOpts.waterDepth;
numSphere       = array.number;

% Water parameters
rho = 1025;             % kg/m^3, water density
g = 9.80665;            % m/s^2
numApprox = 4;          % number of equations

% Range of frequencies used to calculate array performance
w_array = siteOpts.waveFreqs;
numFreqs = length(w_array);

% Wave structure
wave.waterDensity	= rho;
wave.angle          = siteOpts.location.waveAngle;
numWaveAngle        = length(wave.angle);

Hs_hist  = siteOpts.location.Hs;
Tp_hist  = siteOpts.location.Tp;
N_hist   = siteOpts.location.N_hist;

% PTO specification
I3 = eye(3);

% Parameters for each sphere according to its radius
CtArray     = zeros(numSphere*3);
KtArray     = zeros(numSphere*3);
J           = zeros(numSphere*3);
M           = zeros(numSphere*3);               % Array mass matrix
dPTOarray   = zeros(numSphere*3, 1);
kPTOarray   = zeros(numSphere*3, 1);

for ii = 1:numSphere

    idx = 3*(ii-1)+(1:1:3);     % for 3 modes

    % Settings for the sphere of this radius
    a = array.radius(ii);
       
    %% Detect the tether angle for this body
    alpha	= buoy.tetherAngle;
    V       = buoy.volume;
    mass	= buoy.mass;

    % Net tether force
    Ft = (rho*g*V - mass*g);

    % Initial tether length
    s0n = (oceanDepth - subDepthT - a - a*cos(alpha))/cos(alpha);

    % Initial tention in a leg
    t0 = Ft/(3*cos(alpha));
    gam0 = t0/s0n;
    
    % Mass matrix
    M(idx, idx) = mass*I3;

    for count_waveAngle = 1:numWaveAngle
    
        beta = wave.angle(count_waveAngle);
        
        % Identification of vectors according to Scruggs (2013)
        % Unit vectors along tethers (from ocean bottom to the attachment point)
        es01 = [-sin(alpha)*cos(beta);         sin(alpha)*sin(beta);        cos(alpha)];
        es02 = [ sin(alpha)*sin(pi/6+beta);    sin(alpha)*cos(pi/6+beta);	cos(alpha)];
        es03 = [ sin(alpha)*sin(pi/6-beta);   -sin(alpha)*cos(pi/6-beta);	cos(alpha)];

        % Kt and Ct from Scruggs (2013)
        % Eq.11
        Gt1 = (-I3*es01);
        Gt2 = (-I3*es02);
        Gt3 = (-I3*es03);

        J(idx, idx, count_waveAngle) = [es01'; es02'; es03'];          % 07/01/2016

        for count_w = 1:length(w_array) % 06/04/2016

            kPTO = buoy.kPTO;  % 04/06/2018
            dPTO = buoy.dPTO;  % 04/06/2018

            % Eq.13
            Ct1 = dPTO*(Gt1*Gt1');
            Ct2 = dPTO*(Gt2*Gt2');
            Ct3 = dPTO*(Gt3*Gt3');

            % Eq.14
            Kt1 = ((kPTO-gam0)*(Gt1*Gt1') + gam0*I3);
            Kt2 = ((kPTO-gam0)*(Gt2*Gt2') + gam0*I3);
            Kt3 = ((kPTO-gam0)*(Gt3*Gt3') + gam0*I3);

            Ct = Ct1 + Ct2 + Ct3;
            Kt = Kt1 + Kt2 + Kt3;

            % PTO impedance for the array
            CtArray(idx, idx, count_w, count_waveAngle) = Ct;
            KtArray(idx, idx, count_w, count_waveAngle) = Kt;
            dPTOarray(idx, count_w) = dPTO*ones(3,1);           % 07/01/2016
            kPTOarray(idx, count_w) = kPTO*ones(3,1);           % 22/02/2016

        end
    end
end

Parray          = zeros(numFreqs, numWaveAngle);
ParrayBuoy      = zeros(numSphere, numFreqs, numWaveAngle);
Parray_hist     = zeros(length(Hs_hist), length(Tp_hist), numWaveAngle);
ParrayBuoy_hist = zeros(numSphere, length(Hs_hist), length(Tp_hist), numWaveAngle);
ParrayBuoy_AAP  = zeros(numSphere, 1);
Pisolated_hist  = zeros(length(Hs_hist), length(Tp_hist), numWaveAngle);

% parpool();

%% Calculate power function over the range of wave angles P(w, beta)
parfor count_w = 1:numFreqs
    
    w = w_array(count_w);
    K = w^2/g;

    % Wave power (infinite water depth)
    Pw = rho*g^2/(4*w);
    
    % Array hydrodynamics
    [A, B, Xw] = arraySubmergedSphereParfor(array, wave, w, K, numApprox, 1);
    
    for count_waveAngle = 1:numWaveAngle
        
        X = Xw(:, count_waveAngle);
        
        % PTO impedance
        Zpto = (CtArray(:,:,count_w, count_waveAngle) - 1i*KtArray(:,:,count_w, count_waveAngle)/w);

        % Buoy impedance
        Zbuoy = 1i*(M + A)*w + B;

        % Velocity vector
        Uarray = (Zbuoy + Zpto)^(-1)*X;

        % Total absorbed power of the array
        Parray(count_w, count_waveAngle) = real(1/4*(Uarray'*X + X'*Uarray) - 1/2*Uarray'*B*Uarray);

        %% To calculate power for each buoy individually (07/01/2016)
        ddelta_s = J(:,:,count_waveAngle)*Uarray;
        
        PowerTether = reshape(abs(ddelta_s).^2.*dPTOarray(:, count_w)/2, [3 numSphere]);
        
        ParrayBuoy(:, count_w, count_waveAngle) = sum(PowerTether, 1)';

    end
end

% poolobj = gcp('nocreate');
% delete(poolobj);

%% Calculate power matrix
for count_Hs = 1:length(Hs_hist)
    
    Hs = Hs_hist(count_Hs);
    
    for count_Tp = 1:length(Tp_hist)
        
        Tp = Tp_hist(count_Tp);
        S = spectrum_PMw(Hs, Tp, w_array);
        
        for count_waveAngle = 1:numWaveAngle
            
            Parray_hist(count_Hs, count_Tp, count_waveAngle)        = trapz(w_array, 2*S.data.*Parray(:, count_waveAngle));
            Pisolated_hist(count_Hs, count_Tp, count_waveAngle)     = trapz(w_array, 2*S.data.*buoy.PowerIsolated);
            
            for count_buoy = 1:numSphere
            
                ParrayBuoy_hist(count_buoy, count_Hs, count_Tp, count_waveAngle)	= trapz(w_array, 2*S.data.*ParrayBuoy(count_buoy, :, count_waveAngle)');
                
            end
            
        end
    end
end

%% Calculate averaged annual power output
P               = Parray_hist.*N_hist;
Parray_AAP      = sum(P(:))./(sum(N_hist(:)));

Pi              = Pisolated_hist.*N_hist;
Pisolated_AAP   = sum(Pi(:))./(sum(N_hist(:)));

for count_buoy = 1:numSphere
    P = squeeze(ParrayBuoy_hist(count_buoy, :, :, :)).*N_hist;
    ParrayBuoy_AAP(count_buoy)	= sum(P(:))./(sum(N_hist(:)));
end

q_AAP = Parray_AAP/(numSphere*Pisolated_AAP);