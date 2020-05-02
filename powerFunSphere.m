function [Power, c, ceq, U] = powerFunSphere(x, sphere)

% 06/04/2016 - only translational motion for ComSci

kt = x(1);
ct = x(2);

% at - angle of tethers to the vertical
% kt - stiffness of the PTO
% ct - damping of the PTO

% Parameters of the system
a = sphere.radius;
A = sphere.addedMass;
B = sphere.damping;
X = sphere.excitationForce;
w = sphere.waveFrequency;
k = sphere.waveNumber;
f = sphere.submergenceDepth;
M = sphere.massMatrix;
h = sphere.oceanDepth;
Ft = sphere.tetherPretention;
at = sphere.tetherAngle;
try
    beta = sphere.waveAngle;
catch
    beta = 0;
end

try
    Fdrift = sphere.driftForce;
catch
    Fdrift = 0;
end

if beta ~= 0 && Fdrift ~= 0
    error('Change algorithm to incorporate beta and drift at the same time!!!')
end

I3 = eye(3);

% Identification of vectors according to Scruggs (2013)
% Unit vectors along tethers (from ocean bottom to the attachment point)
es01 = [-sin(at)*cos(beta);         sin(at)*sin(beta);      cos(at)];
es02 = [ sin(at)*sin(pi/6+beta);    sin(at)*cos(pi/6+beta);	cos(at)];
es03 = [ sin(at)*sin(pi/6-beta);   -sin(at)*cos(pi/6-beta);	cos(at)];
       
s0n = (h - f - a*cos(at))/cos(at);

% Initial tention in a leg - function of an inclination angle
t0(1,1) = Ft/(3*cos(at)) - 2*Fdrift/(3*sin(at));
t0(2,1) = Ft/(3*cos(at)) +   Fdrift/(3*sin(at));
t0(3,1) = t0(2,1);

gam0 = t0/s0n;

% Buoy impedance matrix
Zb = 1i*(M+A)*w + B;

% Kt and Ct from Scruggs (2013)
% Eq.11
Gt1 = (-I3*es01);
Gt2 = (-I3*es02);
Gt3 = (-I3*es03);

% Eq.13
Ct1 = ct*(Gt1*Gt1');
Ct2 = ct*(Gt2*Gt2');
Ct3 = ct*(Gt3*Gt3');

% Eq.14
Kt1 = ((kt-gam0(1))*(Gt1*Gt1') + gam0(1)*I3);
Kt2 = ((kt-gam0(2))*(Gt2*Gt2') + gam0(2)*I3);
Kt3 = ((kt-gam0(3))*(Gt3*Gt3') + gam0(3)*I3);

Ct = Ct1 + Ct2 + Ct3;
Kt = Kt1 + Kt2 + Kt3;

% PTO impedance
Zu = Ct - 1i*Kt/w; 

% Total impedance of the system (PTO + buoy)
Zs = Zb + Zu;

U = Zs^(-1)*X;
S = U/(1i*w);

delta_s = [es01'; es02'; es03']*S;

g1 = 3;

Power = real(1/4*(U'*X + X'*U) - 1/2*U'*B*U);
c = abs(delta_s) - g1;
ceq = [];
