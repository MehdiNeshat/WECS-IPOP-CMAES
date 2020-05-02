% 04/06/2018 - chenged to calculate the spectrum using given frequencies


function Spec = spectrum_PMw(Hs, Tp, w)

S = 5*pi^4*Hs^2/Tp^4*1./w.^5.*exp(-20*pi^4/Tp^4*1./w.^4);

% Spectrum structure
Spec.type = 'Pierson-Moskowitz';
Spec.data = S';	% structure that exports the determined spectra
Spec.Tp = Tp;   % export-structure of peak period
Spec.Hs = Hs;	% export-structure of significant wave-height
Spec.freq = w';
Spec.frequnits = 'rad/s';
Spec.densunits = 'm^2 /(rad/s)';
Spec.nfreq = length(w);