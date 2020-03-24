% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Function    : InConstDefine
%   Last edited  : 1 November, 2018 - SW
%   Description  : called by Input.m
%                   Define constants before the main loop
%   Outstanding issues : 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

edgeType = mu./(2*pi*(1-nu));
uConst = 1./(2*pi*(1-nu));
ffric = friction.*tau_FR.*bG;

T = 273; % Temperature in K

inf = 1.0E12;
tol = 1E-8; % tolerance
damping = 0.5;
reps = 2*min(bG);
Le = 6*min(bG); % annihilation distance

% Thermal activation parameters
xfreq = 1E11; % Derby frequency (Hz)
xboltzmann = 1.381E-23; % Boltzman constant (J/K)
l_obs = 1/sqrt(0.01); % 3/(100*b*obsDensity);

xalpha = xfreq/l_obs*bG.*exp(-deltaF./(xboltzmann*T)); % constant
xbeta = (l_obs*1E-6*bG.*bG.*gammazero)/(xboltzmann*T); % tau should use TPa