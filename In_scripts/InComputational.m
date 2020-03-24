% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Function     : InComputational
%   Last edited  : 1 November, 2018 - SW
%   Description  : called by Input.m
%                    Input computational parameters
%   Outstanding issues :                 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arraySize = 1000;       % max. number of dislocations
iMax = 10000;           % max. number of increments
nrppd = 1000;           % no. of increment at which save the dislocation structure.
plotdis = 0;            % plot dislocation structure during simulation? 0 - no; 1- yes
pointSize = 12;         % marker size of plot
dis_size = 0.02;        % dislocation plot size 
dt0= 0.5E-9;            % small time step (s)
% f_inc = 0.5E-9*2000/(Udot/((bbox(2,2)-bbox(1,2)))); % time step increase speed in large time step loop
mObsJump = 1;           % use to select obstacle escaping method. 1 - non thermal, 2 - thermally activated.
fname = ['Results'];