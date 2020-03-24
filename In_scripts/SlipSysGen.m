% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Function    : SlipSysGen
%   Last edited  : 1 November, 2018 - SW
%   Description  : called by Input.m
%                     Generate slip system angles for each grain
%   Outstanding issues :
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lambda, R] = SlipSysGen(nSystems,ngr,orientation,dlambda)

% Define slip systems of each grain with slip shat and normal nhat
lambda = zeros(max(nSystems),ngr);

% rotation matrix for coordinate transform of each slip system
R = zeros(2,2,max(nSystems),ngr);

for ng = 1:ngr
    for isys = 1:nSystems(ng)
        lambda(isys,ng) = orientation(ng)+dlambda(ng)*(isys-1);
        while lambda(isys,ng)<0; lambda(isys,ng) = lambda(isys,ng)+pi; end
        while lambda(isys,ng)>=pi; lambda(isys,ng) = lambda(isys,ng)-pi; end
        R(1:2,1:2,isys,ng) = [cos(lambda(isys,ng)) sin(lambda(isys,ng));...
                            -sin(lambda(isys,ng)) cos(lambda(isys,ng))]; 
    end
end