% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Function    : remove dislocations due to annihilation or escape
%   Description : call by disescape.m and Main.m
%                 flag = 0: annihilation
%                 flag = 1: Left surface
%                 flag = 2: Top surface
%                 flag = 3: Right surface
%                 flag = 4: Bottom surface
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,plane,pinned,irmbound,eta_obs,...
    vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,ngsourceOut,planeOut,nout] = ...
    removeDis(flag,disindex,b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,plane,...
    pinned,irmbound,eta_obs,vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,ngsourceOut,planeOut,nout)

inf = 1.0E12;
% flag==0 %remove dislocations from the list due to annihilation
if flag==2 % flag = 2, Top Crystal Boundary
    bOut = [bOut; b(disindex==1)];
    xdisOut = [xdisOut; xdis(disindex==1)+inf*cos(alpha(disindex==1))];
    ydisOut = [ydisOut; ydis(disindex==1)+inf*sin(alpha(disindex==1))];
    typeOut = [typeOut; type(disindex==1)];
    alphaOut = [alphaOut; alpha(disindex==1)];
    ngsourceOut = [ngsourceOut; ngsource(disindex==1)];
    planeOut = [planeOut; plane(disindex==1)];
    nout = nout + sum(disindex);
elseif flag==3 % flag = 3, Right Crystal Boundary
    bOut = [bOut; b(disindex==1)];
    xdisOut = [xdisOut; xdis(disindex==1)+inf*cos(alpha(disindex==1)).*sign(cos(alpha(disindex==1)))];
    ydisOut = [ydisOut; ydis(disindex==1)+inf*sin(alpha(disindex==1)).*sign(cos(alpha(disindex==1)))];
    typeOut = [typeOut; type(disindex==1)];
    alphaOut = [alphaOut; alpha(disindex==1)];
    ngsourceOut = [ngsourceOut; ngsource(disindex==1)];    
    planeOut = [planeOut; plane(disindex==1)];
    nout = nout + sum(disindex);
elseif flag==4 % flag = 4, Bottom Crystal Boundary
    bOut = [bOut; b(disindex==1)];
    xdisOut = [xdisOut; xdis(disindex==1)-inf*cos(alpha(disindex==1))];
    ydisOut = [ydisOut; ydis(disindex==1)-inf*sin(alpha(disindex==1))];
    typeOut = [typeOut; type(disindex==1)];
    alphaOut = [alphaOut; alpha(disindex==1)];
    ngsourceOut = [ngsourceOut; ngsource(disindex==1)];
    planeOut = [planeOut; plane(disindex==1)];
    nout = nout + sum(disindex);
elseif flag==1 % flag = 1, Left Crystal Boundary
    bOut = [bOut; b(disindex==1)];
    xdisOut = [xdisOut; xdis(disindex==1)-inf*cos(alpha(disindex==1)).*sign(cos(alpha(disindex==1)))];
    ydisOut = [ydisOut; ydis(disindex==1)-inf*sin(alpha(disindex==1)).*sign(cos(alpha(disindex==1)))];
    typeOut = [typeOut; type(disindex==1)];
    alphaOut = [alphaOut; alpha(disindex==1)];
    ngsourceOut = [ngsourceOut; ngsource(disindex==1)];
    planeOut = [planeOut; plane(disindex==1)];
    nout = nout + sum(disindex);
end

ndis = ndis-sum(disindex);

z = zeros(sum(disindex),1);
% update dislocation data structure
b(disindex==1) = [];            b = [b; z];
xdis(disindex==1) = [];         xdis = [xdis; z];
ydis(disindex==1) = [];         ydis = [ydis; z];
type(disindex==1) = [];         type = [type; z];
alpha(disindex==1) = [];        alpha = [alpha; z];
ngsource(disindex==1) = [];     ngsource = [ngsource; z];
source(disindex==1) = [];       source = [source; z];
rdis(disindex==1) = [];         rdis = [rdis; z];
plane(disindex==1) = [];        plane = [plane; z];
pinned(disindex==1) = [];       pinned = [pinned; z];
irmbound(disindex==1) = [];     irmbound = [irmbound; z];
vdispre(disindex==1) = [];      vdispre = [vdispre; z];
eta_obs(disindex==1) = [];      eta_obs = [eta_obs; z];