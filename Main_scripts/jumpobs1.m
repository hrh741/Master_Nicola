% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Function    : dislocation jump trough obstacle method 1
%   Description : call by Main.m
%                 Non-thermal activated process
%                 If dislocation is already pinned(pinned>0), set velocity
%                 to zeros;
%                 if resolved shear stress exceed obstacle strength, pass
%                 it and reset pinned = 0
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vdis,pinned] = jumpobs1(activeplanes,plane,ngsource,planeObs,irmbound,...
    vdis,dt,rdis,pinned,taudis,tau_obs,rObs,ndis)

for p = 1:length(activeplanes)
    [m,ngp] = find(plane(1:ndis)==activeplanes(p));
    ndisp = length(m);
    ng = unique(ngsource(m));
    [mO,ngp] = find(planeObs==activeplanes(p));
    nObsp = length(mO);
    for i = 1:ndisp
        % if dislocation was pinned at GB, do not update position, otherwise:
        if irmbound(m(i))==0
            dxi = vdis(m(i))*dt;
            rd = rdis(m(i));
            if pinned(m(i))>0 && abs(taudis(m(i)))<tau_obs(ngsource(m(i)))
                vdis(m(i)) = 0;
            elseif pinned(m(i))>0 && abs(taudis(m(i)))>=tau_obs(ngsource(m(i)))
                % if resolved shear stress exceed obstacle strength-release
                jold = pinned(m(i));
                pinned(m(i)) = 0;
                for j = 1:nObsp
                    if j~=jold
                        dist = rObs(mO(j),ng)-rd;
                        % if dislocation was moved beyond an obstacle, pin it
                        if dxi/dist>=1
                            vdis(m(i)) = dist/dt;
                            pinned(m(i)) = j;
                        end
                    end
                end
            else
                for j = 1:nObsp
                    dist = rObs(mO(j),ng)-rd;
                    % if dislocation was moved beyond an obstacle, pin it
                    if dxi/dist>=1
                        vdis(m(i)) = dist/dt;
                        pinned(m(i)) = j;
                    end
                end
            end
        end
    end
end