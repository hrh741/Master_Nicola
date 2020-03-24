% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Function    : dislocation jump trough obstacle method 2
%   Description : call by Main.m
%                 Thermal activated process
%                 If dislocation is already pinned(pinned>0), set velocity
%                 to zeros;
%                 if the time being pinned exceed a critical value, pass
%                 it and reset pinned = 0
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vdis,pinned,eta_obs] = jumpobs2(activeplanes,plane,ngsource,planeObs,irmbound,...
    vdis,dt,rdis,pinned,taudis,rObs,xalpha,xbeta,eta_obs,ndis,tau_FR,reps)

for p = 1:length(activeplanes) % each slip plane
    [m,ngp] = find(plane(1:ndis)==activeplanes(p));
    ndisp = length(m);
    ng = unique(ngsource(m));
    [mO,ngp] = find(planeObs==activeplanes(p));
    nObsp = length(mO);
    for i = 1:ndisp % each dislocation
        % if dislocation was pinned at GB, do not update position, otherwise:
        if irmbound(m(i))==0
            dxi = vdis(m(i))*dt;
            rd = rdis(m(i));
            if pinned(m(i))>0 && abs(taudis(m(i)))>tau_FR(ng)
                if eta_obs(m(i))<1
                    vdis(m(i)) = 0;
                    xGAMMA = xalpha(ng)*sinh((abs(taudis(m(i))))*xbeta(ng));
                    tobs = 1/xGAMMA;
                    eta_obs(m(i)) = eta_obs(m(i))+dt/tobs;
                    if eta_obs(m(i))>1; eta_obs(m(i))=1; end
                else % release from current obstacle and check the other obstacles
                    jold = pinned(m(i));
                    pinned(m(i)) = 0;
                    eta_obs(m(i)) = 0;
                    for j = 1:nObsp
                        if j~=jold
                            dist = rObs(mO(j),ng)-rd;
                            % if dislocation was moved beyond an obstacle, pin it
                            if dxi/dist>=1
                                vdis(m(i)) = dist/dt;
                                dxi = vdis(m(i))*dt;
                                pinned(m(i)) = j;
                            end
                        end
                    end
                    if pinned(m(i))>0 && sum(pinned(m)==pinned(m(i)))>1
                        vdis(m(i)) = vdis(m(i))-sign(vdis(m(i)))*reps/dt;
                        pinned(m(i)) = 0;
                    end
                end
            elseif pinned(m(i))>0 && abs(taudis(m(i)))<=tau_FR(ng)
                vdis(m(i)) = 0;
                eta_obs(m(i)) = 0;
            else
                eta_obs(m(i)) = 0;
                for j = 1:nObsp
                    dist = rObs(mO(j),ng)-rd;
                    % if dislocation was moved beyond an obstacle, pin it
                    if dxi/dist>=1
                        vdis(m(i)) = dist/dt;
                        dxi = vdis(m(i))*dt;
                        pinned(m(i)) = j;
                    end
                end
                if pinned(m(i))>0 && sum(pinned(m)==pinned(m(i)))>1
                    vdis(m(i)) = vdis(m(i))-sign(vdis(m(i)))*reps/dt;
                    pinned(m(i)) = 0;
                end
            end
        end
    end
end