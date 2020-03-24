% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Function    : Deal with Grain Boundaries
%   Description : call by Main.m
%                 if dislocation has moved beyond a grain boundary, keep
%                 it there and set irmbound = 1 (rdis<_p_rend(1)) or 
%                 irmbound = 2 (rdis>_p_rend(2))
%   
%   GB energy has to be input for given problem
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [irmbound,vdis] = dealGB(irmbound,ndis,p_rend,vdis,rdis,dt,plane,reps)
for i = 1:ndis
    if irmbound(i)==0
        rtrial = rdis(i) + vdis(i)*dt;
        if rtrial<=p_rend(plane(i),1)
            if any(irmbound(plane==plane(i))==1)
                vdis(i) = (p_rend(plane(i),1)+reps-rdis(i))/dt;
                if vdis(i) == 0
                    vdis(i) = -reps/dt;
                end
            else
                irmbound(i) = 1;
                vdis(i) = (p_rend(plane(i),1)-rdis(i))/dt;
            end
        elseif rtrial>=p_rend(plane(i),2)
            if any(irmbound(plane==plane(i))==2)
                vdis(i) = (p_rend(plane(i),2)-reps-rdis(i))/dt;
                if vdis(i) == 0
                    vdis(i) = reps/dt;
                end
            else
                irmbound(i) = 2;
                vdis(i) = (p_rend(plane(i),2)-rdis(i))/dt;
            end
        end
    end
end
