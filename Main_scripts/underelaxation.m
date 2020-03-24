% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Function    : Determine Underrelaxation Factor c for each dislocation
%   Description : call by Main.m
%                 If the velocity has changed direction use 0.5 for c
%                 If dislocation has same type as either of neighboring 
%                 dislocations, using an additional factor of 0.75 for each
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vdis] = underelaxation(ndis,activeplanes,vdis,vdispre,type,plane)
c = ones(ndis,1);
for p = 1:length(activeplanes)
    
    m = find(plane(1:ndis)==activeplanes(p));
    ndisp = length(m);
    for i = 1:ndisp
        idis = m(i);
        if vdis(idis)*vdispre(idis)<0 
            c(idis) = 0.5*c(idis);
            if i>1
                jdis = m(i-1);
                if (type(jdis)*type(idis))>0
                    c(idis) = 0.75*c(idis);
                end
            end
            if i<ndisp
                jdis = m(i+1);
                if (type(jdis)*type(idis))>0
                    c(idis) = 0.75*c(idis);
                end
            end
        end
    end
end
vdis(1:ndis) = vdis(1:ndis).*c(1:ndis);