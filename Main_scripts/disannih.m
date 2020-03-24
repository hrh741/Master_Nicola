
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Function    : Dislocation annihilation
%   Description : call by Main.m
%                 flag dislocations if oppsite sign on the same slip plane
%                 within critical distance Le
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [destroy] = disannih(plane,ndis,rdis,type,Le)
destroy = zeros(ndis,1);
activeplanes = unique(plane(1:ndis));
activeplanes(activeplanes<0) = [];
for p = 1:length(activeplanes)
    m = find(plane(1:ndis)==activeplanes(p));
    ndisp = length(m);
    [r,k] = sort(rdis(m));
    typep = type(m(k));
    if ndisp>=2
        dtrial = zeros(ndisp,1);
        for i = 2:ndisp
            j = i-1;
            if (r(i)-r(j)<=Le) && (typep(i)+typep(j)==0) && (dtrial(j)==0)
                dtrial(j) = 1;
                dtrial(i) = 1; 
            end
        end
        destroy(m(k)) = dtrial;
    end
end
