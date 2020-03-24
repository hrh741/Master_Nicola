% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Function    : AddNodes
%   Last edited : 6 November, 2018 - SW
%   Description : called by InMeshGen.m
%                     Add nodes on the edge of each polygon
%   Outstanding issues : 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Outnodes] = AddNodes(Npolygon,minelesize)

Outnodes = Npolygon(1,:);

for iline = 1:size(Npolygon,1)-1
    % size(Npolygon,1) = number of rows of Npolygon i.e. 1
    
    twonodes = Npolygon(iline:iline+1,:);
    L = sqrt(sum((twonodes(1,:)-twonodes(2,:)).^2));
    Nnodes = ceil(L/minelesize)-1; % ceil round towards plus infinity

    for inode = 1:Nnodes
        newnode = [twonodes(1,:)+(twonodes(2,:)-twonodes(1,:))*inode/(Nnodes+1)];
        Outnodes = [Outnodes;newnode];
    end
    Outnodes = [Outnodes;twonodes(2,:)];
end
% inpolygon(x,y,pv,qv)
[in,on] = inpolygon(Outnodes(:,1),Outnodes(:,2),Npolygon(:,1),Npolygon(:,2));
if any(on~=1); plot(Outnodes(find(on~=1),1),Outnodes(find(on~=1),2), 'g*'); error('Add node error'); end

%   figure; hold on; axis equal;
%   plot(Npolygon(:,1),Npolygon(:,2),'bo');
%   plot(Outnodes(:,1),Outnodes(:,2),'rx');
%   hold off
