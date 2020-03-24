% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Function    : findbelong
%   Description : call by Main.m
%                 find the dislocations belong to which element
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [belong] = findbelong(ndis,xdis,ydis,plane,p,t,elegrainindex,ngL,flag,gridbelong)

m2 = find(gridbelong==ngL);
ndis = length(m2);
xdis = xdis(m2);
ydis = ydis(m2);

belong = zeros(ndis,1);
nele = size(t,1);
for i = 1:nele
    if isempty(elegrainindex) || isempty(ngL)
        in = inpolygon(xdis,ydis,p(t(i,[1:3,1]),1),p(t(i,[1:3,1]),2));
        belong(in) = i;
    else
        if elegrainindex(i)==ngL
            in = inpolygon(xdis,ydis,p(t(i,[1:3,1]),1),p(t(i,[1:3,1]),2));
            belong(in) = i;
        end
    end
end

while any(belong==0)
    if flag==1
         n = find(belong==0);
        error('error in findbelong.m for grid point');
       
        closest=n+1;
        if any(n==length(belong))
            closest(end)=closest(end)-2;
        end
        belong(n) = belong(closest);
    elseif flag==2
        disp('error in findbelong.m for dislocation');
        disp('dislocation index:')
        n = find(belong==0);
        
        %     disp(num2str(n));
        
        %     disp('finding element again')
        %     for i = 1:nele
        %         in = inpolygon(xdis(n),ydis(n),p(t(i,[1:3,1]),1),p(t(i,[1:3,1]),2));
        %         if in == 1
        %             belong(n) = i;
        %         end
        %     end
        
        disp('Assigning closest dislocation,s belong')
        for i=1:length(n)
            dis_sameplane = find(plane(1:ndis)==plane(n(i)));
            dis_sameplane = dis_sameplane(dis_sameplane~=n(i));
            dist2 = (xdis(dis_sameplane)-xdis(n(i))).^2 + (ydis(dis_sameplane)-ydis(n(i))).^2; % find distances from problem dislocation
            closest = dis_sameplane(dist2 == min(dist2)); % find closest dislocation
            belong(n(i)) = belong(closest);
        end
    end
    %     error('!')
end