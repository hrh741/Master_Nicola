% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Function     : InSets
%   Last edited  : 2 November, 2018 - SW
%   Description  : called by Meshgen.m
%                   Define sets of nodes
%                    BCs 1,2 strain controlled, BCs 3 stress controlled
%   Outstanding issues : 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Stop = [];
Sbot = [];
Sleft = [];
Sright = [];
Svoid = [];
Svoid = nodeSurfaceVoid;


for i = 1:size(p,1)
    if ismember(i,nodeSurfacef) % return 1 and 0
        if (p(i,1)<bbox(2,1)-1e-8 && p(i,1)>0+1e-8 && p(i,2)>(bbox(2,2)-1e-8))
            % left && right && top 
            Stop = [Stop; i];
        end
        if (p(i,1)<bbox(2,1)-1e-8 && p(i,1)>0+1e-8 && p(i,2)<(bbox(1,2)+1e-8))
            Sbot = [Sbot; i];
        end
        if (p(i,1)<=bbox(1,1)+1e-8) 
            Sleft = [Sleft; i];
        end
        if (p(i,1)>=bbox(2,1)-1e-8) 
            Sright = [Sright; i];
        end
    end
%     if ismember(i,nodeSurfaceVoid)
%                 Svoid = [Svoid; i];
%     end

end

% first internal node
[m,n] = sort(p(Sleft,2)); % sort in ascending order
internalNode = Sleft(n(round(length(n)/2)));
% % second internal node
% [m,n] = sort(p(Sright,2));
% internalNode = [internalNode; Sright(n(round(length(n)/2)))];

%% plot
figure; hold on;

for ng = 1:ngr
    xnodes = polynode{ng,1}; % position of grain nodes
    fill(xnodes(:,1),xnodes(:,2),rand(1,3));
    holeind=find([holes{:,2}]'==ng);
    for i_h=1:length(holeind)
        Hn = holes{holeind(i_h),1}; % hole nodes
        fill(Hn(:,1),Hn(:,2),rand(1,3));
        axis equal
    end
end
for i = 1:size(nc,1)
    plot(p(nc(i,[1:3,1]),1),p(nc(i,[1:3,1]),2),'k-','LineWidth',1);
    % text(mean(p(t(i,[1:3,1]),1)),mean(p(t(i,[1:3,1]),2)),num2str(i));
end

plot(p(Stop,1),p(Stop,2),'rv');
plot(p(Sbot,1),p(Sbot,2),'b^');
plot(p(Sleft,1),p(Sleft,2),'g>');
plot(p(Sright,1),p(Sright,2),'m<');
plot(p(Svoid,1),p(Svoid,2),'kx')

plot(p(internalNode,1),p(internalNode,2),'k*');

axis equal; axis off;
ax=axis;axis(ax*1.001);

hold off

%%

if BCs==1 % strain controlled right surface
    gammaMixed2X = [Sleft]; % 2*gammaMixed2-1 = x-fixed & 2*gammaMixed2 = y-free
    gammaMixed2X(gammaMixed2X==internalNode) = []; % move internal node from gammaMixed2 to gammau
    gammauX = internalNode;
    gammatX = [Stop;Sbot;Svoid;Sright];    
    
    gammaMixed1 = []; % 2*gammaMixed1 = y-fixed & 2*gammaMixed1-1 = x-free % USED WHEN BOTTOM SURFACE FIXED AND TOP COMPRESSED
%    gammaMixed1(gammaMixed2==internalNode) = []; % move internal node from gammaMixed1 to gammau
    gammaMixed2 = [Sleft;Sright]; % 2*gammaMixed2-1 = x-fixed & 2*gammaMixed2 = y-free
    gammaMixed2(gammaMixed2==internalNode) = []; % move internal node from gammaMixed2 to gammau
    gammau = internalNode; 
    gammat = [Stop;Sbot;Svoid];
    
    fixedDofs = [2*gammaMixed1; 2*gammaMixed2-1; 2*gammau-1; 2*gammau];
    freeDofs =  [2*gammaMixed1-1; 2*gammaMixed2; 2*gammat-1; 2*gammat];
    
    fixedDofsX = [2*gammaMixed2X-1; 2*gammauX-1; 2*gammauX];
    freeDofsX =  [2*gammaMixed2X; 2*gammatX-1; 2*gammatX];
    
    
    if length([fixedDofs;freeDofs])>length(unique([fixedDofs;freeDofs]))
        error('DOF error in boundary conditions');
    end
    if length([fixedDofsX;freeDofsX])>length(unique([fixedDofsX;freeDofsX]))
        error('DOF error in boundary conditions');
    end
    

    
    
elseif BCs==2 % strain controlled top surface
    gammaMixed1 = [Stop;Sbot]; % y fixed: 2*gammaMixed1 = fixed & x free: 2*gammaMixed1-1 = free
    gammaMixed1(gammaMixed1==internalNode) = []; % move internal node from gammaMixed1 to gammau
    gammaMixed2 = []; %  2*gammaMixed2-1 = x-fixed & 2*gammaMixed2 = y-free
%     gammaMixed2(gammaMixed2==internalNode) = []; % move internal node from gammaMixed2 to gammau
    gammau = internalNode;
    gammat = [Sleft;Sright;Svoid];
    
    fixedDofs = [2*gammaMixed1; 2*gammaMixed2-1; 2*gammau-1; 2*gammau];
    freeDofs =  [2*gammaMixed1-1; 2*gammaMixed2; 2*gammat-1; 2*gammat];
    
    if length([fixedDofs;freeDofs])>length(unique([fixedDofs;freeDofs]))
        error('DOF error in boundary conditions');
    end
    
elseif BCs==3 % stress controlled right surf
    gammaMixed1 = [];
    gammaMixed2 = [Sleft]; % x fixed: 2*gammaMixed2-1 = fixed & y free: 2*gammaMixed2 = free
    gammaMixed2(gammaMixed2==internalNode) = []; % move internal node from gammaMixed2 to gammau
    gammau = internalNode;
    gammat = [Stop;Sbot;Sright;Svoid];
    
    fixedDofs = [2*gammaMixed1; 2*gammaMixed2-1; 2*gammau-1; 2*gammau];
    freeDofs =  [2*gammaMixed1-1; 2*gammaMixed2; 2*gammat-1; 2*gammat];
    
    if length([fixedDofs;freeDofs])>length(unique([fixedDofs;freeDofs]))
        error('DOF error in boundary conditions');
    end
    
elseif BCs==4 % other - here defined as all boundaries being traction free
    gammaMixed1 = [];
    gammaMixed2 = [];
    gammau = internalNode;
    gammat = [Stop;Sbot;Sleft;Sright;Svoid];
    gammat(gammat==internalNode) = []; % move internal node from gammat to gammau

    
    fixedDofs = [2*gammaMixed1; 2*gammaMixed2-1; 2*gammau-1; 2*gammau];
    freeDofs =  [2*gammaMixed1-1; 2*gammaMixed2; 2*gammat-1; 2*gammat];
    
    if length([fixedDofs;freeDofs])>length(unique([fixedDofs;freeDofs]))
        error('DOF error in boundary conditions');
    end
end