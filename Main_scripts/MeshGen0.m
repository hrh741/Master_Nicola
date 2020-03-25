% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Function     : MeshGen
%   Last edited  : 2 November, 2018 - SW
%   Description  : called by Main
%                   Mesh each grain individually, using linear triangular elements
%                   combine into global mesh, define neighbouring
%                   nodes for analytical dislocation forces on nodes
%   Outstanding issues :
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% flag=0;
% pause('on')
% figure; hold on;
Npolygon = []; % coordinates of predefined boundary "nodes" for each grain
for ng = 1:ngr
    xnodes = polynode{ng,1}; % grain nodes position
    % plot(xnodes(:,1), xnodes(:,2), 'k-.')
    bnodes = AddNodes(xnodes,minelesize);
        holeind=find([holes{:,2}]'==ng);
        for i_h=1:length(holeind)
            Hn = holes{holeind(i_h),1}; % hole nodes
            bnodes = [bnodes; AddNodes(Hn,minelesize)];
        end
    % plot(bnodes(:,1), bnodes(:,2), 'bo')
    Npolygon{ng} = bnodes;
    bnodes=[];
end
%%

cd ./utilities/distmesh/;
p = [];     % nodes positions (x,y)
nc = [];    % nodes connectivity
grainelenum = zeros(ngr,1);  % number of elements in each grain (4x1)
elegrainindex = []; % use to index which grain the element belongs to

for ng = 1:ngr
    pv = Npolygon{ng}; % Boundary nodes coordinates of grain
    xnodes = polynode{ng,1};  % Grain coordinates
    fprintf('Mesh grain: %s in process ...\n', ...
        num2str(ng));
    
    figure;
    % p1 = Mesh element nodes position, t1 = Local node number
    [p1,t1] = distmesh2d_SW0(minelesize,elesize,bbox,pv,xnodes,holes,2,ng); % polynode{ng,2}
    Nold = size(p,1);
    p = [p; p1]; % p contains global coordinated of nodes 
    nc = [nc; t1+Nold]; % nodes connectivity
    grainelenum(ng,1) = size(t1,1); % Number of elements in a grain 4x1
    if polynode{ng,2}==2; saveas(gcf,['../../../',directoryOut,'BetaMesh',num2str(elesize),'.fig']); end
    pause
    close all;
end

% find which grain each element belongs to
ng = 1;
for iele = 1:size(nc,1) % 1: number of elements 
    if iele > sum(grainelenum(1:ng,1)) % grainelenum(1,1)
        ng = ng+1;
    end
    elegrainindex(iele) = ng;
end

cd ../../;

%% remove common nodes on the edge
p_temp = round(p*1e10)*1e-10; % temp node position used to find common nodes
dltrow = [];
for ng = 1:ngr
    for i = 1:size(Npolygon{ng},1) % Npolygon = position of boundary nodes 
        Npolygon_temp = round(Npolygon{ng}*1e10)*1e-10;
        nrow = find(ismember(p_temp,Npolygon_temp(i,:),'rows'));
 % ismember(A,S,'rows')= returns 1 when rows of A are also in S

        k2 = [];
        for k = 1:length(nrow)
            if sum(ismember(dltrow,nrow(k)))
                k2 = [k2; k];
            end
        end
        nrow(k2) = [];
        if length(nrow)>1
            dltrow = [dltrow; nrow(2:end)];
            for j = 2:length(nrow)
                nc(ismember(nc,nrow(j))) = nrow(1);
            end
        end
    end
end

p(dltrow,:) = [];
k = nc*0;
for i = 1:length(dltrow)
    k(nc>dltrow(i)) = k(nc>dltrow(i))+1;
end
nc = nc-k;

% Check if extra nodes added to grain boundary
figure;hold on;
for ng = 1:ngr
    xnodes = polynode{ng,1}; % Grain nodes position
    fill(xnodes(:,1),xnodes(:,2),rand(1,3)); % rand(1,3)
    holeind=find([holes{:,2}]'==ng);
    for i_h=1:length(holeind) % 1:number of holes
        Hn = holes{holeind(i_h),1}; % hole nodes
        fill(Hn(:,1),Hn(:,2),rand(1,3));
    end
end
for i = 1:size(nc,1) % number of rows in nc, nc=nodes connectivity 
    plot(p(nc(i,[1:3,1]),1),p(nc(i,[1:3,1]),2),'k-','LineWidth',1);
    %     text(mean(p(nc(i,[1:3,1]),1)),mean(p(nc(i,[1:3,1]),2)),num2str(i));
end
axis equal; axis off;
ax=axis;axis(ax*1.001);

nodeonGB = [];
for ng = 1:ngr
    xnodes = polynode{ng,1};
    fd=dpolySW(p,xnodes);
    holeind=find([holes{:,2}]'==ng);
    for i_h=1:length(holeind)
        Hn = holes{holeind(i_h),1}; % hole nodes
        xnodes=[xnodes; NaN NaN; Hn];
        fd = ddiff(fd,dpolySW(p,Hn));
    end
    nodeonGB = find(abs(fd)<=1e-8);
    xyval = setdiff(p(nodeonGB,:),Npolygon{ng},'rows');
    plot(xyval(:,1),xyval(:,2), 'b.', 'MarkerSize', 10);
    pause
end
close gcf;

%% find nodes on the boundaries of the whole sample
pv = samplenode; % samplenode = bbox
fd=dpolySW(p,pv);
nodeSurface = find(abs(fd)<=1e-8);
nodeSurfaceVoid = [];
nodeSurfaceVoid1= [];
nodeSurfacef = [];

    %  finding nodes on the surface of Void
    for ng = 1:ngr
        holeind=find([holes{:,2}]'==ng);
        for i_h=1:length(holeind)
            Hn = holes{holeind(i_h),1}; % hole nodes
            fd1 = dpolySW(p,Hn);
            nodeSurfaceVoid1 = find(abs(fd1)<=1e-8);
            nodeSurfaceVoid = [nodeSurfaceVoid; nodeSurfaceVoid1];
        end
    end

nodeSurfacef = [nodeSurface; nodeSurfaceVoid]; % just for plotting
% Surface nodes of sample and voids
nodeSurfacef = nodeSurface; % incorporating voids

% Define Sets
InSets;
% saveas(gcf,['../',directoryOut,'Mesh',num2str(elesize),'.fig']);
% pause;

%% Find nodes on the boundaries of each grain, and define node A and B for analytical tildanodes
% nodeonGB anticlockwise node ordering

nodeonGB = []; % nodes on the whole boundaries of each grain ordered anticlockwise
% nodeonGB = nodeoninterface+nodeonsample % including holes
nodeingrain = [];
normalnodeGB = [];
nodeoninterface = []; % nodes on the interfaces between grains
nodeonsample = []; % nodes on the boundaries of the whole sample for each grain

for ng=1:ngr+size(holes,1) % ngr + number of rows of the holes
%     if flag==1
        if ng<ngr+1; xnodes = polynode{ng,1}; else xnodes = holes{ng-ngr,1}; end
%     else
%         if ng<ngr+1; xnodes = polynode{ng,1}; end
%     end
    
    % SW June 2017
    fd=dpolySW(p,xnodes); % excluding holes for now
    nodeonGB{ng,1} = find(abs(fd)<=1e-8);
    
%     if flag==1
        holeind=find([holes{:,2}]'==ng);
        for i_h=1:length(holeind) % including holes
            Hn = holes{holeind(i_h),1}; % hole nodes
            fd = ddiff(fd,dpolySW(p,Hn));
        end
%     end
    % remove small value numerical error caused by distmesh2d.m
    nodeingrain{ng,1} = find(fd<1e-8);
    
    % sort nodeonGB
    acw_order = convhull(p(nodeonGB{ng,1}(:),1),p(nodeonGB{ng,1}(:),2));
    nodeonGB{ng,1} = nodeonGB{ng,1}(acw_order(1:end-1));
    nodeoninterface{ng,1} = nodeonGB{ng,1}(~ismember(nodeonGB{ng,1},nodeSurface),1);
    nodeonsample{ng,1} = nodeonGB{ng,1}(ismember(nodeonGB{ng,1},nodeSurface),1);
    
    gngrain = nodeonGB{ng,1};
    gngraincompile = zeros(length(gngrain),4); % gnit, gnA, gnB, boundangle
    for i=1:length(gngrain)
        gnit = gngrain(i);
        if i==1
            gn1 = gngrain(end); gn2 = gngrain(i+1);
        elseif i==length(gngrain)
            gn1 = gngrain(i-1); gn2 = gngrain(1);
        else
            gn1 = gngrain(i-1);
            gn2 = gngrain(i+1);
        end
        
        if all(gnit ~= Svoid)
            
            % Deal with corner nodes % ADD HOLE CORNERS
            triarea = abs(0.5*((p(gnit,1)-p(gn2,1))*(p(gn1,2)-p(gnit,2))-(p(gnit,1)-p(gn1,1))*(p(gn2,2)-p(gnit,2)))); % check if gn, gn1 and gn2 collinear
            %        if rank([p(gn1,:)-p(gnit,:); p(gn2,:)-p(gnit,:)])>1 % not collinear
            if triarea>1e-10 % not collinear
                if any([gammaMixed1;gammaMixed2;gammau]==gnit) % if on any fixed dof node, consider node to be on fixed dof boundary
                    if any([gammaMixed1;gammaMixed2;gammau]==gn1)
                        gn2 = gnit;
                    elseif any([gammaMixed1;gammaMixed2;gammau]==gn2)
                        gn1 = gnit;
                    else; error('not collinear, on [gammaMixed1;gammaMixed2;gammau] but something wrong');
                    end
                elseif any(gammat==gnit) % if on gammat
                    if any(cell2mat(nodeoninterface)==gn1) % if on gammaB, move to and consider part of gammaT
                        gn1 = gnit;
                    elseif any(cell2mat(nodeoninterface)==gn2)
                        gn2 = gnit;
                    elseif BCs==3 %
                        if any(Sright==gn1) % on upper right corner, consider node as part of right bound where load applied
                            gn2 = gnit;
                        elseif any(Sright==gn2) % on lower right corner, consider node as part of right bound where load applied
                            gn1 = gnit;
                        else % on void surface
                            gn1=gnit;
                        end
                    else; error('not collinear, on gammat but something wrong');
                    end
                elseif any(cell2mat(nodeoninterface)==gn1) && any(cell2mat(nodeoninterface)==gn2) % if all on gammaB
                    gn1 = gnit;
                end
                
            end
        end
        normaltemp = [(p(gn2,1)-p(gn1,1)),(p(gn2,2)-p(gn1,2))];
        normaltemp = normaltemp/sqrt(normaltemp(1)^2+normaltemp(2)^2);        
        normaltemp = [normaltemp(2) -normaltemp(1)];
       
        
        if ng>ngr; normaltemp = -normaltemp; end % flip direction of normal if ng>ngr (i.e. for holes) 
        % out of the material i.e. inside of void
       
            
        normalnodeGB{ng,1}(i,:) = normaltemp; % normal to boundary unit vector
        
        if p(gn1,2) == p(gn2,2) % if horizontal boundary
            if p(gn1,1) > p(gn2,1)
                gnA = gn2; gnB = gn1;
            else; gnA = gn1; gnB = gn2;
            end
        else % inclined boundary
            if p(gn1,2) > p(gn2,2)
                gnA = gn2; gnB = gn1;
            else; gnA = gn1; gnB = gn2;
            end
        end
        
        if gnit == gnA && gnit == gnB
            error('gnit is both gnA and gnB')
        end
        
        % find the angle of boundary, such that normal once aligned is always [0; 1]
        % rotate normal at node by 90 anticlockwise
        normaltemp = normalnodeGB{ng,1}(i,:); % this is unit vector still
        normaltemp = [normaltemp(2) -normaltemp(1)];
        angletemp = acos(dot([1 0],normaltemp/norm(normaltemp))); % angle between this vector and horizontal
        if normaltemp(2) < 0
            angletemp = 2*pi - angletemp;
        end
        
        gngraincompile(i,:) = [gnit; gnA; gnB; angletemp];
                
    end
        if ng<=ngr
            nodeonGB{ng,1} = gngraincompile;
        else % add hole info to its relevant grain
            grainind = holes{ng-ngr,2}; % which grain, hole belongs to
            nodeonGB{grainind,1}     = [nodeonGB{grainind,1}; gngraincompile];
            normalnodeGB{grainind,1} = [normalnodeGB{grainind,1}; normalnodeGB{ng,1}];
            nodeonsample{grainind,1} = [nodeonsample{grainind,1}; nodeonsample{ng,1}];
        end  
%     nodeonsample{ng,1} = nodeonGB{ng,1}(ismember(nodeonGB{ng,1}(:,1),nodeSurface),1);
end
% remove info for nodes in holes

nodeonGB(ngr+1:end,:)=[];
normalnodeGB(ngr+1:end,:)=[];
nodeoninterface(ngr+1:end,:)=[];
nodeingrain(ngr+1:end,:) = [];
nodeonsample(ngr+1:end,:)=[];
% now nodeonGB{ng,1}(1,2,3,4) = globalnodenumber, gnA, gnB, grainboundaryangle

%% plot nodeonGB and normalnodeGB
% figure; hold on;
% for ng=1:ngr
%     for i=1:length(nodeonGB{ng,1})
%         plot(p(nodeonGB{ng,1}(i),1),p(nodeonGB{ng,1}(i),2),'ro'); 
%         %plot(p(nodeonGB{ng,1}(i,2:3),1),p(nodeonGB{ng,1}(i,2:3),2),'bo');
%         %plot(p(nodeonGB{ng}(i,1),1)+[0,normalnodeGB{ng,1}(i,1)],p(nodeonGB{ng}(i,1),2)+[0,normalnodeGB{ng,1}(i,2)],'k-');
%        %pause
%     end
% %      %pause
% end

%%
% plot nodes on GB
% for ng = 1:ngr
%     %     ng=2;
%     figure; hold on;
%     for i = 1:ngr
%         xnodes = polynode{i,1};
%         plot(xnodes(:,1),xnodes(:,2),'k-','LineWidth',2);
%     end
%     xnodes = polynode{ng};
%     fill(xnodes(:,1),xnodes(:,2),rand(1,3));
% %      for i = 1:size(nc,1)
% %          plot(p(nc(i,[1:3,1]),1),p(nc(i,[1:3,1]),2),'k-','LineWidth',1);
% %      end
%     
%     plot(p(nodeonGB{ng}(:,1),1),p(nodeonGB{ng}(:,1),2),'ro','MarkerSize',10)
%     
%     for i = 1:length(nodeonGB{ng})
%         plot(p(nodeonGB{ng}(i,1),1)+[0,normalnodeGB{ng,1}(i,1)],p(nodeonGB{ng}(i,1),2)+[0,normalnodeGB{ng,1}(i,2)],'m-');
%         ptplot = plot(p(nodeonGB{ng}(i,2:3),1), p(nodeonGB{ng}(i,2:3),2), 'r*');
%         text(p(nodeonGB{ng}(i,1),1),p(nodeonGB{ng}(i,1),2), num2str(180*nodeonGB{ng}(i,4)/pi), 'FontSize', 14);
%         delete(ptplot)
%         pause
%     end
%     plot(p(nodeoninterface{ng},1),p(nodeoninterface{ng},2),'b*')
%     plot(p(nodeonsample{ng},1),p(nodeonsample{ng},2),'g+')
%     plot(p(nodeingrain{ng},1),p(nodeingrain{ng},2),'y.')
%     axis equal; axis off;
%     ax=axis;axis(ax*1.001);
%     hold off
%     pause
% %     close gcf
% end

%% plot
figure;hold on;

for ng = 1:ngr
    xnodes = polynode{ng,1};
    fill(xnodes(:,1),xnodes(:,2),rand(1,3));
%     holeind=find([holes{:,2}]'==ng);
%     for i_h=1:length(holeind)
%         Hn = holes{holeind(i_h),1}; % hole nodes
%         fill(Hn(:,1),Hn(:,2),rand(1,3));
%     end
%     pause
end
for i = 1:size(nc,1)
    plot(p(nc(i,[1:3,1]),1),p(nc(i,[1:3,1]),2),'k-','LineWidth',1);
%     text(mean(p(nc(i,[1:3,1]),1)),mean(p(nc(i,[1:3,1]),2)),num2str(i));
end

axis equal; axis off;
ax=axis;axis(ax*1.001);

saveas(figure,['../',directoryOut,'Mesh.fig'])
hold off

fprintf('Generate %d elements and %d nodes.\n', ...
    size(nc,1),size(p,1));
fprintf(logFID, 'Generate %d elements and %d nodes.\n', ...
    size(nc,1),size(p,1));
fprintf('InMeshGen Done!\n')


%%
%
