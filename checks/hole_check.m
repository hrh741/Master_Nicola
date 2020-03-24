% hole nodes
jH=[];
for iH=1:size(holes{1,1},1)
    [j1H,k] = find(holes{1,1}(iH,:)==p);
    jH = [jH; j1H];
end
% indH = unique(jH);
pindH = p(jH,:); % for the sake of plotting


% because 4 nodes of grains are also included in hole
% Grain nodes
jG=[];
for iG=1:size(gnodes{1,1},1)
    [j1G,k] = find(gnodes{1,1}(iG,:)==p); % gnodes contain the grain boundary nodes
    jG = [jG; j1G];
end
% indG = unique(jG);
pindG = p(jG,:); % for the sake of plotting


Hnode_no = setdiff(jH, jG);
pHnode_no = p(Hnode_no,:); % for the sake of plotting

%%
figure; hold on;
for ng = 1:size(holes,1)
    xnode = holes{ng,1};
    for inh=1:size(xnode,1)
        plot(xnode(inh,1),xnode(inh,2),'ro','LineWidth',1);
        pause
    end
%     plot(xnode(:,1),xnode(:,2),'b.','LineWidth',1);
end

%%
figure; hold on;
for inh=1:size(Svoid,1)
    plot(p(Svoid(inh),1),p(Svoid(inh),2),'ro','LineWidth',1);
    pause
end

%%
figure; hold on;
for ng = 1:size(holes,1)
    xnode = holes{ng,1};
    plot(xnode(:,1),xnode(:,2),'b.','LineWidth',1);
end
plot(p(Svoid,1),p(Svoid,2),'ro','LineWidth',1);

%%
figure; hold on; axis equal;
scatter(pindG(:,1), pindG(:,2), 'kx')
scatter(pindH(:,1), pindH(:,2), 'ko')
scatter(pHnode_no(:,1), pHnode_no(:,2))
hold off;


%%
figure; hold on; axis equal;
plot(pindG(:,1), pindG(:,2))
plot(pindH(:,1), pindH(:,2))
plot(pHnode_no(:,1), pHnode_no(:,2))
hold off;
