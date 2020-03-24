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