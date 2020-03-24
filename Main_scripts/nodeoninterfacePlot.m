figure; hold on;
nInterface = []; nInterfacep = [];
nInterfacep1 = [];

col = ['r','g','b','k'];
for j = 1:size(nodeoninterface,1)
nInterface = nodeoninterface{j,1};
nInterfacep = p(nInterface,:);
plot(nInterfacep(:,1),nInterfacep(:,2),'MarkerEdgeColor', col(j), 'Marker','o', 'LineStyle','none')
axis([-1 5 -1 3])
title('nodeoninterface')
pause


end