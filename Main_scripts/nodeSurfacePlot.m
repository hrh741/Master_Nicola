nSurface = []; nSurfacep = [];
nSurface1 = [];

for i = 1:size(nodeSurface,1)
nSurface = nodeSurface(i,1);
nSurface1 = p(nSurface,:);
nSurfacep = [nSurfacep; nSurface1];
plot(nSurfacep(:,1),nSurfacep(:,2), 'ro')
axis([-1 5 -1 3])
title('nodeSurface')
% pause


end