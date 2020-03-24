figure; 
   hold on;
for ng = 1:ngr
    plot(p(nodeonsample{ng,1},1),p(nodeonsample{ng,1},2),'g+')
    axis([-1 5 -1 3])
    hold off
    pause
    
end



%%
figure; hold on;
nSample = []; nSamplep = [];
col = ['r','g','b','k'];
for i = 1:size(nodeonsample,1)
nSample = nodeonsample{i,1};
nSamplep = p(nSample,:);
plot(nSamplep(:,1),nSamplep(:,2), 'MarkerEdgeColor', col(i), 'Marker','o', 'LineStyle','none')
axis([-1 5 -1 3])
title('nodeonsample')
pause

end