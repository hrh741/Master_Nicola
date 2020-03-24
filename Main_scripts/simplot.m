% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Function    : Plot subroutine
%   Description : call by Main.m
%                 plot dislocation structure and stress strain curve
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function simplot(FEND,UEND,increment,t,xdis,ydis,ndis,type,alpha,source,...
%     xSource,ySource,sourcecount,xObs,yObs,obscount,bbox,polynode,...
%     ngr,dis_size,pointSize)
%%

figure(1); clf; hold on;

for ng = 1:ngr
    Npolygon1 = polynode{ng};
    plot(Npolygon1(:,1),Npolygon1(:,2),'k-','Linewidth',2);
%     holeind=find([holes{:,2}]'==ng);
%     bnodes=[];
%     for i_h=1:length(holeind)
%         Hn = holes{holeind(i_h),1}; % hole nodes
%         bnodes = [bnodes; Hn];
%         plot(bnodes(:,1), bnodes(:,2), 'k-','Linewidth',2)
%     end
    
end

for ng = 1:ngr
    plot(xSource(1:sourcecount(ng),ng),ySource(1:sourcecount(ng),ng),'y.','MarkerSize',pointSize)
%     plot(xObs(1:obscount(ng),ng),yObs(1:obscount(ng),ng),'c.','MarkerSize',pointSize)
end

plot(xdis(type==1),ydis(type==1), 'r.', 'Markersize', 8);
plot(xdis(type==-1),ydis(type==-1), 'b.', 'Markersize', 8);
plot(xdis(source==0),ydis(source==0), 'c.', 'Markersize', 8); % debris
plot(xdis(source==-1),ydis(source==-1), 'y.', 'Markersize', 8); % re-emmission
plot(xdis(source==-2),ydis(source==-2), 'g.', 'Markersize', 8); % transmitted


% for idis = 1:ndis
%     if source(idis)==0 % debris
%         plot(xdis(idis)+[-1,1]*dis_size/2*cos(alpha(idis)),ydis(idis)+[-1,1]*dis_size/2*sin(alpha(idis)),'b-','LineWidth',2)
%         plot(xdis(idis)+[0,dis_size*0.8*cos(alpha(idis)+pi/2*type(idis))],ydis(idis)+[0,dis_size*0.8*sin(alpha(idis)+pi/2*type(idis))],'b-','LineWidth',2)
%     elseif source(idis)==-1 % re-emmission
%         plot(xdis(idis)+[-1,1]*dis_size/2*cos(alpha(idis)),ydis(idis)+[-1,1]*dis_size/2*sin(alpha(idis)),'r-','LineWidth',2)
%         plot(xdis(idis)+[0,dis_size*0.8*cos(alpha(idis)+pi/2*type(idis))],ydis(idis)+[0,dis_size*0.8*sin(alpha(idis)+pi/2*type(idis))],'r-','LineWidth',2)
%     else % normal
%         plot(xdis(idis)+[-1,1]*dis_size/2*cos(alpha(idis)),ydis(idis)+[-1,1]*dis_size/2*sin(alpha(idis)),'k-','LineWidth',2)
%         plot(xdis(idis)+[0,dis_size*0.8*cos(alpha(idis)+pi/2*type(idis))],ydis(idis)+[0,dis_size*0.8*sin(alpha(idis)+pi/2*type(idis))],'k-','LineWidth',2)
%     end
% end

axis equal; axis off;
ax = axis; axis(ax*1.001);
drawnow;
hold off

%%

figure(2); clf; hold on;
set(gcf,'Color','w'); box on;
xlabel('Strain, \epsilon [ - ]','FontSize', 12);
ylabel('Stress, \sigma [MPa]','FontSize', 12);
% ylabel('\sigma/\mu','FontSize', 12);
set(gca,'LineWidth',1.5,'FontSize',12,'Box','on','Color','w');
plot(UEND(1:increment)/bbox(2,1),1e6*FEND(1:increment)/bbox(2,2),'-','LineWidth',2);
% plot(UEND(1:increment)/bbox(2,1),FEND(1:increment)/(mu*bbox(2,2)),'r-','LineWidth',2);
hold off
% figure(3); clf; hold on;
% set(gcf,'Color','w'); box on;
% xlabel('Time [s]','FontSize', 12);
% ylabel('Stress, \sigma [MPa]','FontSize', 12);
% set(gca,'LineWidth',1.5,'FontSize',12,'Box','on','Color','w');
% plot(TINC(1:increment),1e6*FEND(1:increment)/bbox(2,2),'k-','LineWidth',2);
% hold off

drawnow;