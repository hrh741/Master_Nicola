% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Function    : dislocation escape from free surfaces
%   Description : call by Main.m
%   Outstanding issues: Dislocations escaping into holes not defined.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tol2=1e-4;
% noutb4 = nout;

% %% exit from left surface
% disnodes = [xdis(1:ndis) ydis(1:ndis)];
% d_min = dpolySW(disnodes,samplenode);
% distemp = (abs(d_min)<=tol & irmbound(1:ndis)>0);
%
% idis = distemp & ydis(1:ndis)>0+tol & ydis(1:ndis)<bbox(2,2)-tol & xdis(1:ndis)<bbox(1,1)+tol;
% if sum(idis)>0
%     [b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,plane,pinned,irmbound,eta_obs,...
%         vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,ngsourceOut,planeOut,nout] = ...
%         removeDis(1,idis,b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,...
%         plane,pinned,irmbound,eta_obs,vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,...
%         ngsourceOut,planeOut,nout);
% end
%

% %% exit from right surface
% disnodes = [xdis(1:ndis) ydis(1:ndis)];
% d_min = dpolySW(disnodes,samplenode);
% distemp = (abs(d_min)<=tol & irmbound(1:ndis)>0);
%
% idis = distemp & ydis(1:ndis)>0+tol & ydis(1:ndis)<bbox(2,2)-tol & xdis(1:ndis)>bbox(2,1)-tol;
% if sum(idis)>0
%     [b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,plane,pinned,irmbound,eta_obs,...
%         vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,ngsourceOut,planeOut,nout] = ...
%         removeDis(3,idis,b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,...
%         plane,pinned,irmbound,eta_obs,vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,...
%         ngsourceOut,planeOut,nout);
% end


%% exit from top surface
disnodes = [xdis(1:ndis) ydis(1:ndis)];
d_min = dpolySW(disnodes,samplenode); % dpolySW(p,pv)
distemp = (abs(d_min)<=tol & irmbound(1:ndis)>0);

idis = distemp & ydis(1:ndis)>bbox(2,2)-tol;
if sum(idis)>0
    [b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,plane,pinned,irmbound,eta_obs,...
        vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,ngsourceOut,planeOut,nout] = ...
        removeDis(2,idis,b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,...
        plane,pinned,irmbound,eta_obs,vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,...
        ngsourceOut,planeOut,nout);
end

%% exit from bottom surface
disnodes = [xdis(1:ndis) ydis(1:ndis)];
d_min = dpolySW(disnodes,samplenode); % dpolySW(p,pv)
distemp = (abs(d_min)<=tol & irmbound(1:ndis)>0);
idis = distemp & ydis(1:ndis)<0+tol;
if sum(idis)>0
    [b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,plane,pinned,irmbound,eta_obs,...
        vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,ngsourceOut,planeOut,nout] = ...
        removeDis(4,idis,b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,...
        plane,pinned,irmbound,eta_obs,vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,...
        ngsourceOut,planeOut,nout);
end

%%  exit from void surface
% if holetype==0 || holetype==2 % Circular and Elliptical Void
%     % above mean
%     distemp = zeros(ndis,1);
%     disnodes = [xdis(1:ndis) ydis(1:ndis)];
%     for i_h=1:size(holes)
%         Hn = holes{i_h,1}; % hole nodes
%         d_min=dpolySW(disnodes,Hn);
%         distemp1 = (abs(d_min)<=tol & irmbound(1:ndis)>0);
%         distemp = (distemp==1 | distemp1==1);
%         distemp = distemp;
% %         idis = distemp & ydis(1:ndis)>mean([bbox(1,2),bbox(2,2)]); % +tol;
%         idis = distemp & ydis(1:ndis)>mean([max(Hn(:,2)),min(Hn(:,2))]);
%         if sum(idis)>0
%             [b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,plane,pinned,irmbound,eta_obs,...
%                 vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,ngsourceOut,planeOut,nout] = ...
%                 removeDis(4,idis,b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,...
%                 plane,pinned,irmbound,eta_obs,vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,...
%                 ngsourceOut,planeOut,nout);
%         end
%     end
%     
%     % below mean
%     distemp = zeros(ndis,1);
%     disnodes = [xdis(1:ndis) ydis(1:ndis)];
%     for i_h=1:size(holes)
%         Hn = holes{i_h,1}; % hole nodes
%         d_min=dpolySW(disnodes,Hn);
%         distemp1 = (abs(d_min)<=tol & irmbound(1:ndis)>0);
%         distemp = (distemp==1 | distemp1==1);
%         distemp = distemp;
% %         idis = distemp & ydis(1:ndis)<mean([bbox(1,2),bbox(2,2)]); %-tol;
%         idis = distemp & ydis(1:ndis)<mean([max(Hn(:,2)),min(Hn(:,2))]);
%         if sum(idis)>0
%             [b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,plane,pinned,irmbound,eta_obs,...
%                 vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,ngsourceOut,planeOut,nout] = ...
%                 removeDis(2,idis,b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,...
%                 plane,pinned,irmbound,eta_obs,vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,...
%                 ngsourceOut,planeOut,nout);
%         end
%     end
%     
%     % ==mean & xdis>0
%     distemp = zeros(ndis,1);
%     disnodes = [xdis(1:ndis) ydis(1:ndis)];
%     for i_h=1:size(holes)
%         Hn = holes{i_h,1}; % hole nodes
%         d_min=dpolySW(disnodes,Hn);
%         distemp1 = (abs(d_min)<=tol & irmbound(1:ndis)>0);
%         distemp = (distemp==1 | distemp1==1);
%         distemp = distemp;
%         
% %         idis = distemp & ydis(1:ndis)==mean([bbox(1,2),bbox(2,2)]) & xdis(1:ndis)>mean([bbox(1,1),bbox(2,1)]);
%         idis = distemp & ydis(1:ndis)==mean([max(Hn(:,2)),min(Hn(:,2))]) & xdis(1:ndis)>mean([max(Hn(:,1)),min(Hn(:,1))]);
%         if sum(idis)>0
%             [b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,plane,pinned,irmbound,eta_obs,...
%                 vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,ngsourceOut,planeOut,nout] = ...
%                 removeDis(3,idis,b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,...
%                 plane,pinned,irmbound,eta_obs,vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,...
%                 ngsourceOut,planeOut,nout);
%         end
%     end
%     
%     % ==mean & xdis<0
%     distemp = zeros(ndis,1);
%     disnodes = [xdis(1:ndis) ydis(1:ndis)];
%     for i_h=1:size(holes)
%         Hn = holes{i_h,1}; % hole nodes
%         d_min=dpolySW(disnodes,Hn);
%         distemp1 = (abs(d_min)<=tol & irmbound(1:ndis)>0);
%         distemp = (distemp==1 | distemp1==1);
%         distemp = distemp;
%         
% %         idis = distemp & ydis(1:ndis)==mean([bbox(1,2),bbox(2,2)]) & xdis(1:ndis)<mean([bbox(1,1),bbox(2,1)]);
%         idis = distemp & ydis(1:ndis)==mean([max(Hn(:,2)),min(Hn(:,2))]) & xdis(1:ndis)<mean([max(Hn(:,1)),min(Hn(:,1))]);
%         if sum(idis)>0
%             [b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,plane,pinned,irmbound,eta_obs,...
%                 vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,ngsourceOut,planeOut,nout] = ...
%                 removeDis(1,idis,b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,...
%                 plane,pinned,irmbound,eta_obs,vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,...
%                 ngsourceOut,planeOut,nout);
%         end
%     end
% end
% 
% if holetype==1 % Rectangular Void
%     disnodes = [xdis(1:ndis) ydis(1:ndis)];
%     distemp = zeros(ndis,1);
%     for i_h=1:size(holes)
%         Hn = holes{i_h,1}; % hole nodes
%         d_min=dpolySW(disnodes,Hn);
%         distemp1 = (abs(d_min)<=tol & irmbound(1:ndis)>0);
%         distemp = (distemp==1 | distemp1==1);
%         distemp = distemp;
%         
%         % exit from left surface of void
%         idis = distemp & ydis(1:ndis)>min(Hn(:,2))+tol & ydis(1:ndis)<max(Hn(:,2))-tol & xdis(1:ndis)<min(Hn(:,1))+tol;
%         if sum(idis)>0
%             [b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,plane,pinned,irmbound,eta_obs,...
%                 vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,ngsourceOut,planeOut,nout] = ...
%                 removeDis(3,idis,b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,...
%                 plane,pinned,irmbound,eta_obs,vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,...
%                 ngsourceOut,planeOut,nout);
%         end
%     end
%     
%     distemp = zeros(ndis,1);
%     disnodes = [xdis(1:ndis) ydis(1:ndis)];
%     for i_h=1:size(holes)
%         Hn = holes{i_h,1}; % hole nodes
%         d_min=dpolySW(disnodes,Hn);
%         distemp1 = (abs(d_min)<=tol & irmbound(1:ndis)>0);
%         distemp = (distemp==1 | distemp1==1);
%         distemp = distemp;
%         
%         % exit from right surface of void
%         idis = distemp & ydis(1:ndis)>min(Hn(:,2))+tol & ydis(1:ndis)<max(Hn(:,2))-tol & xdis(1:ndis)>max(Hn(:,1))-tol;
%         if sum(idis)>0
%             [b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,plane,pinned,irmbound,eta_obs,...
%                 vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,ngsourceOut,planeOut,nout] = ...
%                 removeDis(1,idis,b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,...
%                 plane,pinned,irmbound,eta_obs,vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,...
%                 ngsourceOut,planeOut,nout);
%         end
%     end
%     
%     distemp = zeros(ndis,1);
%     disnodes = [xdis(1:ndis) ydis(1:ndis)];
%     for i_h=1:size(holes)
%         Hn = holes{i_h,1}; % hole nodes
%         d_min=dpolySW(disnodes,Hn);
%         distemp1 = (abs(d_min)<=tol & irmbound(1:ndis)>0);
%         distemp = (distemp==1 | distemp1==1);
%         distemp = distemp;
%         
%         % exit from top surface of void
%         idis = distemp & ydis(1:ndis)>max(Hn(:,2))-tol;
%         if sum(idis)>0
%             [b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,plane,pinned,irmbound,eta_obs,...
%                 vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,ngsourceOut,planeOut,nout] = ...
%                 removeDis(4,idis,b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,...
%                 plane,pinned,irmbound,eta_obs,vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,...
%                 ngsourceOut,planeOut,nout);
%         end
%     end
%     
%     distemp = zeros(ndis,1);
%     disnodes = [xdis(1:ndis) ydis(1:ndis)];
%     for i_h=1:size(holes)
%         Hn = holes{i_h,1}; % hole nodes
%         d_min=dpolySW(disnodes,Hn);
%         distemp1 = (abs(d_min)<=tol & irmbound(1:ndis)>0);
%         distemp = (distemp==1 | distemp1==1);
%         distemp = distemp;
%         
%         % exit from bottom surface of void
%         idis = distemp & ydis(1:ndis)<min(Hn(:,2))+tol;
%         if sum(idis)>0
%             [b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,plane,pinned,irmbound,eta_obs,...
%                 vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,ngsourceOut,planeOut,nout] = ...
%                 removeDis(2,idis,b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,...
%                 plane,pinned,irmbound,eta_obs,vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,...
%                 ngsourceOut,planeOut,nout);
%         end
%     end
% end

% % c++ stores arrays with row-major ordering
% % node_slipC created in applyBC01, before dispMex
