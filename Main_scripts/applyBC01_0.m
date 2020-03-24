% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Function    : Boundary condition subroutine
%   Description : call by Main.m
%                 apply boundary conditions corresponding to InBC01.m
%
%                 Type 01:
%                   displacement controlled loading with constant strain
%                   rate, apply on the right surface of the sample and the
%                   left surface is constained.
%                   % Exx = +constant
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DD subproblem
% tic
uhat = zeros(2*(mno),ngr);
rDD  = zeros(2*(mno),1);
if ndis>0 || nout>0
    for ng = 1:ngr
        u_node = [gammaMixed1_grain{ng}; gammaMixed2_grain{ng}; gammau_grain{ng}]; % left and right boundary nodes
        
        utilda = displacementMex(int32(u_node),int32(nSystems),p(:,1),p(:,2),lambdaC,...
            uConst(ng),nu(ng),b(ngsource==ng),mno,length(u_node),sum(ngsource==ng),...
            sum(ngsourceOut==ng),xdis(ngsource==ng),ydis(ngsource==ng),alpha(ngsource==ng),...
            int32(type(ngsource==ng)),xdisOut(ngsourceOut==ng),ydisOut(ngsourceOut==ng),...
            alphaOut(ngsourceOut==ng),int32(typeOut(ngsourceOut==ng)),bOut(ngsourceOut==ng));
        
        %     ftilda = tildaforce(nodeonGB{ng},normalnodeGB{ng},intptxy{ng},intptw{ng},...
        %         intptbelong{ng},sum(ngsource==ng),xdis(ngsource==ng),ydis(ngsource==ng),type(ngsource==ng),...
        %         alpha(ngsource==ng),edgeType(ng),mno,b(ngsource==ng));
        
        ftilda = tildanodes(nodeonGB{ng,1},xdis(ngsource==ng),ydis(ngsource==ng),type(ngsource==ng),...
            alpha(ngsource==ng),b(ngsource==ng),p,edgeType(ng),mno,nSystems(ng),lambda(:,ng));
        
        fhat = zeros(2*(mno),1);
        fhat(freeDofs_grain{ng,1})  = 0 - ftilda(freeDofs_grain{ng,1});  % zero traction on sample zero traction
        uhat(fixedDofs_grain{ng,1},ng) = 0 - utilda(fixedDofs_grain{ng,1}); % zero displacement on grain interfaces
        
        fhat = fhat-kg_grain{ng}(:,fixedDofs_grain{ng,1})*uhat(fixedDofs_grain{ng,1},ng);
        fhat(fixedDofs_grain{ng,1}) = bcwt_grain(ng,1)*uhat(fixedDofs_grain{ng,1},ng);
        uhat(allDofs_grain{ng,1},ng) = (Umatrix_grain{ng,1}\(Lmatrix_grain{ng,1}\fhat(allDofs_grain{ng,1})));
        rhat = kg_grain{ng,1}*uhat(:,ng); % CHECK
        rDD  = rDD+rhat+ftilda;
    end
end
% toc
% disp('DD problem done!')
%% FE elasticity problem
% tic
uFE  = zeros(2*(mno),1);
fFE  = zeros(2*(mno),1);
u    = zeros(2*(mno),1);
fa   = zeros(2*(mno),1);

u(2*Sleft-1) = 0;
u(2*Sright-1) = U;
u(2*internalNode) = 0;
uFE(fixedDofs) = u(fixedDofs);
fFE(freeDofs) = fa(freeDofs);
fFE(interfaceDofs) = -rDD(interfaceDofs);
fFE = fFE - kg(:,fixedDofs)*uFE(fixedDofs);
fFE(fixedDofs) = bcwt*uFE(fixedDofs);
% uFE = Umatrix\(Lmatrix\fFE);
uFE(p_rcm) = Umatrix\(Lmatrix\fFE(p_rcm)); % With reverse Cuthill-Mckee applied
rFE = kg*uFE;
% toc
% disp('FE peoblem done!')
%% force/dispacement at the loading boundary
fend = sum(rFE(2*Sright-1)+rDD(2*Sright-1));
uend = mean(uFE(2*Sright-1));