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
    uhat = zeros(2*(mno),ngr);
    rDD  = zeros(2*(mno),1);

if (ndis+nout)>0
    for ng = 1:ngr
        
        u_node = [gammaMixed1_grain{ng}; gammaMixed2_grain{ng}; gammau_grain{ng}];
        % nout included
        utilda = displacementMex(int32(u_node),int32(nSystems),p(:,1),p(:,2),lambdaC,...
            uConst(ng),nu(ng),b(ngsource==ng),mno,length(u_node),sum(ngsource==ng),...
            sum(ngsourceOut==ng),xdis(ngsource==ng),ydis(ngsource==ng),alpha(ngsource==ng),... % sum(ngsourceOut==ng)
            int32(type(ngsource==ng)),xdisOut(ngsourceOut==ng),ydisOut(ngsourceOut==ng),...
            alphaOut(ngsourceOut==ng),int32(typeOut(ngsourceOut==ng)),bOut(ngsourceOut==ng));       
        
        ftilda = tildanodes(nodeonGB{ng,1},xdis(ngsource==ng),ydis(ngsource==ng),type(ngsource==ng),...
            alpha(ngsource==ng),b(ngsource==ng),p,edgeType(ng),mno,nSystems(ng),lambda(:,ng));
            
        fhat = zeros(2*(mno),1);
        fhat(freeDofs_grain{ng,1})  = 0 - ftilda(freeDofs_grain{ng,1});     % zero traction on sample zero traction
        uhat(fixedDofs_grain{ng,1},ng) = 0 - utilda(fixedDofs_grain{ng,1}); % zero displacement on grain interfaces
        
        fhat = fhat-kg_grain{ng}(:,fixedDofs_grain{ng,1})*uhat(fixedDofs_grain{ng,1},ng); % Matrix reduction
        fhat(fixedDofs_grain{ng,1}) = bcwt_grain(ng,1)*uhat(fixedDofs_grain{ng,1},ng); % Replacement of fhat values at fixedDofs
        uhat(allDofs_grain{ng,1},ng) = (Umatrix_grain{ng,1}\(Lmatrix_grain{ng,1}\fhat(allDofs_grain{ng,1}))); % It will gives us back u=0 at fixed BC
        rhat = kg_grain{ng,1}*uhat(:,ng); % new fhat
        rDD  = rDD+rhat+ftilda; % fhat + ftilda = Traction, which is applied in FE elasticity problem (-rDD)
    end
end


 %% displacement x due to slipping
% uend_stepx = 0;
% for iout = 1:nout
%     
%     islpOut = planeOut(iout); % plane number
%     slangleOut = alphaOut(iout);  
%     
%     % distance pf each node from out_dislocation
%     dx = 0 - xdisOut(iout); % consider nodes to be at origin
%     dy = 0 - ydisOut(iout);
%     r2 = dx*dx+dy*dy;
%     dx_temp = [dx dy]*[ cos(slangleOut) sin(slangleOut)]';
%     dy_temp = [dx dy]*[-sin(slangleOut) cos(slangleOut)]';
%     
%     r2inv = 1/r2;
%     
%     % 2*dx_temp because only considering abovenodes, unlike in dispMex
%     ux_temp = (0.5*dx_temp*dy_temp*r2inv-(1-nu)*atan(dx_temp/dy_temp))*typeOut(iout)*bOut(iout)*uConst;
%     uy_temp = 0.0;
%     ux = cos(slangleOut)*ux_temp - sin(slangleOut)*uy_temp;
%     uy = sin(slangleOut)*ux_temp + cos(slangleOut)*uy_temp;
%     
%     uend_stepx = uend_stepx+ux;
% end
% 
%    uend_stepx= abs(uend_stepx);
%    
%  %
% % FE elasticity problem
% uFE  = zeros(2*(mno),1);
% fFE  = zeros(2*(mno),1);
% u    = zeros(2*(mno),1);
% fa   = zeros(2*(mno),1);
% 
% u(2*Sleft-1) = 0;
% % u(2*Sright-1) = U;
% u(2*Sright-1) = U-2*uend_stepx;
% u(2*internalNode) = 0;
% uFE(fixedDofs) = u(fixedDofs);
% fFE(freeDofs) = fa(freeDofs);
% fFE(interfaceDofs) = -rDD(interfaceDofs);
% fFE = fFE - kg(:,fixedDofs)*uFE(fixedDofs);
% fFE(fixedDofs) = bcwt*uFE(fixedDofs);
% % uFE = Umatrix\(Lmatrix\fFE);
% uFE(p_rcm) = Umatrix\(Lmatrix\fFE(p_rcm)); % With reverse Cuthill-Mckee applied
% rFE = kg*uFE;
% 
%  
% % force/dispacement at the loading boundary
% fend = sum(rFE(2*Sright-1)+rDD(2*Sright-1));
% % uend = mean(uFE(2*Sright-1));
% uend = mean(uFE(2*Sright-1))+2*uend_stepx;


%% XFEM within each increment   
% XFEM subproblem
fenr_invhat = zeros(2*mno,1);
% [fenr_invhat, stress_invhat] = XFEM_displacement(planeOut, nout, nc, p, typeOut, slangle, );
XFEM_displacement;

uenr  = zeros(2*(mno),1);
fenr  = zeros(2*(mno),1);
fenr = fenr - fenr_invhat;

uenr(p_rcmX) = UmatrixX\(LmatrixX\fenr(p_rcmX)); % With reverse Cuthill-Mckee applied
fenr = kg*uenr;

% FE subproblem
uFE  = zeros(2*(mno),1);
fFE  = zeros(2*(mno),1);
u    = zeros(2*(mno),1);
fa   = zeros(2*(mno),1);

u(2*Sleft-1) = 0;
u(2*Sright-1) = U - 2*(uenr(2*Sright-1));

u(2*internalNode) = 0;
uFE(fixedDofs) = u(fixedDofs);
fFE(freeDofs) = fa(freeDofs);
fFE(interfaceDofs) = -rDD(interfaceDofs); % Negative traction on GrainInterfaces(FE)

fFE = fFE - kg(:,fixedDofs)*uFE(fixedDofs);
fFE(fixedDofs) = bcwt*uFE(fixedDofs);  

% uFE = Umatrix\(Lmatrix\fFE);
uFE(p_rcm) = Umatrix\(Lmatrix\fFE(p_rcm)); % With reverse Cuthill-Mckee applied
rFE = kg*uFE;

%%
fend = sum(rFE(2*Sright-1)+rDD(2*Sright-1));% + bcwt*fenr(2*Sright-1));
uend = mean(uFE(2*Sright-1)+2*(uenr(2*Sright-1)));

