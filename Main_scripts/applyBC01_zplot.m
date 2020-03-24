
%% DD subproblem
    uhat = zeros(2*(mno),ngr);
    rDD  = zeros(2*(mno),1);
    
if (ndis+nout)>0
    for ng = 1:ngr
    
    u_node = [gammaMixed1_grain{ng}; gammaMixed2_grain{ng}; gammau_grain{ng}];
    % nout excluded
    utilda = displacementMexZplot(int32(u_node),int32(nSystems),p(:,1),p(:,2),lambdaC,...
        uConst(ng),nu(ng),b(ngsource==ng),mno,length(u_node),sum(ngsource==ng),...
        sum(ngsourceOut==ng),xdis(ngsource==ng),ydis(ngsource==ng),alpha(ngsource==ng),... % sum(ngsourceOut==ng)
        int32(type(ngsource==ng)),xdisOut(ngsourceOut==ng),ydisOut(ngsourceOut==ng),...
        alphaOut(ngsourceOut==ng),int32(typeOut(ngsourceOut==ng)),bOut(ngsourceOut==ng));
   
    ftilda = tildanodes(nodeonGB{ng,1},xdis(ngsource==ng),ydis(ngsource==ng),type(ngsource==ng),...
        alpha(ngsource==ng),b(ngsource==ng),p,edgeType(ng),mno,nSystems(ng),lambda(:,ng));

    fhat = zeros(2*(mno),1);
    % Shishvan 2010 i.e. hat = DD - tilda 
    fhat(freeDofs_grain{ng,1})  = 0 - ftilda(freeDofs_grain{ng,1});     % zero traction on sample zero traction (for DD problem) 
    uhat(fixedDofs_grain{ng,1},ng) = 0 - utilda(fixedDofs_grain{ng,1}); % zero displacement on grain interfaces (for DD problem)
    
    fhat = fhat-kg_grain{ng}(:,fixedDofs_grain{ng,1})*uhat(fixedDofs_grain{ng,1},ng); % find traction at fixedDofs
    fhat(fixedDofs_grain{ng,1}) = bcwt_grain(ng,1)*uhat(fixedDofs_grain{ng,1},ng);
    uhat(allDofs_grain{ng,1},ng) = (Umatrix_grain{ng,1}\(Lmatrix_grain{ng,1}\fhat(allDofs_grain{ng,1})));
    rhat = kg_grain{ng,1}*uhat(:,ng); % CHECK
    rDD  = rDD + rhat + ftilda; % Total traction of DD sub-problem (Shishvan 2010)
    end
end

%% FE elasticity problem
uFE  = zeros(2*(mno),1);
fFE  = zeros(2*(mno),1);
u    = zeros(2*(mno),1);
fa   = zeros(2*(mno),1);

u(2*Sleft-1) = 0;
u(2*Sright-1) = U;
u(2*internalNode) = 0;
uFE(fixedDofs) = u(fixedDofs);
fFE(freeDofs) = fa(freeDofs);
fFE(interfaceDofs) = -rDD(interfaceDofs); % Negative traction on GrainInterfaces(FE)

fFE = fFE - kg(:,fixedDofs)*uFE(fixedDofs);
fFE(fixedDofs) = bcwt*uFE(fixedDofs);

% uFE = Umatrix\(Lmatrix\fFE);
uFE(p_rcm) = Umatrix\(Lmatrix\fFE(p_rcm)); % With reverse Cuthill-Mckee applied
rFE = kg*uFE;


%% Enr step problem
fenr_invhat = zeros(2*mno,1);
% fenr_invhat = XFEM_displacement(planeOut, nout, nc, p, typeOut, slangle, );
XFEM_displacement;

uenr  = zeros(2*(mno),1);
fenr  = zeros(2*(mno),1);
fenr = fenr - fenr_invhat;

fenr = fenr - kg(:,fixedDofs)*uenr(fixedDofs);
fenr(fixedDofs) = bcwt*uenr(fixedDofs);

uenr(p_rcm) = Umatrix\(Lmatrix\fenr(p_rcm)); % With reverse Cuthill-Mckee applied
fenr = kg*uenr;


%% force/dispacement at the loading boundary
fend = sum(rFE(2*Sright-1)+rDD(2*Sright-1));
uend = mean(uFE(2*Sright-1));
