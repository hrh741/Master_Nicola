% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Function     : FEdefine
%   Last edited  : 6 November, 2018 - SW
%   Description  : called by Main.m
%                  define Finite Element related variables - for linear
%                  triangular elements
%   Outstanding issues :
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [mno, mel, D, kg, Lmatrix, Umatrix, bcwt, kg_grain, Lmatrix_grain, Umatrix_grain, ...
    bcwt_grain, gammaMixed1_grain, gammaMixed2_grain, gammau_grain, gammat_grain, ...
    fixedDofs_grain, freeDofs_grain, allDofs_grain, freeDofs, interfaceDofs, B, p_rcm, LmatrixX, UmatrixX, p_rcmX] ...      
    = FEdefine(logFID, p, nc, ngr, nu, e, mu, nodeonsample, nodeingrain, nodeoninterface, ...
    gammaMixed1, gammaMixed2, gammau, gammat, fixedDofs, fixedDofsX, freeDofs, elegrainindex)


tic;
mno = size(p,1);    % no. of nodes
mel = size(nc,1);   % no. of elements

D = zeros(3,3,ngr); % Material Constant used in Hooke's law Stress = D*strain
for ng = 1:ngr
    la = nu(ng)*e(ng)/((1+nu(ng))*(1-2*nu(ng)));
    
    % plane strain
    D(1,1,ng) = la + 2*mu(ng);
    D(1,2,ng) = la;
    D(1,3,ng) = 0;
    D(2,2,ng) = la + 2*mu(ng);
    D(2,3,ng) = 0;
    D(3,3,ng) = mu(ng);
    
    % plane stress
    % D(1,1,ng) = e(ng)/(1-nu(ng)^2);
    % D(1,2,ng) = e(ng)*nu(ng)/(1-nu(ng)^2);
    % D(1,3,ng) = 0;
    % D(2,2,ng) = D(1,1);
    % D(2,3,ng) = 0;
    % D(3,3,ng) = e(ng)/2/(1+nu(ng));
    
    D(2,1,ng) = D(1,2,ng);
    D(3,1,ng) = D(1,3,ng);
    D(3,2,ng) = D(2,3,ng);
end

% 1 integration point of each element
z = zeros(1,2);
z(1,1:2) = [1/3, 1/3]; % Integration point for triangular element

ns = zeros(3,2);
% derivative of shape function(i.e.Na(s1,s2))  w.r.t. sj evaluated at
% integration point q  ns(q,a,j) = (dNa/dsj)s=zq
ns(1,1) = 1;
ns(1,2) = 0;
ns(2,1) = 0;
ns(2,2) = 1;
ns(3,1) = -1;
ns(3,2) = -1;

% calculations within loop for integration point
J = zeros(2,2); % Jacobian Matrix for triangular elements
detJ = 0; % determinant of Jacobian
ke = zeros(6,6,mel); % local stiffness matrix 6x6x1300
w = 1/2; % weight of integration point
for iele = 1:mel % all elements
    for i = 1:2 
        for j = 1:2 
            J(i,j) = ns(1:3,j)'*p(nc(iele,1:3),i); % 10-18
% Jacobian = 'Derivatives of shape functions w.r.t parent domain' * 'position of element(x1,x2,x3)'
        end
    end
    detJ = det(J);
    invJ = inv(J);
    
    nx = zeros(3,2); % derivative of shape functions w.r.t global x,y
    B = zeros(3,6); % (1:3,# dof/element) calculated within loop for each integration point
    for a = 1:3
        for j = 1:2
            for i = 1:2
                nx(a,j) = nx(a,j)+ns(a,i)*invJ(i,j); % 10-18 right top
            end
        end
        B(1,(a-1)*2+1) = nx(a,1);
        B(1,(a-1)*2+2) = 0;
        B(2,(a-1)*2+1) = 0;
        B(2,(a-1)*2+2) = nx(a,2);
        B(3,(a-1)*2+1) = nx(a,2);
        B(3,(a-1)*2+2) = nx(a,1);
    end
    
    ke(:,:,iele) = ke(:,:,iele)+w*detJ*B'*D(:,:,elegrainindex(iele))*B; % 11-7
end

fprintf('Define Local Stiffness Matrix \n')
t = toc
fprintf(logFID,'Define local stiffness matrix ke ... Time = %f \n', t);


% ensure ke is symmetric eg remove any very small entries due to numerical
% error.
tic
kg = sparse(mno*2,mno*2);
for iele = 1:mel
    for a = 1:3
        for b = 1:3
            gna = nc(iele,a);
            gnb = nc(iele,b);
            for i = 1:2
                for j = 1:2
                    kg(2*(gna-1)+i,2*(gnb-1)+j) = kg(2*(gna-1)+i,2*(gnb-1)+j)+ke(2*(a-1)+i,2*(b-1)+j,iele);
                end
            end
        end
    end
end
fprintf('Define Global Stiffness Matrix \n')
t = toc
fprintf(logFID,'Define global stiffness matrix kg ... Time = %f \n', t);


fprintf(logFID,'reformatting K ...\n');
% reformat kg due to applied displacements
K = kg;
bcwt = mean(diag(kg)); % = trace(K)/length(K)
bcwt = full(bcwt);
for n = 1:length(fixedDofs)
    i = fixedDofs(n);
    K(:,i) = 0;
    K(i,:) = 0;
    K(i,i) = bcwt;
end
KX = kg;                                        %
for nX = 1:length(fixedDofsX)                    %
    iX = fixedDofsX(nX);                         %
    KX(:,iX) = 0;                               %
    KX(iX,:) = 0;                               %
    KX(iX,iX) = bcwt;                           %
end             

% fixedDof_XFEM = [2*Sleft-1; 2*Sleft];         %
% KX = kg;                                        %
% for nX = 1:length(fixedDof_XFEM)                    %
%     iX = fixedDof_XFEM(nX);                         %
%     KX(:,iX) = 0;                               %
%     KX(iX,:) = 0;                               %
%     KX(iX,iX) = bcwt;                           %
% end                                             %

% apply reverse Cuthill-Mckee reordering
p_rcm = symrcm(K);   p_rcmX = symrcm(KX);                 %
rcm_K = K(p_rcm,p_rcm);   rcm_KX = KX(p_rcmX, p_rcmX);    %
fprintf(logFID,'LU decomposition of K ...\n');
fprintf('LU decomp \n')
tic
[Lmatrix,Umatrix] = lu(rcm_K);  [LmatrixX,UmatrixX] = lu(rcm_KX);    %
t = toc
fprintf(logFID,'LU decomposition time = %f\n',t);

%% stiffness matrix for each grain
tic;
kg_grain = [];
for ng = 1:ngr
    fprintf(logFID,'Define global stiffness matrix kg_grain of grain %d ...\n',ng);
    kg_temp = sparse(mno*2,mno*2);
    for iele = 1:mel % = 1:number of elements
        if elegrainindex(iele)==ng
            for a = 1:3
                for b = 1:3
                    gna = nc(iele,a);
                    gnb = nc(iele,b);
                    for i = 1:2
                        for j = 1:2
                            kg_temp(2*(gna-1)+i,2*(gnb-1)+j) = kg_temp(2*(gna-1)+i,2*(gnb-1)+j)+ke(2*(a-1)+i,2*(b-1)+j,iele);
                            %
                        end
                    end
                end
            end
        end
    end
    kg_grain{ng,1} = kg_temp;
end


bcwt_grain = zeros(ngr,1);
gammaMixed1_grain = cell(ngr,1);
gammaMixed2_grain = cell(ngr,1);
gammau_grain = cell(ngr,1);
gammat_grain = cell(ngr,1);
fixedDofs_grain = cell(ngr,1);
freeDofs_grain = cell(ngr,1);
allDofs_grain = cell(ngr,1);
Lmatrix_grain = cell(ngr,1);
Umatrix_grain = cell(ngr,1);

for ng = 1:ngr
    fprintf(logFID,'reformatting K_grain of grain %d\n',ng);
    K_grain = kg_grain{ng,1};
    
    allDofs_grain{ng,1} = sort([2*nodeingrain{ng,1}-1;2*nodeingrain{ng,1}]);
    
    bcwt_temp = mean(diag(kg_grain{ng}(allDofs_grain{ng,1},allDofs_grain{ng,1}))); % ZEBANG HAS CHANGED HERE! COMPARE!Z: mean(diag(kg_grain{ng}))
    bcwt_temp = full(bcwt_temp); % = trace(K)/length(K)
    bcwt_grain(ng,1) = bcwt_temp;
    
    gammaMixed1_grain{ng,1} = nodeonsample{ng,1}(ismember(nodeonsample{ng,1},gammaMixed1));
    gammaMixed2_grain{ng,1} = nodeonsample{ng,1}(ismember(nodeonsample{ng,1},gammaMixed2));
    gammat_grain{ng,1} = nodeonsample{ng,1}(ismember(nodeonsample{ng,1},gammat));
    gammau_grain{ng,1} = nodeonsample{ng,1}(ismember(nodeonsample{ng,1},gammau));
    gammau_grain{ng,1} = [gammau_grain{ng,1};nodeoninterface{ng,1}];
    
    fixedDofs_grain{ng,1} = [2*gammaMixed1_grain{ng,1}; 2*gammaMixed2_grain{ng,1}-1; 2*gammau_grain{ng,1}-1; 2*gammau_grain{ng,1}];
    freeDofs_grain{ng,1} =  [2*gammaMixed1_grain{ng,1}-1; 2*gammaMixed2_grain{ng,1}; 2*gammat_grain{ng,1}-1; 2*gammat_grain{ng,1}];
    
    if length([fixedDofs_grain{ng,1};freeDofs_grain{ng,1}])>...
            length(unique([fixedDofs_grain{ng,1};freeDofs_grain{ng,1}]))
        error(['DOF error in boundary conditions of grain ', num2str(ng)]);
    end
    
    for n = 1:length(fixedDofs_grain{ng,1})
        i = fixedDofs_grain{ng,1}(n);
        K_grain(:,i) = 0;
        K_grain(i,:) = 0;
        K_grain(i,i) = bcwt_grain(ng);
    end
    Ktemp = K_grain(allDofs_grain{ng,1},allDofs_grain{ng,1}); % so size is not (2*mno, 2*mno) henceforth !
    
    fprintf(logFID,'LU decomposition of K_grain of grain %d\n',ng);
    tic
    [Ltemp,Utemp] = lu(Ktemp);
    Lmatrix_grain{ng,1} = Ltemp;
    Umatrix_grain{ng,1} = Utemp;
    t = toc;
    fprintf(logFID,'LU decomposition time of grain %d is %f\n',ng, t);
end
fprintf('LU decomp for each grain \n')
toc


%% put all interface Dofs into one vector
% (if the interface node is u-controlled, then it is removed from this vector)
allinterfacenode = []; % allinterfacenode = cell2mat(nodeoninterface{:});
for ng = 1:ngr
    allinterfacenode = [allinterfacenode; nodeoninterface{ng}];
end
allinterfacenode = unique(allinterfacenode);
interfaceDofs = [2*allinterfacenode-1; 2*allinterfacenode];
interfaceDofs(ismember(interfaceDofs,fixedDofs)) = [];

%% remove the nodal Dofs from the global freeDofs if it belongs to an interface node
freeDofs = unique([freeDofs;interfaceDofs]);
if length([fixedDofs;freeDofs])>length(unique([fixedDofs;freeDofs]))
    error('DOF error in boundary conditions');
end

fprintf(logFID,'Finite Element define all done! \n\n');



