% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Function    : Tildanodes
%   Last edited : 6 November, 2018 - SW
%   Description : call by ApplyBCXX.m
%                 Calculates nodal force using analytical solution
%
%                    B o--------o
%                     /        /
%         N(s)       /        /
%             \  P  o--------o      Find analytical solution of force on P
%              \   /        /       due to * dislocation.
%               \ /    *   /
%             A  o--------o
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ftilda = tildanodes(nodes,xdis,ydis,type,alpha,bng,xnodes,edgeType,mno,nSystems,lambda)

ftilda = zeros(mno*2,1);
delta = 0.0001*max(bng); % b

for j = 1:length(nodes)
    gn = nodes(j,1); gnA = nodes(j,2); gnB = nodes(j,3); anglextos = nodes(j,4); 
    normalvec = [0; 1];
    RtIR = zeros(2,2);
    
    if anglextos >= pi; anglextos = anglextos-pi; normalvec = [0; -1]; end
    
    Rxy_s = [cos(anglextos) sin(anglextos); -sin(anglextos) cos(anglextos)];
    sB = Rxy_s*xnodes(gnB,:)';
    sP = Rxy_s*xnodes(gn,:)';
    sA = Rxy_s*xnodes(gnA,:)';
    
    if any(abs(sA(2,1) - sP(2,1)) > 1E-5 || abs(sP(2,1) - sB(2,1))> 1E-5) 
        format long;
        % fprintf('Node %i , t value of A %f , t value of P %f /n/n', gn, sA(2,1), sP(2,1));
        sB(2,1) = sP(2,1);
        sA(2,1) = sP(2,1);
        
%         error('ERROR - nodal t values not the same');
    end
    
    for a =1:nSystems % slip system i {1,2,3}
        theta = lambda(a);
        i=find(alpha==theta); % i is dislocation indices on system a
        if size(i)==0; i = zeros(0,1); end
        
        sD = Rxy_s*[xdis(i)'; ydis(i)'];
        sType = type(i);
        b_all = bng(i);
        
        thetas_dis = theta - anglextos; % anticlockwise positive from s to dislocation tilda coord
        Rs_tild = [cos(thetas_dis) sin(thetas_dis); -sin(thetas_dis) cos(thetas_dis)];
        
        if sA(1,1) == sP(1,1); IAP11 = 0; IAP22 = 0; IAP12 = 0;
        else [IAP11, IAP22, IAP12] = anainteg(sD(1,:), sA(1,1), sP(1,1), -(sD(2,:)-sA(2,1)), thetas_dis, sType, b_all, delta);
        end
        if sB(1,1) == sP(1,1); IBP11 = 0; IBP22 = 0; IBP12 = 0;
        else [IBP11, IBP22, IBP12] = anainteg(sD(1,:), sB(1,1), sP(1,1), -(sD(2,:)-sB(2,1)), thetas_dis, sType, b_all, delta);
        end
        IP11 = IAP11 - IBP11; IP22 = IAP22 - IBP22; IP12 = IAP12 - IBP12;
        RtIR = RtIR + Rs_tild'*[IP11 IP12; IP12 IP22]*Rs_tild; % this is sum(R'*I'R) over all dislocations
    end
    force3_ana_s = edgeType.*RtIR*normalvec;
    ftilda(2*gn-1:2*gn) = Rxy_s'*force3_ana_s;
    
    if any(isnan(ftilda(2*gn-1:2*gn)))
        fprintf('Node %i , Node A %i , Node B %i /n/n', gn, gnA, gnB);
        fprintf('sP %f , sA %f , sB %f /n/n', sP(1,1), sA(1,1), sB(1,1));
        disp('Error - ftilda NaN values. To continue press any button...'); pause;
        ftilda(2*gn-1:2*gn) = 0;
    end
end

end