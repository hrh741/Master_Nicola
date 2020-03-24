% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Function    : trial velocity of dislocations 
%   Description : call by Main.m
%                 calculate trial velocity using mobility law
%                 vdis = fPK/B;
%                 fPK = taudis*b;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [vdis,taudis] = vdistrial(ndis,xdis,ydis,plane,type,alpha,ngsource,b,edgeType,...
    ngr,irmbound,uhat,uFE,ncC,p,DC,friction,Bdrag,tau_FR,vcutoff,elegrainindex,nc,mel,dt)

vdis = zeros(ndis,1);
taudis = zeros(ndis,1);
sigma = zeros(2,2);
s_sum = zeros(ndis,4);
sDxy = zeros(ndis,6);

for ng = 1:ngr
    
    m = (ngsource==ng);
    
%     [sD11, sD22, sD12, sD21] = tildaStressMex(xdis(m),ydis(m),int32(type(m)),alpha(m),sum(m),2,...
%         xdis(m),ydis(m),sum(m),edgeType(ng),b(m));
    
    % Using Chakravarthy and Curtin 2010 MSMSE
    [sD11, sD22, sD12, sD21, sD11x, sD11y, sD12x, sD12y, sD22x, sD22y]...
        = tildaStress2Mex(xdis(m),ydis(m),int32(type(m)),alpha(m),sum(m),2,...
        xdis(m),ydis(m),sum(m),edgeType(ng),b(m));

    [belong] = findbelong(sum(m),xdis(m),ydis(m),plane(m),p,nc,elegrainindex,ng,2);
    [sA11, sA22, sA12, sA21] = hatStressMex(uhat(:,ng),int32(ncC),p(:,1),p(:,2),DC,...
        int32(elegrainindex),int32(belong),sum(m),mel);
    
    s_sum(m,1:4) = s_sum(m,1:4)+[sD11, sD22, sD12, sD21]+[sA11, sA22, sA12, sA21];
    sDxy(m,1:6) = [sD11x, sD11y, sD12x, sD12y, sD22x, sD22y];
end
[belong] = findbelong(ndis,xdis,ydis,plane,p,nc,[],[],2);
[sFE11, sFE22, sFE12, sFE21] = hatStressMex(uFE,int32(ncC),p(:,1),p(:,2),DC,...
    int32(elegrainindex),int32(belong),ndis,mel);

s_sum(:,1:4) = s_sum(:,1:4) + [sFE11, sFE22, sFE12, sFE21];

for i = 1:ndis
    sigma = [s_sum(i,[1,3]);s_sum(i,[4,2])]; 
        taudis(i) = [cos(alpha(i)),sin(alpha(i))]*sigma...
            *[cos(alpha(i)+pi/2);sin(alpha(i)+pi/2)];
    if irmbound(i)==0
        
        % Chakravarthy & Curtin (2010) Backward Euler
        sigma_grad = [sDxy(i,1)*cos(alpha(i))+sDxy(i,2)*sin(alpha(i)), sDxy(i,3)*cos(alpha(i))+sDxy(i,4)*sin(alpha(i));... % in direction of slip
            sDxy(i,3)*cos(alpha(i))+sDxy(i,4)*sin(alpha(i)), sDxy(i,5)*cos(alpha(i))+sDxy(i,6)*sin(alpha(i))];
        PK_grad = [cos(alpha(i)),sin(alpha(i))]*sigma_grad...
            *[cos(alpha(i)+pi/2);sin(alpha(i)+pi/2)];
        
        if friction==0
%             vdis(i) = taudis(i)*type(i)*b(i)/Bdrag;
            vdis(i) = taudis(i)*type(i)*b(i)/Bdrag/(1-(PK_grad*type(i)*b(i)*1e-6*dt)/Bdrag);
        else
            if abs(taudis(i))<tau_FR
                vdis(i) = 0;
            else
                vdis(i) = (taudis(i)*type(i)-sign(taudis(i)*type(i))*tau_FR)...
                    *b(i)/Bdrag;
            end
        end
        % cut-off velocity
        if abs(vdis(i))>vcutoff
            vdis(i) = sign(vdis(i))*vcutoff;
        end
    end
end