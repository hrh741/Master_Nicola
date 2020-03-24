% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Function    : dislocation nucleation
%   Description : call by Main.m
%                 calculate stress on each sources 
%                 if tau>=tauSource and Ln>L,add two more dislocations
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigma = zeros(2,2);
s_sum = zeros(sourcecountC,4);

for ng = 1:ngr
    
    m = (ngsource==ng);
    
    [sD11, sD22, sD12, sD21] = tildaStressMex(xdis(m),ydis(m),int32(type(m)),alpha(m),sum(m),1,...
        xSource(1:sourcecount(ng),ng),ySource(1:sourcecount(ng),ng),sourcecount(ng),edgeType(ng),b(m));

    [belong] = findbelong(sourcecount(ng),xSource(1:sourcecount(ng),ng),ySource(1:sourcecount(ng),ng),planeSource(1:sourcecount(ng),ng),p,nc,elegrainindex,ng);
    [sA11, sA22, sA12, sA21] = hatStressMex(uhat(:,ng),int32(ncC),p(:,1),p(:,2),DC,...
        int32(elegrainindex),int32(belong),sourcecount(ng),mel);
    
    s_sum([sum(sourcecount(1:ng-1))+1:sum(sourcecount(1:ng))],1:4) = s_sum([sum(sourcecount(1:ng-1))+1:sum(sourcecount(1:ng))],1:4)+[sD11, sD22, sD12, sD21]+[sA11, sA22, sA12, sA21];
end

[belong] = findbelong(sourcecountC,xSourceC,ySourceC,planeSourceC,p,nc,[],[]);
[sFE11, sFE22, sFE12, sFE21] = hatStressMex(uFE,int32(ncC),p(:,1),p(:,2),DC,...
    int32(elegrainindex),int32(belong),sourcecountC,mel);

s_sum(1:sourcecountC,1:4) = s_sum(1:sourcecountC,1:4) + [sFE11, sFE22, sFE12, sFE21];


for ng = 1:ngr
    for isource = 1:sourcecount(ng)
        i = sum(sourcecount(1:(ng-1)))+isource;
        sigma = [s_sum(i,[1,3]);s_sum(i,[4,2])]; 
        tauisource = [cos(alphaSource(isource,ng)),sin(alphaSource(isource,ng))]*sigma*[cos(alphaSource(isource,ng)+pi/2);sin(alphaSource(isource,ng)+pi/2)];
        if abs(tauisource)>=tauSource(isource,ng) % start timing for nucleation
            tn(isource,ng) = tn(isource,ng)+sign(tauisource)*dt;
            Ln(isource,ng) = Ln(isource,ng)+L(isource,ng)/(tnuc(isource,ng)*tauSource(isource,ng)/tauisource)*dt;
        else % restart the nucleation clock if tauSource drops     
            Ln(isource,ng) = 0;
        end
        if (abs(tauisource)>=tauSource(isource,ng)) && (abs(Ln(isource,ng))>=L(isource,ng))
            delta1 = abs(L(isource,ng));
            delta2 = -abs(L(isource,ng));
            pin_temp1 = 0;
            pin_temp2 = 0;
            % Check whether new dislocation overshoot obstacles
            if obscount(ng)~=0
                m = find(planeObs(1:obscount(ng),ng)==planeSource(isource,ng) &...
                    abs(rObs(1:obscount(ng),ng)-rSource(isource,ng))<L(isource,ng));
                if ~isempty(m)
                    [r,j] = sort(rObs(m,ng)-rSource(isource,ng));
                    if any(r>0)
                        delta1 = min(r(r>0));
                        pin_temp1 = 1;
                    end
                    if any(r<0)
                        delta2 = max(r(r<0));
                        pin_temp2 = 1;
                    end
                end
            end
            
            % Check that dislocation not already present at nucleation site
            m = find(plane(1:ndis)==planeSource(isource,ng) &...
                rdis(1:ndis)<=(rSource(isource,ng)+delta1) &...
                rdis(1:ndis)>=(rSource(isource,ng)+delta2));
            if ~isempty(m)
                [r,j] = sort(rdis(m)-rSource(isource,ng));
                if any(r>0)
                    delta1 = min(r(r>0))/2;
                    pin_temp1 = 0;
                end
                if any(r<0)
                    delta2 = max(r(r<0))/2;
                    pin_temp2 = 0;
                end
            end
            
            ndis = ndis + 1;
            b(ndis) = bG(ng);
            xdis(ndis) = xSource(isource,ng)+delta1*cos(alphaSource(isource,ng)); 
            ydis(ndis) = ySource(isource,ng)+delta1*sin(alphaSource(isource,ng));
            rdis(ndis) = rSource(isource,ng)+delta1;
            plane(ndis) = planeSource(isource,ng);
            alpha(ndis) = alphaSource(isource,ng); 
            pinned(ndis) = pin_temp1;
            irmbound(ndis) = 0;
            source(ndis) = i;
            ngsource(ndis) = ng;
            eta_obs(ndis) = 0;
            if tauisource>0 % nucleate -.+ edge dipole
                type(ndis) = 1;
            else % nucleate +.- edge dipole
                type(ndis) = -1;
            end
            vdispre(ndis,1) = delta1/abs(tn(isource,ng));
            
            ndis = ndis + 1;
            b(ndis) = bG(ng);
            xdis(ndis) = xSource(isource,ng)+delta2*cos(alphaSource(isource,ng)); 
            ydis(ndis) = ySource(isource,ng)+delta2*sin(alphaSource(isource,ng));
            rdis(ndis) = rSource(isource,ng)+delta2;
            plane(ndis) = planeSource(isource,ng);
            alpha(ndis) = alphaSource(isource,ng); 
            pinned(ndis) = pin_temp2;
            irmbound(ndis) = 0;
            source(ndis) = i;
            ngsource(ndis) = ng;
            eta_obs(ndis) = 0;
            if tauisource>0  % nucleate -.+ edge dipole, indexed at ndis+1:ndis+2
                type(ndis) = -1;
            else % nucleate +.- edge dipole
                type(ndis) = 1;
            end
            vdispre(ndis,1) = delta2/abs(tn(isource,ng));
            
            Ln(isource,ng) = 0;
            tn(isource,ng) = 0;
        end
    end
end