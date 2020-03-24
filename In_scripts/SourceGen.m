% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Function    : SourceGen
%   Last edited  : 1 November, 2018 - SW
%   Description : call by Input.m
%                 generate slip planes, FR sources and obstacles
%                   store information relevant for sliptransmission
%   Outstanding issues :
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xSource,ySource,planeSource,alphaSource,rSource,tauSource,LG,sourcecount,tnuc,...
    xObs,yObs,planeObs,alphaObs,rObs,obscount,nslp,p_rend,p_xend,p_yend,p_connect]...
    = SourceGen(bG,polynode,holes,bbox,ngr,nSystems,sourceDensity,obsDensity,edgeType,lambda,tau_FR,sdev_FR,...
    tnuc0,directoryOut,samplenode,logFID,obsRand,ObsSpac)

testsourcegen = 1;
singlesource  = 0;

% if exist(oldSource)
%     fprintf(logFID,'Load old source/obstacle structure from :\n');
%     fprintf(logFID,[oldSource,'\n']);
%     load(oldSource,...
%     'xSource','ySource','planeSource','alphaSource','rSource','tauSource',...
%     'LG','sourcecount','tnuc','xObs','yObs','planeObs','alphaObs','rObs',...
%     'obscount','nslp','p_rend','p_xend','p_yend','p_connect','bG',...
%     'polynode','bbox','ngr','nSystems','sourceDensity','obsDensity','edgeType',...
%     'lambda','tau_FR','sdev_FR','tnuc0','samplenode','oldSource','obsRand','ObsSpac');
% else
    
ObsDummy = cell(ngr,3);
%% generate slip planes
fprintf(logFID,'Total grain number: %d \n',...
                                ngr);
tic;
xlimit = bbox(:,1)';
ylimit = bbox(:,2)';

for ng = 1:ngr
    d = 200*bG(ng); % distance between slip planes i.e. 200 for Nicola
    fprintf(logFID,'Define slip planes in Grain %d \n',...
                                    ng);
    Gbn = polynode{ng,1};
    
    holeind=find([holes{:,2}]'==ng);
    
    for isys = 1:nSystems(ng)
        j = 0; % use to count planes
        if lambda(isys,ng)==0
            x = xlimit;
            y = ylimit(1)*[1,1]+d;
            while y(1)<ylimit(2)
                [xi,yi] = polyxpoly(x,y,Gbn(:,1),Gbn(:,2));
                if ~isempty(xi)
                    
                    for i_h=1:length(holeind)
                        Hn = holes{holeind(i_h),1}; % hole nodes
                        [x_ih, y_ih] = polyxpoly(x,y,Hn(:,1),Hn(:,2));
                        xi = [xi; x_ih]; yi=[yi; y_ih]; % add pts of intersection with hole
                    end
                    
                    [xi, indexxi] = sort(xi); % sort by x
                    yi = yi(indexxi);
                    xq = (xi(2:end)-xi(1:end-1))/2+xi(1:end-1); % average value, excluding 1st+last points average
                    yq = (yi(2:end)-yi(1:end-1))/2+yi(1:end-1);
                    [in1,on1] = inpolygon(xq,yq,Gbn(:,1),Gbn(:,2)); % in grain
                    q = find(in1==1);
                    
                    out_t = [];
                    for i_h=1:length(holeind)
                        Hn = holes{holeind(i_h),1}; % hole nodes
                        [out1,on2] = inpolygon(xq,yq,Hn(:,1),Hn(:,2)); % In HOLE
                        q_o = find(out1==1); 
                        out_t = [out_t; q_o];                         
                    end
                    
                    for k = 1:sum(in1) % =length(q)
                        if ~any(out_t==q(k))
                            j = j+1;
                            qk_next = q(k)+1;
                            ObsDummy{ng,isys}(j,1:6) = [xi(q(k)),yi(q(k)),xi(qk_next),yi(qk_next),sqrt((xi(qk_next)-xi(q(k)))^2+(yi(qk_next)-yi(q(k)))^2),0];
                        end
                    end
                end
                y = y+d;
            end
        elseif lambda(isys,ng)==pi/2
            y = ylimit;
            x = xlimit(1)*[1,1]+d;
            while x(1)<xlimit(2)
                [xi,yi] = polyxpoly(x,y,Gbn(:,1),Gbn(:,2));
                if ~isempty(xi)
                    
                    for i_h=1:length(holeind)
                        Hn = holes{holeind(i_h),1}; % hole nodes
                        [x_ih, y_ih] = polyxpoly(x,y,Hn(:,1),Hn(:,2));
                        xi = [xi; x_ih]; yi=[yi; y_ih]; % add pts of intersection with hole
                    end

                    [yi, indexyi] = sort(yi);
                    xi = xi(indexyi);
                    xq = (xi(2:end)-xi(1:end-1))/2+xi(1:end-1); % average value, excluding 1st point
                    yq = (yi(2:end)-yi(1:end-1))/2+yi(1:end-1);
                    [in1,on1] = inpolygon(xq,yq,Gbn(:,1),Gbn(:,2));
                    q = find(in1==1);
                    
                    out_t = [];
                        for i_h=1:length(holeind)
                            Hn = holes{holeind(i_h),1}; % hole nodes
                            [out1,on2] = inpolygon(xq,yq,Hn(:,1),Hn(:,2)); % In HOLE
                            q_o = find(out1==1);
                            out_t = [out_t; q_o];
                        end
                    
                    for k = 1:sum(in1) % =length(q)
                        if ~any(out_t==q(k))
                            j = j+1;
                            qk_next = q(k)+1;
                            ObsDummy{ng,isys}(j,1:6) = [xi(q(k)),yi(q(k)),xi(qk_next),yi(qk_next),sqrt((xi(qk_next)-xi(q(k)))^2+(yi(qk_next)-yi(q(k)))^2),0];
                        end
                        
                    end
                end
                x = x+d;
            end
        elseif lambda(isys,ng)<pi/2
            y = ylimit;
            x = [xlimit(1)-(ylimit(2)-ylimit(1))/tan(lambda(isys,ng)), xlimit(1)];
            while x(1)<xlimit(2)
                [xi,yi] = polyxpoly(x,y,Gbn(:,1),Gbn(:,2));
                if ~isempty(xi)
                    
                        for i_h=1:length(holeind)
                            Hn = holes{holeind(i_h),1}; % hole nodes
                            [x_ih, y_ih] = polyxpoly(x,y,Hn(:,1),Hn(:,2));
                            xi = [xi; x_ih]; yi=[yi; y_ih]; % add pts of intersection with hole
                        end
                    
                    [yi, indexyi] = sort(yi);
                    xi = xi(indexyi);
                    xq = (xi(2:end)-xi(1:end-1))/2+xi(1:end-1); % average value, excluding 1st point
                    yq = (yi(2:end)-yi(1:end-1))/2+yi(1:end-1);
                    [in1,on1] = inpolygon(xq,yq,Gbn(:,1),Gbn(:,2));
                    q = find(in1==1);
                    
                    
                    out_t = [];
                        for i_h=1:length(holeind)
                            Hn = holes{holeind(i_h),1}; % hole nodes
                            [out1,on2] = inpolygon(xq,yq,Hn(:,1),Hn(:,2)); % In HOLE
                            q_o = find(out1==1);
                            out_t = [out_t; q_o];
                        end
                    
                    for k = 1:sum(in1) % =length(q)
                        if ~any(out_t==q(k))
                            j = j+1;
                            qk_next = q(k)+1;
                            ObsDummy{ng,isys}(j,1:6) = [xi(q(k)),yi(q(k)),xi(qk_next),yi(qk_next),sqrt((xi(qk_next)-xi(q(k)))^2+(yi(qk_next)-yi(q(k)))^2),0];
                        end
                        
                    end
                end
                x = x+d/sin(lambda(isys,ng));
            end
        elseif lambda(isys,ng)>pi/2
            y = ylimit;
            x = [xlimit(1), xlimit(1)-(ylimit(2)-ylimit(1))*tan(lambda(isys,ng)-pi/2)];
            while x(2)<xlimit(2)
                [xi,yi] = polyxpoly(x,y,Gbn(:,1),Gbn(:,2));
                if ~isempty(xi)
                       
                        for i_h=1:length(holeind)
                            Hn = holes{holeind(i_h),1}; % hole nodes
                            [x_ih, y_ih] = polyxpoly(x,y,Hn(:,1),Hn(:,2));
                            xi = [xi; x_ih]; yi=[yi; y_ih]; % add pts of intersection with hole
                        end

                    [yi, indexyi] = sort(yi);
                    xi = xi(indexyi);
                    xq = (xi(2:end)-xi(1:end-1))/2+xi(1:end-1); % average value, excluding 1st point
                    yq = (yi(2:end)-yi(1:end-1))/2+yi(1:end-1);
                    [in1,on1] = inpolygon(xq,yq,Gbn(:,1),Gbn(:,2));
                    q = find(in1==1);
                    
                    out_t = [];
                        for i_h=1:length(holeind)
                            Hn = holes{holeind(i_h),1}; % hole nodes
                            [out1,on2] = inpolygon(xq,yq,Hn(:,1),Hn(:,2)); % In HOLE
                            q_o = find(out1==1);
                            out_t = [out_t; q_o];
                        end
                    
                    for k = 1:sum(in1) % =length(q)
                        if ~any(out_t==q(k))
                            j = j+1;
                            qk_next = q(k)+1;
                            ObsDummy{ng,isys}(j,1:6) = [xi(q(k)),yi(q(k)),xi(qk_next),yi(qk_next),sqrt((xi(qk_next)-xi(q(k)))^2+(yi(qk_next)-yi(q(k)))^2),0];
                        end
                        
                    end
                end
                x = x+d/cos(lambda(isys,ng)-pi/2);
            end    
        end
    end
end
toc
fprintf('Slip plane generation done. \n');
fprintf(logFID,'Slip plane generation done. \n\n');
                            
p_rend = [];
p_xend = [];
p_yend = [];
nslp = [];

j = 0;
for ng = 1:ngr
    for isys = 1:nSystems(ng)
        for islp = 1:size(ObsDummy{ng,isys},1)
            j = j+1;
            p_rend(j,1:2) = [0,ObsDummy{ng,isys}(islp,5)];
            p_xend(j,1:2) = [ObsDummy{ng,isys}(islp,1),ObsDummy{ng,isys}(islp,3)];
            p_yend(j,1:2) = [ObsDummy{ng,isys}(islp,2),ObsDummy{ng,isys}(islp,4)];
        end
        nslp(isys,ng) = size(ObsDummy{ng,isys},1);
    end
end

figure(1);clf;hold on
box on;
set(gcf,'Color','w')
set(gca,'FontSize',20,'LineWidth',4)
set(gcf,'units','normalized','outerposition',[0.2 0.1 0.5 0.7])
axis equal;
for ng = 1:ngr
    xnode = polynode{ng,1};
    plot(xnode(:,1),xnode(:,2),'k-','LineWidth',4);
    
    for isys = 1:nSystems(ng)
        for j = 1:size(ObsDummy{ng,isys},1)
            plot(ObsDummy{ng,isys}(j,[1,3]),ObsDummy{ng,isys}(j,[2,4]),'-','LineWidth',2);
            if ObsDummy{ng,isys}(j,6)==1
                plot(ObsDummy{ng,isys}(j,[1,3]),ObsDummy{ng,isys}(j,[2,4]),'ko-','LineWidth',2);
            end
%             text(mean(ObsDummy{ng,isys}(j,[1,3])),mean(ObsDummy{ng,isys}(j,[2,4])),num2str(j))
%             pause(0.1);
        end
%         pause(1)
    end
end

    for ng = 1:size(holes,1)
        xnode = holes{ng,1};
        plot(xnode(:,1),xnode(:,2),'k-','LineWidth',4);
    end

axis equal; axis off;
ax=axis;axis(ax*1.001);
hold off
saveas(figure(1),['../',directoryOut,'SlipPlaneStructure.fig'])

maxplanes = 1e5;
xold = zeros(2,maxplanes,max(nSystems),ngr);
xpts = zeros(2,maxplanes,max(nSystems),ngr);
prevplane1 = 0;
prevplane2 = 0;
sourcecount = zeros(ngr,1);
obscount = zeros(ngr,1);
totallen = 0;

if (~isempty(obsRand) && ~isempty(ObsSpac))
    critical_length = ObsSpac+50*max(bG);
else
    critical_length = 750*max(bG);
end
    
for ng = 1:ngr
    for isys = 1:nSystems(ng)
        for islp = 1:size(ObsDummy{ng,isys},1)
            if ObsDummy{ng,isys}(islp,5)>=critical_length
                totallen = totallen+ObsDummy{ng,isys}(islp,5);
            end
        end
    end
end
totalsource  = 1.05*sourceDensity*polyarea(samplenode(1:end-1,1),samplenode(1:end-1,2)); % *constant added because later on sources close boundaries will be removed
avelen1 = totallen/totalsource;

%% generate sources
tic;
nran = 0;
for ng = 1:ngr
    fprintf(logFID,'Distribute FR sources in Grain %d \n',...
                                    ng);
    xnum = rand(100000,1);
    for isys = 1:nSystems(ng)
        for islp = 1:size(ObsDummy{ng,isys},1)
            
            if ObsDummy{ng,isys}(islp,5)>=critical_length
                nran = nran+1;
                scprob = ObsDummy{ng,isys}(islp,5)/avelen1;
                iprob = floor(scprob); % rounds a number to the next smaller integer
                fprob = scprob-iprob;
                n = iprob;

                if fprob>xnum(nran); n=n+1; end
            else
                n = 0;
            end
            
            if singlesource == 1
                if isys==1 && islp == round(0.5*size(ObsDummy{ng,isys},1)); n=1;
                else; n=0; end
            end
            
            if n~=0
                tau = normrnd(tau_FR(ng),sdev_FR(ng),n,1);
                if singlesource == 1; tau = tau_FR(ng); end % single source 
                L = max(edgeType(ng)*bG(ng)./tau);

                jcond = 1;
                i = 1;
                while i<=n
                    nold = find(xold(:,islp,isys,ng)~=0);
                    xpts(i,islp,isys,ng) = critical_length/2+(ObsDummy{ng,isys}(islp,5)-critical_length)*rand(1);
                    if singlesource == 1; xpts(i,islp,isys,ng) = 0.5*(ObsDummy{ng,isys}(islp,5)); end % single source  
                    if isempty(nold)
                        xold(1,islp,isys,ng) = xpts(i,islp,isys,ng);
                        i = i+1;
                    else
                        if min(abs(xpts(i,islp,isys,ng)-xold(:,islp,isys,ng)))<=L
                            i = i;
                            jcond = jcond+1;
                        else
                            xold(length(nold)+1,islp,isys,ng) = xpts(i,islp,isys,ng);
                            i = i+1;
                        end
                    end
                    if jcond>1000
                        fprintf(logFID,'two points too close, reduce no. of obstacles\n ');
                        fprintf(logFID,'Grain %d, slip system %d, slip plane index %d \n',ng,isys,islp);
                        i=i+1; 
                    end
                end

                if lambda(isys,ng)<=pi/2
                    xSource(sourcecount(ng)+1:sourcecount(ng)+n,ng) = ObsDummy{ng,isys}(islp,1)+xold(1:n,islp,isys,ng)*cos(lambda(isys,ng));
                    ySource(sourcecount(ng)+1:sourcecount(ng)+n,ng) = ObsDummy{ng,isys}(islp,2)+xold(1:n,islp,isys,ng)*sin(lambda(isys,ng));
                elseif lambda(isys,ng)>pi/2
                    xSource(sourcecount(ng)+1:sourcecount(ng)+n,ng) = ObsDummy{ng,isys}(islp,1)-xold(1:n,islp,isys,ng)*cos(pi-lambda(isys,ng));
                    ySource(sourcecount(ng)+1:sourcecount(ng)+n,ng) = ObsDummy{ng,isys}(islp,2)+xold(1:n,islp,isys,ng)*sin(pi-lambda(isys,ng));
                end
                tauSource(sourcecount(ng)+1:sourcecount(ng)+n,ng) = tau;
                LG(sourcecount(ng)+1:sourcecount(ng)+n,ng) = edgeType(ng)*bG(ng)./tauSource(sourcecount(ng)+1:sourcecount(ng)+n,ng)/2;
                planeSource(sourcecount(ng)+1:sourcecount(ng)+n,ng) = prevplane1+islp;
                alphaSource(sourcecount(ng)+1:sourcecount(ng)+n,ng) = lambda(isys,ng);
                rSource(sourcecount(ng)+1:sourcecount(ng)+n,ng) = xold(1:n,islp,isys,ng);
                sourcecount(ng) = sourcecount(ng)+n;
            end
        end
        prevplane1 = prevplane1+size(ObsDummy{ng,isys},1);
    end
end
tnuc = tnuc0/mean(tau_FR)*tauSource;
toc
fprintf('FR sources distribution done. \n');
fprintf(logFID,'FR sources distribution done. \n\n');

%% generate obstacles
tic;
if (~isempty(obsRand) && ~isempty(ObsSpac)) % double-ended obstacles
    fprintf(logFID,'Distribute bouble-ended obstacles. \n');
    obscount = sourcecount*2;
    
    for ng = 1:ngr
        fprintf(logFID,'Distribute obstacles in Grain %d \n',...
                                        ng);
        for isys = 1:nSystems(ng)
            a = find(alphaSource(1:sourcecount(ng),ng) == lambda(isys,ng));
            if lambda(isys,ng)<=pi/2
                xObs(a,ng) = xSource(a,ng)-ObsSpac/2*cos(lambda(isys,ng));
                yObs(a,ng) = ySource(a,ng)-ObsSpac/2*sin(lambda(isys,ng));
                xObs(sourcecount(ng)+a,ng) = xSource(a,ng)+ObsSpac/2*cos(lambda(isys,ng));
                yObs(sourcecount(ng)+a,ng) = ySource(a,ng)+ObsSpac/2*sin(lambda(isys,ng));
            elseif lambda(isys,ng)>pi/2
                xObs(a,ng) = xSource(a,ng)+ObsSpac/2*cos(pi-lambda(isys,ng));
                yObs(a,ng) = ySource(a,ng)-ObsSpac/2*sin(pi-lambda(isys,ng));
                xObs(sourcecount(ng)+a,ng) = xSource(a,ng)-ObsSpac/2*cos(pi-lambda(isys,ng));
                yObs(sourcecount(ng)+a,ng) = ySource(a,ng)+ObsSpac/2*sin(pi-lambda(isys,ng));
            end
            planeObs(a,ng) = planeSource(a,ng);
            planeObs(sourcecount(ng)+a,ng) = planeSource(a,ng);
            alphaObs(a,ng) = alphaSource(a,ng);
            alphaObs(sourcecount(ng)+a,ng) = alphaSource(a,ng);
            rObs(a,ng) = rSource(a,ng)-ObsSpac/2;
            rObs(sourcecount(ng)+a,ng) = rSource(a,ng)+ObsSpac/2;
        end
    end
else
    fprintf(logFID,'Distribute random obstacles. \n');
    nran = 0;
    totalobs = obsDensity*polyarea(samplenode(1:end-1,1),samplenode(1:end-1,2));
    avelen2 = totallen/totalobs;
    for ng = 1:ngr
        fprintf(logFID,'Distribute obstacles in Grain %d \n',...
                                        ng);
        xnum = rand(100000,1);
        for isys = 1:nSystems(ng)
            for islp = 1:size(ObsDummy{ng,isys},1)
                % determine actual number of points on current slip plane each plane has [iprob] 
                % points with a probability  [fprob] of having another one
                if ObsDummy{ng,isys}(islp,5)>=critical_length
                    nran = nran+1;
                    scprob = ObsDummy{ng,isys}(islp,5)/avelen2;
                    iprob = floor(scprob); % rounds a number to the next smaller integer
                    fprob = scprob-iprob;
                    n = iprob;

                    if fprob>xnum(nran); n=n+1; end
                else
                    n = 0;
                end

                if n~=0
                    jcond = 1; 
                    i = 1;
                    while i<=n
                        nold = find(xold(:,islp,isys,ng)~=0);
                        xpts(i,islp,isys,ng) = critical_length/2+(ObsDummy{ng,isys}(islp,5)-critical_length)*rand(1);
                        if isempty(nold)
                            xold(1,islp,isys,ng) = xpts(i,islp,isys,ng);
                            i = i+1;
                        else
                            if min(abs(xpts(i,islp,isys,ng)-xold(:,islp,isys,ng)))<=L/2
                                i = i;
                                jcond = jcond+1;
                            else
                                xold(length(nold)+1,islp,isys,ng) = xpts(i,islp,isys,ng);
                                i = i+1;
                            end
                        end
                        if jcond>1000; 
                            fprintf(logFID,'two points too close, reduce no. of obstacles\n ');
                            fprintf(logFID,'Grain %d, slip system %d, slip plane index %d \n',ng,isys,islp);
                            i=i+1; 
                        end
                    end
                    nold = find(xold(:,islp,isys,ng)~=0);
                    if lambda(isys,ng)<=pi/2
                        xObs(obscount(ng)+1:obscount(ng)+n,ng) = ObsDummy{ng,isys}(islp,1)+xold(length(nold)-n+1:length(nold),islp,isys,ng)*cos(lambda(isys,ng));
                        yObs(obscount(ng)+1:obscount(ng)+n,ng) = ObsDummy{ng,isys}(islp,2)+xold(length(nold)-n+1:length(nold),islp,isys,ng)*sin(lambda(isys,ng));
                    elseif lambda(isys,ng)>pi/2
                        xObs(obscount(ng)+1:obscount(ng)+n,ng) = ObsDummy{ng,isys}(islp,1)-xold(length(nold)-n+1:length(nold),islp,isys,ng)*cos(pi-lambda(isys,ng));
                        yObs(obscount(ng)+1:obscount(ng)+n,ng) = ObsDummy{ng,isys}(islp,2)+xold(length(nold)-n+1:length(nold),islp,isys,ng)*sin(pi-lambda(isys,ng));
                    end
                    planeObs(obscount(ng)+1:obscount(ng)+n,ng) = prevplane2+islp;
                    alphaObs(obscount(ng)+1:obscount(ng)+n,ng) = lambda(isys,ng);
                    rObs(obscount(ng)+1:obscount(ng)+n,ng) = xold(length(nold)-n+1:length(nold),islp,isys,ng);
                    obscount(ng) = obscount(ng)+n;
                end
            end
            prevplane2 = prevplane2+size(ObsDummy{ng,isys},1);
        end
    end
end
toc;
fprintf('obstacles distribution done. \n'); 
fprintf(logFID,'obstacles distribution done. \n\n'); 

%% delete sources too close to boundary
figure(2);clf;hold on
box on;
set(gcf,'Color','w')
set(gca,'FontSize',20,'LineWidth',4)
set(gcf,'units','normalized','outerposition',[0.2 0.1 0.5 0.7])
axis equal;
jslp = 0;
for ng = 1:ngr
    Gbn = polynode{ng};
    plot(Gbn(:,1),Gbn(:,2),'k-','LineWidth',4);
    plot(xSource(1:sourcecount(ng),ng),ySource(1:sourcecount(ng),ng),'b.','MarkerSize',15)
    
end

    for ng = 1:size(holes,1)
        xnode = holes{ng,1};
        plot(xnode(:,1),xnode(:,2),'k-','LineWidth',4);
    end
    axis equal; axis off;
    ax=axis;axis(ax*1.001);

if singlesource == 0
for ng=1:ngr
    xnodes = polynode{ng,1};
    fd=dpolyDD([xSource(1:sourcecount(ng),ng),ySource(1:sourcecount(ng),ng)],xnodes); % excluding holes for now
        holeind=find([holes{:,2}]'==ng);
        for i_h=1:length(holeind)
            Hn = holes{holeind(i_h),1}; % hole nodes
            xnodes=[xnodes; NaN NaN; Hn];
            fd = ddiffDD(fd,dpolyDD([xSource(1:sourcecount(ng),ng),ySource(1:sourcecount(ng),ng)],Hn));
        end
    scatter3(xSource(1:sourcecount(ng),ng),ySource(1:sourcecount(ng),ng),fd);
    %     pause
    %   %  remove small value numerical error caused by distmesh2d.m
    keep_source = find(fd<-0.12);
    toaddvec = length(planeSource)-length(keep_source);
    planeSource(:,ng) = [planeSource(keep_source,ng); zeros(toaddvec,1)];
    rSource(:,ng) = [rSource(keep_source,ng); zeros(toaddvec,1)];
    xSource(:,ng) = [xSource(keep_source,ng); zeros(toaddvec,1)];
    ySource(:,ng) = [ySource(keep_source,ng); zeros(toaddvec,1)];
    alphaSource(:,ng) = [alphaSource(keep_source,ng); zeros(toaddvec,1)];
    tauSource(:,ng) = [tauSource(keep_source,ng); zeros(toaddvec,1)];
    LG(:,ng) = [LG(keep_source,ng); zeros(toaddvec,1)];
    tnuc(:,ng) = [tnuc(keep_source,ng); zeros(toaddvec,1)];
    sourcecount(ng) = length(keep_source);
    plot(xSource(1:sourcecount(ng),ng),ySource(1:sourcecount(ng),ng),'r.','MarkerSize',15)
end
end
saveas(figure(2),['../',directoryOut,'SourceStructure.fig']);

%% Put Sources and Obstacles in right order based on their position
totalplane = 0;
for ng = 1:ngr
    for isys = 1:nSystems(ng)
        totalplane = totalplane+size(ObsDummy{ng,isys},1);
    end
end
for i = 1:totalplane
    n = [];
    ngp = [];
    [n,ngp] = find(planeSource==i);
    ng = unique(ngp);
    [r,j] = sort(rSource(n,ng));
    rSource(n,ng) = r;
    x_temp = []; 
    x_temp = xSource(n,ng);
    xSource(n,ng) = x_temp(j);
    x_temp = []; 
    x_temp = ySource(n,ng);
    ySource(n,ng) = x_temp(j);
    x_temp = []; 
    x_temp = alphaSource(n,ng);
    alphaSource(n,ng) = x_temp(j);
    x_temp = []; 
    x_temp = tauSource(n,ng);
    tauSource(n,ng) = x_temp(j);
    x_temp = []; 
    x_temp = LG(n,ng);
    LG(n,ng) = x_temp(j);
    x_temp = []; 
    x_temp = tnuc(n,ng);
    tnuc(n,ng) = x_temp(j);
    
    n = [];
    ngp = [];
    if sum(obscount(ng))~=0
        [n,ngp] = find(planeObs==i);
        ng = unique(ngp);
        [r,j] = sort(rObs(n,ng));
        rObs(n,ng) = r;
        x_temp = [];
        x_temp = xObs(n,ng);
        xObs(n,ng) = x_temp(j);
        x_temp = [];
        x_temp = yObs(n,ng);
        yObs(n,ng) = x_temp(j);
        x_temp = [];
        x_temp = alphaObs(n,ng);
        alphaObs(n,ng) = x_temp(j);
    end
end
if sum(obscount)==0; rObs = 0; xObs = 0; yObs = 0; alphaObs = 0; planeObs = 1; end

%%
% test sourcege-----------------------------------------------------------
% if testsourcegen==1
%     figure(2);clf;hold on
%     box on;
%     set(gcf,'Color','w')
%     set(gca,'FontSize',20,'LineWidth',4)
%     set(gcf,'units','normalized','outerposition',[0.2 0.1 0.5 0.7])
%     axis equal;
%     jslp = 0;
%     for ng = 1:ngr
%         Gbn = polynode{ng};
%         plot(Gbn(:,1),Gbn(:,2),'k-','LineWidth',4);
% 
%         for isys = 1:nSystems(ng)
%             for j = 1:size(ObsDummy{ng,isys},1)
%                 jslp = jslp + 1;
%                 plot(ObsDummy{ng,isys}(j,[1,3]),ObsDummy{ng,isys}(j,[2,4]),'-','LineWidth',2);
% %                 text(mean(ObsDummy{ng,isys}(j,[1,3])),mean(ObsDummy{ng,isys}(j,[2,4])),num2str(jslp));
%                 if ObsDummy{ng,isys}(j,6)==1
%                     plot(ObsDummy{ng,isys}(j,[1,3]),ObsDummy{ng,isys}(j,[2,4]),'ko-','LineWidth',2);
%                 end
%             end
%         end
%         plot(xSource(1:sourcecount(ng),ng),ySource(1:sourcecount(ng),ng),['b','.'],'MarkerSize',15)
%         plot(xObs(1:obscount(ng),ng),yObs(1:obscount(ng),ng),['r','o'],'MarkerSize',5)
% %         pause
%     end
%     
%     if flag==1
%         for ng = 1:size(holes,1)
%             xnode = holes{ng,1};
%             plot(xnode(:,1),xnode(:,2),'k-','LineWidth',4);
%         end
%         axis equal; axis off;
%         ax=axis;axis(ax*1.001);
%         hold off
%     end
%     saveas(figure(2),['../',directoryOut,'SourceStructure.fig'])
% end
%--------------------------------------------------------------------------

%%
tic; fprintf(logFID,'Find connected plane information:\n');
prevplane1 = 0;
j = 0;
for ng = 1:ngr
    check_dis = (50*bG(ng))^2;
    fprintf(logFID,'Find connected plane information in grain %d \n',...
                                        ng);
    for isys = 1:nSystems(ng)
        for islp = 1:nslp(isys,ng)
            j = j+1;
            prevplane2 = 0;
            p_grain1 = [];
            p_alpha1 = [];
            p_plane1 = [];
            p_pend1 = [];
            p_grain2 = [];
            p_alpha2 = [];
            p_plane2 = [];
            p_pend2 = [];
            jL = 0;
            for ngL = 1:ngr
                
                if ngL ~= ng
                    for isysL = 1:nSystems(ngL)
                                                
                        for islpL = 1:nslp(isysL,ngL)
                            jL = jL+1;
                            if ((p_xend(jL,1)-p_xend(j,1))^2+(p_yend(jL,1)-p_yend(j,1))^2)<check_dis
                                p_plane1 = [p_plane1,jL];
                                p_pend1 = [p_pend1,1];
                                p_grain1 = [p_grain1,ngL];
                                p_alpha1 = [p_alpha1,lambda(isysL,ngL)];
                            elseif ((p_xend(jL,2)-p_xend(j,1))^2+(p_yend(jL,2)-p_yend(j,1))^2)<check_dis
                                p_plane1 = [p_plane1,jL];
                                p_pend1 = [p_pend1,2];
                                p_grain1 = [p_grain1,ngL];
                                p_alpha1 = [p_alpha1,lambda(isysL,ngL)];
                            end
                            if ((p_xend(jL,1)-p_xend(j,2))^2+(p_yend(jL,1)-p_yend(j,2))^2)<check_dis
                                p_plane2 = [p_plane2,jL];
                                p_pend2 = [p_pend2,1];
                                p_grain2 = [p_grain2,ngL];
                                p_alpha2 = [p_alpha2,lambda(isysL,ngL)];
                            elseif ((p_xend(jL,2)-p_xend(j,2))^2+(p_yend(jL,2)-p_yend(j,2))^2)<check_dis
                                p_plane2 = [p_plane2,jL];
                                p_pend2 = [p_pend2,2];
                                p_grain2 = [p_grain2,ngL];
                                p_alpha2 = [p_alpha2,lambda(isysL,ngL)];
                            end
                        end
                        
                    end
                else
                    for isysL = 1:nSystems(ngL)
                        jL = jL+nslp(isysL,ngL);
                    end
                end
            end
            p_connect{j,1} = p_plane1;
            p_connect{j,2} = p_grain1;
            p_connect{j,3} = p_alpha1;
            p_connect{j,4} = p_pend1;
            
            p_connect{j,5} = p_plane2;
            p_connect{j,6} = p_grain2;
            p_connect{j,7} = p_alpha2;
            p_connect{j,8} = p_pend2;
        end
    end
end
fprintf('Find connected plane information done  \n');
t = toc
fprintf(logFID,'Find connected plane information done within %f seconds \n',...
                                        t);

% end
fprintf(logFID,'Sources and obstacles generation done.\n\n');
% save([directoryOut,'source.mat'],...
%     'xSource','ySource','planeSource','alphaSource','rSource','tauSource',...
%     'LG','sourcecount','tnuc','xObs','yObs','planeObs','alphaObs','rObs',...
%     'obscount','nslp','p_rend','p_xend','p_yend','p_connect','bG',...
%     'polynode','bbox','ngr','nSystems','sourceDensity','obsDensity','edgeType',...
%     'lambda','tau_FR','sdev_FR','tnuc0','directoryOut','samplenode','oldSource',...
%     'logFID','obsRand','ObsSpac');

close all;


sfin = sum(sourcecount)/((bbox(2,1)-bbox(1,1))*(bbox(2,2)-bbox(1,2)));
fprintf('Source density in the end is : %f um^-2 \n\n', sfin)
end


function d=dpolyDD(p,pv)

%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

np=size(p,1);
nvs=size(pv,1)-1;

ds=dsegment(p,pv);
d=min(ds,[],2);

d=(-1).^(inpolygon(p(:,1),p(:,2),pv(:,1),pv(:,2))).*d;
d(abs(d)<eps)=0;
end

function d=ddiffDD(d1,d2), d=max(d1,-d2);
end