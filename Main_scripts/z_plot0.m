% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Function    : Plot contour subroutine
%   Description : called after simulations
%                 plot contours of strain field
%                 dislocation structure superposed is also implemented
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SDV = [1:4]; % 1 slip dist, , 2 strain, 3 stress, 4 deformed mesh

switch BCs
    case 1
        applyBC01;
    case 2
        applyBC02;
    otherwise
        error('BCs error in Main.m');
end

%% new grid (finer) for calculation
% % DOES NOT WORK, error in sending value of sum(m2) into findbelong, and then hatstress
% % mx_g = 100; my_g = mx_g; % no of elements in x&y-direction
% % w_g = xcst/mx_g; % elements width
% % h_g = ycst/my_g; % elements height
% % mno_g = (mx_g+1)*(my_g+1);
% % xgrid = zeros(mno_g,2);
% % 
% % for j = 1:mx_g+1 % nodes in x-direction
% %     for i = 1:my_g+1 % nodes in y-direction
% %         gn = j+(i-1)*(mx_g+1); % global node index
% %         xgrid(gn,1) = (j-1)*w_g;
% %         xgrid(gn,2) = ycst-(i-1)*h_g;
% %     end
% % end

% % fd = dpolySW(xgrid,holes{1,1});
% % xgrid = xgrid(fd>=0,:);
% % plot(xgrid(:,1), xgrid(:,2), 'r.');

%% if want to use FE mesh for plotting, uncomment
mno_g = mno;
xgrid = p;


%% calculate displacement and stress component on grid nodes
sigma = zeros(2,2);
sAall = zeros(mno_g,4);
sDall = zeros(mno_g,4);
sFEall = zeros(mno_g,4);
utilda = zeros(2*mno,1);

 % figure; hold on;
 % grid nodes belong to which grain
pbelong = zeros(mno,1);
gridbelong = zeros(mno_g,1);
for ng = 1:ngr
    Npolygon1 = polynode{ng,1}; % outside boundary should be anticlockwise
    if ispolycw(Npolygon1(:,1), Npolygon1(:,2)); Npolygon1 = flipud(Hn); end % check
    holeind=find([holes{:,2}]'==ng);
    for i_h=1:length(holeind)
        Hn = holes{holeind(i_h),1}; % hole nodes were ACW, change to CW using flipud
        if ~ispolycw(Hn(:,1), Hn(:,2)); Hn = flipud(Hn); end
        Npolygon1 = [Npolygon1; NaN NaN; Hn];
    end
    [in,on] = inpolygon(p(:,1),p(:,2),Npolygon1(:,1),Npolygon1(:,2));
    pbelong(in|on) = ng;
    [in,on] = inpolygon(xgrid(:,1),xgrid(:,2),Npolygon1(:,1),Npolygon1(:,2));
    gridbelong(in|on) = ng;
%         plot(Npolygon1(:,1),Npolygon1(:,2),'-');
%         h1 = plot(p(pbelong==ng,1),p(pbelong==ng,2),'m.');
%         % pause;
%         delete(h1);
%         h2 = plot(xgrid(gridbelong==ng,1),xgrid(gridbelong==ng,2),'r.');
%         % pause
%         delete(h2);
end

for ng = 1:ngr
    
    m = (ngsource==ng);
    mp = (pbelong==ng);
    m2 = (gridbelong==ng);
    
    
    [sD11, sD22, sD12, sD21] = tildaStressMex(xdis(m),ydis(m),int32(type(m)),alpha(m),sum(m),1,...
        xgrid(m2,1),xgrid(m2,2),sum(m2),edgeType(ng),b(m));
    
    [belong] = findbelong(sum(m2),xgrid(m2,1),xgrid(m2,2),plane,p,nc,elegrainindex,ng,1);
    [sA11, sA22, sA12, sA21] = hatStressMex(uhat(:,ng),int32(ncC),p(:,1),p(:,2),DC,...
        int32(elegrainindex),int32(belong),sum(m2),mel);
    
    sAall(m2,1:4) = [sA11, sA22, sA12, sA21];
    sDall(m2,1:4) = [sD11, sD22, sD12, sD21];
    
    allnodes=[1:mno]';
    u_node = allnodes(mp);
    utilda_ng = displacementMex(int32(u_node),int32(nSystems),p(:,1),p(:,2),lambdaC,...
        uConst(ng),nu(ng),b(ngsource==ng),mno,sum(mp),sum(ngsource==ng),...
        sum(ngsourceOut==ng),xdis(ngsource==ng),ydis(ngsource==ng),alpha(ngsource==ng),...
        int32(type(ngsource==ng)),xdisOut(ngsourceOut==ng),ydisOut(ngsourceOut==ng),...
        alphaOut(ngsourceOut==ng),int32(typeOut(ngsourceOut==ng)),bOut(ngsourceOut==ng));

    utilda([2*u_node; 2*u_node-1]) = 0; % for nodes on boundary to only keep utilda from one grain
    utilda = utilda + utilda_ng;
    
end

[belong] = findbelong(mno,xgrid(:,1),xgrid(:,2),plane,p,nc,[],[],1);
[sFE11, sFE22, sFE12, sFE21] = hatStressMex(uFE,int32(ncC),p(:,1),p(:,2),DC,...
    int32(elegrainindex),int32(belong),mno_g,mel);
sFEall = [sFE11, sFE22, sFE12, sFE21];

% sEn not physically real elastic stress
% sEn11 = stress_invhat(belong,1);
% sEn22 = stress_invhat(belong,2);
% sEn12 = stress_invhat(belong,3);
% sEn21 = stress_invhat(belong,4);
% sEnall = [sEn11, sEn22, sEn12, sEn21];

% [strain11, strain22, strain12, strain21] = hatStrainMex(uFE+uhat+utilda+uenr,int32(ncC),p(:,1),p(:,2),DC,...
%     int32(elegrainindex),int32(belong),mno_g,mel);


 close all;
%% plot Slip distribution

if any(SDV==1)
    
%     [strain11, strain22, strain12, strain21] = hatStrainMex(utilda,int32(ncC),p(:,1),p(:,2),DC,...
%         int32(elegrainindex),int32(belong),mno_g,mel);
    [strain11, strain22, strain12, strain21] = hatStrainMex(uFE+uhat+utilda+uenr,int32(ncC),p(:,1),p(:,2),DC,...
    int32(elegrainindex),int32(belong),mno_g,mel);
    
    numpixel = 1000;
    plotcomp = zeros(mno_g,1);
    maxslip = 0.05;
    for gn = 1:mno_g
        ng = gridbelong(gn);
        for isys = 1:nSystems
            strain = [strain11(gn), strain12(gn); strain21(gn), strain22(gn)];
            plotcomp(gn) = plotcomp(gn)+abs([cos(lambda(isys,ng)),sin(lambda(isys,ng))]*strain(:,:)*[cos(lambda(isys,ng)+pi/2);sin(lambda(isys,ng)+pi/2)]);
        end
        if plotcomp(gn) > maxslip; plotcomp(gn) = maxslip; end
    end
    
    F = scatteredInterpolant(xgrid(:,1),xgrid(:,2),plotcomp, 'nearest');
    [X, Y] = meshgrid([min(xgrid(:,1)):(max(xgrid(:,1))-min(xgrid(:,1)))/(numpixel-1):max(xgrid(:,1))], [min(xgrid(:,2)):(max(xgrid(:,2))-min(xgrid(:,2)))/(numpixel-1):max(xgrid(:,2))]);
    Z = F(X, Y);
    
    % remove holes from the plot
    for i = 1:size(X,2)
        plotcomp1=zeros(size(X,1),1);
        for i_h=1:size(holes,1)
            Hn = holes{i_h,1}; % hole nodes
            fd=dpolySW([X(:,i) ,Y(:,i)],Hn);
            plotcomp1 = (fd<0);
            Z(plotcomp1,i) = nan;
        end
    end
    
    Slip_Distr = figure('Visible','on');clf;hold on;
    axis equal;
    title('Slip Distribution');
    set(gca,'Visible','off')
    set(gcf,'Color','w')
    pcolor(X,Y,Z);
    shading interp
    colorbar
    %colormap(jet)
    colormap(flipud(gray))
    caxis([0 maxslip])
    
    for ng = 1:ngr
        Npolygon1 = polynode{ng};
        plot(Npolygon1(:,1),Npolygon1(:,2),'k-','Linewidth',2)
    end
    for ng = 1:size(holes,1)
        xnode = holes{ng,1};
        plot(xnode(:,1),xnode(:,2),'-','LineWidth',1);
    end
    % for idis = 1:ndis
    %     plot(xdis(idis)+[-1,1]*dis_size/2*cos(alpha(idis)),ydis(idis)+[-1,1]*dis_size/2*sin(alpha(idis)),'k-','LineWidth',2)
    %     plot(xdis(idis)+[0,dis_size*0.8*cos(alpha(idis)+pi/2*type(idis))],ydis(idis)+[0,dis_size*0.8*sin(alpha(idis)+pi/2*type(idis))],'k-','LineWidth',2)
    % end
    
    hold on;
    plot(xdis(type==1),ydis(type==1), 'r.', 'Markersize', 6);
    plot(xdis(type==-1),ydis(type==-1), 'b.', 'Markersize', 6);
    plot(xdis(source==0),ydis(source==0), 'c.', 'Markersize', 6); % debris
    plot(xdis(source==-1),ydis(source==-1), 'y.', 'Markersize', 6); % re-emmission
    plot(xdis(source==-2),ydis(source==-2), 'g.', 'Markersize', 6); % transmitted
    
    
    axis equal; axis off;
    ax = axis; axis(ax*1.001);
    hold off
    
    % clearvars -except mno_g gridbelong strain xgrid bbox polynode sAall sDall sFEall uhat uFE utilda ncC belong p DC elegrainindex belong mel lambda
    
end

%% plot Effective Strain

if any(SDV==2)
    
    [strain11, strain22, strain12, strain21] = hatStrainMex(uFE+uhat+utilda,int32(ncC),p(:,1),p(:,2),DC,...
        int32(elegrainindex),int32(belong),mno_g,mel);
    
    numpixel = 1000;
    plotcomp = sqrt(2/3*(strain11.^2+strain12.^2+strain22.^2+strain21.^2));
    % plotcomp = srain22;
    % plotcomp = uhat(1:2:end)+utilda(1:2:end)+uFE(1:2:end);
    
    F = scatteredInterpolant(xgrid(:,1),xgrid(:,2),plotcomp);
    [X, Y] = meshgrid([min(xgrid(:,1)):(max(xgrid(:,1))-min(xgrid(:,1)))/(numpixel-1):max(xgrid(:,1))], [min(xgrid(:,2)):(max(xgrid(:,2))-min(xgrid(:,2)))/(numpixel-1):max(xgrid(:,2))]);
    Z = F(X, Y);
 
    % remove holes from plot
    for i = 1:size(X,2)
         plotcomp1=zeros(size(X,1),1);
         for i_h=1:size(holes,1)
             Hn = holes{i_h,1}; % hole nodes
             fd=dpolySW([X(:,i) ,Y(:,i)],Hn);
             plotcomp1 = (fd<0);
             Z(plotcomp1,i) = nan;
         end
     end
    
    hfig = figure('Visible','on');clf;hold on;
    axis equal;
    title('Effective Strain', 'FontSize', 14);
    set(gca,'Visible','off')
    set(gcf,'Color','w')
    pcolor(X,Y,Z);
    shading interp
    figbar=colorbar;
    title(figbar,'MPa', 'FontSize', 14);
    colormap(jet)
        
    for ng = 1:ngr
        Npolygon1 = polynode{ng};
        plot(Npolygon1(:,1),Npolygon1(:,2),'k-','Linewidth',2)
    end
    for ng = 1:size(holes,1)
        xnode = holes{ng,1};
        plot(xnode(:,1),xnode(:,2),'-','LineWidth',1);
        
    end
    % for idis = 1:ndis
    %     plot(xdis(idis)+[-1,1]*dis_size/2*cos(alpha(idis)),ydis(idis)+[-1,1]*dis_size/2*sin(alpha(idis)),'k-','LineWidth',2)
    %     plot(xdis(idis)+[0,dis_size*0.8*cos(alpha(idis)+pi/2*type(idis))],ydis(idis)+[0,dis_size*0.8*sin(alpha(idis)+pi/2*type(idis))],'k-','LineWidth',2)
    % end
    
    axis equal; axis off;
    ax = axis; axis(ax*1.001);
    hold off
    
end

%% plot Stress

if any(SDV==3)
    numpixel = 1000;
    plotcomp = (sAall(:,1)+sDall(:,1)+sFEall(:,1))*1e6; % sigma xx MPa
   
    F = scatteredInterpolant(xgrid(:,1),xgrid(:,2),plotcomp);
    [X, Y] = meshgrid([min(xgrid(:,1)):(max(xgrid(:,1))-min(xgrid(:,1)))/(numpixel-1):max(xgrid(:,1))], [min(xgrid(:,2)):(max(xgrid(:,2))-min(xgrid(:,2)))/(numpixel-1):max(xgrid(:,2))]);
    Z = F(X, Y);
    for i = 1:size(X,2)
         plotcomp1=zeros(size(X,1),1);
         for i_h=1:size(holes,1)
             Hn = holes{i_h,1}; % hole nodes
             fd=dpolySW([X(:,i) ,Y(:,i)],Hn);
             plotcomp1 = (fd<0);
             Z(plotcomp1,i) = nan;
         end
     end
        
    hfig = figure('Visible','on');clf;hold on;
    axis equal;
    title('\sigma_{xx}', 'FontSize', 14);
    set(gca,'Visible','off')
    set(gcf,'Color','w')
    pcolor(X,Y,Z);
    shading interp
    figbar=colorbar;
    title(figbar,'MPa', 'FontSize', 14);
    colormap(jet)
%     caxis([0 1500])
        
    for ng = 1:ngr
        Npolygon1 = polynode{ng};
        plot(Npolygon1(:,1),Npolygon1(:,2),'k-','Linewidth',2)
    end
    for ng = 1:size(holes,1)
        xnode = holes{ng,1};
        plot(xnode(:,1),xnode(:,2),'-','LineWidth',1);
    end
%     for idis = 1:ndis
%         plot(xdis(idis)+[-1,1]*dis_size/2*cos(alpha(idis)),ydis(idis)+[-1,1]*dis_size/2*sin(alpha(idis)),'k-','LineWidth',2)
%         plot(xdis(idis)+[0,dis_size*0.8*cos(alpha(idis)+pi/2*type(idis))],ydis(idis)+[0,dis_size*0.8*sin(alpha(idis)+pi/2*type(idis))],'k-','LineWidth',2)
%     end
    
    axis equal; axis off;
    ax = axis; axis(ax*1.001);
    hold off
    
end


%% 4 plot deformed mesh

% applyBC01_zplot;

if any(SDV==4)
    % Z = zeros((my_g+1),(mx_g+1));
       u = uFE + uhat + utilda + uenr;
    X = p;
    scale=1;
    for gn = 1:mno_g
        X(gn,1) = p(gn,1)+(scale*u(2*gn-1));
        X(gn,2) = p(gn,2)+(scale*u(2*gn));
    end
    Def_mesh = figure('Visible','on');clf;hold on;
    set(gca,'Visible','off')
    set(gcf,'Color','w')
    title('Deformed mesh');
    axis equal;
    for k = 1:size(nc,1)
        plot(X(nc(k,[1:3,1]),1),X(nc(k,[1:3,1]),2),'r-')
    end
     for ng = 1:ngr
        Npolygon1 = polynode{ng};
        plot(Npolygon1(:,1),Npolygon1(:,2),'k-','Linewidth',2)
    end
    for ng = 1:size(holes,1)
        xnode = holes{ng,1};
        plot(xnode(:,1),xnode(:,2),'-','LineWidth',1);
    end
%     for idis = 1:ndis
%         plot(xdis(idis)+[-1,1]*dis_size/2*cos(alpha(idis)),ydis(idis)+[-1,1]*dis_size/2*sin(alpha(idis)),'k-','LineWidth',2)
%         plot(xdis(idis)+[0,dis_size*0.8*cos(alpha(idis)+pi/2*type(idis))],ydis(idis)+[0,dis_size*0.8*sin(alpha(idis)+pi/2*type(idis))],'k-','LineWidth',2)
%     end
end


%%

function d=dpolySW(p,pv)

%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

np=size(p,1);
nvs=size(pv,1)-1;

ds=dsegment(p,pv);
d=min(ds,[],2);

d=(-1).^(inpolygon(p(:,1),p(:,2),pv(:,1),pv(:,2))).*d;
d(abs(d)<eps)=0;
end

function d=ddiff(d1,d2), d=max(d1,-d2);
end