% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Function    : Main
%   Description : For high strain rates
%                 Without thermal activation
%                 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause('off')

logFID = fopen([directoryOut,'Log.txt'],'at');

close all;
dbstop if error

cd ./Main_scripts/;

% Mesh the sample
MeshGen0; 
close all;

%% Slip Planes and Mesh elements Intersection
IntersectionMatrix;


%%
[mno, mel, D, kg, Lmatrix, Umatrix, bcwt, kg_grain, Lmatrix_grain, Umatrix_grain, ...
    bcwt_grain, gammaMixed1_grain, gammaMixed2_grain, gammau_grain, gammat_grain, ...
    fixedDofs_grain, freeDofs_grain, allDofs_grain, freeDofs, interfaceDofs, B, p_rcm, LmatrixX, UmatrixX, p_rcmX]... 
    = FEdefine(logFID, p, nc, ngr, nu, e, mu, nodeonsample, nodeingrain, nodeoninterface, ...
    gammaMixed1, gammaMixed2, gammau, gammat, fixedDofs, fixedDofsX, freeDofs, elegrainindex );  


%%
MatrixInitial;
MatrixReshape;

if BCs==1 || BCs==2
    t = (0.0009*xcst)/Udot; % already apply initial elastic strain
elseif BCs == 3    
    t = (0.75*totalT); % already apply initial load
end

%% Main loop
tic
% zero point at increment = 1;
for increment = 2:total_inc
    
%     if size(xdis,1)~=1000; error('!'); end
    
dt = dt0;
t = t + dt;
U = Udot*t;
F = Fdot*t;
% tic
switch BCs
    case 1
        applyBC01; 
    case 2
        applyBC02; 
    case 3
        applyBC03; 
    otherwise
        error('BCs error in Main.m');
end

%  if size(xdis,1)~=1000; error('!'); end
[vdis,taudis] = vdistrial(ndis,xdis(1:ndis),ydis(1:ndis),plane(1:ndis),type(1:ndis),alpha(1:ndis),...
    ngsource(1:ndis),b(1:ndis),edgeType,ngr,irmbound(1:ndis),uhat,uFE,ncC,p,DC,...
    friction,Bdrag,tau_FR,vcutoff,elegrainindex,nc,mel,dt);

activeplanes = unique(plane(1:ndis));
activeplanes(activeplanes<0) = [];
n_actplanes = length(activeplanes);

switch mObsJump
    case 1
        % Non-thermal activated process
        [vdis,pinned] = jumpobs1(activeplanes,plane,ngsource,planeObs,irmbound,...
            vdis,dt,rdis,pinned,taudis,tau_obs,rObs,ndis);
    case 2
        % Thermal activated process
        [vdis,pinned,eta_obs] = jumpobs2(activeplanes,plane,ngsource,planeObs,irmbound,...
            vdis,dt,rdis,pinned,taudis,rObs,xalpha,xbeta,eta_obs,ndis,tau_FR,reps);
    otherwise
        error('mObsJump is not specified in Main.m');
end

[vdis] = underelaxation(ndis,activeplanes,vdis,vdispre,type,plane);

[irmbound,vdis] = dealGB(irmbound,ndis,p_rend,vdis,rdis,dt,plane,reps);

[vdis, irmbound, pinned, eta_obs] = reorderingMex(int32(activeplanes), ...
    int32(plane(1:ndis)), int32(pinned(1:ndis)), int32(irmbound(1:ndis)), ndis,...
    rdis(1:ndis), vdis(1:ndis), eta_obs(1:ndis), reps, dt, damping, n_actplanes);

% initially update dislocation position, without checking numerical error
% if ndis>0
    vdispre = vdis;
    rdis(1:ndis) = rdis(1:ndis) + vdis(1:ndis)*dt;
    xdis(1:ndis) = xdis(1:ndis) + vdis(1:ndis)*dt.*cos(alpha(1:ndis));
    ydis(1:ndis) = ydis(1:ndis) + vdis(1:ndis)*dt.*sin(alpha(1:ndis));
% end

%  if size(xdis,1)~=1000; error('!'); end

disescape00;

%  if size(xdis,1)~=1000; error('!'); end

[destroy] = disannih(plane,ndis,rdis,type,Le);
if any(destroy==1)
    [b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,plane,pinned,irmbound,eta_obs,...
        vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,ngsourceOut,planeOut,nout] = ...
        removeDis(0,destroy,b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,...
        plane,pinned,irmbound,eta_obs,vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,...
        ngsourceOut,planeOut,nout);
end

%  if size(xdis,1)~=1000; error('!'); end

% now check numerical error and update dislocation position
    xdistrial = xdis(1:ndis); 
    ydistrial = ydis(1:ndis); 
    rdistrial = rdis(1:ndis);
    vdis = vdispre;
    
    for ng = 1:ngr
        Npolygon1 = polynode{ng,1}; % outside boundary should be anticlockwise
        if ispolycw(Npolygon1(:,1), Npolygon1(:,2)); Npolygon1 = flipud(Hn); end % check
        holeind=find([holes{:,2}]'==ng);
        for i_h=1:length(holeind)
            Hn = holes{holeind(i_h),1}; % hole nodes were ACW, change to CW
            if ~ispolycw(Hn(:,1), Hn(:,2)); Hn = flipud(Hn); end
            Npolygon1 = [Npolygon1; NaN NaN; Hn];
        end
        m1 = find(ngsource==ng);
        [in,on] = inpolygon(xdistrial(m1),ydistrial(m1),Npolygon1(:,1),Npolygon1(:,2));
        m = find(~(in|on));
		count_numerror = 0;
        while ~isempty(m)             
            count_numerror = count_numerror + 1;
            if any(vdis(m1(m))==0)  && all(source(m1(m))>0) ; error('vdis'); end
            %             vdis(m1(m)) = vdis(m1(m))*0.99999999;
            if count_numerror<2500;  vdis(m1(m)) = vdis(m1(m))*0.999;
            else;  vdis(m1(m)) = vdis(m1(m))*0.97;
                %             else; vdis(m1(m)) = 0;
            end
            xdistrial(m1(m)) = xdis(m1(m)) - vdis(m1(m))*dt.*cos(alpha(m1(m)));
            ydistrial(m1(m)) = ydis(m1(m)) - vdis(m1(m))*dt.*sin(alpha(m1(m)));
            rdistrial(m1(m)) = rdis(m1(m)) - vdis(m1(m))*dt;
            [in,on] = inpolygon(xdistrial(m1),ydistrial(m1),Npolygon1(:,1),Npolygon1(:,2));
            m = find(~(in|on));
            if count_numerror>1000; disp('in here'); break; end
        end
    end
    vdispre = vdis;
    rdis(1:ndis) = rdistrial;
    xdis(1:ndis) = xdistrial;
    ydis(1:ndis) = ydistrial;

%  if size(xdis,1)~=1000; error('!'); end


if any(irmbound~=0) && GBtrans == 1       
    [destroy,ndis,xdis,ydis,rdis,type,alpha,ngsource,b,irmbound,pinned,eta_obs,...
        source,vdis,vdispre,plane] ...
    = distrans(ndis,xdis,ydis,rdis,type,...
    alpha,ngsource,b,irmbound,...
    p_connect,bG,mu,p_xend,p_yend,p_rend,pinned,...
    eta_obs,source,vdis,vdispre,plane,reps,...
    logFID,taudis);

%  if size(xdis,1)~=1000; error('!'); end

    if any(destroy==1)
%         vdispre(destroy==1)=[]; vdispre = [vdispre; zeros(sum(destroy),1)];
        [b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,plane,pinned,irmbound,eta_obs,...
            vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,ngsourceOut,planeOut,nout] = ...
            removeDis(0,destroy,b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,...
            plane,pinned,irmbound,eta_obs,vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,...
            ngsourceOut,planeOut,nout);
    end
end

%  if size(xdis,1)~=1000; error('!'); end

disnuc;

%  if size(xdis,1)~=1000; error('!'); end

[destroy] = disannih(plane,ndis,rdis,type,Le);
if any(destroy==1)
    [b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,plane,pinned,irmbound,eta_obs,...
        vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,ngsourceOut,planeOut,nout] = ...
        removeDis(0,destroy,b,xdis,ydis,type,alpha,ngsource,ndis,source,rdis,...
        plane,pinned,irmbound,eta_obs,vdispre,bOut,xdisOut,ydisOut,typeOut,alphaOut,...
        ngsourceOut,planeOut,nout);
end

%  if size(xdis,1)~=1000; error('!'); end
savedata;

if mod(increment,100)==0
    fprintf(logFID,'>> Work done:%%%f, Increment:%d, ndis:%d, nout:%d, t:%e, dt:%e.\n',t/totalT*100,increment,ndis,nout,t,dt);
    fprintf('>> Work done:%6.2f, dt: %e, ndis: %i, nout %i \n',t/totalT*100,dt,ndis,nout);
    simplot
end

end

tt = toc

save(fullfile(['../',directoryOut,fname,'.mat']));
simplot;

saveas(figure(1),['../',directoryOut,fname,'_dislocations.fig'])
saveas(figure(2),['../',directoryOut,fname,'_stress-strain.fig'])
% saveas(figure(3),['../',directoryOut,fname,'_stress-time.fig'])
fprintf(logFID,'Time elapse of main loop is: %f s\n',tt);
fprintf(logFID,'Simulation end: %s \n',datestr(now));
fprintf(logFID,'All done! \n');
fclose all;
dbclear all;


cd ../
