% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Function    : save dislocation data at current increment
%   Description : call by Main.m
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fend = sum(rFE(2*Sright-1)+rDD(2*Sright-1));
% uend = mean(uFE(2*Sright-1));

UEND(increment) = uend;
FEND(increment) = fend;
NDIS(increment) = ndis;
DT(increment) = dt;
TINC(increment) = t;
   
if increment == 1 || increment == 2
    savedata2;
elseif mod(increment,nrppd)==0 && ndis>0    
    simplot;
    save(['../',directoryOut,'increment_',num2str(increment),'.mat'],'uhat','rhat',...
        'U','b','xdis','ydis','type','alpha','ngsource','ndis','source','rdis',...
        'plane','pinned','irmbound','eta_obs','vdispre','bOut','xdisOut',...
        'ydisOut','typeOut','alphaOut','planeOut','ngsourceOut','nout','increment','t',...
        'FEND','UEND','NDIS','dt','uFE','rFE','rDD','DT','TINC')
end


