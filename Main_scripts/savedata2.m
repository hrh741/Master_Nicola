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
   
save(fullfile(['../',directoryOut,'increment_',num2str(increment),'.mat']));
simplot;
saveas(figure(1),['../',directoryOut,fname,'_dislocations.fig'])
saveas(figure(2),['../',directoryOut,fname,'_stress-strain.fig'])
saveas(figure(3),['../',directoryOut,fname,'_stress-time.fig'])

