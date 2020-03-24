% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Function     : MatrixReshape
%   Last edited  : 6 November, 2018 - SW
%   Description  : called by Main.m
%                   Linearise matirces for mexing (C++) - matrices to vector
%   Outstanding issues : 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reshape(x, numel(x),1)
RC = reshape(R,numel(R),1);
lambdaC = reshape(lambda,numel(lambda),1);
ncC = reshape(nc,numel(nc),1);
DC = reshape(D,numel(D),1);


xSourceC = [];
ySourceC = [];
planeSourceC = [];
for ng = 1:ngr
    xSourceC     = [xSourceC;xSource(1:sourcecount(ng),ng)];
    ySourceC     = [ySourceC;ySource(1:sourcecount(ng),ng)];
    planeSourceC = [planeSourceC;planeSource(1:sourcecount(ng),ng)];
end
sourcecountC = sum(sourcecount);