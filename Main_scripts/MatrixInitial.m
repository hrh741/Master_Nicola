% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Function     : MatrixInitial
%   Last edited  : 6 November, 2018 - SW
%   Description  : called by Main.m
%                    Initialize matrices
%   Outstanding issues : 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tn = 0*xSource;
Ln = 0*xSource;
t = 0;
told = 0;
uhat = zeros(2*mno,1);
uFE = zeros(2*mno,1);
rhat = zeros(2*mno,1);
u = zeros(2*mno,1);   
fa = zeros(2*mno,1); 
U = 0;
F = 0;
increment = 0;
saved_count=0;

ndis = 0;
ndispre = 0;
deltardismax = 0;
b = zeros(arraySize,1);
xdis = zeros(arraySize,1);
ydis = zeros(arraySize,1);
type = zeros(arraySize,1);
plane = zeros(arraySize,1);
alpha = zeros(arraySize,1);
source = zeros(arraySize,1);
ngsource = zeros(arraySize,1);
rdis = zeros(arraySize,1);
pinned = zeros(arraySize,1);
remove = zeros(arraySize,1);
irmbound = zeros(arraySize,1);
eta_obs = zeros(arraySize,1);
vdis = zeros(arraySize,1);
vdispre = zeros(arraySize,1);

nout = 0;
bOut = [];
xdisOut = [];
ydisOut = []; 
typeOut = []; 
alphaOut = [];
ngsourceOut = [];
planeOut = [];
node_slip = [];

UEND = zeros(iMax,1);
FEND = zeros(iMax,1);
TINC = zeros(iMax,1);
NDIS = zeros(iMax,1);
DT = zeros(iMax,1);

if any(BCs==[1,2])
    total_inc = ceil(Umax/Udot/dt0);
elseif BCs==3
    total_inc = ceil(Fmax/Fdot/dt0);
end