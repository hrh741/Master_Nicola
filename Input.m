%
% Last update 7th July, 2019 - SW
% Outstanding issues: Holes/voids not included InSets, MeshGen, FEdefine as traction free boundary (freedofs)
%
%
%
%-------------------------------%

fclose all; close all; clear all;
% Name the output directory
directoryOut = sprintf('./Nicola try9 %s/', datestr(now,'dd-mm-yyyy HH-MM'));
pause('on')
% Load old structure
pre_source = 1; % 0 - No; 1 - Yes

if pre_source==0
%% Crystal/Grain input
xcst = 12; % length of specimen (microns)
ycst = 1; % width of specimen (microns)
ngrx = 8; % number of grains along x-direction
ngry = 1; % number of grains along y-direction
ngr = ngrx*ngry; % number of total grains
% for bimaterial
% matdef = [2;3]; % which grains are second material type

%% Finite element triangulation parameter
elesize = 0.125; % max dimesnion of meshing triangles 
minelesize = 0.05; % min element size 

cd ./In_scripts/; InGraingen; cd ../;

%% Material parameters
% deltaF and gammazero are related to thermally activated processes
% Initiallization 
material{1} = 'Ti6242 alpha'; material{2} = 'Ti6242 beta';
nu=zeros(ngr,1); e=zeros(ngr,1); bG=zeros(ngr,1); tau_FR=zeros(ngr,1); deltaF=zeros(ngr,1); gammazero=zeros(ngr,1); % bG = Burger Vectore
dlambda = zeros(ngr,1); orientation = zeros(ngr,1); % dlambda is the angle between the slip systems
% units: nu(Poisson's ratio)=constant ; e=TPa; bG=um; tau_FR=TPa; deltaF=J;
% gammazero(ng) ; dlambda=rad; orientation=rad

a=0; b=59;
for ng=1:ngr
    nu(ng) = 0.34; e(ng) = 0.11; bG(ng) = 0.25E-3; tau_FR(ng) = 1E-4; deltaF(ng) = 9.01E-20; gammazero(ng) = 9E-06;
    dlambda(ng) = 60/180.0*pi; orientation(ng) = ((b-a).*rand(1,1) + a)/180.0*pi;
    
%     if polynode{ng,2}==1 % first material
%         nu(ng) = 0.34; e(ng) = 0.11; bG(ng) = 0.25E-3; tau_FR(ng) = 1E-4; deltaF(ng) = 9.01E-20; gammazero(ng) = 9E-06;
%         dlambda(ng) = 60/180.0*pi; orientation(ng) = ((b-a).*rand(1,1) + a)/180.0*pi;
%     elseif polynode{ng,2}==2 % second material
%         nu(ng) = 0.34; e(ng) = 0.11; bG(ng) = 0.25E-3; tau_FR(ng) = 1E-4; deltaF(ng) = 9.01E-20; gammazero(ng) = 9E-06;
%         dlambda(ng) = 60/180.0*pi; orientation(ng) = ((b-a).*rand(1,1) + a)/180.0*pi;
%     elseif polynode{ng,2}==3 % second material
%         nu(ng) = 0.34; e(ng) = 0.11; bG(ng) = 0.25E-3; tau_FR(ng) = 1E-4; deltaF(ng) = 9.01E-20; gammazero(ng) = 9E-06;
%         dlambda(ng) = 60/180.0*pi; orientation(ng) = 15/180.0*pi;
%     elseif polynode{ng,2}==4 % second material
%         nu(ng) = 0.34; e(ng) = 0.11; bG(ng) = 0.25E-3; tau_FR(ng) = 1E-4; deltaF(ng) = 9.01E-20; gammazero(ng) = 9E-06;
%         dlambda(ng) = 60/180.0*pi; orientation(ng) = 0/180.0*pi;
%     elseif polynode{ng,2}==5 % second material
%         nu(ng) = 0.34; e(ng) = 0.11; bG(ng) = 0.25E-3; tau_FR(ng) = 1E-4; deltaF(ng) = 9.01E-20; gammazero(ng) = 9E-06;
%         dlambda(ng) = 60/180.0*pi; orientation(ng) = 30/180.0*pi;
%     elseif polynode{ng,2}==6 % second material
%         nu(ng) = 0.34; e(ng) = 0.11; bG(ng) = 0.25E-3; tau_FR(ng) = 1E-4; deltaF(ng) = 9.01E-20; gammazero(ng) = 9E-06;
%         dlambda(ng) = 60/180.0*pi; orientation(ng) = ((b-a).*rand(1,1) + a)/180.0*pi;
%     elseif polynode{ng,2}==7 % second material
%         nu(ng) = 0.34; e(ng) = 0.11; bG(ng) = 0.25E-3; tau_FR(ng) = 1E-4; deltaF(ng) = 9.01E-20; gammazero(ng) = 9E-06;
%         dlambda(ng) = 60/180.0*pi; orientation(ng) = ((b-a).*rand(1,1) + a)/180.0*pi;
%     elseif polynode{ng,2}==8 % second material
%         nu(ng) = 0.34; e(ng) = 0.11; bG(ng) = 0.25E-3; tau_FR(ng) = 1E-4; deltaF(ng) = 9.01E-20; gammazero(ng) = 9E-06;
%         dlambda(ng) = 60/180.0*pi; orientation(ng) = ((b-a).*rand(1,1) + a)/180.0*pi;
%     else
%         error('Grain material property desciption');
%     end
end

mu = e./(2*(1+nu)); % Shear modulus (TPa)
sourceDensity = 15;     % Sources density (/micron^2)
sdev_FR = 0.2*tau_FR.*ones(ngr,1);   % standard deviation of Gaussian distribution(TPa)
obsDensity = 0;        % Obstacles density (/micron^2)
tau_obs = 3*tau_FR;     % strength of all Obstacles (TPa)
nSystems = 3*ones(ngr,1); % number of slip systems of each grain
Bdrag = 1E-16;          % Drag coefficient, (TPa.s)
tnuc0 = 10.0E-9;        % nucleation time (s)
vcutoff = 20E6;         % Cut-off velocity (microns/s)
friction = 0;           % Inertia of the crystal against dislocation movement. 
                        % Does friction exit? 0 - no; 1 - yes
GBtrans = 0;            % 0 - hard GB; 1 - penetrable GB


%% Define the applied boundary conditions
BCs = 1; % fixed left boundary, disp applied to right boundary
         % (1) Strain control with constant strain rate right surface; 
         % (2) Strain control with constant strain rate top surface;
         % (3) Stress control with constant loading rate;
         % (4) Traction free all surfaces.
Umax = 0.006*(bbox(2,1)-bbox(1,1)); % maximum applied displacement (microns)
Udot = 2500*(bbox(2,1)-bbox(1,1)); % applied displacement rate s^-1
totalT = Umax/Udot; % total loading time in seconds
Fdot = 0; % not used in displacement-controlled loading

% BCs = 2;
% Umax = 0.01*(bbox(2,2)-bbox(1,2)); % maximum applied displacement (microns)
% Udot = 2000*(bbox(2,2)-bbox(1,2)); % applied displacement rate s^-1
% totalT = Umax/Udot; % total loading time in second
% Fdot = 0; % not used in displacement-controlled loading

% BCs = 3;
% Umax = []; % maximum applied displacement (microns)
% Udot = []; % applied displacement rate;
% Appliedstress = 700; % MPa
% Max applied force needs to be calculated according generated mesh
% Fmax = Appliedstress*1e-6*(bbox(2,2)-bbox(1,2))/17; % maximum applied Force (N)
% Fdot = 0.5; % applied force rate - reach Fmax at Xs;
% totalT = Fmax/Fdot; % total loading time in second

%% Other input functions
cd ./In_scripts/; 
InComputational;        % Computational parameters
InConstDefine;          % Define constants such annilihation dist, edgeconst, etc
cd ../;
setLogFile;             % make output directory

%% Source generation
cd ./In_scripts/; 
[lambda, R] = SlipSysGen(nSystems,ngr,orientation,dlambda);

[xSource,ySource,planeSource,alphaSource,rSource,tauSource,L,sourcecount,tnuc,...
    xObs,yObs,planeObs,alphaObs,rObs,obscount,nslp,p_rend,p_xend,p_yend,p_connect]...
    = SourceGen(bG,polynode,holes,bbox,ngr,nSystems,sourceDensity,obsDensity,edgeType,lambda,tau_FR,sdev_FR,...
    tnuc0,directoryOut,samplenode,logFID,[],[]);
cd ../;


elseif pre_source==1
    disp('Loading old source structure. Done!');
    load('source.mat');
    directoryOut = sprintf('./Nicola try2 New %s/', datestr(now,'dd-mm-yyyy HH-MM'));
    setLogFile;

    
% Change values of variables here
%     elesize = 0.125; % max dimension of meshing triangles
%     nrppd = 5000; % no. of increment at which save the dislocation structure.
%    GBtrans = 1; % 0 - hard GB; 1 - penetrable GB
%     Umax = 0.02*(bbox(2,1)-bbox(1,1)); % maximum applied displacement (microns)
    
end


%%
save(fullfile(directoryOut,'source.mat'));

