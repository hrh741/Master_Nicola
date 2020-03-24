% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Function     : InGraingen
%   Last edited  : 1 November, 2018 - SW
%   Description  : called by Input.m
%                     Generate grain vertices for rectangular grains
%   Outstanding issues : Holes defined
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

polynode = cell(ngr,2); % 1. nodal positions, 23 material type
holes = cell(1,2); % 1. vertices 2. Grain number hole belongs

% define nodes of each grain (i.e. 5)
% Ist column contains x-components and 2nd contains y-components
outbox = [0 xcst xcst 0 0; 0 0 ycst ycst 0]'; 
samplenode = outbox;
% Total size of the sample which contain grains, 2x2 matrix
% co-ordinates of bottom-left corner and upper-right corner
bbox = [min(samplenode(:,1)),min(samplenode(:,2)); max(samplenode(:,1)),max(samplenode(:,2))]; 

xgrain = xcst / ngrx;    % Length of each grain in x-direction
ygrain = ycst / ngry;    % Length of each grain in y-direction

% if flag==1
holetype = 2; % 0. Circular 1. Rectangular 2. Elliptical

    if holetype == 0
    % Circular hole
    radC = 0.125*xgrain; % radius of circle
    thetastep = (0.5*elesize/radC);
    
    elseif holetype == 1
    % Rectangular grain
    holeedge_x = 0.5*xgrain; % Length of hole in x-direction
    holeedge_y = 0.5*ygrain; % Length of hole in y-direction

    elseif holetype == 2
    % Elliptical Hole
    ellip_x = 0.25*xgrain;  % Major axis
    ellip_y = 0.125*ygrain; % Minor axis
    thetastep = 0.5*minelesize;
    Periellip = 2*pi*sqrt(((ellip_x)^2+(ellip_y)^2)/2); % Perimeter of an Ellipse
    npoints = round(1*(Periellip/minelesize)); % no of points on Ellipse
    
    end
% end
% Define each grain 
for ngy = 1:ngry
for ngx = 1:ngrx

% 1+(ngry-1)*ngrx   ...      ...     ngry*ngrx
% .                 .                .
% 1+ngrx            2+ngrx   ...     2*ngrx
% 1                 2        ...     ngrx

ng = ngx+(ngy-1)*ngrx; % Grain Number 

% 'Grain' vertices position
xleft  = (ngx-1)*xgrain;
xright = ngx*xgrain;
ylower = (ngy-1)*ygrain;
yupper = ngy*ygrain;

grainvertices = [xleft ylower; xright ylower; xright yupper; xleft yupper; xleft ylower];
polynode{ng,1} = grainvertices; % defined anticlockwise 
polynode{ng,2} = 1; % Material type

% if flag==1 
    nholesng = 1; % number of holes in one grain
    % Defining center of hole
    holecenter = [mean([xleft,xright]) mean([ylower,yupper])];

    if holetype == 1
    % Defining a Rectangular hole
    xleft  = holecenter(1,1)-0.25*holeedge_x;
    xright = holecenter(1,1)+0.25*holeedge_x;
    yupper = holecenter(1,2)+0.25*holeedge_y;
    ylower = holecenter(1,2)-0.25*holeedge_y;
    holevertices = [xleft ylower; xright ylower; xright yupper; xleft yupper; xleft ylower];
    % boundary nodes of the hole
    holes{ng,1} = holevertices; % for rectangular hole

    elseif holetype == 2
    % Defining an ellipse
    theta = 0:thetastep:2*pi;
    x1 = ellip_x*cos(theta)+holecenter(1,1);
    y1 = ellip_y*sin(theta)+holecenter(1,2);   
    x1(end) = x1(1); y1(end) = y1(1);
    pt = interparc(npoints,x1,y1,'csape'); % to create equi-distance nodes
    % boundary nodes of the hole
    holes{ng,1} = [pt(:,1) pt(:,2)]; 

    elseif holetype == 0
    % Defining a circle
    theta = 0:thetastep:2*pi;
    x1 = radC*cos(theta)+holecenter(1,1);
    y1 = radC*sin(theta)+holecenter(1,2);
    x1(end) = x1(1); y1(end) = y1(1);
    % boundary nodes of the hole
    holes{ng,1} = [x1; y1]'; % for elliptical and circular hole

    end

    holes{ng,2} = ng; % number of grain associated with the hole, saved in 2nd column
% end

% define different material
if exist('matdef','var'); if any(ng == matdef); polynode{ng,2} = 2; end; end

end
end

clearvars xgrain ygrain ngx ngy grainvertices

%% plot 
figure; clf; hold on
for ng = 1:ngr
    xnode = polynode{ng};
    plot(xnode(:,1),xnode(:,2),'b-','LineWidth',1.5);
    text(mean(xnode(:,1)),mean(xnode(:,2)),num2str(ng));
end

    for ng = 1:size(holes,1)
        xnode = holes{ng,1};
        plot(xnode(:,1),xnode(:,2),'r-o','LineWidth',1);
    end
    
plot(outbox(:,1),outbox(:,2),'k-','LineWidth',1.5)
axis equal; axis off;
ax=axis;axis(ax*1.001);
hold off
