% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Function     : InGraingen
%   Last edited  : 1 November, 2018 - SW
%   Description  : called by Input.m
%                     Generate grain vertices for rectangular grains
%   Outstanding issues : Holes not defined
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

polynode = cell(ngr,2); % 1 nodal positions, 2 material type

nholes = 0; % HOLES NOT DEFINED
holes = cell(nholes,2); 

% define nodes of each grain (i.e. 5)
% Ist column contains x-components and 2nd contains y-components
outbox = [0 xcst xcst 0 0; 0 0 ycst ycst 0]'; 
samplenode = outbox;
% Total size of the sample which contain grains, 2x2 matrix
% co-ordinates of bottom-left corner and upper-right corner
bbox = [min(samplenode(:,1)),min(samplenode(:,2)); max(samplenode(:,1)),max(samplenode(:,2))]; 

xgrain = xcst / ngrx; % Length of each grain in x-direction
ygrain = ycst / ngry; % Length of each grain in y-direction

% Define each grain 
for ngy = 1:ngry
for ngx = 1:ngrx

% 1+(ngry-1)*ngrx   ...      ...     ngry*ngrx
% .                 .                .
% 1+ngrx            2+ngrx   ...     2*ngrx
% 1                 2        ...     ngrx

ng = ngx+(ngy-1)*ngrx; % Grain Number 

xleft = (ngx-1)*xgrain;
xright = ngx*xgrain;
ylower = (ngy-1)*ygrain;
yupper = ngy*ygrain;

grainvertices = [xleft ylower; xright ylower; xright yupper; xleft yupper; xleft ylower];
polynode{ng,1} = grainvertices; % defined anticlockwise 
polynode{ng,2} = 1; 

% define different material
if exist('matdef','var'); if any(ng == matdef); polynode{ng,2} = 2; end; end

end
end

clearvars xgrain ygrain ngx ngy grainvertices

%% plot 
figure; clf; hold on
for ng = 1:ngr
    xnode = polynode{ng};
    plot(xnode(:,1),xnode(:,2),'r-','LineWidth',1.5);
    text(mean(xnode(:,1)),mean(xnode(:,2)),num2str(ng));
end
% for ng = 1:size(holes,1)
%     xnode = holes{ng,1};
%     plot(xnode(:,1),xnode(:,2),'-','LineWidth',1);
% end
plot(outbox(:,1),outbox(:,2),'k-','LineWidth',1.5)
axis equal; axis off;
ax=axis;axis(ax*1.001);
hold off
