function [p,t]=distmesh2d_SW0(h0,h1,bbox,pfix,Npolygon1,holes,nonuni,ng,varargin)
%DISTMESH2D 2-D Mesh Generator using Distance Functions.
%   [P,T]=DISTMESH2D(FD,FH,H0,BBOX,PFIX,FPARAMS)
%
%      P:         Node positions (Nx2)
%      T:         Triangle indices (NTx3)
%      FD:        Distance function d(x,y)
%      FH:        Scaled edge length function h(x,y)
%      H0:        Initial edge length
%      BBOX:      Bounding box [xmin,ymin; xmax,ymax]
%      PFIX:      Fixed node positions (NFIXx2)
%      FPARAMS:   Additional parameters passed to FD and FH
%
%   Example: (Uniform Mesh on Unit Circle)
%      fd=@(p) sqrt(sum(p.^2,2))-1;
%      [p,t]=distmesh2d(fd,@huniform,0.2,[-1,-1;1,1],[]);
%
%   Example: (Rectangle with circular hole, refined at circle boundary)
%       fd=@(p) ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.5));
%       fh=@(p) 0.05+0.3*dcircle(p,0,0,0.5);
%       [p,t]=distmesh2d(fd,fh,0.05,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1]);
%
%   Example: (Polygon)
%      pv=[-0.4 -0.5;0.4 -0.2;0.4 -0.7;1.5 -0.4;0.9 0.1;
%          1.6 0.8;0.5 0.5;0.2 1;0.1 0.4;-0.7 0.7;-0.4 -0.5];
%      [p,t]=distmesh2d(@dpoly,@huniform,0.1,[-1,-1; 2,1],pv,pv);
%
%   Example: (Ellipse)
%      fd=@(p) p(:,1).^2/2^2+p(:,2).^2/1^2-1;
%      [p,t]=distmesh2d(fd,@huniform,0.2,[-2,-1;2,1],[]);
%
%   Example: (Square, with size function point and line sources)
%      fd=@(p) drectangle(p,0,1,0,1);
%      fh=@(p) min(min(0.01+0.3*abs(dcircle(p,0,0,0)), ...
%                   0.025+0.3*abs(dpoly(p,[0.3,0.7; 0.7,0.5]))),0.15);
%      [p,t]=distmesh2d(fd,fh,0.01,[0,0;1,1],[0,0;1,0;0,1;1,1]);
%
%   Example: (NACA0012 airfoil)  
%      hlead=0.01; htrail=0.04; hmax=2; circx=2; circr=4;
%      a=.12/.2*[0.2969,-0.1260,-0.3516,0.2843,-0.1036];
%
%      fd=@(p) ddiff(dcircle(p,circx,0,circr),(abs(p(:,2))-polyval([a(5:-1:2),0],p(:,1))).^2-a(1)^2*p(:,1));
%      fh=@(p) min(min(hlead+0.3*dcircle(p,0,0,0),htrail+0.3*dcircle(p,1,0,0)),hmax);
%
%      fixx=1-htrail*cumsum(1.3.^(0:4)');
%      fixy=a(1)*sqrt(fixx)+polyval([a(5:-1:2),0],fixx);
%      fix=[[circx+[-1,1,0,0]*circr; 0,0,circr*[-1,1]]'; 0,0; 1,0; fixx,fixy; fixx,-fixy];
%      box=[circx-circr,-circr; circx+circr,circr];
%      h0=min([hlead,htrail,hmax]);
%
%      [p,t]=distmesh2d(fd,fh,h0,box,fix);

%
%   See also: MESHDEMO2D, DISTMESHND, DELAUNAYN, TRIMESH.

%   distmesh2d.m v1.1
%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.
    
dptol=.001; ttol=.1; Fscale=1.2; deltat=.2; geps=.001*h0; deps=sqrt(eps)*h0; % deps=sqrt(eps)*h0;
densityctrlfreq=30;

% 1. Create initial distribution in bounding box (equilateral triangles)
[x,y]=meshgrid(bbox(1,1):h0:bbox(2,1),bbox(1,2):h0*sqrt(3)/2:bbox(2,2));
x(2:2:end,:)=x(2:2:end,:)+h0/2;                      % Shift even rows
p=[x(:),y(:)];                                       % List of node coordinates

% SW - define fd
holeind=find([holes{:,2}]'==ng);
fd=dpolySW(p,Npolygon1);
for i_h=1:length(holeind)
    Hn = holes{holeind(i_h),1}; % hole nodes
    fd = ddiffSW(fd,dpolySW(p,Hn));
end

% 2. Remove points outside the region, apply the rejection method
p=p(fd<geps,:);                 % Keep only d<0 points



if nonuni==1
    fd=dpolySW(p,Npolygon1);
    for i_h=1:length(holeind)
        Hn = holes{holeind(i_h),1}; % hole nodes
        fd = ddiffSW(fd,dpolySW(p,Hn));
    end
    
     fh = min(h1,max(h0/3+0.3*abs(fd),h0)); % dpolySW(p,Npolygon1)
%      fh = min(0.25,max(h0/3+0.3*abs(fd),h0)); 
%      max(h0/2+0.3*abs(dpolySW(p,Npolygon1)),h0);
else
    fh = ones(size(p,1),1);
end

r0=1./fh; % feval(fh,p,varargin{:}).^2;              % Probability to keep point
p=p(rand(size(p,1),1)<r0./max(r0),:);                % Rejection method
if ~isempty(pfix), p=setdiff(p,pfix,'rows'); end     % Remove duplicated nodes
pfix=unique(pfix,'rows'); nfix=size(pfix,1);
p=[pfix; p];                                         % Prepend fix points
N=size(p,1);                                         % Number of points N

count=0;
pold=inf;                                            % For first iteration
clf,view(2),axis equal,axis off
while 1
  count=count+1;
  % 3. Retriangulation by the Delaunay algorithm
  if max(sqrt(sum((p-pold).^2,2))/h0)>ttol           % Any large movement?
    pold=p;                                          % Save current positions
    t=delaunayn(p);                                  % List of triangles
    
    pmid=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;    % Compute centroids
    fd=dpolySW(pmid,Npolygon1);
    for i_h=1:length(holeind)
        Hn = holes{holeind(i_h),1}; % hole nodes
        fd = ddiffSW(fd,dpolySW(pmid,Hn));
    end
    
    t=t(fd<-geps,:);         % Keep interior triangles
    
    % 4. Describe each bar by a unique pair of nodes
    bars=[t(:,[1,2]);t(:,[1,3]);t(:,[2,3])]; %3900x2 % Interior bars duplicated
    bars=unique(sort(bars,2),'rows'); %1998x2        % Bars as node pairs
    % 'sort(bars,2)' = sort each of the row in the ascending order
    
    
    % 5. Graphical output of the current mesh
    cla,patch('vertices',p,'faces',t,'edgecol','k','facecol',[.8,.9,1]);
    drawnow;
  end

  % 6. Move mesh points based on bar lengths L and forces F
  barvec=p(bars(:,1),:)-p(bars(:,2),:);              % List of bar vectors
  L=sqrt(sum(barvec.^2,2)); % sum of rows            % L = Bar lengths
 
  
  
  
  if nonuni==1
      
       fd=dpolySW((p(bars(:,1),:)+p(bars(:,2),:))/2,Npolygon1);
    for i_h=1:length(holeind)
        Hn = holes{holeind(i_h),1}; % hole nodes
        fd = ddiffSW(fd,dpolySW((p(bars(:,1),:)+p(bars(:,2),:))/2,Hn));
    end
    
        hbars = min(h1,max(h0/3+0.3*abs(fd),h0));
  else
      %Generate a column vector of same size as position of bars
      hbars = ones(size((p(bars(:,1),:)+p(bars(:,2),:))/2,1),1);
  end
%   hbars=feval(fh,(p(bars(:,1),:)+p(bars(:,2),:))/2,varargin{:});
  L0=hbars*Fscale*sqrt(sum(L.^2)/sum(hbars.^2));     % L0 = Desired lengths
  
  % Density control - remove points that are too close
  if mod(count,densityctrlfreq)==0 && any(L0>2*L)
      p(setdiff(reshape(bars(L0>2*L,:),[],1),1:nfix),:)=[];
      N=size(p,1); pold=inf;
      continue;
  end
  
  F=max(L0-L,0);                                     % Bar forces (scalars)
  Fvec=F./L*[1,1].*barvec;                           % Bar forces (x,y components)
  Ftot=full(sparse(bars(:,[1,1,2,2]),ones(size(F))*[1,2,1,2],[Fvec,-Fvec],N,2));
  Ftot(1:size(pfix,1),:)=0;                          % Force = 0 at fixed points
  p=p+deltat*Ftot;                                   % Update node positions

  % 7. Bring outside points back to the boundary
  d=dpolySW(p,Npolygon1);
  for i_h=1:length(holeind)
      Hn = holes{holeind(i_h),1}; % position hole nodes
      d = ddiffSW(d,dpolySW(p,Hn));
  end
    ix=d>0;                 % Find points outside (d>0), -ve and zero ar inside the polygon

  pnum = [p(ix,1)+deps,p(ix,2)]; % access position(x,y) of points outside
  fd=dpolySW(pnum,Npolygon1);
    for i_h=1:length(holeind)
        Hn = holes{holeind(i_h),1}; % hole nodes
        fd = ddiffSW(fd,dpolySW(pnum,Hn));
    end
  dgradx=(fd-d(ix))/deps; % Numerical
  
  pgrad=[p(ix,1),p(ix,2)+deps];
  fd=dpolySW(pgrad,Npolygon1);
  for i_h=1:length(holeind)
      Hn = holes{holeind(i_h),1}; % hole nodes
      fd = ddiffSW(fd,dpolySW(pgrad,Hn));
  end
  dgrady=(fd-d(ix))/deps; %    gradient
  dgrad2=dgradx.^2+dgrady.^2;
  p(ix,:)=p(ix,:)-[d(ix).*dgradx./dgrad2,d(ix).*dgrady./dgrad2];    % Project
 
    % 8. Termination criterion: All interior nodes move less than dptol (scaled)
  if max(sqrt(sum(deltat*Ftot(d<-geps,:).^2,2))/h0)<dptol, break; end
  
end

% Clean up and plot final mesh
[p,t]=fixmesh(p,t);
simpplot(p,t)

end

function d=dpolySW(p,pv) % p=nodes positions, pv=boundary nodes of the grain

%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

np=size(p,1);
nvs=size(pv,1)-1;

ds=dsegment(p,pv);
d=min(ds,[],2); % access minimum value from each row

d=(-1).^(inpolygon(p(:,1),p(:,2),pv(:,1),pv(:,2))).*d;
d(abs(d)<100*eps)=0;
end


function d=ddiffSW(d1,d2), d=max(d1,-d2);
end
