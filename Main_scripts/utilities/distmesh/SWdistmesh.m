
cd ./utilities/distmesh/; 
p = [];
nc = [];
grainelenum = [];   % number of elements in each grain
elegrainindex = []; % use to index which grain is the element belongs

for ng = 1:ngr
    pv = Npolygon{ng};
    Npolygon1 = polynode{ng,1};
    fprintf('Mesh grain: %s in process ...\n', ...
                        num2str(ng));
    
                    
    if polyarea(pv(:,1),pv(:,2)) < elesize1^2; elesize = elesize1/4; else elesize = elesize1; end 
    
    holeind=find([holes{:,2}]'==ng);
    if isempty(holeind)
        fd=@dpoly; 
    else
%         hole number
        Hn = holes{holesind(i_h),1}
        fd=@(p) 
        dgrain = dpoly(p,Npolygon1); 
        dhole{i} = dpoly(p,Npolygon1
    end

     fd=@(p) ddiff(dcircle(p,circx,0,circr),(abs(p(:,2))-polyval([a(5:-1:2),0],p(:,1))).^2-a(1)^2*p(:,1));

    

    
    figure;
    [p1,t1] = distmesh2d(@dpoly,@huniform,elesize,bbox,pv,Npolygon1);
    Nold = size(p,1);
    p = [p; p1];
    nc = [nc; t1+Nold];
    grainelenum(ng,1) = size(t1,1);
    close all
end