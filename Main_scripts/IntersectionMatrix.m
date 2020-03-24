%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slip planes and Mesh elements intersection
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
    slp_ele_indG = zeros(size(nc,1), sum(nslp,'all'));
    slp_aboveNode_ind = zeros(size(p,1), sum(nslp,'all'));

for islp = 1:sum(nslp,'all') % all slip planes
    aff_ele_ind = zeros(size(nc,1),1);
    if ~isequal(p_xend(islp,:),(p_yend(islp,:))) % plane should not be a point
    slplane = [p_xend(islp,:)' p_yend(islp,:)']; % [x1 y1; x2 y2]
    
    slangle = atan((slplane(2,2)-slplane(1,2) ) / (slplane(2,1)-slplane(1,1))); % calculate angle
    if slangle<0;   slangle = slangle+pi; end
    if slangle>=pi; slangle = slangle-pi; end
    
    min_x = min(slplane(:,1));  max_x = max(slplane(:,1));
    min_y = min(slplane(:,2));  max_y = max(slplane(:,2));
    
    if slangle == pi/2   % vertical
    above_box = [min_x min_y-elesize;...
                max_x max_y+elesize;...
                max_x-1.2*elesize max_y+elesize;...
                min_x-1.2*elesize min_y-elesize;...
                min_x min_y-elesize ];   
    elseif  slangle == 0   % horizontal 
    above_box = [min_x-1.2*elesize min_y;...
                max_x+elesize max_y;...
                max_x+elesize max_y+1.2*elesize;...
                min_x-1.2*elesize min_y+1.2*elesize;...
                min_x-1.2*elesize min_y ];               
    elseif slangle<pi/2   % acute
    above_box = [  min_x-elesize*cos(slangle) min_y-elesize*sin(slangle);...
                max_x+elesize*cos(slangle) max_y+elesize*sin(slangle);...
                max_x+elesize*cos(slangle)-1.2*elesize max_y+elesize*sin(slangle);...
                min_x-elesize*cos(slangle)-1.2*elesize min_y-elesize*sin(slangle);...
                min_x-elesize*cos(slangle) min_y-elesize*sin(slangle) ];        
    elseif slangle>pi/2   % obtuse
    above_box = [  max_x-elesize*cos(slangle) min_y-elesize*sin(slangle);...
                min_x+elesize*cos(slangle) max_y+elesize*sin(slangle);...
                min_x+elesize*cos(slangle)-1.2*elesize max_y+elesize*sin(slangle);...
                max_x-elesize*cos(slangle)-1.2*elesize min_y-elesize*sin(slangle);...
                max_x-elesize*cos(slangle) min_y-elesize*sin(slangle) ];        
    end   
    
    % obtuse should be
%         above_box = [  max_x+elesize*cos(slangle)             min_y-elesize*sin(slangle);...
%                    min_x-1.2*elesize*cos(slangle)             max_y+1.2*elesize*sin(slangle);...            
%                    min_x-1.2*elesize*cos(slangle)-1.2*elesize max_y+1.2*elesize*sin(slangle);...  
%                    max_x+elesize*cos(slangle)-1.2*elesize     min_y-elesize*sin(slangle);...
%                    max_x+elesize*cos(slangle)                 min_y-elesize*sin(slangle)  ];        
          
    fd = dpolySW(p,above_box); % finding nodes inside above_box
    above_ind = find(fd<0); % node numbers   
    
    nc_affected = ismember(nc, above_ind); % binary output with 3 columns
    affected_nodesTemp = nc(nc_affected==1); % = above_ind but with repetitions
%     affected_ele = nc(any(nc_affected,2),:); % affected elements(i.e.
%     rows of nc) with node numbers
%     plot(p(affected_nodes,1), p(affected_nodes,2), 'y*')
%     xi = [];
%     yi = [];
%     xele = [];
%     yele = [];

    aff_ele_ind = any(nc_affected,2); % indices of affected elements (binary 1-column)
    
 for iaff = 1:size(aff_ele_ind,1) % for neighbouring nodes
     if aff_ele_ind(iaff)==1 % elements which are in above_box
        node_no = nc(iaff,:); % 1x3 node numbers
        ele_pos = p(node_no',:); % 3x2 nodal positions
        [xi_temp,yi_temp] = polyxpoly(slplane(:,1), slplane(:,2), ele_pos(:,1), ele_pos(:,2));
        if isempty(xi_temp) || isempty(yi_temp) % nodes not intersecting
            aff_ele_ind(iaff) = 0;              
        end
     end
 end
    slp_ele_indG(:,islp) = aff_ele_ind;
    slp_aboveNode_ind(above_ind,islp) = 1;
    end
    
end

 
% figure; hold on
% plot(above_box(:,1),above_box(:,2), 'r-' , 'Linewidth', 1)
% plot(slplane(:,1), slplane(:,2), 'y-', 'Linewidth', 2)
% plot(p(above_ind,1), p(above_ind,2), 'y*')
% realaff_el = find(slp_ele_indG(:,islp)==1);
% plot(p(nc(realaff_el,:),1),p(nc(realaff_el,:),2),'*b')


fprintf(logFID,'IntersectionMatrix done! \n\n');
fprintf('IntersectionMatrix done! \n')
toc 
