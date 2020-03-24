%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% New neighbouring calculations

aboveNodesMatrix = cell(sum(nslp,'all'),1); 
stress_invhat = zeros(size(nc,1),3);
for iout = 1:nout
    
    islpOut = planeOut(iout); % plane number
    slp_ele_indOut = find(slp_ele_indG(:,islpOut)==1); % element indices intersecting slip plane
    slp_aboveNode_indTemp = find(slp_aboveNode_ind(:,islpOut)==1);  % neighboring node indices above sloip plane
    slplaneOut = [p_xend(islpOut,:)' p_yend(islpOut,:)']; % [x1 y1; x2 y2]
    slangleOut = alphaOut(iout);  
    
    % distance of each node from out_dislocation
    dx = 0 - xdisOut(iout); % consider nodes to be at origin
    dy = 0 - ydisOut(iout);
    r2 = dx*dx+dy*dy;
    dx_temp = [dx dy]*[ cos(slangleOut) sin(slangleOut)]';
    dy_temp = [dx dy]*[-sin(slangleOut) cos(slangleOut)]';
    
    r2inv = 1/r2;
    
    % 2*dx_temp because only considering abovenodes, unlike in dispMex
    ux_temp = (0.5*dx_temp*dy_temp*r2inv-(1-nu(ngsourceOut(iout))*atan(dx_temp/dy_temp))*typeOut(iout)*bOut(iout)*uConst(ngsourceOut(iout)));
%     ux_temp = (0.5*dx_temp*dy_temp*r2inv-(1-nu*atan(dx_temp/dy_temp))*typeOut(iout)*bOut(iout)*uConst(ngsourceOut(iout)));
    uy_temp = 0.0;
    ux = cos(slangleOut)*ux_temp - sin(slangleOut)*uy_temp;
    uy = sin(slangleOut)*ux_temp + cos(slangleOut)*uy_temp;
    
    aboveNodes = []; 

    force_invhatTemp = zeros(6,1);
    for i = 1:size(slp_ele_indOut,1)
        ind = slp_ele_indOut(i);
        eleNode_no = nc(ind,:); % nodes of element
        ele_pos = p(eleNode_no,:); % position of nodes
        x1 = ele_pos(1,1);
        x2 = ele_pos(2,1);
        x3 = ele_pos(3,1);
        y1 = ele_pos(1,2);
        y2 = ele_pos(2,2);
        y3 = ele_pos(3,2);
        
        Atri = 0.5*( (x2*y3-x3*y2) + (x3*y1-x1*y3) + (x1*y2-x2*y1) );
        
        % B_ele = derivative of standard shape functions
        B_eleSt = (1/(2*Atri)) * [ y2-y3    0    y3-y1    0    y1-y2    0
                                     0    x3-x2    0    x1-x3    0    x2-x1
                                   x3-x2  y2-y3  x1-x3  y3-y1  x2-x1  y1-y2 ];
                               
                               
        above_indTemp = ismember(eleNode_no, slp_aboveNode_indTemp);  % Check if that node exist in abovebox
        above_ind = find( above_indTemp==1 ); % either [1] [2] [3] viceversa or  [1;2], [2;3], [1;3] viceversa
        below_ind = find( above_indTemp~=1 );
        
        aboveNodesTemp = eleNode_no(above_ind)'; % node number of nodes which are in above_box

        if size(above_ind,1)==1 % one node in the above box
            
            B_enr = B_eleSt; 
            B_enr(:,[2*below_ind-1 2*below_ind]) =  0; % 3x2
                                 
            
        elseif size(above_ind,1)==2 % two nodes in above box
            % order is important            
            if all(ismember([1 2], above_ind))
                B_enr = (1/(2*Atri)) * [ y2-y1   0   y2-y1   0    0  0
                                           0   x1-x2   0   x1-x2  0  0
                                         x1-x2 y2-y1 x1-x2 y2-y1  0  0 ];
                
            elseif all(ismember([2 3], above_ind))
                B_enr = (1/(2*Atri)) * [ 0  0   y3-y2    0    y3-y2    0
                                         0  0    0    x2-x3    0    x2-x3
                                         0  0  x2-x3  y3-y2  x2-x3  y3-y2 ];
                
            elseif all(ismember([1 3], above_ind))
                B_enr = (1/(2*Atri)) * [ y1-y3    0    0 0   y1-y3    0
                                           0    x3-x1  0 0    0    x3-x1
                                         x3-x1  y1-y3  0 0  x3-x1  y1-y3 ];
                
            end
                        
        elseif size(above_ind,1)==3 || size(below_ind,1)==3
            error('3 nodes are in above_box ')
        end
        u_enr =  [ux uy ux uy ux uy]';
        u_enr([2*below_ind-1 2*below_ind]) = 0;
         
        strain_invhatTemp =  -(B_enr * u_enr);
        stress_invhatTemp = (D(:,:,ngsourceOut(iout)) * strain_invhatTemp); % [xx; yy; xy] 3x1
%         stress_invhatTemp = (D(:,:,ngsourceOut(iout)) * strain_invhatTemp(:,ngsourceOut(iout)));
        stress_invhat(ind,:) = stress_invhat(ind,:) + stress_invhatTemp'; % element wise 3x1
%         stress_invhat(eleNode_no) = stress_invhat(eleNode_no) + stress_invhatTemp;
        force_invhatTemp = Atri * B_eleSt' * stress_invhatTemp;
        
        fenr_invhat(2*eleNode_no-1) = fenr_invhat(2*eleNode_no-1) + force_invhatTemp(1:2:end);
        fenr_invhat(2*eleNode_no) = fenr_invhat(2*eleNode_no) + force_invhatTemp(2:2:end);
        
        aboveNodes = unique([aboveNodes; aboveNodesTemp]);
        aboveNodesMatrix{islpOut,1} = aboveNodes;

    end
    
end
stress_invhat(:,4) = stress_invhat(:,3);


%%

%  figure; hold on
% plot(above_box(:,1),above_box(:,2), 'r-' , 'Linewidth', 1)
% plot(slplaneOut(:,1), slplaneOut(:,2), 'y-', 'Linewidth', 2)
% 
%     plot([x1 x2 x3 x1],[y1 y2 y3 y1],'k-','LineWidth',1);
%   
% 
% %%
% 
%     plot(p([eleNode_no,eleNode_no(1)],1),p([eleNode_no,eleNode_no(1)],2),'b-','LineWidth',1);
% 
%     %%
%     
%     figure; hold on
%    
%         plot3(p([eleNode_no,eleNode_no(1)],1),p([eleNode_no,eleNode_no(1)],2),[1 0 1 1],'b-','LineWidth',1);
% 
%                 plot3(p([eleNode_no,eleNode_no(1)],1),p([eleNode_no,eleNode_no(1)],2),[0 0 0 0],'k--','LineWidth',1);
% xtry = p([eleNode_no],1);
% ytry = p([eleNode_no],2);
% ztry = [1 0 1];

% plot3(ans(:,1), ans(:,2), [1 0 1 1], 'k-', 'Linewidth', 2)
% hold on
% plot3(ans(:,1), ans(:,2), [0 0 0 0], 'r-', 'Linewidth', 2)