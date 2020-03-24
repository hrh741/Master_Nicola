% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Function    : dislocation transmission through GBs
%   Description : call by Main.m
%                     p_connect
% {                       |                       }
% {   1     2     3    4  |   5     6     7    8  }
% {                       |                       }
% { plane grain alpha end | plane grain alpha end }
% {                       |                       }
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [destroy,ndis,xdis,ydis,rdis,type,alpha,ngsource,b,irmbound,...
    pinned,eta_obs,source,vdis,vdispre,plane] ...
    = distrans(ndis,xdis,ydis,rdis,type,alpha,ngsource,b,irmbound,...
    p_connect,bG,mu,p_xend,p_yend,p_rend,pinned,eta_obs,...
    source,vdis,vdispre,plane,reps,logFID,taudis)

destroy = zeros(ndis,1);
xalpha = 1;

for i = 1:ndis
    if irmbound(i)==1
        if taudis(i)*type(i)>0
            irmbound(i) = 0;
            rdis(i) = rdis(i)+1e-19;
        elseif ~isempty(p_connect{plane(i),1}) && taudis(i)*type(i)<0
            % calculate tau_pass >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        % Step 1 : determine the out-going slip plane
            % Criterion 1 : find minimum delta_b
            db = [];
            dtheta = [];
            typetemp = [];
            vecb = b(i)*type(i)*[cos(alpha(i)),sin(alpha(i))];
            for j = 1:length(p_connect{plane(i),1})
                if dot(vecb,[cos(p_connect{plane(i),3}(j)),sin(p_connect{plane(i),3}(j))])>0
                    typetemp = [typetemp;1];
                else
                    typetemp = [typetemp;-1];
                end
                db(j) = norm(bG(p_connect{plane(i),2}(j))*typetemp(j)*[cos(p_connect{plane(i),3}(j)),sin(p_connect{plane(i),3}(j))]-vecb);
                dtheta(j) = abs(p_connect{plane(i),3}(j)-alpha(i));
            end
            m = min(db);
            n = find(db==m);
            % Criterion 2 : if more than one slip plane give the minimum
            % delta_b, find minimum distance between two intersections
            if length(n)>1
                dintsctn = [];
                for j = 1:length(n)
                    dintsctn(j) = pdist2([p_xend(plane(i),1),p_yend(plane(i),1)],...
                        [p_xend(p_connect{plane(i),1}(n(j)),p_connect{plane(i),4}(n(j))),...
                        p_yend(p_connect{plane(i),1}(n(j)),p_connect{plane(i),4}(n(j)))]);
                end
                [mm, nn] = min(dintsctn); % if more than one minimum value, automatically chose the first one.
            else
                nn = 1;
            end
        % Step 2 : calculate tau_pass
            j_index = n(nn);
            E_gb = energyGB(dtheta(j_index),3);
            b_find = p_connect{plane(i),2}(j_index);
            tau_pass = (E_gb*bG(b_find)*1e-9 + xalpha*mu(b_find)*(db(j_index))^2)/bG(b_find)/bG(b_find);
        % Step 3 : determine whether dislocation penetration happened
            if taudis(i)*type(i) <= -tau_pass && ~any(plane==p_connect{plane(i),1}(j_index) & irmbound==p_connect{plane(i),4}(j_index))
            % slip transfer happened
            fprintf(logFID,'slip transfer happened on slip plane: %d in grain: %d .\n',plane(i),ngsource(i));  
                planeold = plane(i);
                alphaold = alpha(i);
                grainold = ngsource(i);
                endold = irmbound(i);
                b(i) = bG(b_find);
                xdis(i) = p_xend(p_connect{plane(i),1}(j_index),p_connect{plane(i),4}(j_index));
                ydis(i) = p_yend(p_connect{plane(i),1}(j_index),p_connect{plane(i),4}(j_index));
                rdis(i) = p_rend(p_connect{plane(i),1}(j_index),p_connect{plane(i),4}(j_index));
                alpha(i) = p_connect{plane(i),3}(j_index);
                ngsource(i) = p_connect{plane(i),2}(j_index);
                source(i) = -2; % transmitted dis flag
                type(i) = typetemp(j_index);
                % pinned(i) = pinned(i); eta_obs(i) = eta_obs(i);
                % vdispre(i) = vdispre(i); source(i) = source(i);
                irmbound(i) = p_connect{plane(i),4}(j_index);
                plane(i) = p_connect{plane(i),1}(j_index); % plane has to update last
                
                % deal with dislocation debris
                xdistemp = (xdis(i)+p_xend(planeold,1))/2;
                ydistemp = (ydis(i)+p_yend(planeold,1))/2;
                mm = find(xdis(1:ndis)>xdistemp-reps & xdis(1:ndis)<xdistemp+reps & ydis(1:ndis)>ydistemp-reps & ydis(1:ndis)<ydistemp+reps & irmbound(1:ndis)<0);
                if length(mm)>1; error('mm'); end
                db_vec = vecb - b(i)*type(i)*[cos(alpha(i)),sin(alpha(i))];
                if ~isempty(mm) % add to old debris
                    db_vec = db_vec+b(mm)*type(mm)*[cos(alpha(mm)),sin(alpha(mm))];
                    alpha(mm) = atan(db_vec(2)/db_vec(1));
                    if alpha(mm)<0; alpha(mm) = alpha(mm)+pi; end
                    b(mm) = norm(db_vec);
                    if b(mm)==0
                        destroy(mm) = 1;
                    end
                    if db_vec(1)>0
                        type(mm) = 1;
                    elseif db_vec(1)<0
                        type(mm) = -1;
                    else
                        if db_vec(2)>0
                            type(mm) = 1;
                        else
                            type(mm) = -1;
                        end
                    end
                    if b(mm)>=max(bG)
                        % dislocation re-emission
                        planetemp = [planeold,p_connect{planeold,1}];
                        alphatemp = [alphaold,p_connect{planeold,3}];
                        graintemp = [grainold,p_connect{planeold,2}];
                        endtemp = [endold,p_connect{planeold,4}];
                        p_test = []; Edelta = []; newdb = []; newout = [];
                        for iplane = 1:length(planetemp)
                            btemp = bG(graintemp(iplane));
                            vectprojtemp = dot(db_vec,[cos(alphatemp(iplane)),sin(alphatemp(iplane))]);
                            if abs(vectprojtemp)>=btemp
                                p_test = [p_test; iplane];
                                vecttemp = db_vec-btemp*[cos(alphatemp(iplane)),sin(alphatemp(iplane))]*sign(vectprojtemp);
                                newdb = [newdb; vecttemp];
                                newout = [newout; btemp*[cos(alphatemp(iplane)),sin(alphatemp(iplane))]*sign(vectprojtemp)];
                            end
                        end
                        for ip_test = 1:length(p_test)
                            btemp = norm(newout(ip_test,:));
                            Edelta(ip_test) = xalpha*mu(-ngsource(mm))*btemp*btemp+xalpha*mu(-ngsource(mm))*(norm(newdb(ip_test,:)))^2;
                            dis_pout = find(plane(1:ndis)==planetemp(p_test(ip_test)));
%                             for idis_pout = 1:length(dis_pout)
%                                 Edelta(ip_test) = Edelta(ip_test) + edgeType(-ngsource(mm))*btemp*btemp*log(max(max(bbox))/pdist([xdis(dis_pout(idis_pout)),ydis(dis_pout(idis_pout));xdis(mm),ydis(mm)]));
%                             end
                            if xalpha*mu(-ngsource(mm))*b(mm)*b(mm) >= Edelta(ip_test) % if more than one slip plane meets the criterion, only chose the first one
                                % deal with old debris
                                b(mm) = norm(newdb(ip_test,:));
                                alpha(mm) = atan(newdb(ip_test,2)/newdb(ip_test,1));
                                if alpha(mm)<0; alpha(mm) = alpha(mm)+pi; end
                                if b(mm)==0
                                    destroy(mm) = 1;
                                end
                                if newdb(ip_test,1)>0
                                    type(mm) = 1;
                                elseif newdb(ip_test,1)<0
                                    type(mm) = -1;
                                else
                                    if newdb(ip_test,2)>0
                                        type(mm) = 1;
                                    else
                                        type(mm) = -1;
                                    end
                                end
                                
                                % add the re-emmision dislocation
                                ndis = ndis + 1;
                                b(ndis) = btemp;
                                plane(ndis) = planetemp(p_test(ip_test));
                                xdis(ndis) = p_xend(plane(ndis),endtemp(p_test(ip_test)));
                                ydis(ndis) = p_yend(plane(ndis),endtemp(p_test(ip_test)));
                                rdis(ndis) = p_rend(plane(ndis),endtemp(p_test(ip_test))); % check no affect
                                alpha(ndis) = atan(newout(ip_test,2)/newout(ip_test,1));
                                if alpha(ndis)<0; alpha(ndis) = alpha(ndis)+pi; end
                                pinned(ndis) = 0;
                                irmbound(ndis) = 0;
                                source(ndis) = -1;
                                ngsource(ndis) = graintemp(p_test(ip_test));
                                eta_obs(ndis) = 0;
								vdis(ndis,1) = 0;
                                vdispre(ndis,1) = 0;
                                if newout(ip_test,1)>0
                                    type(ndis) = 1;
                                elseif newout(ip_test,1)<0
                                    type(ndis) = -1;
                                else
                                    if newout(ip_test,2)>0
                                        type(ndis) = 1;
                                    else
                                        type(ndis) = -1;
                                    end
                                end
                                fprintf(logFID,'dislocation re-emission happened to slip plane: %d .\n',plane(ndis));  
                                break;
                            end
                        end
                    end    
                else % create new debris
                    ndis = ndis + 1;
                    b(ndis) = db(j_index);
                    if b(ndis)==0
                        destroy(ndis) = 1;
                    end
                    xdis(ndis) = xdistemp;
                    ydis(ndis) = ydistemp;
                    rdis(ndis) = 0; % check no affect
                    plane(ndis) = -planeold;
                    alpha(ndis) = atan(db_vec(2)/db_vec(1));
                    if alpha(ndis)<0; alpha(ndis) = alpha(ndis)+pi; end
                    pinned(ndis) = 0;
                    irmbound(ndis) = -1;
                    source(ndis) = 0;
                    ngsource(ndis) = -ngsource(i);
                    eta_obs(ndis) = 0;
					vdis(ndis,1) = 0;
                    vdispre(ndis,1) = 0;
                    if db_vec(1)>0
                        type(ndis) = 1;
                    elseif db_vec(1)<0
                        type(ndis) = -1;
                    else
                        if db_vec(2)>0
                            type(ndis) = 1;
                        else
                            type(ndis) = -1;
                        end
                    end
                end
            end
        end
    elseif irmbound(i)==2
        if taudis(i)*type(i)<0
            irmbound(i) = 0;
            rdis(i) = rdis(i)-1e-19;
        elseif ~isempty(p_connect{plane(i),5}) && taudis(i)*type(i)>0
            % calculate tau_pass >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        % Step 1 : determine the out-going slip plane
            % Criterion 1 : find minimum delta_b
            db = [];
            dtheta = [];
            typetemp = [];
            vecb = b(i)*type(i)*[cos(alpha(i)),sin(alpha(i))];
            for j = 1:length(p_connect{plane(i),5})
                if dot(vecb,[cos(p_connect{plane(i),7}(j)),sin(p_connect{plane(i),7}(j))])>0
                    typetemp = [typetemp;1];
                else
                    typetemp = [typetemp;-1];
                end
                db(j) = norm(bG(p_connect{plane(i),6}(j))*type(i)*[cos(p_connect{plane(i),7}(j)),sin(p_connect{plane(i),7}(j))]-vecb);
                dtheta(j) = abs(p_connect{plane(i),7}(j)-alpha(i));
            end
            m = min(db);
            n = find(db==m);
            % Criterion 2 : if more than one slip plane give the minimum
            % delta_b, find minimum distance between two intersections
            if length(n)>1
                dintsctn = [];
                for j = 1:length(n)
                    dintsctn(j) = pdist2([p_xend(plane(i),2),p_yend(plane(i),2)],...
                        [p_xend(p_connect{plane(i),5}(n(j)),p_connect{plane(i),8}(n(j))),...
                        p_yend(p_connect{plane(i),5}(n(j)),p_connect{plane(i),8}(n(j)))]);
                end
                [mm, nn] = min(dintsctn); % if more than one minimum value, automatically chose the first one.
            else
                nn = 1;
            end
        % Step 2 : calculate tau_pass
            j_index = n(nn);
            E_gb = energyGB(dtheta(j_index),3);
            b_find = p_connect{plane(i),6}(j_index);
            tau_pass = (E_gb*bG(b_find)*1e-9 + xalpha*mu(b_find)*(db(j_index))^2)/bG(b_find)/bG(b_find);
        % Step 3 : determine whether dislocation penetration happened
            if taudis(i)*type(i) >= tau_pass && ~any(plane==p_connect{plane(i),5}(j_index) & irmbound==p_connect{plane(i),8}(j_index))
            % slip transfer happened
            fprintf(logFID,'slip transfer happened on slip plane: %d in grain: %d .\n',plane(i),ngsource(i));  
                planeold = plane(i);
                alphaold = alpha(i);
                grainold = ngsource(i);
                endold = irmbound(i);
                b(i) = bG(b_find);
                xdis(i) = p_xend(p_connect{plane(i),5}(j_index),p_connect{plane(i),8}(j_index));
                ydis(i) = p_yend(p_connect{plane(i),5}(j_index),p_connect{plane(i),8}(j_index));
                rdis(i) = p_rend(p_connect{plane(i),5}(j_index),p_connect{plane(i),8}(j_index));
                alpha(i) = p_connect{plane(i),7}(j_index);
                ngsource(i) = p_connect{plane(i),6}(j_index);
                source(i) = -2; % transmitted dis flag
                type(i) = typetemp(j_index);
                % pinned(i) = pinned(i); eta_obs(i) = eta_obs(i);
                % vdispre(i) = vdispre(i); source(i) = source(i);
                irmbound(i) = p_connect{plane(i),8}(j_index);
                plane(i) = p_connect{plane(i),5}(j_index); % plane has to update last
                
                % deal with dislocation debris
                xdistemp = (xdis(i)+p_xend(planeold,2))/2;
                ydistemp = (ydis(i)+p_yend(planeold,2))/2;
                mm = find(xdis(1:ndis)>xdistemp-reps & xdis(1:ndis)<xdistemp+reps & ydis(1:ndis)>ydistemp-reps & ydis(1:ndis)<ydistemp+reps & irmbound(1:ndis)<0);
                if length(mm)>1; error('mm'); end
                db_vec = vecb - b(i)*type(i)*[cos(alpha(i)),sin(alpha(i))];
                if ~isempty(mm) % add to old debris
                    db_vec = db_vec+b(mm)*type(mm)*[cos(alpha(mm)),sin(alpha(mm))];
                    alpha(mm) = atan(db_vec(2)/db_vec(1));
                    if alpha(mm)<0; alpha(mm) = alpha(mm)+pi; end
                    b(mm) = norm(db_vec);
                    if b(mm)==0
                        destroy(mm) = 1;
                    end
                    if db_vec(1)>0
                        type(mm) = 1;
                    elseif db_vec(1)<0
                        type(mm) = -1;
                    else
                        if db_vec(2)>0
                            type(mm) = 1;
                        else
                            type(mm) = -1;
                        end
                    end
                    if b(mm)>=max(bG)
                        % dislocation re-emission
                        planetemp = [planeold,p_connect{planeold,5}];
                        alphatemp = [alphaold,p_connect{planeold,7}];
                        graintemp = [grainold,p_connect{planeold,6}];
                        endtemp = [endold,p_connect{planeold,8}];
                        p_test = []; Edelta = []; newdb = []; newout = [];
                        for iplane = 1:length(planetemp)
                            btemp = bG(graintemp(iplane));
                            vectprojtemp = dot(db_vec,[cos(alphatemp(iplane)),sin(alphatemp(iplane))]);
                            if abs(vectprojtemp)>=btemp
                                p_test = [p_test; iplane];
                                vecttemp = db_vec-btemp*[cos(alphatemp(iplane)),sin(alphatemp(iplane))]*sign(vectprojtemp);
                                newdb = [newdb; vecttemp];
                                newout = [newout; btemp*[cos(alphatemp(iplane)),sin(alphatemp(iplane))]*sign(vectprojtemp)];
                            end
                        end
                        for ip_test = 1:length(p_test)
                            btemp = norm(newout(ip_test,:));
                            Edelta(ip_test) = xalpha*mu(-ngsource(mm))*btemp*btemp+xalpha*mu(-ngsource(mm))*(norm(newdb(ip_test,:)))^2;
                            dis_pout = find(plane(1:ndis)==planetemp(p_test(ip_test)));
%                             for idis_pout = 1:length(dis_pout)
%                                 Edelta(ip_test) = Edelta(ip_test) + edgeType(-ngsource(mm))*btemp*btemp*log(max(max(bbox))/pdist([xdis(dis_pout(idis_pout)),ydis(dis_pout(idis_pout));xdis(mm),ydis(mm)]));
%                             end
                            if xalpha*mu(-ngsource(mm))*b(mm)*b(mm) >= Edelta(ip_test) % if more than one slip plane meets the criterion, only chose the first one
                                % deal with old debris
                                b(mm) = norm(newdb(ip_test,:));
                                alpha(mm) = atan(newdb(ip_test,2)/newdb(ip_test,1));
                                if alpha(mm)<0; alpha(mm) = alpha(mm)+pi; end
                                if b(mm)==0
                                    destroy(mm) = 1;
                                end
                                if newdb(ip_test,1)>0
                                    type(mm) = 1;
                                elseif newdb(ip_test,1)<0
                                    type(mm) = -1;
                                else
                                    if newdb(ip_test,2)>0
                                        type(mm) = 1;
                                    else
                                        type(mm) = -1;
                                    end
                                end
                                
                                % add the re-emmision dislocation
                                ndis = ndis + 1;
                                b(ndis) = btemp;
                                plane(ndis) = planetemp(p_test(ip_test));
                                xdis(ndis) = p_xend(plane(ndis),endtemp(p_test(ip_test)));
                                ydis(ndis) = p_yend(plane(ndis),endtemp(p_test(ip_test)));
                                rdis(ndis) = p_rend(plane(ndis),endtemp(p_test(ip_test))); % check no affect
                                alpha(ndis) = atan(newout(ip_test,2)/newout(ip_test,1));
                                if alpha(ndis)<0; alpha(ndis) = alpha(ndis)+pi; end
                                pinned(ndis) = 0;
                                irmbound(ndis) = 0;
                                source(ndis) = -1;
                                ngsource(ndis) = graintemp(p_test(ip_test));
                                eta_obs(ndis) = 0;
								vdis(ndis,1) = 0;
                                vdispre(ndis,1) = 0;
                                if newout(ip_test,1)>0
                                    type(ndis) = 1;
                                elseif newout(ip_test,1)<0
                                    type(ndis) = -1;
                                else
                                    if newout(ip_test,2)>0
                                        type(ndis) = 1;
                                    else
                                        type(ndis) = -1;
                                    end
                                end

                                fprintf(logFID,'dislocation re-emission happened to slip plane: %d .\n',plane(ndis));
                                break;
                            end
                        end
                    end    
                else % create new debris
                    ndis = ndis + 1;
                    b(ndis) = db(j_index);
                    if b(ndis)==0
                        destroy(ndis) = 1;
                    end
                    xdis(ndis) = xdistemp;
                    ydis(ndis) = ydistemp;
                    rdis(ndis) = 0; % check no affect
                    plane(ndis) = -planeold;
                    alpha(ndis) = atan(db_vec(2)/db_vec(1));
                    if alpha(ndis)<0; alpha(ndis) = alpha(ndis)+pi; end
                    pinned(ndis) = 0;
                    irmbound(ndis) = -2;
                    source(ndis) = 0;
                    ngsource(ndis) = -ngsource(i);
                    eta_obs(ndis) = 0;
					vdis(ndis,1) = 0;
                    vdispre(ndis,1) = 0;
                    if db_vec(1)>0
                        type(ndis) = 1;
                    elseif db_vec(1)<0
                        type(ndis) = -1;
                    else
                        if db_vec(2)>0
                            type(ndis) = 1;
                        else
                            type(ndis) = -1;
                        end
                    end
                end
            end
        end
    end
end

function E_gb = energyGB(alpha,hardg)
% while alpha>=pi/3
%     alpha = alpha-pi/3;
% end
% alpha = alpha/pi*180;
% if hardg == 3 % no hard grain involved
%     if alpha<=46.4
%         E_gb = -0.5548*(alpha^2)+41.977*alpha;
%     else
%         E_gb = -56.912*alpha+3414.7;
%     end
% else
%     if alpha<20
%         E_gb = 39*alpha+100;
%     elseif alpha>=20 && alpha<40
%         E_gb = 880;
%     else
%         E_gb = -39*alpha+2440; 
%     end
% end
E_gb = 1550; % mJ m^-2

% function [Lnorm, proj, theta] = vectorproj(vector1, vector2)
% % project vector1 to vector2
% % vector2 should be unit vector
% proj = abs(dot(vector1,vector2));
% Lnorm = norm(vector1-proj*vector2);
% theta = atan(Lnorm/proj); % 0<=theta<=90 

