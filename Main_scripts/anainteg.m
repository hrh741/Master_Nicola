   
    function [I11, I22, I12] = anainteg(sdisall, sa, sp, t_deltaall, thetas_tild, type, b_all, delta)
    
    I11 = 0; I22 = 0; I12 = 0;
    
    
    for disno = 1:length(sdisall)
        
        sd = sdisall(1,disno);
        distype = type(disno);
        t_delta = t_deltaall(disno); 
        disb = b_all(disno);
        
        % if too close
        r2 = (sp-sd)^2 + t_delta^2;
        if r2<(delta*delta)
            if r2 == 0
                sd=sd+delta; t_delta = 0;
            else
                rinv = delta/sqrt(r2);
                sd = (sp-sd)*rinv + sd;
                t_delta = t_delta*rinv;
            end
        end
        
        I11_AtoP_1 = 2*cos(thetas_tild)*cos(thetas_tild)*( cos(thetas_tild)*(sd-sa) + 3*sin(thetas_tild)*t_delta )* ( atan( -(sp-sd)/t_delta ) - atan( -(sa-sd)/t_delta ) ) ;
        
        I11_AtoP_2 = 0.5*( 3*sin(thetas_tild) + sin(3*thetas_tild) ) * ( sp-sa ) ;
        
        I11_AtoP_3 = -t_delta* ( ( sd*sd + t_delta*t_delta + sp*sa - sd*(sp+sa) )*cos(3*thetas_tild) - t_delta*(sp-sa)*sin(3*thetas_tild) ) ...
            / ( (sp-sd)*(sp-sd) + t_delta*t_delta ) ...
            + t_delta* ( ( sd*sd + t_delta*t_delta + sa*sa - sd*(sa+sa) )*cos(3*thetas_tild) + 0 ) ... % t_delta*(sa-sa)*sin(3*thetas_tild)
            / ( (sa-sd)*(sa-sd) + t_delta*t_delta ) ;
        
        I11_AtoP_4 = 0.25* ( -3*t_delta*(cos(thetas_tild) + cos(3*thetas_tild)) + (sd-sa)*(3*sin(thetas_tild) + sin(3*thetas_tild)) ) ...
            * ( log( (sp-sd)*(sp-sd) + t_delta*t_delta ) - log( (sa-sd)*(sa-sd) + t_delta*t_delta ) );
        
        I11 = I11 + disb*distype*(I11_AtoP_1+I11_AtoP_2+I11_AtoP_3+I11_AtoP_4) / (sp-sa); % summed over dislocations on one slip system
        
        % ---
        
        I22_AtoP_1 = sin(thetas_tild)*( -t_delta*(1+3*cos(2*thetas_tild)) + (sd - sa)*sin(2*thetas_tild) ) * ( atan( -(sp-sd)/t_delta ) - atan( -(sa-sd)/t_delta ) ) ;
        
        I22_AtoP_2 = 0.25* ( -t_delta*(cos(thetas_tild)-3*cos(3*thetas_tild)) + (sd-sa)*(sin(thetas_tild)-sin(3*thetas_tild)) ) ...
            * ( log( (sp-sd)*(sp-sd) + t_delta*t_delta ) - log( (sa-sd)*(sa-sd) + t_delta*t_delta ) ) ;
        
        I22_AtoP_3 = 0.5*( sin(thetas_tild)-sin(3*thetas_tild) ) * (sp-sa) ;
        
        I22_AtoP_4 = t_delta* ( ( sd*sd + t_delta*t_delta + sp*sa - sd*(sp+sa) )*cos(3*thetas_tild) - t_delta*(sp-sa)*sin(3*thetas_tild) ) ...
            / ( (sp-sd)*(sp-sd) + t_delta*t_delta ) ...
            - t_delta* ( ( sd*sd + t_delta*t_delta + sa*sa - sd*(sa+sa) )*cos(3*thetas_tild) - 0 ) ... % t_delta*(sa-sa)*sin(3*thetas_tild)
            / ( (sa-sd)*(sa-sd) + t_delta*t_delta ) ;
        
        I22 = I22 + disb*distype*(I22_AtoP_1+I22_AtoP_2+I22_AtoP_3+I22_AtoP_4) / (sp-sa);
        
        % ---
               
        I12_AtoP_1 = 0.5* ( cos(thetas_tild) + cos(3*thetas_tild) ) * ( sp-sa );
        
        I12_AtoP_2 = cos(thetas_tild)* ( t_delta*(3*cos(2*thetas_tild)-1) - (sd - sa)*sin(2*thetas_tild) )...
            * ( atan( -(sp-sd)/t_delta ) - atan( -(sa-sd)/t_delta ) ) ;
        
        I12_AtoP_3 = t_delta* ( ( sd*sd + t_delta*t_delta + sp*sa - sd*(sp+sa) )*sin(3*thetas_tild) + t_delta*(sp-sa)*cos(3*thetas_tild) ) ...
            / ( (sp-sd)*(sp-sd) + t_delta*t_delta ) ...
            - t_delta* ( ( sd*sd + t_delta*t_delta + sa*sa - sd*(sa+sa) )*sin(3*thetas_tild) + 0 ) ... % t_delta*(sa-sa)*cos(3*thetas_tild)
            / ( (sa-sd)*(sa-sd) + t_delta*t_delta ) ;
        
        I12_AtoP_4 = 0.25* ( (sd-sa)*(cos(thetas_tild)+cos(3*thetas_tild)) + t_delta*(sin(thetas_tild) + 3*sin(3*thetas_tild)) ) ...
            * ( log( (sp-sd)*(sp-sd) + t_delta*t_delta ) - log( (sa-sd)*(sa-sd) + t_delta*t_delta ) ) ;
        
        I12 = I12 + disb*distype*(I12_AtoP_1+I12_AtoP_2+I12_AtoP_3+I12_AtoP_4) / (sp-sa);
        
        if any(isnan(I11) || isnan(I12) || isnan(I22))
            disno
            sa
            sp
            sd
            t_delta
            disp('error'); pause;
        end
        
    end     
        
    end
