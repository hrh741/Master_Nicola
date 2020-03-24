nout = 1000;


for i=3:nout
    if rem(i,2)==0
        oneortwo = 2;
    else
        oneortwo = 1;
    end
        alphaOut(i) = alphaOut(oneortwo);
        bOut(i) = bOut(oneortwo);
        xdisOut(i) = xdisOut(oneortwo);
        ydisOut(i) = ydisOut(oneortwo);
        typeOut(i) = typeOut(oneortwo);
        ngsourceOut(i) = ngsourceOut(oneortwo);
        planeOut(i) = planeOut(oneortwo);
       
%         node_slip(:,i)=node_slip(:,1);
end

%%
z_plot0