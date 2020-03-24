function d=dpolySW(p,pv) % p=nodes postions, pv=bbox 

%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

np=size(p,1); 
nvs=size(pv,1)-1; 

ds=dsegment(p,pv);
d=min(ds,[],2);

d=(-1).^(inpolygon(p(:,1),p(:,2),pv(:,1),pv(:,2))).*d;
d(abs(d)<eps)=0;
end

function d=ddiffSW(d1,d2), d=max(d1,-d2);
end