function [temp0,T_fact00]=get_microbe_physiology_Tcurve(Tref0,xpar0)
%
%Description
%Calculate the physiology temperature curve, this assumes a key enzyme
%controls all temperature dependent physiological response.
%Author: Jinyun Tang

%reference trait parameters
if exist('xpar0')
    xpar0 = xpar0;
else
    xpar0=[249.544170969785   5341.422691388677   5.617549086429];
end

temp0=(220:340);
T_fact00=fact_murphy(temp0,xpar0);
T_fact0=fact_murphy(Tref0,xpar0);
T_fact00=T_fact00./T_fact0;


end