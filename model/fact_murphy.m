function fact=fact_murphy(temp,enz_par)

%enz_par(:,1)
%enz_par(:,2)
%enz_par(:,3)

x3=-46.0+30.*(1-1.54.*enz_par(:,1).^(-0.268)).*enz_par(:,3);
fact=1./(1+exp(-enz_par(:,1).*(enz_par(:,2)-...
    18.1.*temp+x3.*(temp-373.6-temp.*log(temp./385.2)))./(8.314.*temp)));

end