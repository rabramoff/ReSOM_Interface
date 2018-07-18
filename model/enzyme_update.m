function [enz_flux,enz_decay]=enzyme_update(x,pro_enz_matrix,vid,par,enzyme_adsorp_flux,adsenz_decay_flux)
%usage
%update extracellular enzyme that is potentially active
%
%input arguments
%x: state variables
%pro_enz_matrix: enzyme production rates, nenzymes x nmicrobes
%vid: id structure
%par: parameter structure
%enzyme_adsorp_flux: enzyme adsorption flux by mineral, nsurfaces x nenzymes
%adsenz_decay_flux: nenzymes x nsurfaces

nenzymes=length(vid.enzymes);
enz_decay=zeros(nenzymes,1);

for kk = 1 : nenzymes
    enz_decay(kk)=par(kk).decay_ee*x(vid.enzymes(kk));
end
%the net enzyme flux
enz_flux=sum(pro_enz_matrix,2)-enz_decay-sum(enzyme_adsorp_flux,1)'+sum(adsenz_decay_flux,2); % nenzymes x 1
end