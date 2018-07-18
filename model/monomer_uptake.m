function [monomer_uptake_matrix,monomer_adsorp_flux,enzprostress]=monomer_uptake(monomers, ...
    microbes, surfaces, par_mic, par_surface)
%usage
%
%do microbe and surface competition of monomers
%
%input argument
%monomers: vector of monomers
%microbes: vector of microbes
%surface: vector of mineral surfaces
%par_mic: par_micameter structure
%
%return variable
%

%%%%%%%%%%%


%obtain the number of monomers
nmonomers=length(monomers);

%obtain the number of microbes
nmicrobes=length(microbes);

%compute the amount of effective adsorption surfaces
eff_surfaces=[microbes.*[par_mic.zb],surfaces];

%number of effective surfaces competing for DOC
nsurfaces = length(surfaces);

%obtain the K_matrix
K_matrix=clac_Kmat_Microbe_monomer_surface(par_mic, par_surface, nmonomers, nmicrobes, nsurfaces);

%obtain the Vmax_matrix
Vmax_matrix=clac_Vmax_microbe_monomer_surface(par_mic, par_surface, nmonomers, nmicrobes, nsurfaces);

%obtain the monomer complex
monomer_complex_matrix=ECAcomplex(K_matrix,monomers,eff_surfaces);

%monomer2enzstress is a matrix of size, enz x monomers
[enzprostress_comp,enzprostress]=cal_enz_pro_stress(monomer_complex_matrix(1:nmonomers,1:nmicrobes),...
    eff_surfaces(1:nmicrobes));

%form the monomer uptake matrix which will be used to obtain monomer flux
%into the microbial cells. nmonomers x nmicrobes
monomer_uptake_matrix=Vmax_matrix(:,1:nmicrobes).*monomer_complex_matrix(:,1:nmicrobes);

% get monomer adsorption flux by mineral, nmonomers x nsurfaces
monomer_adsorp_flux=Vmax_matrix(:,nmicrobes+1:nmicrobes+nsurfaces).*monomer_complex_matrix(:,nmicrobes+1:nmicrobes+nsurfaces);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function K_matrix=clac_Kmat_Microbe_monomer_surface(par_mic, par_surface, nmonomers, nmicrobes, nsurfaces)
%
%usage
%compute the affinity K for doc binding
%The K matrix is of the size, #substrates x #surfaces
%here the substrate surfaces include both monomer surfaces and mineral
%surfaces
%K parameter from microbes
K_matrix=zeros(nmonomers, nmicrobes + nsurfaces);
for jj = 1 : nmicrobes
    K_matrix(:,jj)=par_mic(jj).Kaff_monomer';
end
%K parameter for mineral surfaces
for jj = 1 : nsurfaces
    kk=jj+nmicrobes;
    K_matrix(:,kk)=par_surface(jj).Kaff_monomer';
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Vmax_matrix=clac_Vmax_microbe_monomer_surface(par_mic, par_surface, nmonomers, nmicrobes, nsurfaces)
%
%usage
%compute the Vmax par_micameter for enzymatic polymer degradation
%
Vmax_matrix=zeros(nmonomers,nmicrobes + nsurfaces);
%Vmax for microbes
for jj = 1 : nmicrobes
    Vmax_matrix(:,jj)=par_mic(jj).Vmax_micb';
end

%Vmax for mineral surfaces
for jj = 1 : nsurfaces
    kk=jj+nmicrobes;
    Vmax_matrix(:,kk)=par_surface(jj).Vmax_ads_monomer';
end
end
