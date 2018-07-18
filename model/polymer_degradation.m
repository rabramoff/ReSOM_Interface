function [polymerflux,enzymeflux]=polymer_degradation(Polymers, Enzymes, Surfaces, par)
%
%Description
%
%do enzymatic polymer degradation and enzyme adsorption flux, enzymes are subject to mineral surface
%adsorption
%author: Jinyun Tang
%input argements
%Polymers: nx1 array
%Enzymes:  nx1 array
%Surfaces: nx1 array
%
%return variable
%polymerflux: bulk degradation rate for each polymer
%enzymeflux: enzyme adsorption rate by mineral for each enzyme

%%%%%%%%%%%

%local arguments
%K_matrix: (npolymers+nsurfaces)x(nenzymes)
%Vmax_matrix: (npolymers)x(nenzymes), because I have not tried to resolve
%all specific enzymatic degradations, some enzyme can degradate multiple
%polymers

%polymer degradation
%obtain the number of polymers
nPolymers=length(Polymers);

%obtain number of enzymes
nEnzymes=length(Enzymes);

%Form the augmented substrate matrix
substrates=[Polymers,Surfaces];

nSurfaces=length(Surfaces);

%Form the affinity K_matrix
K_matrix=calc_Kmat_Enzyme_polymer_surface(par, nSurfaces, nPolymers, nEnzymes);

%Form the degradation Vmax_matrix
Vmax_matrix=calc_Vmax_Enzyme_polymer_surface(par, nSurfaces, nPolymers, nEnzymes);

%Obtain the substrate-enzyme complex using ECA kinetics
polymer_complex_matrix=ECAcomplex(K_matrix,substrates,Enzymes);

%obtain the polymerflx, npolymers x nenzymes
% polymerflux=(polymer_complex_matrix(1:nPolymers,:).*Vmax_matrix)*ones(nEnzymes,1);
polymerflux=polymer_complex_matrix(1:nPolymers,:).*Vmax_matrix(1:nPolymers,:);

%Obtain enzymeflux, nsurfaces x nenzymes
enzymeflux=polymer_complex_matrix(nPolymers+1:nPolymers+nSurfaces,:).*Vmax_matrix(nPolymers+1:nPolymers+nSurfaces,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function K_matrix=calc_Kmat_Enzyme_polymer_surface(par, nSurfaces, nPolymers, nEnzymes)
%
%usage
%compute the affinity K for enzyme polymer binding
%The K matrix is of the size, #substrates x #enzymes
%here the substrate surfaces include both polymer surfaces and mineral
%surfaces
%
%because of the setup here, each enzyme has affinity parameters for
%different polymers and mineral surfaces, therefore
K_matrix=zeros(nPolymers+nSurfaces,nEnzymes);

for kk = 1 : nEnzymes
    K_matrix(:,kk)=[par(kk).Kaff_ee_polymer';par(kk).Kaff_ee_msurf'];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Vmax_matrix=calc_Vmax_Enzyme_polymer_surface(par, nSurfaces, nPolymers, nEnzymes)
%
%usage
%compute the Vmax parameter for enzymatic polymer degradation
%
Vmax_matrix=zeros(nPolymers+nSurfaces,nEnzymes);

for kk = 1 : nEnzymes
    Vmax_matrix(:,kk)=[par(kk).Vmax_ee';par(kk).Vmax_ads_enzyme'];
end
end