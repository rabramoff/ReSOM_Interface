function [monomer_to_microbe_flux, co2_assim_flux]=monomer_microbe_yield(monomer_uptake_matrix, par_mic, x)
%
%usage
%compute the first step yield rate for each microbe. The first step yield
%may be quite independent from the type of microbes
%author: Jinyun Tang
%input argument
%monomer_uptake_matrix: monomer_microbe binding matrix, nmonomers x
%
%return variable
%monomer_to_Microbe_flux: 
%co2_assim_flux: co2 flux due to monomer assimilation

%get number of monomers
nmonomers=size(monomer_uptake_matrix,1);

%get substrate yield matrix
Yx_matrix=calc_Ymat_microbe_monomer(nmonomers, par_mic);

%get the incoming c flux to microbes
% monomer_to_microbe_flux=ones(1,nmonomers)*(monomer_uptake_matrix.*Yx_matrix); % 1 x nmicrobes
monomer_to_microbe_flux=monomer_uptake_matrix.*Yx_matrix; % nmonomers x nmicrobes

%get the respired co2
%when Yx_matrix is mixed with both carbon and nutrient species, the
%following line should be modified accordingly
% co2_assim_flux=sum(sum(monomer_uptake_matrix))-sum(monomer_to_microbe_flux); % 1 x 1
co2_assim_flux=monomer_uptake_matrix-monomer_to_microbe_flux; % nmonomers x nmicrobes
end

function Yx_matrix=calc_Ymat_microbe_monomer(nmonomers, par_mic)
%
% DESCRIPTION
% calculate the monomer yield rate for different microbes
nmicrobes=length(par_mic);

Yx_matrix=zeros(nmonomers, nmicrobes);

%the following line will be replaced with thermodynamic based calculations
for jj = 1 : nmicrobes
    Yx_matrix(:,jj)=par_mic(jj).Yld_micx_monomer';
end

end