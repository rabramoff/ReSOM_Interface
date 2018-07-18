function dxdt=many_bug_many_substrate(t,x)
%
%DESCRIPTION
%core file for microbial model with many bugs attacking many substrates
%Author: Jinyun Tang, jinyuntang@lbl.gov
%Created Feb, 9, 2014.
%currently, I only deal with carbon substrate, so the first step carbon
%yield is well defined for each substrate by each microbe
%
%All vectors in this code are columnwise

%input arguments
%t: time
%x: state variable vector

global vid;        %id index
global par_mic;    %each par represents a microbe and its enzyme
global par_enz;    %each par represents a enzyme
global par_surface;%each par represents one mineral surface
global par_ss;     %parameter structure for substrate
global input;      %substrate input structure

%obtain number of pools
nmicrobes=length(vid.microbep);
nenzymes=length(vid.enzymes);
nsurfaces=length(vid.surfaces);
nmonomers=length(vid.monomers);
npolymers=length(vid.polymers);

%initialize the trend vector
dxdt=zeros(size(x));

%the substrate-microbe relationship is represented as a matrix with each
%column representing substrates, and each row representing microbes (or
%mineral surface) the matrixces are K parameters, Vmax parameters and 
%substrate yielding rates polymer degradation needs extracellular enzyme,
%which is produced by producers. In this model, mineral surfaces only 
%adsorb enzymes and monomers but do not actively protect them through
%aggregation.

% Obtain polymer degradation flux (by enzyme, npolymers x nenzymes) and enzyme adsorption flux (by mineral, nsurfaces x nenzymes)
[polymer_degrad_flux,enzyme_adsorp_flux]=polymer_degradation(x(vid.polymers), x(vid.enzymes), ...
    x(vid.surfaces), par_enz);

% Obtain monomer uptake flux (by microbe, nmonomers x nmicrobes) and monomer adsorption flux (by mineral, nmonomers x nsurfaces)
[monomer_uptake_matrix,monomer_adsorp_flux,enz_pro_stress_tot]=monomer_uptake(x(vid.monomers), ...
    x(vid.microbep), x(vid.surfaces), par_mic, par_surface);

% turnover of adsorbed enzyme pool
adsenz_decay_flux=zeros(nenzymes,nsurfaces);
for ii=1:nenzymes
    for jj=1:nsurfaces
        adsenz_decay_flux(ii,jj)=par_enz(ii).adsenz_decay(jj)*x(vid.enzymes_ads((ii-1)*nsurfaces+jj));
    end
end

% turnover of adsorbed monomer pool
adsmon_decay_flux=zeros(nmonomers,nsurfaces);
for ii=1:nmonomers
    for jj=1:nsurfaces
        adsmon_decay_flux(ii,jj)=par_surface(jj).adsmon_decay(ii)*x(vid.monomers_ads((ii-1)*nsurfaces+jj));
    end
end

%now assign the effective carbon yield flux to each microbes, nmonomers x nmicrobes
[monomer_to_microbe_flux, assim_co2_flux]=monomer_microbe_yield(monomer_uptake_matrix(:,1:nmicrobes),par_mic,x);

%now do cell level metabolism to find rates of maintenance, population
%growth and enzyme production
[maintenance, growth, pro_enz, cell_death, co2_metab_flux, uptake]=microbe_metabolism(monomer_to_microbe_flux,...
    enz_pro_stress_tot, x, vid, par_mic);

%calculate CUE
for jj = 1 : nmicrobes
    cue(jj) = 1 - (co2_metab_flux(jj)./uptake(jj));
end

%do mass update

%microbial population
for jj = 1 : nmicrobes
    if(plasticity_test(par_mic(jj))) % need change
        %rigid microbe
        dxdt(vid.microbep(jj))=(growth(jj)-...
            cell_death(jj))*x(vid.microbep(jj)); % !! add 13C and 14C stuff for rigid microbe
    else
        dxdt(vid.microbep(jj))=(growth(jj)-cell_death(jj))*x(vid.microbep(jj));     
        dxdt(vid.micc(jj))=sum(monomer_to_microbe_flux(:,jj),1)-(par_mic(jj).kappa_micb...
            -growth(jj)+cell_death(jj))*x(vid.micc(jj));        
    end
end

% enzyme production matrix
pro_enz_bulk=pro_enz.*(x(vid.microbep))'; %bulk enzyme production
pro_enz_matrix=zeros(nenzymes,nmicrobes);
for jj=1:nmicrobes
    pro_enz_matrix(:,jj)=pro_enz_bulk(jj)*par_mic(jj).pro_enz_dist';
end

%bulk microbial mortality
% cell_death_bulk=cell_death.*(x(vid.microbep)+x(vid.micc))';
cell_death_bulk_microbep=cell_death.*(x(vid.microbep))';
cell_death_bulk_micc=cell_death.*(x(vid.micc))';

% update extracellular enzyme
[dxdt(vid.enzymes),enz_decay]=enzyme_update(x,pro_enz_matrix,vid, par_enz,enzyme_adsorp_flux,adsenz_decay_flux);

% update adsorbed enzymes by mineral surfaces
dxdt(vid.enzymes_ads)=reshape(enzyme_adsorp_flux-adsenz_decay_flux',nsurfaces*nenzymes,1);

% update polymers
dxdt(vid.polymers)=input.polymers ...
    -sum(polymer_degrad_flux,2)' ...
    +(par_ss.deadmicrobep2polymers*cell_death_bulk_microbep)' ...
    +(par_ss.deadmicc2polymers*cell_death_bulk_micc)' ...
    +(par_ss.deadenz2polymers*enz_decay)';

% update monomers
dxdt(vid.monomers)=input.monomers ...
    +(par_ss.polymer2monomer*sum(polymer_degrad_flux,2))' ...
    -sum(monomer_uptake_matrix,2)' ...
    +(par_ss.deadmicrobep2monomers*cell_death_bulk_microbep)' ...
    +(par_ss.deadmicc2monomers*cell_death_bulk_micc)' ...
    +(par_ss.deadenz2monomers*enz_decay)' ...
    -sum(monomer_adsorp_flux,2)' ...
    +sum(adsmon_decay_flux,2)';

% update adsorped monomers
dxdt(vid.monomers_ads)=reshape(monomer_adsorp_flux'-adsmon_decay_flux',nsurfaces*nmonomers,1);

% update mineral surfaces (effective binding surfaces)
dxdt(vid.surfaces)= sum(adsmon_decay_flux,1) ...
    + sum(adsenz_decay_flux,1) ...
    -sum(enzyme_adsorp_flux,2)' ...
    -sum(monomer_adsorp_flux,1);

% update co2
dxdt(vid.co2)=sum(assim_co2_flux(:))+sum(co2_metab_flux,2);

 for jj = 1 : nmicrobes
% update cue
         dxdt(vid.cue(jj))=cue(jj);
 end

% do mass balance check
[residual]=mass_balance_check(dxdt,x);

if(abs(residual)>1d-10)
    fprintf('residual=%f\n',residual);
    error('mass balance error');
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [residual]=mass_balance_check(dxdt,x)
%do mass balance check
    global input vid

    idx=[vid.microbep vid.micc vid.monomers vid.monomers_ads vid.polymers vid.enzymes vid.enzymes_ads vid.co2];
    residual=sum(dxdt(idx))-(sum(input.polymers)+sum(input.monomers));
    
end





