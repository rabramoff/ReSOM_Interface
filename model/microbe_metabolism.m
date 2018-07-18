function [maintenance_specific, growth_specific, pro_enz_specific, ...
    death_specific, co2_metab_flux_bulk, uptake]=microbe_metabolism(monomer_to_microbe_flux,...
    enz_pro_stress, x, vid, par_mic)
%usage
%do cell metabolism, maitenance, growth, and enzyme production
%I assume all microbes are rigid withtout reserves 
%reference: Tang and Riley, 2015
%
%density dependent microbial mortality added by Rose Abramoff on 12/19/2016
%following Georgiou et al. 2017
%
%input argument
%monomer_to_microbe_flux: monomer assimilation flux, nmonomers x nmicrobes
%enz_pro_stress:  enzyme production stress, 1 x nmicrobes
%x: state variable
%vid: indice id
%par_mic: par_micameter structure

%return variable
%maintenance: bulk maintenance rate [mol c/ut], ut-> unit of time
%growth: population growth rate     [mol c/ut]
%pro_enz: enzyme production rate    [mol c/ut]
%death: mortality rate              [mol c/ut]
%co2_growth_flux: co2 respiration rate from growth [mol c/ut]
%get number of microbes

nmicrobes=length(par_mic);
%initialize vector of mortality
death_specific=zeros(nmicrobes,1);

%specific maintenance demand
maintenance_demand=calc_microbe_specific_maintenance(par_mic,nmicrobes);

%compute the specific carbon flux to support microbial metabolism
ja=sum(monomer_to_microbe_flux,1);

%metabolite flux
je=zeros(size(ja));
for jj = 1 : nmicrobes
    if(plasticity_test(par_mic(jj)))
        je(jj)=ja(jj); %if true/rigid
    else
        je(jj)=par_mic(jj).kappa_micb*x(vid.micc(jj))/x(vid.microbep(jj)); %if false/plastic
    end
end

%do metabolism budget
[growth_specific,maintenance_specific,pro_enz_specific]=deb_budget(je,maintenance_demand,enz_pro_stress, x, vid, par_mic);

%co2 flux from cell metabolism
co2_metab_flux_bulk=zeros(size(je));    % 1 x nmicrobes
for jj = 1 : nmicrobes
    if(plasticity_test(par_mic(jj)))
        co2_metab_flux_bulk(jj)=ja(jj)-(growth_specific(jj)+pro_enz_specific(jj))*x(vid.microbep(jj));
    else        
        co2_metab_flux_bulk(jj)=je(jj)*x(vid.microbep(jj))-...
            growth_specific(jj)*x(vid.micc(jj))-...
            (growth_specific(jj)+pro_enz_specific(jj))*x(vid.microbep(jj)); 
    end
end

%compute gross uptake for plastic
for jj = 1 : nmicrobes     
        uptake(jj)=je(jj)*x(vid.microbep(jj)); %carbon exported from reserve for use in maintenance, growth, enzymes
end

%compute specific microbial mortality

for kk=1 : nmicrobes
      death_specific(kk)=par_mic(kk).decay_micb0*(1d0+maintenance_demand(kk)/(je(kk)+maintenance_demand(kk))...
          *par_mic(kk).decay_micb1)*x(vid.microbep(kk))^0.5;
end

% disp([num2str(growth_specific(1)) ',' num2str(death_specific(1)) ...
%     ',' num2str(growth_specific(1)-death_specific(1))]);
end


function maintenance_demand=calc_microbe_specific_maintenance(par_mic, nmicrobes)
%
%DESCRIPTION
% calculate the specific microbial maintenance respiration
maintenance_demand=zeros(nmicrobes,1);

for jj = 1 : nmicrobes
    maintenance_demand(jj)=par_mic(jj).mr_micb;

end

end
