%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [growth,maintenance,pro_enz]=deb_budget(je,maintenance_demand,...
    enz_pro_stress, x, vid, par_mic)
%usage
%do microbial metabolism using DEB theory
%
%input arguments
%je: specific reserve export flux, before dilution correction
%maintenance_demand: specific metabolic flux demanded for maintenance 
%enz_pro_stress: enzyme production stress
%x: state variable vector
%vid: variable index structure
%par: microbial parameter structure
%
%return variable
%growth: specific population growth rate
%maintenance: specific maintenance rate
%pro_enz: specific enzyme production rate


nmicrobes=length(maintenance_demand);
maintenance=zeros(nmicrobes,1);
growth=zeros(nmicrobes,1);
pro_enz=zeros(nmicrobes,1);

for kk = 1 : nmicrobes
    if(plasticity_test(par_mic(kk)))
        %rigid microbe
        [growth(kk),maintenance(kk),pro_enz(kk)]=deb_rigid(je(kk),...
            maintenance_demand(kk),enz_pro_stress(kk),par_mic(kk));
    else        
        %plastic microbe
        ec=x(vid.micc(kk))*x(vid.microbep(kk))/(x(vid.microbep(kk)).^2+1d-20);
        [growth(kk),maintenance(kk),pro_enz(kk)]=deb_plastic(je(kk),ec,...
            maintenance_demand(kk),enz_pro_stress(kk),par_mic(kk));
        
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [growth,maintenance,pro_enz]=deb_rigid(je,maintenance_demand, ...
    enz_pro_stress, par_mic)
%usuage
%
%do DEB budget for rigid microbe
%input arguments
%je: specific reserve export flux
%maintenance_demand: specific metabolic flux demanded for maintenance 
%enz_pro_stress: enzyme production stress
%x: state variable vector
%par_mic: microbial parameter structure
%
%return variable
%growth: specific population growth rate
%maintenance: specific maintenance rate
%pro_enz: specific enzyme production rate


%does the specific carbon flux greater than maintenance

if(je>maintenance_demand)

    %actual maintenance rate
    maintenance=maintenance_demand;
    
    %jc flux supports enzyme production and growth        
    jc=je-maintenance_demand;    
    
    %maximum production rate for enzyme       
    pe_max=par_mic.pro_ee1*enz_pro_stress+par_mic.pro_ee0;
    
    %maximum production vector        
    gmax_vec=[par_mic.gmax_micb,pe_max]'; 
         
    %yield rate for different processes    
    
    yld_vec=[par_mic.Yld_micb,par_mic.Yld_ee]';    
    
    %potential production rate    
    
    gp=1./(1./gmax_vec+1./(jc.*yld_vec));    
    
    %scaling factor due to substrate limitation    
    
    scal=min([jc/(sum(gp./yld_vec)),1]);    
    
    %actual population growth rate    
    
    growth=gp(1)*scal;    
    
    %actual enzyme production rate    
    
    pro_enz=gp(2)*scal;
    
else    
    maintenance=je;    
    growth=0;    
    pro_enz=0;    
end 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [growth,maintenance,pro_enz]=deb_plastic(je,ec,maintenance_demand, ...
    enz_pro_stress, par_mic)
%usage
%do DEB budget for plastic microbe
%
%input arguments
%je: specific reserve export flux
%maintenance_demand: specific metabolic flux demanded for maintenance 
%enz_pro_stress: enzyme production stress
%x: state variable vector
%par_mic: microbial parameter structure
%
%return variable
%growth: specific population growth rate
%maintenance: specific maintenance rate
%pro_enz: specific enzyme production rate

pe_max=par_mic.pro_ee1*enz_pro_stress+par_mic.pro_ee0;

%determine if there is growth in population        
isgrw=deb_microbe_init(maintenance_demand,[par_mic.gmax_micb,pe_max]',[],...
        [par_mic.Yld_micb,par_mic.Yld_ee]',je,ec);
      
switch isgrw    
    case -1            
        %no growth        
        %there is no flux for growth, but penalty for mortality        
        growth=0.0;        
        pro_enz=0.0;        
        maintenance=je;            
        
    case 0    
        %limited growth        
        %        
        %potential growth rate > 0        
        %normalize all rates        
        %normalized maximum growth rate        
        gmax_vec=[1,pe_max./par_mic.gmax_micb]';
        %normalized microbial maintenance
        mr_micb=maintenance_demand./par_mic.gmax_micb;        
        %normalized reserve export
        je1=je./par_mic.gmax_micb;        
        %Yld rate vector               
        yld_vec=[par_mic.Yld_micb,par_mic.Yld_ee]';
        %solve for the growth rate
        gp=deb_microbe(mr_micb,gmax_vec,0,yld_vec,je1,ec);
        %population growth        
        growth=gp(1).*par_mic.gmax_micb;       
        %enzyme production        
        pro_enz=gp(2).*par_mic.gmax_micb;        
        %maintenance
        maintenance=maintenance_demand;    
        
    case 1    
        %maximum growth       
        growth=par_mic.gmax_micb;        
        %maximum enzyme production        
        pro_enz=pe_max;        
        %maintenance        
        maintenance=maintenance_demand;
    
end

end
