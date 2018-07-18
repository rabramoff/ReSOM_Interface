function par=par_pro_correct(par, par_msurf, micb, micc,doc, msurf)
%do back calculation of the maximum growth rate and maximum enzyme production rate
%based on steady state pool sizes

if(micc==0.0)
    %doc uptake
    Fc=par.zb.*doc.*par.Vmax_micb./(par.Kaff_monomer+doc...
        +par.zb.*micb+msurf.*par.Kaff_monomer./par_msurf.Kaff_monomer);
    jc=Fc.*par.Yld_micx_monomer-par.mr_micb;
    %maximum enzyme production rate
    par.pro_ee0=1./(1./par.pro_ee-1./(jc.*par.Yld_ee));
    %maximum population growth rate
    par.gmax_micb=1./(1./par.decay_mic-1./(jc.*par.Yld_micb));
else
    %carbon export flux
    jc=(par.kappa_micb-par.decay_mic).*micc./micb-par.mr_micb;
    %maximum enzyme production rate
    par.pro_ee0=1./(1./par.pro_ee-1./(jc.*par.Yld_ee));
    %maximum population growth rate
    par.gmax_micb=1./(1./par.decay_mic-1./(jc.*par.Yld_micb));
end

end