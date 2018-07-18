function [doc,micb,micc,ee,som,rco2,cue_all,cue_ee]=one_box_ss(input, par_mic, par_enz, par_msurf, par_ss, msurf)
%Jinyun Tang, Nov 25, 2013
%steady-state, one-box deb model
%doc : dissolved organic carbon
%micb : microbial population
%micc : microbial reserve
%ee : extracellular enzyme
%som: polymeric soc
%rco2: respiration
%cue_all: cue considering extraenzyme as cost
%cue_ee: doese not consider extraenzyme as cost
if(nargin>=1) 
    input_som = input.polymers;                                            %input rate of somc         (gC/day)
    input_doc = input.monomers;                                            %input rate of doc          (gC/day)
    decay_mic = par_mic.decay_mic;                                         %microbial death rate       1/day
    decay_ee  = par_enz.decay_ee;                                          %enzyme decay rate          1/day
    pro_ee    = par_mic.pro_ee;                                            %enzyme production rate     1/day/(g mic C)
    yield_b   = par_mic.Yld_micb;                                          %growth efficiency of enzyme and microbes  (g mic C/g res C)
    yield_x   = par_mic.Yld_micx_monomer;                                  %assimilation efficiency from doc uptake   (g res C/g DOC C)
    mr        = par_mic.mr_micb;                                           %microbial maintenance rate                (1/day)
    kappa     = par_mic.kappa_micb;                                        %reserve turnover rate                     (1/day)
    z         = par_mic.zb;                                                %scaling factor between transporter and microbial cell biomass
    vmax_doc  = par_mic.Vmax_micb;                                         %maximum doc uptake rate                   (1/day)
    vmax_som  = par_enz.Vmax_ee;                                           %maximum som degradation rate              (1/day)
    kb_doc    = par_mic.Kaff_monomer;                                      %microbial doc affinity                    (g C)
    km_doc    = par_msurf.Kaff_monomer;                                    %adsorption surface doc affinity           (g C)
    ke_som    = par_enz.Kaff_ee_polymer;                                   %enzyme affinity to som                    (g ee c)
    km_ee     = par_enz.Kaff_ee_msurf;                                     %enzyme affinity for adsorptive surface    (g ee c)
    M1        = msurf;                                                     %abundance of mineral surface              (g C surface)
    fee2som   = par_ss.deadenz2polymers;                                   %proportation of degraded exoenzyme into som (g som C/g ee C)
else
    input_som=0.1;                                                         %input rate of somc         (gC/day)
    input_doc=0.01;                                                        %input rate of doc          (gC/day)
    decay_mic=1.5d-3;                                                        %microbial death rate       1/day
    decay_ee=1.e-4;                                                        %enzyme decay rate          1/day
    pro_ee = 1d-5;                                                         %enzyme production rate     1/day/(g mic C)
    yield_b=0.3;                                                           %growth efficiency of enzyme and microbes  (g mic C/g res C)
    yield_x=0.4;                                                           %assimilation efficiency from doc uptake   (g res C/g DOC C)
    mr=1d-3;                                                               %microbial maintenance rate                (1/day)
    kappa=1d-1;                                                            %reserve turnover rate                     (1/day)
    z=0.05;                                                                %scaling factor between transporter and microbial cell biomass
    vmax_doc=5d1;                                                          %maximum doc uptake rate                   (1/day)
    vmax_som=7d0;                                                          %maximum som degradation rate              (1/day)
    kb_doc=1d-2;                                                           %microbial doc affinity                    (g C)
    km_doc=1d2;                                                            %adsorption surface doc affinity           (g C)
    ke_som=5d1;                                                            %enzyme affinity to som                    (g ee c)
    km_ee=1d2;                                                             %enzyme affinity for adsorptive surface    (g ee c)
    M1=1d3;                                                                %abundance of mineral surface              (g C surface)
    fee2som=0.2;                                                           %proportation of degraded exoenzyme into som (g som C/g ee C)
end
decay_b=decay_mic+pro_ee;                                                  %bulk decay rate of microbial biomass   (total equivalent decay rate 1/day)    

a11=input_som+input_doc;

if(kappa>1d5)
    micb=a11.*yield_x./(mr+decay_b.*(1./yield_b-yield_x)).*ones(size(M1));
    micc=0;
    a12=(mr+decay_b./yield_b);
    a13=kb_doc+z.*micb+kb_doc./km_doc.*M1;    
    doc=a12.*a13./(z.*vmax_doc.*yield_x-a12);
    cue_all=decay_mic.*yield_x./a12;


    ee=pro_ee.*micb./decay_ee;
    a15=input_som+decay_mic.*micb+fee2som.*decay_ee.*ee;
    
    a16=ke_som+ee+ke_som./km_ee.*M1;

    som=a15.*a16./(ee.*vmax_som-a15);    

    if(som<0)
        fprintf('a15=%f,a16=%f,ee=%f,vmax_som=%f\n',a15,a16,ee,vmax_som);
    end
    %doc uptake
    Fc=z.*micb.*doc.*vmax_doc./(kb_doc+doc+z.*micb+M1.*kb_doc./km_doc);    
    
    cue_ee=1-a11./Fc;
    cue_ees=decay_b./a12.*yield_x;
    %fprintf('cue_ee=%f,cue_ees=%f,tm1=%f,tm2=%f\n',cue_ee,cue_ees,decay_b,a12);
   
else
    
    aa=(mr+decay_b./yield_b)./(kappa-decay_mic);

    a12=(kappa./yield_x-decay_mic).*aa-decay_b;

    micb=a11./a12.*ones(size(M1));

    ee=pro_ee.*micb./decay_ee;

    micc=aa.*micb;

    a13=kb_doc+z.*micb+kb_doc./km_doc.*M1;

    a14=kappa.*aa./yield_x;

    doc=a14.*a13./(z.*vmax_doc-a14);

    a15=input_som+decay_mic.*micb+fee2som.*decay_ee.*ee;

    a16=ke_som+ee+ke_som./km_ee.*M1;

    som=a15.*a16./(ee.*vmax_som-a15);

    cue_all=(1+1./aa).*decay_mic./kappa.*yield_x;

    cue_ee=1-(a11./micb).*yield_x./(kappa.*aa);

end

rco2=a11./micb;
%cue=1-(a11./micb+pro_ee).*(1-decay_mic./kappa)./aa.*yield_x;



%max(abs(pro_ee.*micb-decay_ee.*ee))
%Fs=ee.*som.*vmax_som./(ke_som.*(1+som./ke_som+ee./ke_som+M1./km_ee));
%Fc=z.*micb.*doc.*vmax_doc./(kb_doc.*(1+doc./kb_doc+z.*micb./kb_doc+M1./km_doc));
%max(abs(yield_x.*Fc-kappa.*micc))

if(any([doc,micb,micc,ee,som]<0))

    print_negval(doc,'doc');
    print_negval(micb,'micb');
    print_negval(micc,'micc');
    print_negval(ee,'ee');
    print_negval(som,'soc');
    
    error('negative state variable encountered in one_box_ss');
end
end

function print_negval(invar,varname)

    if(any(invar<0))
        fprintf('%s\n',varname);
        invar
    end
end

