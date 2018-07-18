function par=set_microbe_par_default(scal, opt)

% Description
%   Setup parameters for microbes
%   Author: Jinyun Tang

global cpar;
global vid;

nmicrobes=length(vid.microbep);
nenzymes=length(vid.enzymes);
nsurfaces=length(vid.surfaces);
nmonomers=length(vid.monomers);
npolymers=length(vid.polymers);

if(nargin==0)
    scal=1;
end
if(nargin<=1)
    opt=1;
end

% - scalars

par.decay_micb0=cpar{6}.*scal;                                           % microbial death rate       (1/day)
par.decay_micb1=0d0;                                                       % scaling parameter to account for mortality change from respiration stress
par.gmax_micb  =cpar{8}.*scal;                                                % growth rate (1/day)
par.pro_ee0    =cpar{7}.*scal;                                                % constintutive enzyme production rate     (1/day)
par.pro_ee1    =0d-6.*scal; %cpar{7}.*scal; %                                 % inductive enzyme production rate  (1/day)
par.Yld_micb   =0.8;                                                       % growth efficiency of enzyme and microbes  (g mic C/g res C)
par.Yld_ee     =0.8;
par.mr_micb    =cpar{3}.*scal;                                        % microbial maintenance rate                (1/day)
par.kappa_micb =cpar{4}.*scal;                                     % reserve turnover rate                     (1/day)
par.zb         =0.05;                                                      % scaling factor between transporter and microbial cell biomass
par.micb0      =1d-4;                                                      % half saturating microbial biomass for mortality (gC/m3), I have used an alternative value 1d-3 for sensitivity test

% - vectors

par.Kaff_monomer =repmat(1d0,1,nmonomers);                                 % 1 x nmonomers, microbial doc affinity                    (g C/m3)
par.Vmax_micb  =cpar{5}.*scal;                                      % 1 x nmonomers, maximum doc uptake rate                   (1/day)
par.Yld_micx_monomer   =repmat(0.5,1,nmonomers);                           % 1 x nmonomers, assimilation efficiency from doc/monomer uptake   (g res C/g DOC C)
par.pro_enz_dist = repmat(1/nenzymes,1,nenzymes);                          % 1 x nenzymes, distribution of bulk enzyme biomass to each enzyme
%par.gmax_micb  =repmat(cpar{8}.*scal,1,nmicrobes);

% the two parameters below are only for analytic solution test
par.pro_ee=cpar{15}.*scal;
par.decay_mic = par.decay_micb0; 
end