function par=set_enzyme_par_default(scal, opt)

% Description:
%   Setup parameters for enzymes
%   Added by XZhu
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

par.decay_ee   =cpar{9}.*scal;                                             % enzyme decay rate          (1/day)

% - vectors

par.Vmax_ee    =cpar{10}.*scal;                             % 1 x npolymers, maximum som degradation rate              (1/day)
par.Kaff_ee_polymer =repmat(20d1,1,npolymers);                              % 1 x npolymers, enzyme affinity to polymer                (g ee c/m3)
par.Kaff_ee_msurf= repmat(5d1,1,nsurfaces);                                 % 1 x nsurfaces, enzyme affinity for adsorptive surface    (g ee c/m3)
par.Vmax_ads_enzyme = cpar{11}.*scal;                     % 1 x nsurfaces, maximum enzyme adsorption rate (1/day)
par.adsenz_decay= cpar{12}.*scal;                          % 1 x nsurfaces, decay rate of adsorbed enzyme (1/day)

end