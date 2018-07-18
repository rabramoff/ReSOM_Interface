function par=set_msurface_par_default( )
%
%Description
%set up parameters for monomer-mineral surface interaction
%author: Jinyun Tang
global cpar;
global vid;
nmicrobes=length(vid.microbep);
nenzymes=length(vid.enzymes);
nsurfaces=length(vid.surfaces);
nmonomers=length(vid.monomers);
npolymers=length(vid.polymers);

par.Kaff_monomer=repmat(25d0,1,nmonomers);                    	% 1 x nmonomers, adsorption surface doc affinity           (g C/m3)
par.Vmax_ads_monomer=cpar{13};      % 1 x nmonomers, maximum monomer adsorption rate (1/day)
par.adsmon_decay=cpar{14};          % 1 x nmonomers, arbitary, need change
par.adspoly_decay=cpar{14}/100;          % 1 x npolymers, arbitary, need change
end