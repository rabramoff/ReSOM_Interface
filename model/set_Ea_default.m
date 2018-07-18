function TEa=set_Ea_default()
%
%set default activation energy for different processes
%created by Jinyun Tang, Jan 27, 2014

global vid;
nmicrobes=length(vid.microbep);
nenzymes=length(vid.enzymes);
nsurfaces=length(vid.surfaces);
nmonomers=length(vid.monomers);
npolymers=length(vid.polymers);

rgas=8.31446;                        %universal gas constant, [J/K/mol]

TEa.Ea_vmax_micb            = repmat(45d3/rgas,nmicrobes,nmonomers);        % nmicrobes x nmonomers, monomer uptake               [K]
TEa.Ea_vmax_ee              = repmat(45d3/rgas,nenzymes,npolymers);         % nenzymes x npolymers, polymer degradation          [K]
TEa.Ea_Kaff_monomer_micb    = repmat(15d3/rgas,nmicrobes,nmonomers);        % nmicrobes x nmonomers, monomer addsorption to microbe, affinity parameter [K]
TEa.Ea_Kaff_ee_polymer_micb = repmat(15d3/rgas,nenzymes,npolymers);         % nenzymes x npolymers, som adsorption to enzyme, affinity parameter [K]
TEa.Ea_mr_micb              = repmat(60d3/rgas,nmicrobes,1);                % nmicrobes x 1, maintenance  [K],  0.625ev, Brown 2004
TEa.Ea_kappa_micb           = repmat(60d3/rgas,nmicrobes,1);                % nmicrobes x 1, reserve export [K], 0.625ev, Brown 2004

TEa.Ea_Kaff_monomer_msurf   = repmat(10d3/rgas,nsurfaces,nmonomers);        % nsurfaces x nmonomers, doc adsorption to mineral    [K]
TEa.Ea_Kaff_ee_msurf        = repmat(10d3/rgas,nenzymes,nsurfaces);        % nenzymes x nsurfaces, enzyme adsorption to mineral [K]
TEa.Ea_Kaff_polymer_polymer   = repmat(10d3/rgas,nsurfaces,nmonomers);        % nsurfaces x nmonomers, polymer adsorption to polymers    [K]
TEa.Ea_Kaff_polymer_msurf   = repmat(10d3/rgas,nsurfaces,nmonomers);        % nsurfaces x nmonomers, polymer adsorption to mineral    [K]

TEa.Ea_Vmax_ads_monomer     = repmat(45d3/rgas,nsurfaces,nmonomers);        % nsurfaces x nmonomers, monomer adsorption [K]
TEa.Ea_Vmax_ads_enzyme      = repmat(45d3/rgas,nenzymes,nsurfaces);         % nenzymes x nsurfaces, enzyme adsorption [K]
TEa.Ea_Vmax_ads_polymer     = repmat(45d3/rgas,nsurfaces,npolymers);        % nsurfaces x npolymers, polymer-polymer adsorption [K]
TEa.Ea_Vmax_ads_polymer_msurf      = repmat(45d3/rgas,npolymers,nsurfaces);         % npolymers x nsurfaces, polymer-mineral adsorption [K]
end