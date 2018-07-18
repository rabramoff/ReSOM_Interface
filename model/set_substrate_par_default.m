function par_ss=set_substrate_par_default()
%
% DESCRIPTION
% set up substrate default parameters

global vid;
nmicrobes=length(vid.microbep);
nenzymes=length(vid.enzymes);
nsurfaces=length(vid.surfaces);
nmonomers=length(vid.monomers);
npolymers=length(vid.polymers);

par_ss.deadmicrobep2polymers=repmat(1/npolymers,npolymers,nmicrobes);    % all dead microbep into polymers
par_ss.deadmicrobep2monomers=zeros(nmonomers,nmicrobes);

par_ss.deadmicc2polymers=zeros(npolymers,nmicrobes);    % all dead micc into monomers
par_ss.deadmicc2monomers=repmat(1/nmonomers,nmonomers,nmicrobes);

par_ss.deadenzpart = 0.2;
par_ss.deadenz2polymers=repmat(par_ss.deadenzpart/npolymers,npolymers,nenzymes);     % 20% of dead enz into polymers
par_ss.deadenz2monomers=repmat((1-par_ss.deadenzpart)/nmonomers,nmonomers,nenzymes);

par_ss.polymer2monomer=repmat(1/nmonomers,nmonomers,npolymers);

end