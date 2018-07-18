function run_INTERFACE(reps,continueOn,surfinit,litterQ, expt, matfldir)

if continueOn == 0
    iofile=[matfldir,'/mbms_noIso_spinup_',num2str(surfinit),'surfinit_',num2str(litterQ),'litterQ.mat'];
else
    iofile=[matfldir,'/mbms_noIso_',num2str(surfinit),'surfinit_',num2str(litterQ),'litterQ_',num2str(expt),'expt.mat'];
end

% Set number of pools
n_polymers=1;   % cellulose-hemicellulose and ligin
n_monomers=1;   % dissolved organic carbon (compounds of varying complexity)
n_enzymes=1;    % use a generic enzyme to minimize model parameters (BG/Phenox activities reflect corresponding degradation rates) 
n_microbep=1;   % microbes can directly uptake DOC 
n_micc=n_microbep;
n_surfaces=1;
n_co2=1;
n_enzymes_ads=n_enzymes*n_surfaces; % consider every surface has all enzymes/monomers (e.g., 3enzymes x 2surfaces network: enz1-sur1,enz1-sur2,enz2-sur1,...,enz3-sur2)
n_monomers_ads=n_monomers*n_surfaces;
n_all=n_microbep+n_micc+n_enzymes+n_polymers+n_monomers+n_surfaces+n_co2+n_enzymes_ads+n_monomers_ads;

% Set key parameters those results in microbe biomass of 3.89% total organic carbon
cpar0=[300/365*.8 300/365*.2 4d-4*57.653333 1d-1*0.536542 2d-2*546.713826 ...
    1.5d-5*875.965050 0.002215393623082 0.102525103775396 1d-5*611.427356 2d-2*120.666087 ...
    .001 .006 .01 .006 5d-6*386.294367]';
global cpar;            % a vector of considered calibrating parameters
cpar=[
    {repmat(cpar0(1)/n_polymers,1,n_polymers)};       % 1 x npolymers, polymer input, 1/day
    {repmat(cpar0(2)/n_monomers,1,n_monomers)};       % 1 x nmonomers, monomer input, 1/day
    {cpar0(3)};       % 1 x nmicrobes, microbial maintenance rate (1/day)
    {cpar0(4)};       % 1 x nmicrobes, reserve turnover rate, 1/day
    {repmat(cpar0(5),1,n_monomers)};       % 1 x nmonomers, maximum doc uptake rate                   (1/day)
    
    {cpar0(6)};       % 1 x nmicrobes, microbial death rate       (1/day)
    {cpar0(7)};       % 1 x nenzymes, maximum enzyme production rate     (1/day) [inverted from analytic solution]
    {cpar0(8)};       % 1 x nmicrobes, maximum microbial growth rate (1/day) [inverted from analytic solution]
    {cpar0(9)};       % 1 x nenzymes, enzyme decay rate          (1/day)
    {repmat(cpar0(10),1,n_polymers)};      % 1 x npolymers, maximum som degradation rate              (1/day)
    
    {repmat(cpar0(11),1,n_surfaces)};      % 1 x nsurfaces, maximum enzyme adsorption rate (1/day)
    {repmat(cpar0(12),1,n_surfaces)};      % 1 x nsurfaces, decay rate of adsorbed enzyme (1/day)
    {repmat(cpar0(13),1,n_monomers)};      % 1 x nmonomers, maximum monomer adsorption rate (1/day)
    {repmat(cpar0(14),1,n_monomers)};      % 1 x nmonomers, decay rate of adsorbed monomer (1/day)
    {cpar0(15)};      % 1 x nmicrobes, enzyme production rate used for analytic solution only (1/day)
    ];

%
global vid;             % id index
global par_mic;         % each par represents a microbe
global par_enz;         % each par represents a enzyme
global par_surface;     % each par represents a mineral surface
global par_ss;          % parameter structure for substrate
global input;           % substrate input structure

% Set id index - C cycling
vid.microbep=1:n_microbep; 
vid.micc=vid.microbep(end)+(1:n_microbep);
vid.surfaces=vid.micc(end)+(1:n_surfaces);
vid.monomers=vid.surfaces(end)+(1:n_monomers);
vid.monomers_ads=vid.monomers(end)+(1:n_monomers_ads);
vid.polymers=vid.monomers_ads(end)+(1:n_polymers);
vid.enzymes=vid.polymers(end)+(1:n_enzymes);
vid.enzymes_ads=vid.enzymes(end)+(1:n_enzymes_ads);
vid.co2=vid.enzymes_ads(end)+(1:n_co2);
vid.cue=vid.co2(end)+(1:n_microbep); %track cue
vid.defactoTurnover=vid.cue(end)+(1:n_polymers);

% Set external input (should change to time series)
input.polymers= 400/365;%cpar{1}.*ones(1,n_polymers); %0.6575 gC/m3 day
input.monomers= 100/365;%cpar{2}.*ones(1,n_monomers); %0.1644 gC/m3 day

% Set initial states
if continueOn == 1
    load([matfldir,'/mbms_noIso_spinup_',num2str(surfinit),'surfinit_',num2str(litterQ),'litterQ.mat'])
    x = YOUT_ctrl(end,:);
else
x(vid.microbep)=repmat(20/n_microbep,1,n_microbep);
x(vid.micc)=repmat(5/n_microbep,1,n_microbep);
x(vid.surfaces)=repmat(surfinit/n_surfaces,1,n_surfaces);
x(vid.monomers)=repmat(20/n_monomers,1,n_monomers);
x(vid.monomers_ads)=reshape(repmat(.2*x(vid.monomers)/n_surfaces,n_surfaces,1),n_surfaces*n_monomers,1);
x(vid.polymers)=repmat(200/n_polymers,1,n_polymers);
x(vid.enzymes)=repmat(.1/n_enzymes,1,n_enzymes);
x(vid.enzymes_ads)=reshape(repmat(.2*x(vid.enzymes)/n_surfaces,n_surfaces,1),n_surfaces*n_enzymes,1);
x(vid.co2)=0;
x(vid.cue)=repmat(0.3,1, n_microbep); %set initial cue
x(vid.defactoTurnover)=repmat(1/n_polymers,1,n_polymers); %set initial defacto Turnover (of SOM/Litter to compare to Bovkin et al. 2012 BG)
end

% Set key parameters (same for each function group of microbe/surface, might change to consider different groups)
par_mic=set_microbe_par_default();          % microbe-related parameters
for i=2:n_microbep 
    par_mic(i)=set_microbe_par_default(); 
end;

par_enz=set_enzyme_par_default();           % enzyme-related parameters
for i=2:n_enzymes
    par_enz(i)=set_enzyme_par_default();
end
par_surface=set_msurface_par_default();     % mineral surface-related parameters
for i=2:n_surfaces
    par_surface(i)=set_msurface_par_default();
end;

par_ss=set_substrate_par_default();         % substrate-related parameters

% % Get analytical solutions for 1B1S1C1E1M model
% [doc,micb,micc,ee,som,rco2,cue_all,cue_ee]=one_box_ss(input, par_mic, par_enz, par_surface, par_ss, x(vid.surfaces));
% fprintf('doc=%f,micb=%f,micc=%f,ee=%f,som=%f,rco2=%f\n',doc,micb,micc,ee,som,rco2);
% 
% % Update maximum microbial growth rate and maximum enzyme production rate,
% %   inverted from analytic solution, to make sure numerical solutions are
% %   same as analytic solution (valid only for 1B1S1C1E1M model)
% par_mic=par_pro_correct(par_mic,par_surface, micb, micc,doc, x(vid.surfaces));

% % Make a copy of those temperature dependent parameters
par_surface_ref=par_surface;
par_mic_ref=par_mic;
par_enz_ref=par_enz;
par_enz_ref.Vmax_ee = litterQ; % max decomp rate 
par_enz_ref.Kaff_ee_polymer = 200; % kaff of decomp (currently default)

% Define activation energies
TEa=set_Ea_default();
TEa.Ea_vmax_ee = 45000/8.31446; %temp sensitivity of decomp (currently default)

% Set up and run the model
dt=40; % days
tend= 365*reps;
kend=fix(tend/dt);
TOUT_ctrl=zeros(kend+1,1);
YOUT_ctrl=zeros(kend+1,length(x)); %TEMP solution
YOUT_ctrl(1,:)=x;
for kk = 1 : kend
    t=(kk-0.5)*dt;
%     % -- Incorporate temperature effects on parameters
   
    % Define physiological temperature response curve
    Tref=280;%293.15;          %reference temperature
        xpar = [249.544170969785   5341.422691388677   5.617549086429;
            290	5353	6.06;
            234	5347	6.18;
            291	5354	6.01;
            395	5338	5.96;
            249	5354	6.14; 
            915	5324	4.6;
            275	5360	5.83;
            236	5318	5.31;
            453	5358	5.19];
    A1=0;         %seasonal amplitude
    A2=0;          %daily amplitude
    w1=2.*pi/365;
    w2=2.*pi;
    temp=Tref+A1.*sin(w1.*t-pi/2)+A2.*sin(w2.*t);
    
        if(continueOn ==1 && t>10*365)  
            if expt == 1 %increase temperature by 2 degrees K
                temp=2+Tref+A1.*sin(w1.*t-pi/2)+A2.*sin(w2.*t);
            else
                if expt == 2 %increase temperature by 5 degrees K
                    temp=5+Tref+A1.*sin(w1.*t-pi/2)+A2.*sin(w2.*t);
                else %increase monomer and polymer input by 30%
                    if expt == 3
                        input.polymers=((400/365)+0.3*(400/365)).*ones(1,n_polymers);
                        input.monomers=((100/365)+0.3*(100/365)).*ones(1,n_monomers);
                    else
                        if expt == 4 %double monomer and polymer input
                           input.polymers=(400/365)*2.*ones(1,n_polymers);
                           input.monomers=(100/365)*2.*ones(1,n_monomers);
                        else
                            if expt == 5 %increase monomer input by 30%
                                input.monomers=((100/365)+0.3*(100/365)).*ones(1,n_monomers);
                            else
                                if expt == 6 %increase temperature by 2 degrees K and monomer and polymer input by 30%
                                    temp=2+Tref+A1.*sin(w1.*t-pi/2)+A2.*sin(w2.*t); 
                                    input.polymers=((400/365)+0.3*(400/365)).*ones(1,n_polymers);
                                    input.monomers=((100/365)+0.3*(100/365)).*ones(1,n_monomers);
                                else
                                    if expt == 7 %increase temperature by 2 degrees K and monomer input by 30%
                                        temp=2+Tref+A1.*sin(w1.*t-pi/2)+A2.*sin(w2.*t); 
                                        input.monomers=((100/365)+0.3*(100/365)).*ones(1,n_monomers);
                                    else
                                        if expt == 8
                                            input.polymers=0;
                                            input.monomers=0;  
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
    fref0=temp/Tref;    % non-equilibrium but not related with enzyme activities
    Tinv=1/temp-1/Tref;
      
    n_biggest = max(n_microbep, n_enzymes);
    fref = zeros(n_biggest,1);
        for i=1:n_biggest
        [temp0,T_fact00]=get_microbe_physiology_Tcurve(Tref,xpar(i,:));   % active enzyme fraction in total enzyme vs temperaure
        T_fact=interp1(temp0,T_fact00,temp);    % get active enzyme fraction at given temperature
        fref(i)=T_fact*(temp/Tref);  % non-equilibrium processes which are also related with enzyme activities
        end
    
    for i=1:n_microbep
        % equilibrium processes - Arrhenius equation
        par_mic(i).Kaff_monomer=par_mic_ref(i).Kaff_monomer.*exp(-TEa.Ea_Kaff_monomer_micb(i,:)*Tinv);
        par_mic(i).mr_micb=par_mic_ref(i).mr_micb.*exp(-TEa.Ea_mr_micb(i)*Tinv);
        
        % non-equilibrium and enzyme processes - Eyring's transition state theory and Murphy's equations
        par_mic(i).Vmax_micb=par_mic_ref(i).Vmax_micb*fref(i).*exp(-TEa.Ea_vmax_micb(i,:)*Tinv); % maxium DOC uptake rate by microbe
        par_mic(i).kappa_micb=par_mic_ref(i).kappa_micb*fref(i).*exp(-TEa.Ea_kappa_micb(i)*Tinv); % reserve turnover rate 
    end
    for i=1:n_enzymes
        % equilibrium processes
        par_enz(i).Kaff_ee_polymer=par_enz_ref(i).Kaff_ee_polymer.*exp(-TEa.Ea_Kaff_ee_polymer_micb(i,:)*Tinv);
        par_enz(i).Kaff_ee_msurf=par_enz_ref(i).Kaff_ee_msurf.*exp(-TEa.Ea_Kaff_ee_msurf(i,:)*Tinv);
        
        % non-equilibrium and enzyme processes
        par_enz(i).Vmax_ee=par_enz_ref(i).Vmax_ee*fref(i).*exp(-TEa.Ea_vmax_ee(i,:)*Tinv); % maximum SOC degradation rate by enzyme

        % non-equilibrium and non-enzyme processes (mineral adsorption processes)
        par_enz(i).Vmax_ads_enzyme=par_enz_ref(i).Vmax_ads_enzyme*fref0.*exp(-TEa.Ea_Vmax_ads_enzyme(i,:)*Tinv);
    end
    for i=1:n_surfaces
        % equilibrium processes
        par_surface(i).Kaff_monomer=par_surface_ref(i).Kaff_monomer.*exp(-TEa.Ea_Kaff_monomer_msurf(i,:)*Tinv);

        % non-equilibrium and non-enzyme processes
        par_surface(i).Vmax_ads_monomer=par_surface_ref(i).Vmax_ads_monomer*fref0.*exp(-TEa.Ea_Vmax_ads_monomer(i,:)*Tinv);
    end

%     % -- Run core program
    
    YOUT_ctrl(kk+1,:)=adptmbbks1(@many_bug_many_substrate,YOUT_ctrl(kk,:),TOUT_ctrl(kk)+dt/2,dt);
    TOUT_ctrl(kk+1)=TOUT_ctrl(kk)+dt;
%     if mod(kk,100)==1 disp(kk); end;
end

save(iofile,'YOUT_ctrl','TOUT_ctrl','vid');

end
