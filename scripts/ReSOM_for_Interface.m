    %Script to set up plots for Sulman et al. model comparison
%Run scripts and modifications by Rose Abramoff
%Core model code by Jinyun Tang
%Preconditions
%runtime 0-50 years
%constant temperature: 293.15 K (20 deg C)
%constant C inputs: 500 gC/m2/year (fraction to soc/doc is same as default)
%depth: 20cm

%Exp 1: +2 K
    %Run 1: low clay, high quality litter
    %Run 2: low clay, low  quality litter
    %Run 3: mid clay, high quality litter
    %Run 4: mid clay, low quality litter
    %Run 5: high clay, high quality litter
    %Run 6: high clay, low quality litter 
%Exp 2: +5 K
    %Run 1: low clay: 5% clay, 50% sand, high quality litter: 1.37% N, 16.6% lignin
    %Run 2: low clay, low  quality litter: 0.73% N, 24.4% lignin
    %Run 3: mid clay: 40% clay, 10% sand, high quality litter
    %Run 4: mid clay, low quality litter
    %Run 5: high clay: 70% clay, 30% sand, high quality litter
    %Run 6: high clay, low quality litte
%Exp 3: +30% monomers and polymers (total input)
    %Run 1: low clay, high quality litter
    %Run 2: low clay, low  quality litter
    %Run 3: mid clay, high quality litter
    %Run 4: mid clay, low quality litter
    %Run 5: high clay, high quality litter
    %Run 6: high clay, low quality litte
%Exp 4: +100% monomers and polymers (total input)
    %Run 1: low clay, high quality litter
    %Run 2: low clay, low  quality litter
    %Run 3: mid clay, high quality litter
    %Run 4: mid clay, low quality litter
    %Run 5: high clay, high quality litter
    %Run 6: high clay, low quality litte
%Exp 5: +30% monomers
    %Run 1: low clay, high quality litter
    %Run 2: low clay, low  quality litter
    %Run 3: mid clay, high quality litter
    %Run 4: mid clay, low quality litter
    %Run 5: high clay, high quality litter
    %Run 6: high clay, low quality litte
%Exp 6: +2 K and +30% monomers and polymers (total input)
    %Run 1: low clay, high quality litter
    %Run 2: low clay, low  quality litter
    %Run 3: mid clay, high quality litter
    %Run 4: mid clay, low quality litter
    %Run 5: high clay, high quality litter
    %Run 6: high clay, low quality litte
%Exp 7: +2 K and +30% monomers
    %Run 1: low clay, high quality litter
    %Run 2: low clay, low  quality litter
    %Run 3: mid clay, high quality litter
    %Run 4: mid clay, low quality litter
    %Run 5: high clay, high quality litter
    %Run 6: high clay, low quality litte
%Exp 8: Remove all inputs
    %Run 1: low clay, high quality litter
    %Run 2: low clay, low  quality litter
    %Run 3: mid clay, high quality litter
    %Run 4: mid clay, low quality litter
    %Run 5: high clay, high quality litter
    %Run 6: high clay, low quality litte
    
    %Start timer
    tic
    %spinup the four runs for 300 years to SS
    continueOn = 0;
    reps = 300;
    expt = 0; %no experiment
    surfinit = [6029.120 7332.271 8445.413]; %using both Jagadamma et al. 2014 & Doetterl et al. 2015
    litterQ = [2.4133-2.4133*0.43 2.4133-2.4133*0.10]; %recalcitrant & labile
        
    %Create output folder
    [status,results]=system('pwd');
    matfldir=['output'];
    system(['mkdir ', matfldir]);
    addpath(pwd);
    addpath([pwd,'/model']);
    addpath([pwd,'/scripts']);
    
    for i = 1:length(surfinit)
        for j = 1:length(litterQ)
 
    run_INTERFACE(reps,continueOn,surfinit(i),litterQ(j),expt,matfldir)
    load([matfldir,'/mbms_noIso_spinup_',num2str(surfinit(i)),'surfinit_',num2str(litterQ(j)),'litterQ.mat'])
 
    %SS test
    100*(YOUT_ctrl(round(reps*365/40),1:8)-YOUT_ctrl(round((reps-50)*365/40),1:8))./YOUT_ctrl(round((reps-50)*365/40),1:8) < 0.5 %each pool should change less than 0.5% in 50 years
        end
    end
    
    %after 10 years, run experiment for another 50
    continueOn = 1;
    reps = 60;
    expt = [0 1 2 3 4 5 6 7 8]; %Exp 1:8
    
    for i = 1:length(surfinit)
        for j = 1:length(litterQ)
            for k = 1:length(expt)

                run_INTERFACE(reps,continueOn,surfinit(i),litterQ(j),expt(k),matfldir)
    
            end
        end
    end

    %End timer
    toc
