function [enz_pro_stress_comp,enz_pro_stress_tot]=cal_enz_pro_stress(monomer_mic_complex_matrix,eff_mic_surface,par_ss)
%usage
%compute enzyme production stress
%
%Input arguments
%monomer_complex_matrix: the monomer binding matrix, nomonmers x nmicrobes
%eff_mic_surface: effective transporter surface to adsorb monomers
%Output
%enz_pro_stress_comp, enzyme production stress for each monomers
%enz_pro_stress_tot, total enzyme production stress

%In this version, I set enzyme production stress to 1 for simplicity, I
%will try a more mechanistic idea later, Jinyun Tang

n_microbes=length(eff_mic_surface);
enz_pro_stress_tot=repmat(1d0,1,n_microbes);

enz_pro_stress_comp=[];

end