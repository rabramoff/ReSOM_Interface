function tf=plasticity_test(par_mic)
%
%DESCRIPTION
%do microbial plasticity test
%return 1 if the microbe is regid
tf=(par_mic.kappa_micb>1d5);

end