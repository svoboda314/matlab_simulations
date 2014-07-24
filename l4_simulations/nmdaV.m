
function ginf=nmdaV(v)
% voltage-dependece of the NMDA-R
% v, membrane potential in mV
% from Gabbiani et al 1994

ginf=1.358./(1+exp(-0.062*v)*1.2/3.57);