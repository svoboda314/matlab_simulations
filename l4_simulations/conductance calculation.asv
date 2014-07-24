
to=10;%membrane time constant of the postsynaptic neuron in ms 
% to = 20 for ss; to = 10 for fs  
ts= 3; % time constant of synaptic current in ms
% ts = 3 for all ampa / gaba-A
psp =1.8;%peak psps amplitude in mV
dV =25;% driving force during measurement in mV
c =1; %specific membrane capacitance (1 uF/cm^2)

f = (to/(to-ts))*((ts/to)^(ts/(to-ts))-(ts/to)^(to/(to-ts)));

cs=-psp*c/(dV*f)
