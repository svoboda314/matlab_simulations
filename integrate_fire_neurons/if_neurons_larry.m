
duration = 2.0; % duration of simulation; in seconds
tstep=0.1; % time steps, in ms
times=0:tstep:(duration*1000); %in ms
tsteps=length(times);



%% pair of neurons inhibiting each other; alpha synapses

nneur=2; % in mV; # of neurons
el = -70; % mV; resting potential
er = -80; % mV; reset after AP
et = -54; % mV; AP threshold  
tau = 20; % ms; neuronal timeconstant
Rm = 100; % megaohm
RmIe= 25; % in mV

esi = -80; % mV; inhibitory reversal potential
dgi=0.0; % ##; change in inhibnition per ap
gi =zeros(nneur, tsteps);
tau_i=10; % in ms

ese = 0; % mV; excitatory reversal potential
dge=0.05; % ##; change in inhibnition per ap
ge =zeros(nneur, tsteps);
tau_e=10; % in ms

vm =zeros(nneur, tsteps)+el;

cm = [0 1; 1 0]; % all to all
neurons=zeros(nneur,1); % neurons that spike 

RmIe=zeros(nneur, tsteps);
RmIe(1, 500:end)=25;
RmIe(2, 600:end)=25;

spiketotal=0;
for t=2:(tsteps-1)    
   vm(:, t)=vm(:, t-1)+(el - vm(:, t-1))*tstep/tau - gi(:, t-1).*(vm(:, t-1)-esi)*tstep/tau ...
        -ge(:, t-1).*(vm(:, t-1)-ese)*tstep/tau+(RmIe(:, t-1))*tstep/tau; %*[0 1]';     
    % detect spikes
    spiking_neurons=find(vm(:, t)>et);
    neurons=zeros(nneur,1); % reset spiking neuron vector
    neurons(spiking_neurons)=1; % set spikeing neurons to 1
        spiketotal=spiketotal+sum(neurons);

    if sum(neurons) > 0
    temp1=(cm*neurons);
%    t
    temp2=((((t+1):tsteps)-t)*tstep/tau).*exp(1-((((t+1):tsteps)-t))*tstep/tau);
    gi(:, (t+1):tsteps)=gi(:, (t+1):tsteps)+...
        dgi*temp1*temp2;
    
    ge(:, (t+1):tsteps)=ge(:, (t+1):tsteps)+...
        dge*temp1*temp2;
    end
    vm(spiking_neurons, t)=er;    
end








