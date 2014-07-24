
% Dayan and Abbott - page 188
% noticed that al very sensitive to initial conditions; shape of the
% synaptic current etc


%%

duration = 2; % duration of simulation; in seconds
tstep=0.1; % time steps, in ms
times=0:tstep:(duration*1000); %in ms
tsteps=length(times);

%% single neuron IF, conductance-based

el = -70; % mV
et = -50; % mV
tau = 5; % ms
Rm = 100; % megaohm
vm =zeros(1, tsteps)+el;


Ie=zeros(1, tsteps);

cstep=5;
nsteps=10;

current_amp=(cstep:cstep:(nsteps*cstep))+200;
%nsteps=length(current_amp);

rate=zeros(1, nsteps);
rate_theory=zeros(1, nsteps);
figure;
hold
for j=1:nsteps
    
Ie(2000:6000)=current_amp(j); % pA

st=[];
for i=2:tsteps
    vm(i) = vm(i-1)-(vm(i-1)-el)*tstep/tau + 0.001*Rm*Ie(i-1)*tstep/tau;
    if vm(i) > et
        vm(i)=el;
        st=[st, i*tstep];
    end
end

rate(j)=1000/mean(diff(st));

%temp=1/(tau*0.001*log(Rm*current_amp(j)/(Rm*current_amp(j)+1000*(el-et))));


rate_theory(j) = 1/(tau*0.001*log(Rm*current_amp(j)/(Rm*current_amp(j)+1000*(el-et))));;

plot(t, vm+j*30)

end
%plot(current_amp, rate)
hold

figure; plot(current_amp, rate,'o')
hold
%rate_theory = 1/(tau*log(Rm*current_amp*1000./(Rm*current_amp*1000+el-et)));
plot(current_amp, rate_theory, 'r')

%% pair of IF neurons inhibiting each other; exponential synapses

nneur=2; % in mV; # of neurons
el = -70; % mV; resting potential
er = -80; % mV; reset after AP
et = -54; % mV; AP threshold  
tau = 20; % ms; neuronal timeconstant
Rm = 100; % megaohm
RmIe= 25; % in mV

esi = -80; % mV; inhibitory reversal potential
dgi=0.05; % ##; change in inhibnition per ap
gi =zeros(nneur, tsteps);
tau_i=10; % in ms

ese = 0; % mV; excitatory reversal potential
dge=0.0; % ##; change in excitation per ap
ge =zeros(nneur, tsteps);
tau_e=10; % in ms

vm =zeros(nneur, tsteps)+el;

cm = [0 1; 1 0]; % all to all
neurons=zeros(nneur,1);

RmIe=zeros(nneur, tsteps);
RmIe(1, 500:end)=25; % current injection neuron 1
RmIe(2, 750:end)=0; % current injection neuron 2

for t=2:tsteps
    vm(:, t)=vm(:, t-1)+(el - vm(:, t-1))*tstep/tau - gi(:, t-1).*(vm(:, t-1)-esi)*tstep/tau ...
        -ge(:, t-1).*(vm(:, t-1)-ese)*tstep/tau+(RmIe(:, t-1))*tstep/tau;    
%    vm(1, t)=vm(1, t)+RmIe*tstep/tau; 
%    vm(2, t)=vm(2, t)+RmIe*tstep/tau; 
    % detect spikes
    spiking_neurons=find(vm(:, t)>et);
    neurons=zeros(nneur,1);
    neurons(spiking_neurons)=1;
    gi(:, t)=gi(:, t-1)*(1-tstep/tau_i)+dgi*cm*neurons;
    ge(:, t)=ge(:, t-1)*(1-tstep/tau_e)+dge*cm*neurons;
    %vm(spiking_neurons-1, t)=0;
    vm(spiking_neurons, t)=er;
    
end

%% pair of neurons inhibiting each other; alpha synapses

nneur=2; % in mV; # of neurons
el = -70; % mV; resting potential
er = -80; % mV; reset after AP
et = -54; % mV; AP threshold  
tau = 20; % ms; neuronal timeconstant
Rm = 100; % megaohm
RmIe= 25; % in mV

esi = -80; % mV; inhibitory reversal potential
dgi=0.05; % ##; change in inhibnition per ap
gi =zeros(nneur, tsteps);
tau_i=10; % in ms

ese = 0; % mV; excitatory reversal potential
dge=0.0; % ##; change in inhibnition per ap
ge =zeros(nneur, tsteps);
tau_e=10; % in ms

vm =zeros(nneur, tsteps)+el;

cm = [0 1; 1 0]; % all to all connections
neurons=zeros(nneur,1); % neurons that spike at time t 

RmIe=zeros(nneur, tsteps);
RmIe(1, 500:end)=25; % current injection neuron 1
RmIe(2, 700:end)=25; % current injection neuron 2

%figure; hold
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
        %    t
    temp1=(cm*neurons);
    temp2=((((t+1):tsteps)-t)*tstep/tau).*exp(1-((((t+1):tsteps)-t))*tstep/tau);
    gi(:, (t+1):tsteps)=gi(:, (t+1):tsteps)+...
        dgi*temp1*temp2;
    
    ge(:, (t+1):tsteps)=ge(:, (t+1):tsteps)+...
        dge*temp1*temp2;
%    plot((t+1):tsteps, temp2)
    end

    vm(spiking_neurons, t)=er;    
end

%% plotting
subplot(2, 1, 1)
title('rmgs=0.05, inhibition')
plot((1:2000)*tstep, vm(:,1:2000)')
subplot(2, 1, 2)
plot(((tsteps-1999):tsteps)*tstep, vm(:,(tsteps-1999):tsteps)')
%title('rmgs=0.05, inhibition')






