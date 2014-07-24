

% single cell receives input
% intracortical connectivity is turned off

% TO BE CHANGED

l4_l4_cp=0; %0.25; % connection probability
l4_l4_dge=0.02; % conductance in fraction change to reversal

l4_fs_cp=0; %0.4;
l4_fs_dge=0.05; 

fs_l4_cp=0; %0.5;
fs_l4_dgi=0.06; % 0.06 best estimate

fs_fs_cp=0; %0.5;
fs_fs_dgi=0.12; 

vpm_l4_cp=1; %0.5; % Bruno & Sakmann
vpm_l4_dge=0.004; % Beierlein; uF/cm^2
vpm_l4_dt=0; % in seconds; delay of VPM-l4 compared to VPM-fs; 'models' the faster rise-time of the conductance onto fs cells (Cruikshank and Connors)

%ratio_fs_l4=2; %4/2.5; % ratio of conductance fs/l4=4 (Cruikshank and Connors); ratio of input impedances 1/3 (Beierlein)
vpm_fs_cp=1; %0.5 %0.75; 
vpm_fs_dge=0.01; %0.1; %(Cruikshank and Connors)

% END TO BE CHANGED

% ****** simulation paramters *****
global net sim
sim.fractional_system_size=0.1;
sim.duration = .5; % duration of simulation; in seconds
sim.tstep=0.1; % time steps, in ms
sim.times=0:sim.tstep:(sim.duration*1000); %in ms
sim.tsteps=length(sim.times);

% ****** presynaptic VPM neurons ****** 
net.vpm.n=5; % round(250*sim.fractional_system_size); % number of VPM neurons
net.vpm.f=15; % mean rate of modulation of VPM neurons; i.e whisking
net.vpm.sr=50; %50; % mean spike rate
net.vpm.spikes=zeros(net.vpm.n, sim.tsteps);% initialize spike array; zero == no spikes; one == spike

net.vpm.sr_reg=200; % for testing

% ****** set-up the L4 NETWORK ****** 
% l4 (i.e. stellate) neurons
net.l4.n=10; %round(1600*sim.fractional_system_size); 
net.l4.th=-50; % ap threshold; in mv
net.l4.rs=-80; % reset after AP 
net.l4.el=-70; % resting potential
net.l4.ee=0;
net.l4.ei=-80;
net.l4.tau=10; %20; % membrane time constant; in ms
net.l4.cap=1; % specific capacitance, umF/cm^2
net.l4.spikes=zeros(net.l4.n, sim.tsteps); % zero == no spikes; one == spike
net.l4.vm=zeros(net.l4.n, sim.tsteps)+net.l4.el; % in vm
net.l4.gi=zeros(net.l4.n, sim.tsteps); % in CHECK
net.l4.ge=zeros(net.l4.n, sim.tsteps); % in CHECK
net.l4.tge=3; % time constant of excitatory conductance; in ms
net.l4.tgi=3; % time constant of inihibitory conductance; in ms

% fs neurons
net.fs.n=10; %round(200*sim.fractional_system_size);
net.fs.th=-50; % in mv LOWER THRESHOLD?
net.fs.rs=-80;
net.fs.el=-70;
net.fs.ee=0;
net.fs.ei=-80;
net.fs.tau=10; % in ms
net.fs.cap=1;% specific capacitance, umF/cm^2
net.fs.spikes=zeros(net.fs.n, sim.tsteps);
net.fs.vm=zeros(net.fs.n, sim.tsteps)+net.fs.el;
net.fs.gi=zeros(net.fs.n, sim.tsteps);
net.fs.ge=zeros(net.fs.n, sim.tsteps);
net.fs.tge=3; % in ms; note: assume that all excitatory conductnaces have the same tau
net.fs.tgi=3; % in ms; note: assume that all excitatory conductnaces have the same tau

% l4 - l4 connections; square connection matruc, but diagonals are zero
net.l4_l4.cp=l4_l4_cp; % connection probability
net.l4_l4.dge=l4_l4_dge*net.l4.tau/(net.l4.tge*net.l4.cap*sim.fractional_system_size); % in fractional units of driving force; i.e 
net.l4_l4.cm=rand(net.l4.n, net.l4.n);
net.l4_l4.cm(find(net.l4_l4.cm < net.l4_l4.cp))=1;
net.l4_l4.cm(find(net.l4_l4.cm < 1))=0;
net.l4_l4.cm(1:(net.l4.n+1):(net.l4.n*net.l4.n))=0;% set diagonal elements to zero using a(1:(n+1):n^2)=0

% l4 - fs connections; ; not square 
net.l4_fs.cp=l4_fs_cp;
net.l4_fs.dge=l4_fs_dge*net.fs.tau/(net.fs.tge*net.fs.cap*sim.fractional_system_size); % in mv; psps per ap
net.l4_fs.cm=rand(net.fs.n, net.l4.n);
net.l4_fs.cm(find(net.l4_fs.cm < net.l4_fs.cp))=1;
net.l4_fs.cm(find(net.l4_fs.cm < 1))=0;

% fs - l4 connections; ; not square
net.fs_l4.cp=fs_l4_cp;
net.fs_l4.dgi=fs_l4_dgi*net.l4.tau/(net.l4.tgi*net.l4.cap*sim.fractional_system_size); % CHECK OUT
net.fs_l4.cm=rand(net.l4.n, net.fs.n);
net.fs_l4.cm(find(net.fs_l4.cm < net.fs_l4.cp))=1;
net.fs_l4.cm(find(net.fs_l4.cm < 1))=0;

% fs - fs connections; square connection matruc, but diagonals are zero
net.fs_fs.cp=fs_fs_cp;
net.fs_fs.dgi=fs_fs_dgi*net.fs.tau/(net.fs.tgi*net.fs.cap*sim.fractional_system_size); % in nS; conductance per ap CHECK OUT
net.fs_fs.cm=rand(net.fs.n, net.fs.n);
net.fs_fs.cm(find(net.fs_fs.cm < net.fs_l4.cp))=1;
net.fs_fs.cm(find(net.fs_fs.cm < 1))=0;
net.fs_fs.cm(1:(net.fs.n+1):(net.fs.n*net.fs.n))=0;

% VPM - l4 connections; not square ...
net.vpm_l4.cp=vpm_l4_cp; % connection probability
net.vpm_l4.dge=vpm_l4_dge*net.l4.tau/(net.l4.tge*net.l4.cap*sim.fractional_system_size); % in  CHECK OUT
net.vpm_l4.dt=round(vpm_l4_dt*1000/sim.tstep); % in samples; models the fact that vpm-l4 synapses are slower than vpm-fs synapses
net.vpm_l4.cm=rand(net.l4.n, net.vpm.n);
net.vpm_l4.cm(find(net.vpm_l4.cm<net.vpm_l4.cp))=1;
net.vpm_l4.cm(find(net.vpm_l4.cm < 1))=0;

% generate VPM - fs connection matrix; not square ...
net.vpm_fs.cp=vpm_fs_cp; % connection probability
net.vpm_fs.dge=vpm_fs_dge*net.fs.tau/(net.fs.tge*net.fs.cap*sim.fractional_system_size); % 
net.vpm_fs.cm=rand(net.fs.n, net.vpm.n);
net.vpm_fs.cm(find(net.vpm_fs.cm<net.vpm_fs.cp))=1;
net.vpm_fs.cm(find(net.vpm_fs.cm < 1))=0;
