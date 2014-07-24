
% TO BE CHANGED

l4_l4_cp=0.25; % connection probability
l4_l4_dge=0.015; % conductance in fraction change to reversal

l4_fs_cp=0.4;
l4_fs_dge=0.05; 

fs_l4_cp=0.5;
fs_l4_dgi=0.1; 

fs_fs_cp=0.5;
fs_fs_dgi=0.1; 

vpm_l4_cp=0.5; % Bruno & Sakmann
vpm_l4_dge=0.075; % Beierlein; 0.05 for one synapse; could be smalle
vpm_l4_dt=0; % in seconds; delay of VPM-l4 compared to VPM-fs; 'models' the faster rise-time of the conductance onto fs cells (Cruikshank and Connors)

ratio_fs_l4=2; %4/2.5; % ratio of conductance fs/l4=4 (Cruikshank and Connors); ratio of input impedances 1/3 (Beierlein)
vpm_fs_cp=0.5; 
vpm_fs_dge=vpm_l4_dge*ratio_fs_l4; %(Cruikshank and Connors)

% END TO BE CHANGED

% ****** simulation paramters *****
sim.fractional_system_size=0.25;
sim.duration = .5; % duration of simulation; in seconds
sim.tstep=0.1; % time steps, in ms
sim.times=0:sim.tstep:(sim.duration*1000); %in ms
sim.tsteps=length(sim.times);

% ****** presynaptic VPM neurons ****** 
net.vpm.n=round(250*sim.fractional_system_size); % number of VPM neurons
net.vpm.f=15; % mean rate of modulation of VPM neurons; i.e whisking
net.vpm.sr=50; %50; % mean spike rate
net.vpm.spikes=zeros(net.vpm.n, sim.tsteps);% initialize spike array; zero == no spikes; one == spike

% ****** set-up the L4 NETWORK ****** 
% l4 (i.e. stellate) neurons
net.l4.n=round(1600*sim.fractional_system_size); 
net.l4.th=-50; % ap threshold; in mv
net.l4.rs=-80; % reset after AP 
net.l4.el=-70; % resting potential
net.l4.ee=0;
net.l4.ei=-80;
net.l4.tau=10; % membrane time constant; in ms
net.l4.spikes=zeros(net.l4.n, sim.tsteps); % zero == no spikes; one == spike
net.l4.vm=zeros(net.l4.n, sim.tsteps)+net.l4.el; % in vm
net.l4.gi=zeros(net.l4.n, sim.tsteps); % in CHECK
net.l4.ge=zeros(net.l4.n, sim.tsteps); % in CHECK
net.l4.tge=5; % time constant of excitatory conductance; in ms
net.l4.tgi=5; % time constant of inihibitory conductance; in ms

% fs neurons
net.fs.n=round(200*sim.fractional_system_size);
net.fs.th=-50; % in mv LOWER THRESHOLD?
net.fs.rs=-80;
net.fs.el=-70;
net.fs.ee=0;
net.fs.ei=-80;
net.fs.tau=5; % in ms
net.fs.spikes=zeros(net.fs.n, sim.tsteps);
net.fs.vm=zeros(net.fs.n, sim.tsteps)+net.fs.el;
net.fs.gi=zeros(net.fs.n, sim.tsteps);
net.fs.ge=zeros(net.fs.n, sim.tsteps);
net.fs.tge=5; % in ms; note: assume that all excitatory conductnaces have the same tau
net.fs.tgi=5; % in ms; note: assume that all excitatory conductnaces have the same tau

% l4 - l4 connections; square connection matruc, but diagonals are zero
net.l4_l4.cp=l4_l4_cp; % connection probability
net.l4_l4.dge=l4_l4_dge; % in fractional units of driving force; i.e 
net.l4_l4.cm=rand(net.l4.n, net.l4.n);
net.l4_l4.cm(find(net.l4_l4.cm < net.l4_l4.cp))=1;
net.l4_l4.cm(find(net.l4_l4.cm < 1))=0;
net.l4_l4.cm(1:(net.l4.n+1):(net.l4.n*net.l4.n))=0;% set diagonal elements to zero using a(1:(n+1):n^2)=0

% l4 - fs connections; ; not square 
net.l4_fs.cp=l4_fs_cp;
net.l4_fs.dge=l4_fs_dge; % in mv; psps per ap
net.l4_fs.cm=rand(net.fs.n, net.l4.n);
net.l4_fs.cm(find(net.l4_fs.cm < net.l4_fs.cp))=1;
net.l4_fs.cm(find(net.l4_fs.cm < 1))=0;

% fs - l4 connections; ; not square
net.fs_l4.cp=fs_l4_cp;
net.fs_l4.dgi=fs_l4_dgi; % CHECK OUT
net.fs_l4.cm=rand(net.l4.n, net.fs.n);
net.fs_l4.cm(find(net.fs_l4.cm < net.fs_l4.cp))=1;
net.fs_l4.cm(find(net.fs_l4.cm < 1))=0;

% fs - fs connections; square connection matruc, but diagonals are zero
net.fs_fs.cp=fs_fs_cp;
net.fs_fs.dgi=fs_fs_dgi; % in nS; conductance per ap CHECK OUT
net.fs_fs.cm=rand(net.fs.n, net.fs.n);
net.fs_fs.cm(find(net.fs_fs.cm < net.fs_l4.cp))=1;
net.fs_fs.cm(find(net.fs_fs.cm < 1))=0;
net.fs_fs.cm(1:(net.fs.n+1):(net.fs.n*net.fs.n))=0;

% VPM - l4 connections; not square ...
net.vpm_l4.cp=vpm_l4_cp; % connection probability
net.vpm_l4.dge=vpm_l4_dge; % in  CHECK OUT
net.vpm_l4.dt=round(vpm_l4_dt*1000/sim.tstep); % in samples; models the fact that vpm-l4 synapses are slower than vpm-fs synapses
net.vpm_l4.cm=rand(net.l4.n, net.vpm.n);
net.vpm_l4.cm(find(net.vpm_l4.cm<net.vpm_l4.cp))=1;
net.vpm_l4.cm(find(net.vpm_l4.cm < 1))=0;

% generate VPM - fs connection matrix; not square ...
net.vpm_fs.cp=vpm_fs_cp; % connection probability
net.vpm_fs.dge=vpm_fs_dge; % in nS
net.vpm_fs.cm=rand(net.fs.n, net.vpm.n);
net.vpm_fs.cm(find(net.vpm_fs.cm<net.vpm_fs.cp))=1;
net.vpm_fs.cm(find(net.vpm_fs.cm < 1))=0;

global net sim
%testFunc
iniVPMspikes;
simulateL4network;


%% plot mean psth
h_f21=figure(21)
set(h_f21, 'Position', [440 75 560 725], 'Color',[1 1 1]);
delete(findobj(h_f21, 'type','axes'));
% 1 ms bins
l4=sum(reshape(mean(net.l4.spikes(:,1:(sim.tsteps-1))), round(1/sim.tstep), round(sim.tstep*sim.tsteps))); 
tt=1:length(l4);
fs=sum(reshape(mean(net.fs.spikes(:,1:(sim.tsteps-1))), round(1/sim.tstep), round(sim.tstep*sim.tsteps))); 
vpm=sum(reshape(mean(net.vpm.spikes(:,1:(sim.tsteps-1))), round(1/sim.tstep), round(sim.tstep*sim.tsteps))); 

xLim=[0 sim.tsteps*sim.tstep];
%yLim=[0 1];
subplot(3,1,1)
set(gca, 'XLim', xLim)
line(tt, vpm, 'Color', 'k')
title('VPM')
subplot(3,1,2)
set(gca, 'XLim', xLim)
line(tt, l4, 'Color', 'k')
ylabel('Spike rate (1/ms)')
title('l4 stellate')
subplot(3,1,3)
set(gca, 'XLim', xLim)
line(tt, fs, 'Color', 'k')
xlabel('Time (ms)')
title('l4 fs')

% text
subplot(3,1,1)
xLim=get(gca, 'XLim');
yLim=get(gca, 'YLim');
text((xLim(2)-xLim(1))*0.8, (yLim(2)-yLim(1)), strcat('v-4-cp=', num2str(vpm_l4_cp)), 'FontSize', 10, 'BackgroundColor', 'w');
text((xLim(2)-xLim(1))*0.8, (yLim(2)-yLim(1))*0.85, strcat('v-4-dge=', num2str(vpm_l4_dge)), 'FontSize', 10, 'BackgroundColor', 'w');
text((xLim(2)-xLim(1))*0.8, (yLim(2)-yLim(1))*0.70, strcat('v-f-cp=', num2str(vpm_fs_cp)), 'FontSize', 10, 'BackgroundColor', 'w');
text((xLim(2)-xLim(1))*0.8, (yLim(2)-yLim(1))*0.55, strcat('v-f-dge=', num2str(vpm_fs_dge)), 'FontSize', 10, 'BackgroundColor', 'w');

subplot(3,1,2)
xLim=get(gca, 'XLim');
yLim=get(gca, 'YLim');
text((xLim(2)-xLim(1))*0.8, (yLim(2)-yLim(1)), strcat('4-4-cp=', num2str(l4_l4_cp)), 'FontSize', 10, 'BackgroundColor', 'w');
text((xLim(2)-xLim(1))*0.8, (yLim(2)-yLim(1))*0.85, strcat('4-4-dge=', num2str(l4_l4_dge)), 'FontSize', 10, 'BackgroundColor', 'w');
text((xLim(2)-xLim(1))*0.8, (yLim(2)-yLim(1))*0.70, strcat('f-4-cp=', num2str(fs_l4_cp)), 'FontSize', 10, 'BackgroundColor', 'w');
text((xLim(2)-xLim(1))*0.8, (yLim(2)-yLim(1))*0.55, strcat('f-4-dgi=', num2str(fs_l4_dgi)), 'FontSize', 10, 'BackgroundColor', 'w');

subplot(3,1,3)
xLim=get(gca, 'XLim');
yLim=get(gca, 'YLim');
text((xLim(2)-xLim(1))*0.8, (yLim(2)-yLim(1)), strcat('f-f-cp=', num2str(fs_fs_cp)), 'FontSize', 10, 'BackgroundColor', 'w');
text((xLim(2)-xLim(1))*0.8, (yLim(2)-yLim(1))*0.85, strcat('f-f-dgi=', num2str(fs_fs_dgi)), 'FontSize', 10, 'BackgroundColor', 'w');


%% plot all rasters
h_f41=figure(41);
set(h_f41, 'Position', [440 75 560 725], 'Color',[1 1 1]);
%set(h_f41, 'Color',[1 1 1]);
delete(findobj(h_f41, 'type','axes'));

xLim=[0 sim.tsteps*sim.tstep];

nn=net.vpm.n; spks=net.vpm.spikes;
subplot(3,1, 1)
set(gca, 'XLim', xLim)
title('VPM')
for i=1:nn
    temp=find(spks(i, :)==1);
    line(temp*sim.tstep, i*ones(1, length(temp)), 'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor','k','MarkerEdgeColor','k', 'LineStyle', 'none')
end

nn=net.l4.n; spks=net.l4.spikes;
subplot(3,1, 2)
set(gca, 'XLim', xLim)
ylabel('Neurons')
title('l4 stellate')
for i=1:nn
    temp=find(spks(i, :)==1);
    line(temp*sim.tstep, i*ones(1, length(temp)), 'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor','k','MarkerEdgeColor','k', 'LineStyle', 'none')
end

nn=net.fs.n; spks=net.fs.spikes;
subplot(3,1, 3)
set(gca, 'XLim', xLim)
xlabel('Time (ms)')
title('layer fs')
for i=1:nn
    temp=find(spks(i, :)==1);
    line(temp*sim.tstep, i*ones(1, length(temp)), 'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor','k','MarkerEdgeColor','k', 'LineStyle', 'none')
end

% text
subplot(3,1,1)
xLim=get(gca, 'XLim');
yLim=get(gca, 'YLim');
text((xLim(2)-xLim(1))*0.8, (yLim(2)-yLim(1)), strcat('v-4-cp=', num2str(vpm_l4_cp)), 'FontSize', 10, 'BackgroundColor', 'w');
text((xLim(2)-xLim(1))*0.8, (yLim(2)-yLim(1))*0.85, strcat('v-4-dge=', num2str(vpm_l4_dge)), 'FontSize', 10, 'BackgroundColor', 'w');
text((xLim(2)-xLim(1))*0.8, (yLim(2)-yLim(1))*0.70, strcat('v-f-cp=', num2str(vpm_fs_cp)), 'FontSize', 10, 'BackgroundColor', 'w');
text((xLim(2)-xLim(1))*0.8, (yLim(2)-yLim(1))*0.55, strcat('v-f-dge=', num2str(vpm_fs_dge)), 'FontSize', 10, 'BackgroundColor', 'w');

subplot(3,1,2)
xLim=get(gca, 'XLim');
yLim=get(gca, 'YLim');
text((xLim(2)-xLim(1))*0.8, (yLim(2)-yLim(1)), strcat('4-4-cp=', num2str(l4_l4_cp)), 'FontSize', 10, 'BackgroundColor', 'w');
text((xLim(2)-xLim(1))*0.8, (yLim(2)-yLim(1))*0.85, strcat('4-4-dge=', num2str(l4_l4_dge)), 'FontSize', 10, 'BackgroundColor', 'w');
text((xLim(2)-xLim(1))*0.8, (yLim(2)-yLim(1))*0.70, strcat('f-4-cp=', num2str(fs_l4_cp)), 'FontSize', 10, 'BackgroundColor', 'w');
text((xLim(2)-xLim(1))*0.8, (yLim(2)-yLim(1))*0.55, strcat('f-4-dgi=', num2str(fs_l4_dgi)), 'FontSize', 10, 'BackgroundColor', 'w');

subplot(3,1,3)
xLim=get(gca, 'XLim');
yLim=get(gca, 'YLim');
text((xLim(2)-xLim(1))*0.8, (yLim(2)-yLim(1)), strcat('f-f-cp=', num2str(fs_fs_cp)), 'FontSize', 10, 'BackgroundColor', 'w');
text((xLim(2)-xLim(1))*0.8, (yLim(2)-yLim(1))*0.85, strcat('f-f-dgi=', num2str(fs_fs_dgi)), 'FontSize', 10, 'BackgroundColor', 'w');


%% PLOT a few cells - subthreshold
h_f11=figure(11);
set(h_f11, 'Color',[1 1 1]);
nn=10
%vm=net.fs.vm; ge=net.fs.vm; gi=net.fs.gi;
vm=net.l4.vm;ge=net.l4.ge; gi=net.l4.gi;

%figure;
subplot(3,1,1)
dd=std(vm(1, :));
for i=1:nn
line(sim.times,vm(i, :)+dd*(i-1))
end
%plot(vm(1, :))

subplot(3,1,2)
dd=std(ge(1, :));
for i=1:nn
line(sim.times,ge(i, :)+dd*(i-1))
end
%plot(ge(1, :))
subplot(3,1,3)
dd=std(gi(1, :));
for i=1:nn
line(sim.times,gi(i, :)+dd*(i-1))
end
%plot(gi(1, :))

%% various stats

vpm_sr = mean(mean(net.vpm.spikes))*sim.tsteps;
fprintf('vpm spike rate = %d\n', vpm_sr);


l4_sr = mean(mean(net.l4.spikes))*sim.tsteps;
fprintf('l4 ss spike rate = %d\n', l4_sr);

fs_sr = mean(mean(net.fs.spikes))*sim.tsteps;
fprintf('fs spike rate = %d\n', fs_sr);

l4_vm=mean(mean(net.l4.vm))
fs_vm=mean(mean(net.fs.vm))

l4_spikes=sum(sum(net.l4.spikes))
fs_spikes=sum(sum(net.fs.spikes))

%% print figure to pdf file

export_fig('fileName.pdf',-pdf');





%% old



%% plot one raster
h_f31=figure(31);
delete(findobj(h_f31, 'type','axes'));
nn=net.vpm.n; spks=net.vpm.spikes;
%nn=net.l4.n; spks=net.l4.spikes;
nn=net.fs.n; spks=net.fs.spikes;

xLim=[0 sim.tsteps*sim.tstep];
%figure;
ha=axes;
set(ha, 'XLim', xLim)

for i=1:nn
    temp=find(spks(i, :)==1);
    line(temp*tstep, i*ones(1, length(temp)), 'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor','k','MarkerEdgeColor','k', 'LineStyle', 'none')
end