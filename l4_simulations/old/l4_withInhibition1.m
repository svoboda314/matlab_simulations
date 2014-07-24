
% ****** general simulation parameters****** 
duration = 1; % dureation of simulation; in seconds
tstep=1; % time steps, in ms
t=0:1:(duration*1000); %in ms
tsteps=length(t);

% ****** presynaptic VPM neurons ****** 
net.vpm.n=100; % number of VPM neurons
net.vpm.f=15; % frequency of VPM neurons
net.vpm.sr=50; % spike rate
net.vpm.spikes=zeros(net.vpm.n, tsteps);

% generate VPM activity ...

for i=1:net.vpm.n
    rate=net.vpm.sr*(1+sin(net.vpm.f*2*pi*t/1000+rand*2*pi));
    net.vpm.st{i}=inhPoisson(rate, duration);
    net.vpm.st{i}=[net.vpm.st{i}', 500, 501]';% add synchonized spikes at t=500 AD HOC - FIX
    net.vpm.spikes(i, setdiff(net.vpm.st{i}, [0]))=1;
end

% todo1 - pull out synchorny

% binsize = 0.01;
% bins=(binsize/2):binsize:(duration-binsize/2);
% hist_bins=zeros(1, length(bins));
% %hist_bins=hist_bins+hist(net.vpm.st{i}, bins)';
% for i=1:net.vpm.n
%     hist_bins=hist_bins+hist(net.vpm.st{i}, bins);
% end
% hist_bins=hist_bins/net.vpm.n;

% ****** set-up the L4 NETWORK ****** 
% l4 (i.e. stellate) neurons
net.l4.n=250; 
net.l4.th=30; % in mv
net.l4.tau=5; % in ms
net.l4.spikes=zeros(net.l4.n, tsteps); % zero == no spikes; one == spike
net.l4.vm=zeros(net.l4.n, tsteps); % in vm
net.l4.g=zeros(net.l4.n, tsteps); % in nS
net.l4.tau_g=50; % in ms

% fs neurons
net.fs.n=50;
net.fs.th=20; % in mv LOWER THRESHOLD?
net.fs.tau=5; % in ms
net.fs.spikes=zeros(net.fs.n, tsteps);
net.fs.vm=zeros(net.fs.n, tsteps);
net.fs.g=zeros(net.fs.n, tsteps);
net.fs.tau_g=50; % in ms

% l4 - l4 connections; square connection matruc, but diagonals are zero
net.l4_l4.cp=0.5; % connection probability
net.l4_l4.psp=.1;
net.l4_l4.cm=rand(net.l4.n, net.l4.n);
net.l4_l4.cm(find(net.l4_l4.cm < net.l4_l4.cp))=1;
net.l4_l4.cm(find(net.l4_l4.cm < 1))=0;
net.l4_l4.cm(1:(net.l4.n+1):(net.l4.n*net.l4.n))=0; % set diagonal elements to zero using a(1:(n+1):n^2)=0

% l4 - fs connections; ; not square 
net.l4_fs.cp=0;
net.l4_fs.psp=1; % in mv; psps per ap
net.l4_fs.cm=rand(net.fs.n, net.l4.n);
net.l4_fs.cm(find(net.l4_fs.cm < net.l4_fs.cp))=1;
net.l4_fs.cm(find(net.l4_fs.cm < 1))=0;

% fs - l4 connections; ; not square
net.fs_l4.cp=.5;
net.fs_l4.s=0.01; % in nS; conductance per ap CHECK OUT
net.fs_l4.cm=rand(net.l4.n, net.fs.n);
net.fs_l4.cm(find(net.fs_l4.cm < net.fs_l4.cp))=1;
net.fs_l4.cm(find(net.fs_l4.cm < 1))=0;

% fs - fs connections; square connection matruc, but diagonals are zero
net.fs_fs.cp=0;
net.fs_fs.s=0.05; % in nS; conductance per ap CHECK OUT
net.fs_fs.cm=rand(net.fs.n, net.fs.n);
net.fs_fs.cm(find(net.fs_fs.cm < net.fs_l4.cp))=1;
net.fs_fs.cm(find(net.fs_fs.cm < 1))=0;
net.fs_fs.cm(1:(net.fs.n+1):(net.fs.n*net.fs.n))=0;

% VPM - l4 connections; not square ...
net.vpm_l4.cp=0.5; % connection probability
net.vpm_l4.psp=.5; % psp in mv
net.vpm_l4.cm=rand(net.l4.n, net.vpm.n);
net.vpm_l4.cm(find(net.vpm_l4.cm<net.vpm_l4.cp))=1;
net.vpm_l4.cm(find(net.vpm_l4.cm < 1))=0;

% generate VPM - fs connection matrix; not square ...
net.vpm_fs.cp=0.5; % connection probability
net.vpm_fs.psp=1; % psp in mv
net.vpm_fs.cm=rand(net.fs.n, net.vpm.n);
net.vpm_fs.cm(find(net.vpm_fs.cm<net.vpm_fs.cp))=1;
net.vpm_fs.cm(find(net.vpm_fs.cm < 1))=0;

% net.l4.vm=zeros(net.l4.n, tsteps);
% net.l4.g=zeros(net.l4.n, tsteps);
% %vm_l4=zeros(net.l4.n, tsteps);
% 
% net.fs.vm=zeros(net.fs.n, tsteps);
% net.fs.g=zeros(net.fs.n, tsteps);
%vm_fs=zeros(net.fs.n, tsteps);

%vpm_spikes=zeros(1, net.vpm.n);
%l4_spikes=zeros(1, net.l4.n);

% TO DO
% fs spikes
% fs input from VPM and L4
% inhbition onto L4 neurons
% etc

for i=2:tsteps
    % detect L4 spikes
    temp_l4=find(net.l4.vm(:, i-1) > net.l4.th); 
    net.l4.spikes(temp_l4, i-1)=1; % track spikes
    % detect fs spikes
    temp_fs=find(net.fs.vm(:, i-1) > net.fs.th); 
    net.fs.spikes(temp_fs, i-1)=1; % track spikes
    % compute inhibition onto l4 cells
    net.l4.g(:, i)=net.l4.g(:, i-1)*exp(-tstep/net.l4.tau_g)+net.fs_l4.s*net.fs_l4.cm*net.fs.spikes(:, i-1);
    % compute inhibition onto fs cells
    net.fs.g(:, i)=net.fs.g(:, i-1)*exp(-tstep/net.fs.tau_g)+net.fs_fs.s*net.fs_fs.cm*net.fs.spikes(:, i-1);
%     for j=1:net.vpm.n
%         vpm_spikes(j)=length(find(net.vpm.st{j}==i)); % spikes in VPM at time=i
%     end
    % update l4 voltage
    net.l4.vm(:, i)=net.l4.vm(:, i-1)*exp(-tstep/net.l4.tau)+net.vpm_l4.psp*net.vpm_l4.cm*net.vpm.spikes(:, i-1)+  ...% vpm_spikes'+ ...
        net.l4_l4.psp*net.l4_l4.cm*net.l4.spikes(:, i-1)-...
        net.l4.g(:, i-1).*net.l4.vm(:, i-1); %
    net.l4.vm(temp_l4, i)=0;       % reset to 0   
    %update fs voltage
    net.fs.vm(:, i)=net.fs.vm(:, i-1)*exp(-tstep/net.fs.tau)+net.vpm_fs.psp*net.vpm_fs.cm*net.vpm.spikes(:, i-1)+ ... %vpm_spikes'+ ...
        net.l4_fs.psp*net.l4_fs.cm*net.l4.spikes(:, i-1)-...
        net.fs.g(:, i-1).*net.fs.vm(:, i-1);
    net.fs.vm(temp_fs, i)=0;       % reset to 0   

%    net.L4.spikes(temp_l4, i)=1; % track spikes
end

% ****** Plots ****** 
% raster
xLim=[0 tsteps];
%xLim=[495 550];
hf1=figure('Units', 'normalized', 'Position', [0.5 0.1 0.4 0.8]);
har=axes('Units', 'normalized','Position', [0.1 0.6 0.8 0.35]);
set(har, 'XLim', xLim)

for i=1:net.l4.n
    temp=find(net.l4.spikes(i, :)==1);
    line(temp, i*ones(1, length(temp)), 'Marker', 'o', 'LineStyle', 'none')
end

for i=1:net.fs.n
    temp=find(net.fs.spikes(i, :)==1);
    line(temp, i*ones(1, length(temp)), 'Marker', 'o', 'LineStyle', 'none', 'Color', 'r')
end

% for i=1:net.fs.n
%      temp=find(net.fs.spikes(i, :)==1);
%      line(temp, i*ones(1, length(temp)), 'Marker', 'o', 'LineStyle', 'none')
%  end

% psth
hapsth=axes('Units', 'normalized','Position', [0.1 0.33 0.8 0.21]);
set(hapsth, 'XLim', xLim)
line(t, sum(net.l4.spikes))
line(t, sum(net.fs.spikes), 'Color', 'r')
line(t, sum(net.vpm.spikes), 'Color', 'g')

% mean vmL4
hapsvm_l4=axes('Units', 'normalized','Position', [0.1 0.05 0.8 0.21]);
set(hapsvm_l4, 'XLim', xLim)
line(t, mean(net.l4.vm))

% a few vmL4 traces
hf2=figure('Units', 'normalized', 'Position', [0.5 0.1 0.4 0.8]);
havm_l4=axes('Units', 'normalized','Position', [0.1 .1 0.8 0.8]);
set(havm_l4, 'XLim', xLim)

% i=1;
% line(t,net.l4.vm(i, :)+15*(i-1))
% line(temp, i*ones(1, length(temp)), 'Marker', 'o', 'LineStyle', 'none')
% temp=find(net.l4.spikes(i, :)==1);
% line(temp,   -100+i*ones(1, length(temp)), 'Marker', 'o', 'LineStyle', 'none')
% temp=find(net.fs.spikes(i, :)==1);
% line(temp, -200+i*ones(1, length(temp)), 'Marker', 'o', 'LineStyle', 'none', 'Color', 'r')


for i=1:10
line(t,net.l4.vm(i, :)+15*(i-1))
end

l4_vm=mean(mean(net.l4.vm))
fs_vm=mean(mean(net.fs.vm))

l4_spikes=sum(sum(net.l4.spikes))
fs_spikes=sum(sum(net.fs.spikes))

