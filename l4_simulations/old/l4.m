
% ****** general simulation parameters****** 
duration = 1; % in seconds
tstep=1; %in ms
t=0:1:(duration*1000); %in ms
tsteps=length(t);

% ****** presynaptic VPM neurons ****** 
net.vpm.n=50; % number of VPM neurons
net.vpm.f=15; % frequency of VPM neurons
net.vpm.sr=50; % spike rate
net.vpm_L4.cp=0.5; % connection probability

% ****** generate VPM activity ...****** 

for i=1:net.vpm.n
    rate=net.vpm.sr*(1+sin(net.vpm.f*2*pi*t/1000+rand*2*pi));
    net.vpm.st{i}=inhPoisson(rate, duration);
    net.vpm.st{i}=[net.vpm.st{i}', 500, 501]';
end

% add synchonized spikes at t=500

binsize = 0.01;
bins=(binsize/2):binsize:(duration-binsize/2);
hist_bins=zeros(1, length(bins));
%hist_bins=hist_bins+hist(net.vpm.st{i}, bins)';
for i=1:net.vpm.n
    hist_bins=hist_bins+hist(net.vpm.st{i}, bins);
end
hist_bins=hist_bins/net.vpm.n;

% ****** L4 NETWORK ****** 
% postsynaptic L4 neurons
net.l4.n=250; 
net.L4.th=20; % in mv
net.L4.tau=10; % in ms
net.L4.psp=2; % in mv
net.L4.spikes=zeros(net.l4.n, tsteps); % zero == no spikes; one == spike
%net.l4.vm=zeros(net.l4.n, 1); % in mv
net.L4_L4.cp=0.5;
net.L4_L4.psp=1;

% generate VPM - L4 connection matrix; not square ...
net.vpm_L4.cm=rand(net.l4.n, net.vpm.n);
net.vpm_L4.cm(find(net.vpm_L4.cm>net.vpm_L4.cp))=1;
net.vpm_L4.cm(find(net.vpm_L4.cm < 1))=0;
net.vpm_L4.psp=1;

% generate L4 - L4 connection matrix; square, but diagonals are zero
net.L4_L4.cm=rand(net.l4.n, net.l4.n);
net.L4_L4.cm(find(net.L4_L4.cm > net.L4_L4.cp))=1;
net.L4_L4.cm(find(net.L4_L4.cm < 1))=0;
net.L4_L4.cm(1:(net.l4.n+1):(net.l4.n*net.l4.n))=0; % set diagonal elements to zero using a(1:(n+1):n^2)=0

%sum(sum(net.L4_L4.cm))
%sum(diag(net.L4_L4.cm)) % needs to be zero

vm=zeros(net.l4.n, tsteps);
spikes=zeros(net.l4.n, tsteps);
vpm_spikes=zeros(1, net.vpm.n);
l4_spikes=zeros(1, net.l4.n);

for i=2:tsteps
    temp=find(vm(:, i-1) > net.L4.th); % detect spikes
%    vm(temp, i)=0;       % reset to 0   
    net.L4.spikes(temp, i-1)=1; % track spikes
    for j=1:net.vpm.n
        vpm_spikes(j)=length(find(net.vpm.st{j}==i)); % spikes in VPM at time=i
    end
    vm(:, i)=vm(:, i-1)*exp(-tstep/net.L4.tau)+net.vpm_L4.psp*net.vpm_L4.cm*vpm_spikes'+ ...
        net.L4_L4.psp*net.L4_L4.cm*net.L4.spikes(:, i);
%    temp=find(vm(:, i) > net.L4.th); % detect spikes
    vm(temp, i)=0;       % reset to 0   
%    net.L4.spikes(temp, i)=1; % track spikes
end

% ****** Plots ****** 
% raster
hf1=figure('Units', 'normalized', 'Position', [0.5 0.1 0.4 0.8])
har=axes('Units', 'normalized','Position', [0.1 0.6 0.8 0.35])

for i=1:net.l4.n
    temp=find(net.L4.spikes(i, :)==1);
    line(temp, i*ones(1, length(temp)), 'Marker', 'o', 'LineStyle', 'none')
end

% psth
hapsth=axes('Units', 'normalized','Position', [0.1 0.33 0.8 0.21])
line(t, sum(net.L4.spikes))

% mean vm
hapsvm=axes('Units', 'normalized','Position', [0.1 0.05 0.8 0.21])
line(t, mean(vm))

% a few vm traces
hf2=figure('Units', 'normalized', 'Position', [0.5 0.1 0.4 0.8])
havm=axes('Units', 'normalized','Position', [0.1 .1 0.8 0.8])

for i=1:10
line(t,vm(i, :)+15*(i-1))
end

