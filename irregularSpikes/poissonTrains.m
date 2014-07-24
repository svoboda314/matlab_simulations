
%% Fano factor

isi=1;
T=10;

mu = T/isi;

spikeCounts=random('poiss', mu,10000, 1);
meanSK=mean(spikeCounts);
stdSK=std(spikeCounts);
ff=stdSK^2/meanSK

%% With spike times
isi=0.25;
isiS= random('exp', isi,10000,1);
%spikeT=cumsum(y);
%plot(spikeT, ones(length(spikeT),1), 'o')
x=0.05:.1:1.95;
ct=hist(isiS,x);
plot(x,ct, '-o')
cv=std(isiS)/mean(isiS)

%tend=max(spike(T));

%% fitting exponential distributions
[curve, goodness] = fit( x', ct', 'exp1' );


%% plotting a rasster of Poisson points
trials = 10;
isi=0.2;
nPoints=100;

figure;
hold
for i=1:trials
    isiS= random('exp', isi, nPoints ,1);
    spikeT=cumsum(isiS);
    plot(spikeT, ones(1, nPoints)*i, 'o')
end

%% modulated poisson trains

time = 1; %in seconds
tstep = 0.001 % in seconds
npts=time/tstep; 
tpts=0:tstep:time;
sr= sin(2*pi*tpts * 15)+0.5;







%% psth of rate-modulated spike trains

time = 0.1; %in seconds
tstep = 0.001 % in seconds
npts=time/tstep+1; 
tpts=0:tstep:time;
freq= 15; % in Hz
mean_sr= 100*(.1*sin(2*pi*tpts * freq)+0.5);

% sr1=zeros(1,100);
% for i=1:npts
%     sr1(i)=random('Poisson', mean_sr(i));
% end
sr=random('Poisson', mean_sr);


% neurons that are modulated at freq
% frequency is freq
% phase phi is random [0 pi]
% amplitude modulation is random [0 1]

n_f = 10; %# of neurons
sr_f=zeros(n_f, npts);
touch_r = 100;

figure
nplotRows=n_f+2;
subplot(nplotRows, 2, 1)

for i=1:n_f
    mean_sr= 100*(rand*sin(2*pi*tpts * freq+ 2*rand*pi)+1);
    sr_f(i, :)=random('Poisson', mean_sr);
    sr_f(i, 55) = sr_f(i, 55) + touch_r * rand;
    subplot(nplotRows, 2, (i*2-1))
    bar(sr_f(i, :));
    set(gca, 'YLim', [0 300], 'XTickLabel', [], 'YTickLabel', []);
end
subplot(nplotRows, 2, (i*2+1))
bar(mean(sr_f));
set(gca, 'YLim', [0 300], 'XTickLabel', [], 'YTickLabel', []);
% neurons that are modulated at 2xfreq
% frequency is 2xfreq
% phase phi is random [0 pi]
% amplitude modulation is random [0 1]

n_2f = 10; %# of neurons
sr_2f=zeros(n_2f, npts);

for i=1:n_2f
    mean_sr= 100*(rand*sin(2*pi*tpts * 2*freq+ 2*rand*pi)+1);
    sr_2f(i, :)=random('Poisson', mean_sr);
    sr_2f(i, 55) = sr_2f(i, 55) + touch_r * rand;
    subplot(nplotRows, 2, (i*2))
    bar(sr_2f(i, :));
    set(gca, 'YLim', [0 300], 'XTickLabel', [], 'YTickLabel', []);
end
subplot(nplotRows, 2, (i*2+2))
%bar(mean(sr_2f));
plot(0:100, sin((0:100)*15*2*pi));
set(gca, 'YLim', [-1 1], 'XTickLabel', [], 'YTickLabel', []);

sr_grand= (mean(sr_2f)+mean(sr_f))/2;

%% rate-modulated Poisson train

% first a homogeneous Poisson train
t_step = 0.0001; % time step
t_len = 10; % length of trace
len = t_len/t_step; % # of smaples
rate = 1; % spike rate

temp=rand(1,len); 

%ind=(find(temp < tstep*rate*0.1)); 
spikes=zeros(1, len);
%hist(temp)
ind=(find(temp < rate*t_step));
spike_times = ind*t_step;
isi=diff(spike_times);
%hist(isi)
spikes(ind)=1;
plot(spike_times, ones(1, length(spike_times)), 'o', 'Markersize', 10)
%set(gco, 'Markersize', 10)
hold
plot([-1 0 0 10 10 11], [0 0 .1 .1 0 0], 'k')
set(gco, 'Linewidth', 1.5)
set(gca, 'YLim', [-.03 1.1])
hold off

%plot(spikes, 'o')
sum(spikes)

% now an inhomogeneous process 

x=[-1];
y=[0];
for i=1: length(spike_times)
    x=[x, spike_times(i)-.001, spike_times(i)-.001, spike_times(i), spike_times(i)];
    y=[y, 0, 0.9, 0.9, 0];
end
x=[x, 11];
y=[y, 0];

figure
plot(spike_times, ones(1, length(spike_times)), 'o', 'Markersize', 10)
%set(gco, 'Markersize', 10)
hold
plot(x, y, 'k')
set(gco, 'Linewidth', 1.5)
set(gca, 'YLim', [-.03 1.1])
hold off



%% try inhPoisson, spike thinning method
duration = 10;
t=0:0.001:duration;
rate=100*sin(t*pi/duration);
duration=10+0.001;

sp_times=inhPoisson(rate, duration);

% compute spike count in 100 ms bins
binsize = 0.1;
binnums=round(duration/binsize);
sp_rate=zeros(binnums,1);

tbins=(binsize/2):0.1:(duration-binsize/2);

n=hist(sp_times, tbins);

figure;
%tbins=(binsize/2):0.1:(duration-binsize/2);
bar(tbins, n)
set(gca, 'Xlim', [0 10])

%% one Poisson neuron onto one IF neuron

duration = 10;
tstep=1; %in ms
t=0:0.001:duration;
tsteps=length(t);
sp_rate=100; % in Hz
rate=ones(1, tsteps)*sp_rate;
duration=10+0.001;
epsp = 0.1; % in Vth units
tau = 100; % in seconds

sp_times=inhPoisson(rate, duration);
sp_times=round(sp_times*1000);

vm=zeros(1, tsteps);
sp_out=[];
for i=2:tsteps
    dv=vm(i-1)*(1-exp(-tstep/tau));
    vm(i)=vm(i-1)-dv;
    if ismember(i-1, sp_times)
        vm(i)=vm(i)+epsp;
    end
    if vm(i) > 1
    sp_out=[sp_out, i];
    vm(i)=0;
    end
end

plot(vm)
hold
plot(sp_out, ones(1, length(sp_out)), 'o')

%% n Poisson neuron onto one IF neuron

n_neur = 20;

duration = 10;
tstep=1; %in ms
t=0:1:(duration*1000); %in ms
tsteps=length(t);
sp_rate=40; % in Hz
rate=ones(1, tsteps)*sp_rate;
duration=duration+0.001; %in seconds
epsp = 0.05; % in Vth units
tau = 1000; % in milliseconds

sp_times = cell(1,n_neur); 
for i=1:n_neur
temp=inhPoisson(rate, duration);
temp=round(temp*1000);
sp_times{i}=temp; % in ms
end

vm=zeros(1, tsteps);
sp_out=[];
for i=2:tsteps
    dv=vm(i-1)*(1-exp(-tstep/tau));
    vm(i)=vm(i-1)-dv;
    for j=1:n_neur
    if ismember(i-1, sp_times{j})
        vm(i)=vm(i)+epsp;
    end
    end
    if vm(i) > 1
    sp_out=[sp_out, i];
    vm(i)=0;
    end
end

plot(vm)
hold
plot(sp_out, ones(1, length(sp_out)), 'o')
hold
figure;
hist(diff(sp_out), 5:10:500)

cv= std(diff(sp_out))/mean(diff(sp_out))


spike_rate = length(sp_out)/duration