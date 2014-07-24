
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
end
subplot(nplotRows, 2, (i*2+1))
bar(mean(sr_f));

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
end
subplot(nplotRows, 2, (i*2+2))
bar(mean(sr_2f));

sr_grand= (mean(sr_2f)+mean(sr_f))/2;

