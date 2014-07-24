
function sp_times  = inhPoisson(rate, duration)
% inhomogeneous Poisson train computed by spike thinning
% usage: sp_times = inhPoisston(rate, duration)
% time varying rate, in Hz; specified every ms; length duration*rate; should be less than
%   1000 Hz
% duration in seconds
% output in ms
% ks 11/17/2013

max_r = max(rate); % in Hz
dt = 0.0001; % in seconds
temp=rand(round(duration/dt), 1); %
spikeTimes = find(temp < dt*max_r)*dt; %spike times at max rate

while((round(spikeTimes(end)*1000)+1) > length(rate))
    spikeTimes=spikeTimes(1:(end-1));
end

for i=1:length(spikeTimes)
    inst_rate=rate(round(spikeTimes(i)*1000)+1);
    if rand > inst_rate/max_r; spikeTimes(i)=nan; end % spike thinning
end

sp_times=round(spikeTimes(~isnan(spikeTimes))*1000);
%sp_train=zeros(1,duration*1000);
%sp_train(sp_times)=1; % FIX -- 0 index
