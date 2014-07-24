
function sp_times  = inhPoisson(rate)
% inhomogeneous Poisson train computed by spike thinning
% usage: sp_times = inhPoisston(rate, duration)
% time varying rate, in Hz; specified every ms; max < 1000 
% length of spike train is length(rate)*ms
%
% output spike times, in ms
% ks 11/17/2013

max_r = max(rate); % in Hz

%dt = 0.0001; % in seconds
dt= 0.1; %in ms

%temp=rand(round(duration/dt), 1); %

temp=rand(length(rate)*10, 1); % random % for each 0.1 ms

spikeTimes = round(find(temp < dt*max_r/1000)*dt); %spike times at max rate
spikeTimes=spikeTimes(find(spikeTimes > 0)); 

%while((round(spikeTimes(end)*1000)+1) > length(rate))
%    spikeTimes=spikeTimes(1:(end-1));
%end

for i=1:length(spikeTimes)
    inst_rate=rate(spikeTimes(i));
%    inst_rate=rate(round(spikeTimes(i)*1000)+1);
    if rand > inst_rate/max_r; spikeTimes(i)=nan; end % spike thinning
end

sp_times=round(spikeTimes(~isnan(spikeTimes)));
%sp_train=zeros(1,duration*1000);
%sp_train(sp_times)=1; % FIX -- 0 index