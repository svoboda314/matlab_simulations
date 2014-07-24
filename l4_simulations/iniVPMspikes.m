function out=iniVPMspikes

global net sim

% generate VPM activity ...
tt=0:(sim.duration*1000); % in ms 
touch_time=300;
spike_prob_on_touch=0.8;
A=1;
for i=1:net.vpm.n
    %i
    rate=net.vpm.sr*(1+A*sin(net.vpm.f*2*pi*tt/1000+rand*2*pi)); % CHANGE MODULATION DEPTH
    net.vpm.st{i}=inhPoisson(rate);
    if rand < spike_prob_on_touch
        net.vpm.st{i}=[net.vpm.st{i}', touch_time+round(rand)]';% add synchonized spikes at t=touch_time AD HOC - FIX - Jitter on the order of 5 ms
    end
    net.vpm.spikes(i, (1/sim.tstep)*setdiff(net.vpm.st{i}, [0]))=1;
end

out=1;