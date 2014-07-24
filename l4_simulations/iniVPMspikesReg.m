function out=iniVPMspikesReg

global net sim
%spike_rate=5;

nn=round(1000/(net.vpm.sr_reg*sim.tstep));
st=nn:nn:sim.tsteps;

net.vpm.spikes(:, st)=1;

out=1;