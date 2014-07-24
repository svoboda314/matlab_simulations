
function out=simulateL4network
% Basic dynamics is as in, for example, Dayan & Abbott pg 
%
% tau * dV/dt = El -V - g rm (V - Er); 
% g is a specific membrane conductance with reversal potential Er; separate
% terms are added for different sources of synaptic input 
% rm is the specific membrane resistivity due to leak
% 
% the conductance is parametrized as gi and ge, which are potential 
% changes normalized to the driving force
% 
% Assumptions:
% single compartment IF neurons
% random connectivity across neuron types
% 
% all excitatory / inhibitory conductances impinging on a particular type of neuron
% share dynamics


global net sim

for t=(2+net.vpm_l4.dt):sim.tsteps % 011213 net.vpm_l4.dt because of possibility that L4 input is delayed
    % UPDATE net.l4.vm
    net.l4.vm(:,t)=net.l4.vm(:,t-1)+(net.l4.el-net.l4.vm(:,t-1))*sim.tstep/net.l4.tau - ...
        net.l4.gi(:,t-1).*(net.l4.vm(:,t-1) - net.l4.ei)*sim.tstep/net.l4.tau - ...
        net.l4.ge(:,t-1).*(net.l4.vm(:,t-1) - net.l4.ee)*sim.tstep/net.l4.tau; 
    % UPDATE net.l4.gi
    net.l4.gi(:,t)=net.l4.gi(:,t-1)*(1-sim.tstep/net.l4.tgi)+net.fs_l4.dgi*net.fs_l4.cm*net.fs.spikes(:, t-1);
    % UPDATE net.l4.ge
    net.l4.ge(:,t)=net.l4.ge(:,t-1)*(1-sim.tstep/net.l4.tge)+net.l4_l4.dge*net.l4_l4.cm*net.l4.spikes(:, t-1)+...
        net.vpm_l4.dge*net.vpm_l4.cm*net.vpm.spikes(:, t-1-net.vpm_l4.dt);% 011213
    % UPDATE net.fs.vm
    net.fs.vm(:,t)=net.fs.vm(:,t-1)+(net.fs.el-net.fs.vm(:,t-1))*sim.tstep/net.fs.tau - ...
        net.fs.gi(:,t-1).*(net.fs.vm(:,t-1) - net.fs.ei)*sim.tstep/net.fs.tau - ...
        net.fs.ge(:,t-1).*(net.fs.vm(:,t-1) - net.fs.ee)*sim.tstep/net.fs.tau; 
    % UPDATE net.fs.gi
    net.fs.gi(:,t)=net.fs.gi(:,t-1)*(1-sim.tstep/net.fs.tgi)+net.fs_fs.dgi*net.fs_fs.cm*net.fs.spikes(:, t-1);
    % UPDATE net.fs.ge
    net.fs.ge(:,t)=net.fs.ge(:,t-1)*(1-sim.tstep/net.fs.tge)+net.l4_fs.dge*net.l4_fs.cm*net.l4.spikes(:, t-1)+...
        net.vpm_fs.dge*net.vpm_fs.cm*net.vpm.spikes(:, t-1);
    % RESET net.l4.vm
    temp_l4=find(net.l4.vm(:, t-1) > net.l4.th); % detect L4 spikes
    net.l4.spikes(temp_l4, t)=1;
    net.l4.vm(temp_l4,t)=net.l4.rs;
    % RESET net.fs.vm
    temp_fs=find(net.fs.vm(:, t-1) > net.fs.th); % detect L4 spikes
    net.fs.spikes(temp_fs, t)=1;
    net.fs.vm(temp_fs,t)=net.fs.rs;
end

out=1;