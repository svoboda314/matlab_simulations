
global net sim
%testFunc
iniVPMspikesReg %regular spikes
% iniVPMspikes;
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
%set(h_f11, 'Color',[1 1 1]);
set(h_f11, 'Position', [440 75 560 725], 'Color',[1 1 1]);
delete(findobj(h_f11, 'type','axes'));

nn=10;
%vm=net.fs.vm; ge=net.fs.vm; gi=net.fs.gi;
vm=net.l4.vm;ge=net.l4.ge; gi=net.l4.gi;

%figure;
subplot(3,1,1)
title('l4 ss Vm')
ylabel('Vm (mV)')
dd=std(vm(1, :));
for i=1:nn
line(sim.times,vm(i, :)+dd*(i-1))
end
%plot(vm(1, :))

subplot(3,1,2)
title('ge onto l4 ss')
dd=std(ge(1, :));
for i=1:nn
line(sim.times,ge(i, :)+dd*(i-1))
end
%plot(ge(1, :))
subplot(3,1,3)
xlabel('Time (ms)')
title('gi onto l4 ss')
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