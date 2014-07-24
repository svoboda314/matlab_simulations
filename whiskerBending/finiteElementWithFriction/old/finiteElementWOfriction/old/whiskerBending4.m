
% Whisker would point to a range of angles if it weren't deflected by an
% object 
% 
% Deflections are computed for n_iteration positions
%
% parameters are initialized in the script whiskerSimParams
%
% KS 111409

whiskerSimParams; % initialize parameters regarding whisker shape and simulation
wh.sim.x_pole=0.020; % x position of pole

n_iterations=9;
max_base_angle=wh.sim.base_angle;
% initialize arrays for 
% n_iterations whisker trajectories
trajectories=zeros(n_iterations, 2, wh.sim.n_segments+1);
axial_force=zeros(1, n_iterations);
moment=zeros(1, n_iterations);
% n_iterations axial forces
% n_iterations moments

for i=1:n_iterations
    wh.sim.base_angle = i*max_base_angle/n_iterations;    
    vibr=whiskerDeflectionByObject(wh);
    trajectories(i, 1, :)=vibr.prop.x;
    trajectories(i, 2, :)=vibr.prop.y;
    axial_force(i)=vibr.prop.moment;
    moment(i)=vibr.prop.axial_force;
end

figure; 
axes('XLim', [0 1000*wh.prop.L], 'YLim', 1000*[-0.5*wh.prop.L 0.5*wh.prop.L])
xlabel('x-distance (m)')
ylabel('y-distance (m)')
title_str=strcat('max angle= ', num2str(max_base_angle), '; x-pole = ', num2str(wh.sim.x_pole))
title(title_str)
line (1000*squeeze(trajectories(:, 1, :))',-1000*squeeze(trajectories(:, 2, :))', 'Color', 'k');

figure; 
subplot(2,2,1)
xlabel('whisker angle at follicle (degrees)')
ylabel('axial_force (N)')
line((1:n_iterations)*max_base_angle/n_iterations, axial_force)
title_str=strcat('x-pole = ', num2str(wh.sim.x_pole));
title(title_str)

subplot(2,2,2)
xlabel('whisker angle at follicle (degrees)')
ylabel('moment (Nm)')
line((1:n_iterations)*max_base_angle/n_iterations, moment)
title_str=strcat('x-pole = ', num2str(wh.sim.x_pole));
title(title_str)

subplot(2,2,3)
ylabel('axial_force (N)')
xlabel('moment (Nm)')
line(moment, axial_force)
title_str=strcat('x-pole = ', num2str(wh.sim.x_pole));
title(title_str)