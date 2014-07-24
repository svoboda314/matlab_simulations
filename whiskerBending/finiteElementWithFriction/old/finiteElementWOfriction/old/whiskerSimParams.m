
% Whisker would point to angle base-angle if it weren't deflected by an
% object located at (x_pole, y_pole)
% 
% The deflection is computed iterativley
%
% r - vector from current length element and point of force application
% beta - angle between vertical and r
% force_angle - angle between force and vertical axis

% --------------setting whisker parameters and initialize other params

wh.prop.conical=1;
wh.prop.ym=3.5e9; % young's modulus in pascal
wh.prop.base_radius=1.00e-4; % whisker radius at the base, in meters
wh.prop.L=0.060; % whisker length in meters
wh.sim.n_segments=1000; % # of segments for finite element analysis
wh.sim.segment=wh.prop.L/wh.sim.n_segments; % in meters; total length is n_segments * segment

% Depending on whisker geometry, compute moment of inertia
if ~ wh.prop.conical
    I(1:wh.sim.n_segments)=(pi/4)*wh.prop.base_radius^4; %uniform
else
    I(1:wh.sim.n_segments)=(pi/4)*(wh.prop.base_radius*(1-(1:wh.sim.n_segments)*wh.sim.segment/wh.prop.L)).^4;  %conical
end

wh.prop.I=I;
wh.sim.base_angle = 45; % angle in degrees at the base of the whisker

wh.sim.x_pole=0.020; % x position of the object, in meters
wh.sim.y_pole=0.00; % y position of the object, in meters
wh.sim.delta_force=1/300000; % force increment per iteration in Newtons
